
"""
Allocate data from generation_fuel_eia923 table to generator level.

Net electricity generation and fuel consumption are reported in mutiple ways in the EIA
923. The generation_fuel_eia923 table reports both generation and fuel consumption, and
breaks them down by plant, prime mover, and fuel. In parallel, the generation_eia923
table reports generation by generator, and the boiler_fuel_eia923 table reports fuel
consumption by boiler.

The generation_fuel_eia923 table is more complete, but the generation_eia923 +
boiler_fuel_eia923 tables are more granular. The generation_eia923 table includes only
~55% of the total MWhs reported in the generation_fuel_eia923 table.

This module estimates the net electricity generation and fuel consumption attributable
to individual generators based on the more expansive reporting of the data in the
generation_fuel_eia923 table. The main coordinating function here is
:func:`pudl.analysis.allocate_net_gen.allocate_gen_fuel_by_gen`.

The algorithm we're using assumes:

* The generation_eia923 table is the authoritative source of information about
  how much generation is attributable to an individual generator, if it reports
  in that table.
* The generation_fuel_eia923 table is the authoritative source of information
  about how much generation and fuel consumption is attributable to an entire
  plant.
* The generators_eia860 table provides an exhaustive list of all generators
  whose generation is being reported in the generation_fuel_eia923 table.

We allocate the net generation reported in the generation_fuel_eia923 table on the basis
of plant, prime mover, and fuel type among the generators in each plant that have
matching fuel types. Generation is allocated proportional to reported generation if it's
available, and proportional to each generator's capacity if generation is not available.

In more detail: within each year of data, we split the plants into three groups:

* Plants where ALL generators report in the more granular generation_eia923
  table.
* Plants where NONE of the generators report in the generation_eia923 table.
* Plants where only SOME of the generators report in the generation_eia923
  table.

In plant-years where ALL generators report more granular generation, the total net
generation reported in the generation_fuel_eia923 table is allocated in proportion to
the generation each generator reported in the generation_eia923 table. We do this
instead of using net_generation_mwh from generation_eia923 because there are some small
discrepancies between the total amounts of generation reported in these two tables.

In plant-years where NONE of the generators report more granular generation, we create a
generator record for each associated fuel type. Those records are merged with the
generation_fuel_eia923 table on plant, prime mover code, and fuel type. Each group of
plant, prime mover, and fuel will have some amount of reported net generation associated
with it, and one or more generators. The net generation is allocated among the
generators within the group in proportion to their capacity. Then the allocated net
generation is summed up by generator.

In the hybrid case, where only SOME of of a plant's generators report the more granular
generation data, we use a combination of the two allocation methods described above.
First, the total generation reported across a plant in the generation_fuel_eia923 table
is allocated between the two categories of generators (those that report fine-grained
generation, and those that don't) in direct proportion to the fraction of the plant's
generation which is reported in the generation_eia923 table, relative to the total
generation reported in the generation_fuel_eia923 table.

Note that this methology does not distinguish between primary and secondary
energy_sources for generators. It associates portions of net generation to each
generators in the same plant do not report detailed generation, have the same
prime_mover_code, and use the same fuels, but have very different capacity factors in
reality, this methodology will allocate generation such that they end up with very
similar capacity factors. We imagine this is an uncommon scenario.

This methodology has several potential flaws and drawbacks. Because there is no
indicator of what portion of the energy_source_codes, we associate the net generation
equally among them. In effect, if a plant had multiple generators with the same
prime_mover_code but opposite primary and secondary fuels (eg. gen 1 has a primary fuel
of 'NG' and secondary fuel of 'DFO', while gen 2 has a primary fuel of 'DFO' and a
secondary fuel of 'NG'), the methodology associates the generation_fuel_eia923 records
similarly across these two generators. However, the allocated net generation will still
be porporational to each generator's net generation (if it's reported) or capacity (if
generation is not reported).

"""

import logging
import warnings

# Useful high-level external modules.
import numpy as np
import pandas as pd
import sqlalchemy as sa

import src.load_data as load_data
import src.data_cleaning as data_cleaning

#import pudl.helpers
#from pudl.metadata.fields import apply_pudl_dtypes
from typing import List

logger = logging.getLogger(__name__)

IDX_GENS = ["plant_id_eia", "generator_id", "report_date"]
"""Id columns for generators."""

IDX_PM_FUEL = ["plant_id_eia", "prime_mover_code", "energy_source_code", "report_date"]
"""Id columns for plant, prime mover & fuel type records."""

IDX_FUEL = ["report_date", "plant_id_eia", "energy_source_code"]

DATA_COLS = ["net_generation_mwh", "fuel_consumed_mmbtu","fuel_consumed_for_electricity_mmbtu"]
"""Data columns from generation_fuel_eia923 that are being allocated."""


def allocate_gen_fuel_by_gen(year):
    """
    Allocate gen fuel data columns to generators.

    The generation_fuel_eia923 table includes net generation and fuel
    consumption data at the plant/fuel type/prime mover level. The most
    granular level of plants that PUDL typically uses is at the plant/generator
    level. This method converts the generation_fuel_eia923 table to the level
    of plant/generators.

    Args:
        year: the single year (YYYY) selected for analysis

    Returns:
        pandas.DataFrame: table with columns ``IDX_GENS`` and ``DATA_COLS``.
        The ``DATA_COLS`` will be scaled to the level of the ``IDX_GENS``.

    """
    # LOAD THE DATA
    ###############

    pudl_db = 'sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite'
    pudl_engine = sa.create_engine(pudl_db)
    

    # extract all of the tables from pudl_out early in the process and select
    # only the columns we need. this is for speed and clarity.
    
    # gf contains the more complete generation and fuel data at the plant prime mover level
    gf = load_data.load_pudl_table("generation_fuel_eia923", year=year).loc[
        :, IDX_PM_FUEL + ["net_generation_mwh", "fuel_consumed_mmbtu","fuel_consumed_for_electricity_mmbtu"]
    ].pipe(apply_dtype)
    # load the nuclear gf data
    gf_nuc = load_data.load_pudl_table("generation_fuel_nuclear_eia923", year=year).loc[
        :, IDX_PM_FUEL + ["net_generation_mwh", "fuel_consumed_mmbtu","fuel_consumed_for_electricity_mmbtu"]
    ].pipe(apply_dtype)
    # concat the nuclear data with the main gf dataframe
    gf = pd.concat([gf,gf_nuc], axis=0)
    # remove non-grid connected plants
    gf = data_cleaning.remove_non_grid_connected_plants(gf)

    # gen contrains more granular generation data at the generator level for a subset of generators
    gen = (
        load_data.load_pudl_table("generation_eia923", year=year).loc[:, IDX_GENS + ["net_generation_mwh"]]
        # removes 4 records with NaN generator_id as of pudl v0.5
        .dropna(subset=IDX_GENS)
    ).pipe(apply_dtype)
    # remove non-grid connected plants
    gen = data_cleaning.remove_non_grid_connected_plants(gen)

    # gens contains a complete list of all generators
    gens = load_data.load_pudl_table("generators_eia860", year=year).loc[
        :,
        IDX_GENS
        + [
            "capacity_mw",
            "operational_status",
            "current_planned_operating_date",
            "retirement_date",
        ]
        + list(load_data.load_pudl_table("generators_eia860", year=year).filter(like="energy_source_code")),
    ]
    # remove non-grid connected plants
    gens = data_cleaning.remove_non_grid_connected_plants(gens)
    # get a list of fuel types for later
    gen_primary_fuel = gens.copy()[['plant_id_eia','generator_id','energy_source_code_1']]
    # add the prime mover code to the gens df from generators entity
    gens = gens.merge(load_data.load_pudl_table("generators_entity_eia").loc[:,["plant_id_eia", "generator_id","prime_mover_code"]],
                        how='left', 
                        on=["plant_id_eia", "generator_id"]).pipe(apply_dtype)
    # add records for each month of the year
    gens = data_cleaning.create_monthly_gens_records(gens, year)
    # remove retired generator months
    #gens = remove_nonoperating_generators(gens)
    # TODO: some of the prime mover codes in generators_entity_eia860 are incorrect, so these need to be manually fixed for now
    gens = manually_fix_prime_movers(gens)

    # the gen table is missing some generator ids. Let's fill this using the gens table, leaving a missing value for net generation
    gen = gen.merge(gens[['plant_id_eia','generator_id','report_date']], how='outer', on=['plant_id_eia','generator_id','report_date'])

    
    

    # do the allocation! 
    ####################
    # (this function coordinates the bulk of the work in this module)
    gen_pm_fuel = allocate_gen_fuel_by_gen_pm_fuel(gf, gen, gens, year)

    # aggregate the gen/pm/fuel records back to generator records
    gen_allocated = agg_by_generator(gen_pm_fuel)
    _test_gen_fuel_allocation(gen, gen_allocated)

    # DO FINAL CLEANUP
    ##################

    # add balancing authority and state data 
    gen_allocated = data_cleaning.assign_ba_code_to_plant(gen_allocated, year)
    
    # add primary fuel data for each generator
    gen_allocated = gen_allocated.merge(gen_primary_fuel, how='left', on=['plant_id_eia','generator_id'])

    # fill any missing co2 values with zeros where appropriate
    # for any clean fuels, replace missing co2 values with 0
    clean_fuels = ['SUN','MWH','WND', 'WAT','WH','PUR','NUC']
    gen_allocated.loc[gen_allocated['energy_source_code_1'].isin(clean_fuels), 'co2_mass_tons'] = gen_allocated.loc[gen_allocated['energy_source_code_1'].isin(clean_fuels), 'co2_mass_tons'].fillna(0)

    # if there are any missing co2 records where both net generation and fuel content are zero, fill with zero
    gen_allocated.loc[(gen_allocated['net_generation_mwh'] == 0) & (gen_allocated['fuel_consumed_mmbtu'] == 0) & (gen_allocated['co2_mass_tons'].isna()), 'co2_mass_tons'] = 0

    # rename the energy source code column
    gen_allocated = gen_allocated.rename(columns={'energy_source_code_1':'energy_source_code'})

    return gen_allocated

def manually_fix_prime_movers(df):
    """
    See https://github.com/catalyst-cooperative/pudl/issues/1585
    """

    # fix cogeneration issues
    plant_ids = [2465, 50150, 54268, 54410, 54262]
    for id in plant_ids:
        df.loc[(df['plant_id_eia'] == id) & (df['prime_mover_code'] == 'CT'), 'prime_mover_code'] = 'GT'
        df.loc[(df['plant_id_eia'] == id) & (df['prime_mover_code'] == 'CA'), 'prime_mover_code'] = 'ST'

    # fix missing code
    df.loc[(df['plant_id_eia'] == 50489), 'prime_mover_code'] = 'GT'
    df.loc[(df['plant_id_eia'] == 50489) & (df['generator_id'] == 'C3'), 'prime_mover_code'] = 'ST'

    df.loc[(df['plant_id_eia'] == 10884), 'prime_mover_code'] = 'GT'

    df.loc[(df['plant_id_eia'] == 58946), 'prime_mover_code'] = 'IC'

    df.loc[(df['plant_id_eia'] == 60610), 'prime_mover_code'] = 'OT'

    df.loc[(df['plant_id_eia'] == 7854), 'prime_mover_code'] = 'IC'

    df.loc[(df['plant_id_eia'] == 50628) & (df['generator_id'] == 'GEN1'), 'prime_mover_code'] = 'GT'
    df.loc[(df['plant_id_eia'] == 50628) & (df['generator_id'] == 'GEN2'), 'prime_mover_code'] = 'ST'
    df.loc[(df['plant_id_eia'] == 50628) & (df['generator_id'] == 'GEN3'), 'prime_mover_code'] = 'ST'

    df.loc[(df['plant_id_eia'] == 55088) & ~(df['generator_id'] == 'ST1'), 'prime_mover_code'] = 'GT'

    df.loc[(df['plant_id_eia'] == 52168) & (df['generator_id'] == 'GEN2'), 'prime_mover_code'] = 'ST'
    df.loc[(df['plant_id_eia'] == 52168) & (df['generator_id'] == 'GEN3'), 'prime_mover_code'] = 'GT'
    df.loc[(df['plant_id_eia'] == 52168) & (df['generator_id'] == 'GEN4'), 'prime_mover_code'] = 'GT'

    df.loc[(df['plant_id_eia'] == 55096) & (df['generator_id'] == 'GT'), 'prime_mover_code'] = 'GT'
    df.loc[(df['plant_id_eia'] == 55096) & (df['generator_id'] == 'ST'), 'prime_mover_code'] = 'ST'

    df.loc[(df['plant_id_eia'] == 7887) & (df['generator_id'] == '4'), 'prime_mover_code'] = 'GT'

    df.loc[(df['plant_id_eia'] == 6474) & (df['generator_id'] == 'GT1'), 'prime_mover_code'] = 'GT'
    df.loc[(df['plant_id_eia'] == 6474) & (df['generator_id'] == 'GT2'), 'prime_mover_code'] = 'GT'

    return df


def allocate_gen_fuel_by_gen_pm_fuel(gf, gen, gens, year, drop_interim_cols=True):
    """
    Proportionally allocate net gen from gen_fuel table to generators.

    Two main steps here:
     * associate `generation_fuel_eia923` table data w/ generators
     * allocate `generation_fuel_eia923` table data proportionally

     The association process happens via `associate_generator_tables()`.

     The allocation process (via `calc_allocation_fraction()`) entails
     generating a fraction for each record within a ``IDX_PM_FUEL`` group. We
     have two data points for generating this ratio: the net generation in the
     generation_eia923 table and the capacity from the generators_eia860 table.
     The end result is a `frac` column which is unique for each
     generator/prime_mover/fuel record and is used to allocate the associated
     net generation from the `generation_fuel_eia923` table.

     Args:
        gf (pandas.DataFrame): generator_fuel_eia923 table with columns:
            ``IDX_PM_FUEL`` and `net_generation_mwh` and `fuel_consumed_mmbtu`.
        gen (pandas.DataFrame): generation_eia923 table with columns:
            ``IDX_GENS`` and `net_generation_mwh`.
        gens (pandas.DataFrame): generators_eia860 table with cols:
            ``IDX_GENS``, `capacity_mw`, `prime_mover_code`,
            and all of the `energy_source_code` columns
        drop_interim_cols (boolean): True/False flag for dropping interim
            columns which are used to generate the `net_generation_mwh` column
            (they are mostly the `frac` column and  net generataion reported in
            the original generation_eia923 and generation_fuel_eia923 tables)
            that are useful for debugging. Default is True, which will drop
            the columns.

    Returns:
        pandas.DataFrame
    """
    gen_assoc = associate_generator_tables(gf=gf, gen=gen, gens=gens)

    # Generate a fraction to use to allocate net generation by.
    # These two methods create a column called `frac`, which will be a fraction
    # to allocate net generation from the gf table for each `IDX_PM_FUEL` group
    gen_pm_fuel = prep_alloction_fraction(gen_assoc)
    gen_pm_fuel_frac = calc_allocation_fraction(gen_pm_fuel)
    
    # calculate emissions by generator based on fuel type
    # drop any rows where all values are missing
    gen_pm_fuel_frac = gen_pm_fuel_frac.dropna(how='all', axis=0)
    gen_pm_fuel_frac = data_cleaning.calculate_co2_from_heat_content(gen_pm_fuel_frac, year)

    # do the allocating-ing!
    gen_pm_fuel_frac = (
        gen_pm_fuel_frac.assign(
            # we could x.net_generation_mwh_g_tbl.fillna here if we wanted to
            # take the net gen
            net_generation_mwh=lambda x: x.net_generation_mwh_gf_tbl * x.frac,
            # let's preserve the gf version of fuel consumption (it didn't show
            # up in the tables we pulled together in associate_generator_tables()).
            # TODO: THIS IS A HACK! We need to generate a proper fraction for
            # allocating fuel consumption based on the boiler_fuel_eia923 tbl
            fuel_consumed_mmbtu_gf_tbl=lambda x: x.fuel_consumed_mmbtu,
            fuel_consumed_mmbtu=lambda x: x.fuel_consumed_mmbtu * x.frac,
            fuel_consumed_for_electricity_mmbtu_gf_tbl=lambda x: x.fuel_consumed_for_electricity_mmbtu,
            fuel_consumed_for_electricity_mmbtu=lambda x: x.fuel_consumed_for_electricity_mmbtu * x.frac,
            # NOTE: allocating emissions by fraction of net generation
            # TODO: allocate by fraction of heat input
            co2_mass_tons=lambda x: x.co2_mass_tons * x.frac,
        )
        .pipe(apply_dtype)
        .dropna(how="all")
        .pipe(_test_gen_pm_fuel_output, gf=gf, gen=gen)
    )
    if drop_interim_cols:
        gen_pm_fuel_frac = gen_pm_fuel_frac[
            IDX_PM_FUEL
            + [
                "generator_id",
                "energy_source_code_num",
                "net_generation_mwh",
                "fuel_consumed_mmbtu",
                "fuel_consumed_for_electricity_mmbtu",
                "co2_mass_tons",
            ]
        ]
    return gen_pm_fuel_frac


def agg_by_generator(gen_pm_fuel):
    """
    Aggreate the allocated gen fuel data to the generator level.

    Args:
        gen_pm_fuel (pandas.DataFrame): result of
            `allocate_gen_fuel_by_gen_pm_fuel()`
    """
    data_cols = ["net_generation_mwh", "fuel_consumed_mmbtu","fuel_consumed_for_electricity_mmbtu","co2_mass_tons"]
    gen = (
        gen_pm_fuel.groupby(by=IDX_GENS)[data_cols]
        .sum(min_count=1)
        .reset_index()
        .pipe(apply_dtype)
    )

    return gen


def stack_generators(
    gens, cat_col="energy_source_code_num", stacked_col="energy_source_code"
):
    """
    Stack the generator table with a set of columns.

    Args:
        gens (pandas.DataFrame): generators_eia860 table with cols: ``IDX_GENS``
            and all of the `energy_source_code` columns
        cat_col (string): name of category column which will end up having the
            column names of cols_to_stack
        stacked_col (string): name of column which will end up with the stacked
            data from cols_to_stack

    Returns:
        pandas.DataFrame: a dataframe with these columns: idx_stack, cat_col,
        stacked_col

    """
    esc = list(gens.filter(like="energy_source_code"))
    gens_stack_prep = (
        pd.DataFrame(gens.set_index(IDX_GENS)[esc].stack(level=0))
        .reset_index()
        .rename(columns={"level_3": cat_col, 0: stacked_col})
        .pipe(apply_dtype)
    )

    # merge the stacked df back onto the gens table
    # we first drop the cols_to_stack so we don't duplicate data
    gens_stack = pd.merge(
        gens.drop(columns=esc), gens_stack_prep, on=IDX_GENS, how="outer"
    )
    return gens_stack


def associate_generator_tables(gf, gen, gens):
    """
    Associate the three tables needed to assign net gen to generators.

    Args:
        gf (pandas.DataFrame): generator_fuel_eia923 table with columns:
            ``IDX_PM_FUEL`` and `net_generation_mwh` and `fuel_consumed_mmbtu`.
        gen (pandas.DataFrame): generation_eia923 table with columns:
            ``IDX_GENS`` and `net_generation_mwh`.
        gens (pandas.DataFrame): generators_eia860 table with cols: ``IDX_GENS``
            and all of the `energy_source_code` columns

    TODO: Convert these groupby/merges into transforms.
    """
    stack_gens = stack_generators(
        gens, cat_col="energy_source_code_num", stacked_col="energy_source_code"
    )

    gen_assoc = (
        pd.merge(stack_gens, gen, on=IDX_GENS, how="outer")
        .pipe(remove_nonoperating_generators)
        .merge(
            gf.groupby(by=IDX_PM_FUEL, as_index=False).sum(min_count=1),
            on=IDX_PM_FUEL,
            suffixes=("_g_tbl", "_gf_tbl"),
            how="outer",)
            
    )

    # remove any months for proposed generators where no data is available
    gen_assoc = gen_assoc[~((gen_assoc.operational_status == "proposed") & gen_assoc.net_generation_mwh_g_tbl.isna() &  (gen_assoc.net_generation_mwh_gf_tbl.isna()))]

    gen_assoc = (
        pd.merge(
            gen_assoc,
            gen_assoc.groupby(by=IDX_FUEL)[["capacity_mw", "net_generation_mwh_g_tbl"]]
            .sum(min_count=1)
            .add_suffix("_fuel")
            .reset_index(),
            on=IDX_FUEL,
        )
        .pipe(apply_dtype)
        .pipe(_associate_unconnected_records)
        .pipe(_associate_energy_source_only, gf=gf)
        .pipe(data_cleaning.update_energy_source_codes)
    )
    return gen_assoc


def remove_nonoperating_generators(gen_assoc):
    """
    Remove the retired generators.

    We don't want to associate net generation to generators that are retired
    (or proposed! or any other `operational_status` besides `existing`).

    We do want to keep the generators that retire mid-year and have generator
    specific data from the generation_eia923 table. Removing the generators
    that retire mid-report year and don't report to the generation_eia923 table
    is not exactly a great assumption. For now, we are removing them. We should
    employ a strategy that allocates only a portion of the generation to them
    based on their operational months (or by doing the allocation on a monthly
    basis).

    Args:
        gen_assoc (pandas.DataFrame): table of generators with stacked fuel
            types and broadcasted net generation data from the
            generation_eia923 and generation_fuel_eia923 tables. Output of
            `associate_generator_tables()`.
    """
    existing = gen_assoc.loc[(gen_assoc.operational_status == "existing")]
    # keep the gens for each month until they retire, if they have any data to report in that month
    retiring = gen_assoc.loc[
        (gen_assoc.operational_status == "retired")
        & ((gen_assoc.report_date <= gen_assoc.retirement_date))
        #| (gen_assoc.net_generation_mwh_g_tbl.notnull())
        #| (gen_assoc.net_generation_mwh_gf_tbl.notnull()))
    ]
    new = gen_assoc.loc[(gen_assoc.operational_status == "proposed")]

    gen_assoc_removed = pd.concat([existing, retiring, new])
    return gen_assoc_removed


def _associate_unconnected_records(eia_generators_merged):
    """
    Associate unassociated gen_fuel table records on idx_pm.

    There are a subset of generation_fuel_eia923 records which do not
    merge onto the stacked generator table on ``IDX_PM_FUEL``. These records
    generally don't match with the set of prime movers and fuel types in the
    stacked generator table. In this method, we associate those straggler,
    unconnected records by merging these records with the stacked generators on
    the prime mover only.

    Args:
        eia_generators_merged (pandas.DataFrame)

    """
    # we're associating on the plant/pm level... but we only want to associated
    # these unassocaited records w/ the primary fuel type from stack_generators
    # so we're going to merge on energy_source_code_num and
    idx_pm = [
        "plant_id_eia",
        "prime_mover_code",
        "energy_source_code_num",
        "report_date",
    ]
    # we're going to only associate these unconnected fuel records w/
    # the primary fuel so we don't have to deal w/ double counting
    connected_mask = eia_generators_merged.generator_id.notnull()
    eia_generators_connected = eia_generators_merged[connected_mask]
    eia_generators_unconnected = (
        eia_generators_merged[~connected_mask]
        .rename(columns={"energy_source_code": "energy_source_unconnected"})
        .assign(energy_source_code_num="energy_source_code_1")
        .groupby(by=idx_pm)
        .sum(min_count=1)
        .reset_index()
    )
    eia_generators = (
        pd.merge(
            eia_generators_connected,
            eia_generators_unconnected[
                idx_pm + ["net_generation_mwh_gf_tbl", "fuel_consumed_mmbtu","fuel_consumed_for_electricity_mmbtu"]
            ],
            on=idx_pm,
            suffixes=("", "_unconnected"),
            how="left",
        )
        .assign(
            # we want the main and the unconnected net gen to be added together
            # but sometimes there is no main net gen and sometimes there is no
            # unconnected net gen
            net_generation_mwh_gf_tbl=lambda x: np.where(
                x.net_generation_mwh_gf_tbl.notnull()
                | x.net_generation_mwh_gf_tbl_unconnected.notnull(),
                x.net_generation_mwh_gf_tbl.fillna(0)
                + x.net_generation_mwh_gf_tbl_unconnected.fillna(0),
                np.nan,
            ),
            fuel_consumed_mmbtu=lambda x: np.where(
                x.fuel_consumed_mmbtu.notnull()
                | x.fuel_consumed_mmbtu_unconnected.notnull(),
                x.fuel_consumed_mmbtu.fillna(0)
                + x.fuel_consumed_mmbtu_unconnected.fillna(0),
                np.nan,
            ),
            fuel_consumed_for_electricity_mmbtu=lambda x: np.where(
                x.fuel_consumed_for_electricity_mmbtu.notnull()
                | x.fuel_consumed_for_electricity_mmbtu_unconnected.notnull(),
                x.fuel_consumed_for_electricity_mmbtu.fillna(0)
                + x.fuel_consumed_for_electricity_mmbtu_unconnected.fillna(0),
                np.nan,
            ),
        )  # we no longer need these _unconnected columns
        .drop(
            columns=[
                "net_generation_mwh_gf_tbl_unconnected",
                "fuel_consumed_mmbtu_unconnected",
                "fuel_consumed_for_electricity_mmbtu_unconnected",
            ]
        )
    )
    return eia_generators


def _associate_energy_source_only(gen_assoc, gf):
    """
    Associate the records w/o prime movers with fuel cost.

    The 2001 and 2002 generation fuel table does not include any prime mover
    codes. Because of this, we need to associated these records via their fuel
    types.

    Note: 2001 and 2002 eia years are not currently integrated into PUDL.
    """
    # first fine the gf records that have no PM.
    gf_grouped = gf.groupby(by=IDX_PM_FUEL, dropna=False).sum(min_count=1).reset_index()
    gf_missing_pm = (
        gf_grouped[gf_grouped[IDX_PM_FUEL].isnull().any(axis=1)]
        .drop(columns=["prime_mover_code"])
        .set_index(IDX_FUEL)
        .add_suffix("_fuel")
        .reset_index()
        .astype({"plant_id_eia": pd.Int64Dtype()})
    )

    gen_assoc = pd.merge(
        gen_assoc,
        gf_missing_pm.pipe(apply_dtype),
        how="outer",
        on=IDX_FUEL,
        indicator=True,
    )

    gen_assoc = _associate_energy_source_only_wo_matching_energy_source(gen_assoc)

    if gf_missing_pm.empty:
        logger.info("No records found with fuel-only records. This is expected.")
    else:
        logger.info(
            f"{len(gf_missing_pm)/len(gen_assoc):.02%} records w/o prime movers now"
            f" associated for: {gf_missing_pm.report_date.dt.year.unique()}"
        )
    return gen_assoc


def _associate_energy_source_only_wo_matching_energy_source(gen_assoc):
    """
    Associate the missing-pm records that don't have matching fuel types.

    There are some generation fuel table records which don't associate with
    any of the energy_source_code's reported in for the generators. For these
    records, we need to take a step back and associate these records with the
    full plant.
    """
    idx_plant = ["plant_id_eia", "report_date"]
    gen_assoc = pd.merge(
        gen_assoc,
        gen_assoc.groupby(by=idx_plant, dropna=False)[["capacity_mw"]]
        .sum(min_count=1)
        .add_suffix("_plant")
        .reset_index(),
        on=idx_plant,
        how="left",
    )

    gen_assoc_w_unassociated = (
        pd.merge(
            gen_assoc[(gen_assoc._merge != "right_only") | (gen_assoc._merge.isnull())],
            (
                gen_assoc[gen_assoc._merge == "right_only"]
                .groupby(idx_plant)[
                    ["net_generation_mwh_fuel", "fuel_consumed_mmbtu_fuel","fuel_consumed_for_electricity_mmbtu_fuel"]
                ]
                .sum(min_count=1)
            ),
            on=idx_plant,
            how="left",
            suffixes=("", "_missing_pm"),
        )
        .assign(
            net_generation_mwh_gf_tbl=lambda x: x.net_generation_mwh_gf_tbl.fillna(
                x.net_generation_mwh_fuel  # TODO: what is this?
                + x.net_generation_mwh_fuel_missing_pm.fillna(0)
            ),
            fuel_consumed_mmbtu=lambda x: x.fuel_consumed_mmbtu.fillna(
                x.fuel_consumed_mmbtu_fuel
                + x.fuel_consumed_mmbtu_fuel_missing_pm.fillna(0)
            ),
            fuel_consumed_for_electricity_mmbtu=lambda x: x.fuel_consumed_for_electricity_mmbtu.fillna(
                x.fuel_consumed_for_electricity_mmbtu_fuel
                + x.fuel_consumed_for_electricity_mmbtu_fuel_missing_pm.fillna(0)
            ),
        )
        .drop(columns=["_merge"])
    )
    return gen_assoc_w_unassociated


def prep_alloction_fraction(gen_assoc):
    """
    Make flags and aggregations to prepare for the `calc_allocation_ratios()`.

    In `calc_allocation_ratios()`, we will break the generators out into four
    types - see `calc_allocation_ratios()` docs for details. This function adds
    flags for splitting the generators. It also adds

    """
    # flag whether the generator exists in the
    # generation table (this will be used later on)
    # for calculating ratios to use to allocate net generation
    gen_assoc = gen_assoc.assign(
        in_g_tbl=lambda x: np.where(x.net_generation_mwh_g_tbl.notnull(), True, False)
    )

    gens_gb = gen_assoc.groupby(by=IDX_PM_FUEL, dropna=False)
    # get the total values for the merge group
    # we would use on groupby here with agg but it is much slower
    # so we're gb-ing twice w/ a merge
    # gens_gb.agg({'net_generation_mwh_g_tbl': lambda x: x.sum(min_count=1),
    #              'capacity_mw': lambda x: x.sum(min_count=1),
    #              'in_g_tbl': 'all'},)
    gen_pm_fuel = (
        gen_assoc.merge(  # flag if all generators exist in the generators_eia860 tbl
            gens_gb[["in_g_tbl"]].all().reset_index(),
            on=IDX_PM_FUEL,
            suffixes=("", "_all"),
        )
        .merge(  # flag if some generators exist in the generators_eia860 tbl
            gens_gb[["in_g_tbl"]].any().reset_index(),
            on=IDX_PM_FUEL,
            suffixes=("", "_any"),
        )
        # Net generation and capacity are both proxies that can be used
        # to allocate the generation which only shows up in generation_fuel
        # Sum them up across the whole plant-prime-fuel group so we can tell
        # what fraction of the total capacity each generator is.
        .merge(
            (
                gens_gb[["net_generation_mwh_g_tbl", "capacity_mw"]]
                .sum(min_count=1)
                .add_suffix("_pm_fuel")
                .reset_index()
            ),
            on=IDX_PM_FUEL,
        )
        .assign(
            # fill in the missing generation with zeros (this will help ensure
            # the calculations to run the fractions in `calc_allocation_ratios`
            # can be consistent)
            net_generation_mwh_g_tbl=lambda x: x.net_generation_mwh_g_tbl.fillna(0)
        )
    )
    # Add a column that indicates how much capacity comes from generators that
    # report in the generation table, and how much comes only from generators
    # that show up in the generation_fuel table.
    gen_pm_fuel = pd.merge(
        gen_pm_fuel,
        gen_pm_fuel.groupby(by=IDX_PM_FUEL + ["in_g_tbl"], dropna=False)[
            ["capacity_mw"]
        ]
        .sum(min_count=1)
        .add_suffix("_in_g_tbl_group")
        .reset_index(),
        on=IDX_PM_FUEL + ["in_g_tbl"],
    )
    return gen_pm_fuel


def calc_allocation_fraction(gen_pm_fuel, drop_interim_cols=True):
    """
    Make `frac` column to allocate net gen from the generation fuel table.

    There are three main types of generators:
      * "all gen": generators of plants which fully report to the
        generators_eia860 table.
      * "some gen": generators of plants which partially report to the
        generators_eia860 table.
      * "gf only": generators of plants which do not report at all to the
        generators_eia860 table.
      * "no pm": generators that have missing prime movers.

    Each different type of generator needs to be treated slightly differently,
    but all will end up with a `frac` column that can be used to allocate
    the `net_generation_mwh_gf_tbl`.

    Args:
        gen_pm_fuel (pandas.DataFrame): output of `prep_alloction_fraction()`.
        drop_interim_cols (boolean): True/False flag for dropping interim
            columns which are used to generate the `frac` column (they are
            mostly interim frac columns and totals of net generataion from
            various groupings of generators) that are useful for debugging.
            Default is False.

    """
    # break out the table into these four different generator types.
    no_pm_mask = gen_pm_fuel.net_generation_mwh_fuel_missing_pm.notnull()
    no_pm = gen_pm_fuel[no_pm_mask]
    all_gen = gen_pm_fuel.loc[gen_pm_fuel.in_g_tbl_all & ~no_pm_mask]
    some_gen = gen_pm_fuel.loc[
        gen_pm_fuel.in_g_tbl_any & ~gen_pm_fuel.in_g_tbl_all & ~no_pm_mask
    ]
    gf_only = gen_pm_fuel.loc[~gen_pm_fuel.in_g_tbl_any & ~no_pm_mask]

    logger.info(
        "Ratio calc types: \n"
        f"   All gens w/in generation table:  {len(all_gen)}#, {all_gen.capacity_mw.sum():.2} MW\n"
        f"   Some gens w/in generation table: {len(some_gen)}#, {some_gen.capacity_mw.sum():.2} MW\n"
        f"   No gens w/in generation table:   {len(gf_only)}#, {gf_only.capacity_mw.sum():.2} MW\n"
        f"   GF table records have no PM:     {len(no_pm)}#"
    )
    if len(gen_pm_fuel) != len(all_gen) + len(some_gen) + len(gf_only) + len(no_pm):
        raise AssertionError(
            "Error in splitting the gens between records showing up fully, "
            "partially, or not at all in the generation table."
        )

    # In the case where we have all of teh generation from the generation
    # table, we still allocate, because the generation reported in these two
    # tables don't always match perfectly
    all_gen = all_gen.assign(
        frac_net_gen=lambda x: x.net_generation_mwh_g_tbl
        / x.net_generation_mwh_g_tbl_pm_fuel)
    # if any values are both zero, replace the nan with 1
    all_gen.loc[(all_gen.net_generation_mwh_g_tbl == 0) & (all_gen.net_generation_mwh_g_tbl_pm_fuel == 0), 'frac_net_gen'] = 1
    all_gen = all_gen.assign(frac=lambda x: x.frac_net_gen)

    # _ = _test_frac(all_gen)

    # a brief explaination of the equations below
    # input definitions:
    #   ng == net generation from the generation table (by generator)
    #   ngf == net generation from the generation fuel table (summed by PM/Fuel)
    #   ngt == total net generation from the generation table (summed by PM/Fuel)
    #
    # y = ngt / ngf (fraction of generation reporting in the generation table)
    # z = ng * ngt (fraction of generation from generation table by generator)
    # g = y * z  (fraction of generation reporting in generation table by generator - frac_gen)

    some_gen = some_gen.assign(
        # fraction of the generation that should go to the generators that
        # report in the generation table
        frac_from_g_tbl=lambda x: x.net_generation_mwh_g_tbl_pm_fuel
        / x.net_generation_mwh_gf_tbl,
        # for records within these mix groups that do have net gen in the
        # generation table..
        frac_net_gen=lambda x: x.net_generation_mwh_g_tbl
        / x.net_generation_mwh_g_tbl_pm_fuel)  # generator based net gen from gen table
    
    # if any values are both zero, replace the nan with 1
    # TODO: make sure this is not messing up any calculations
    some_gen.loc[(some_gen.net_generation_mwh_g_tbl_pm_fuel == 0) & (some_gen.net_generation_mwh_gf_tbl == 0), 'frac_from_g_tbl'] = 1
    some_gen.loc[(some_gen.net_generation_mwh_g_tbl == 0) & (some_gen.net_generation_mwh_g_tbl_pm_fuel == 0), 'frac_net_gen'] = 1
    
    some_gen = some_gen.assign(    
        frac_gen=lambda x: x.frac_net_gen * x.frac_from_g_tbl,
        # fraction of generation that does not show up in the generation table
        frac_missing_from_g_tbl=lambda x: 1 - x.frac_from_g_tbl,
        capacity_mw_missing_from_g_tbl=lambda x: np.where(x.in_g_tbl, 0, x.capacity_mw),
        frac_cap=lambda x: x.frac_missing_from_g_tbl
        * (x.capacity_mw_missing_from_g_tbl / x.capacity_mw_in_g_tbl_group),
        # the real deal
        # this could aslo be `x.frac_gen + x.frac_cap` because the frac_gen
        # should be 0 for any generator that does not have net gen in the g_tbl
        # and frac_cap should be 0 for any generator that has net gen in the
        # g_tbl.
        frac=lambda x: np.where(x.in_g_tbl, x.frac_gen, x.frac_cap),
    )
    # _ = _test_frac(some_gen)

    # Calculate what fraction of the total capacity is associated with each of
    # the generators in the grouping.
    gf_only = gf_only.assign(
        frac_cap=lambda x: x.capacity_mw / x.capacity_mw_pm_fuel,
        frac=lambda x: x.frac_cap,
    )
    # _ = _test_frac(gf_only)

    no_pm = no_pm.assign(
        # ratio for the records with a missing prime mover that are
        # assocaited at the plant fuel level
        frac_net_gen_fuel=lambda x: x.net_generation_mwh_gf_tbl
        / x.net_generation_mwh_g_tbl_fuel,
        frac_cap_fuel=lambda x: x.capacity_mw / x.capacity_mw_fuel,
        frac=lambda x: np.where(
            x.frac_net_gen_fuel.notnull() | x.frac_net_gen_fuel != 0,
            x.frac_net_gen_fuel,
            x.frac_cap_fuel,
        ),
    )

    # squish all of these methods back together.
    gen_pm_fuel_ratio = pd.concat([all_gen, some_gen, gf_only, no_pm])
    # null out the inf's
    gen_pm_fuel_ratio.loc[abs(gen_pm_fuel_ratio.frac) == np.inf] = np.NaN
    _ = _test_frac(gen_pm_fuel_ratio)

    # drop all of the columns we needed to get to the `frac` column
    if drop_interim_cols:
        gen_pm_fuel_ratio = gen_pm_fuel_ratio[
            IDX_PM_FUEL
            + [
                "generator_id",
                "energy_source_code_num",
                "frac",
                "net_generation_mwh_gf_tbl",
                "net_generation_mwh_g_tbl",
                "capacity_mw",
                "fuel_consumed_mmbtu",
                "fuel_consumed_for_electricity_mmbtu"
            ]
        ]
    return gen_pm_fuel_ratio


def _test_frac(gen_pm_fuel):
    # test! Check if each of the IDX_PM_FUEL groups frac's add up to 1
    ratio_test_pm_fuel = (
        gen_pm_fuel.groupby(IDX_PM_FUEL)[["frac", "net_generation_mwh_g_tbl","net_generation_mwh_gf_tbl"]]
        .sum(min_count=1)
        .reset_index()
    )

    ratio_test_fuel = (
        gen_pm_fuel.groupby(IDX_FUEL)[["frac", "net_generation_mwh_fuel"]]
        .sum(min_count=1)
        .reset_index()
    )

    frac_test = pd.merge(
        ratio_test_pm_fuel, ratio_test_fuel, on=IDX_FUEL, suffixes=("", "_fuel")
    ).assign(
        frac_pm_fuel=lambda x: x.frac,
        frac=lambda x: np.where(
            x.frac_pm_fuel.notnull(),
            x.frac_pm_fuel,
            x.frac_fuel,
        ),
    )

    frac_test_bad = frac_test[~np.isclose(frac_test.frac, 1) & frac_test.frac.notnull()]
    if not frac_test_bad.empty:
        # raise AssertionError(
        warnings.warn(
            f"Ooopsies. You got {len(frac_test_bad)} records where the "
            "'frac' column isn't adding up to 1 for each 'IDX_PM_FUEL' "
            "group. Check 'calc_allocation_fraction()'\n"
            f"{frac_test_bad}"
        )

    # export this to CSV so that we can investigate
    frac_test_bad.to_csv('../data/temp/frac_test_bad.csv')
    return frac_test_bad


def _test_gen_pm_fuel_output(gen_pm_fuel, gf, gen):
    # this is just for testing/debugging
    def calc_net_gen_diff(gen_pm_fuel, idx):
        gen_pm_fuel_test = pd.merge(
            gen_pm_fuel,
            gen_pm_fuel.groupby(by=idx)[["net_generation_mwh"]]
            .sum(min_count=1)
            .add_suffix("_test")
            .reset_index(),
            on=idx,
            how="outer",
        ).assign(
            net_generation_mwh_diff=lambda x: x.net_generation_mwh_gf_tbl
            - x.net_generation_mwh_test
        )
        return gen_pm_fuel_test

    # make different totals and calc differences for two different indexs
    gen_pm_fuel_test = calc_net_gen_diff(gen_pm_fuel, idx=IDX_PM_FUEL)
    gen_fuel_test = calc_net_gen_diff(gen_pm_fuel, idx=IDX_FUEL)

    gen_pm_fuel_test = gen_pm_fuel_test.assign(
        net_generation_mwh_test=lambda x: x.net_generation_mwh_test.fillna(
            gen_fuel_test.net_generation_mwh_test
        ),
        net_generation_mwh_diff=lambda x: x.net_generation_mwh_diff.fillna(
            gen_fuel_test.net_generation_mwh_diff
        ),
    )

    bad_diff = gen_pm_fuel_test[
        (~np.isclose(gen_pm_fuel_test.net_generation_mwh_diff, 0))
        & (gen_pm_fuel_test.net_generation_mwh_diff.notnull())
    ]
    logger.info(
        f"{len(bad_diff)/len(gen_pm_fuel):.03%} of records have are partially "
        "off from their 'IDX_PM_FUEL' group"
    )
    no_cap_gen = gen_pm_fuel_test[
        (gen_pm_fuel_test.capacity_mw.isnull())
        & (gen_pm_fuel_test.net_generation_mwh_g_tbl.isnull())
    ]
    if len(no_cap_gen) > 15:
        logger.info(f"Warning: {len(no_cap_gen)} records have no capacity or net gen")
    # remove the junk/corrective plants
    fuel_net_gen = gf[gf.plant_id_eia != "99999"].net_generation_mwh.sum()
    fuel_consumed = gf[gf.plant_id_eia != "99999"].fuel_consumed_mmbtu.sum()
    fuel_consumed_for_electricity = gf[gf.plant_id_eia != "99999"].fuel_consumed_for_electricity_mmbtu.sum()
    logger.info(
        "gen v fuel table net gen diff:      "
        f"{(gen.net_generation_mwh.sum())/fuel_net_gen:.1%}"
    )
    logger.info(
        "new v fuel table net gen diff:      "
        f"{(gen_pm_fuel_test.net_generation_mwh.sum())/fuel_net_gen:.1%}"
    )
    logger.info(
        "new v fuel table fuel (mmbtu) diff: "
        f"{(gen_pm_fuel_test.fuel_consumed_mmbtu.sum())/fuel_consumed:.1%}"
    )
    logger.info(
        "new v fuel table fuel for electricty (mmbtu) diff: "
        f"{(gen_pm_fuel_test.fuel_consumed_for_electricity_mmbtu.sum())/fuel_consumed_for_electricity:.1%}"
    )

    gen_pm_fuel_test = gen_pm_fuel_test.drop(
        columns=["net_generation_mwh_test", "net_generation_mwh_diff"]
    )
    return gen_pm_fuel_test


def _test_gen_fuel_allocation(gen, gen_allocated, ratio=0.05):
    gens_test = pd.merge(
        gen_allocated, gen, on=IDX_GENS, suffixes=("_new", "_og")
    ).assign(
        net_generation_new_v_og=lambda x: x.net_generation_mwh_new
        / x.net_generation_mwh_og
    )

    os_ratios = gens_test[
        (~gens_test.net_generation_new_v_og.between((1 - ratio), (1 + ratio)))
        & (gens_test.net_generation_new_v_og.notnull())
    ]
    os_ratio = len(os_ratios) / len(gens_test)
    logger.info(
        f"{os_ratio:.2%} of generator records are more than {ratio:.0%} off from the net generation table"
    )
    if ratio == 0.05 and os_ratio > 0.15:
        warnings.warn(
            f"Many generator records that have allocated net gen more than {ratio:.0%}"
        )



def clean_merge_asof(
    left: pd.DataFrame,
    right: pd.DataFrame,
    left_on: str = "report_date",
    right_on: str = "report_date",
    by: List[str] = [],
) -> pd.DataFrame:
    """
    Merge two dataframes having different ``report_date`` frequencies.
    We often need to bring together data which is reported on a monthly basis, and
    entity attributes that are reported on an annual basis.  The
    :func:`pandas.merge_asof` is designed to do this, but requires that dataframes are
    sorted by the merge keys (``left_on``, ``right_on``, and ``by`` here). We also need
    to make sure that all merge keys have identical data types in the two dataframes
    (e.g. ``plant_id_eia`` needs to be a nullable integer in both dataframes, not a
    python int in one, and a nullable :func:`pandas.Int64Dtype` in the other).  Note
    that :func:`pandas.merge_asof` performs a left merge, so the higher frequency
    dataframe **must** be the left dataframe.
    We also force both ``left_on`` and ``right_on`` to be a Datetime using
    :func:`pandas.to_datetime` to allow merging dataframes having integer years with
    those having datetime columns.
    Because :func:`pandas.merge_asof` searches backwards for the first matching date,
    this function only works if the less granular dataframe uses the convention of
    reporting the first date in the time period for which it reports. E.g. annual
    dataframes need to have January 1st as the date. This is what happens by defualt if
    only a year or year-month are provided to :func:`pandas.to_datetime` as strings.
    Args:
        left: The higher frequency "data" dataframe. Typically monthly in our use
            cases. E.g. ``generation_eia923``. Must contain ``report_date`` and any
            columns specified in the ``by`` argument.
        right: The lower frequency "attribute" dataframe. Typically annual in our uses
            cases. E.g. ``generators_eia860``. Must contain ``report_date`` and any
            columns specified in the ``by`` argument.
        left_on: Column in ``left`` to merge on using merge_asof. Default is
            ``report_date``. Must be convertible to a Datetime using
            :func:`pandas.to_datetime`
        right_on: Column in ``right`` to merge on using :func:`pd.merge_asof`.  Default
            is ``report_date``. Must be convertible to a Datetime using
            :func:`pandas.to_datetime`
        by: Columns to merge on in addition to ``report_date``. Typically ID columns
            like ``plant_id_eia``, ``generator_id`` or ``boiler_id``.
    Returns:
        Merged contents of left and right input dataframes.  Will be sorted by
        ``left_on`` and any columns specified in ``by``. See documentation for
        :func:`pandas.merge_asof` to understand how this kind of merge works.
    Raises:
        ValueError: if ``left_on`` or ``right_on`` columns are missing from their
            respective input dataframes.
        ValueError: if any of the labels referenced in ``by`` are missing from either
            the left or right dataframes.
    """
    # Make sure we've got all the required inputs...
    if left_on not in left.columns:
        raise ValueError(f"Left dataframe has no column {left_on}.")
    if right_on not in right.columns:
        raise ValueError(f"Right dataframe has no {right_on}.")
    missing_left_cols = [col for col in by if col not in left.columns]
    if missing_left_cols:
        raise ValueError(f"Left dataframe is missing {missing_left_cols}.")
    missing_right_cols = [col for col in by if col not in right.columns]
    if missing_right_cols:
        raise ValueError(f"Left dataframe is missing {missing_right_cols}.")

    def cleanup(df, on, by):
        df = apply_dtype(df)
        df.loc[:, on] = pd.to_datetime(df[on])
        df = df.sort_values([on] + by)
        return df

    return pd.merge_asof(
        cleanup(df=left, on=left_on, by=by),
        cleanup(df=right, on=right_on, by=by),
        left_on=left_on,
        right_on=right_on,
        by=by,
        tolerance=pd.Timedelta("365 days"),  # Should never match across years.
    )

def apply_dtype(df):
    "converts columns to standardized dtype based on name"

    # convert all columns with "date" in name to datetime dtype
    df = df.astype({col: "datetime64[ns]" for col in df.columns if "date" in col})

    # convert all columns with "plant_id" in name to integer dtype
    df = df.astype({col: "Int64" for col in df.columns if "plant_id" in col})

    return df
