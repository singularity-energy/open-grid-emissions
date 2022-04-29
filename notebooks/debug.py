#%%
%reload_ext autoreload
%autoreload 2

import sys
sys.path.append('../../hourly-egrid/')

# Useful high-level external modules.
import numpy as np
import pandas as pd
import sqlalchemy as sa

import pudl

year = 2020

pudl_db = 'sqlite:///../data/pudl/pudl_data/sqlite/pudl.sqlite'
pudl_engine = sa.create_engine(pudl_db)

pudl_out = pudl.output.pudltabl.PudlTabl(
    pudl_engine,
    freq='MS',
    start_date=f'{year}-01-01',
    end_date=f'{year}-12-31'
)
#%%
test = pudl.analysis.allocate_net_gen.allocate_gen_fuel_by_generator_energy_source(pudl_out)


# %%

IDX_PM_ESC = ["report_date", "plant_id_eia", "energy_source_code", "prime_mover_code"]
pudl_out.gf_eia923().loc[:, IDX_PM_ESC + ["net_generation_mwh", "fuel_consumed_mmbtu","fuel_consumed_for_electricity_mmbtu"]]
# %%
