# Grid Gross Loss: Current Method

*Working document. It describes the current calculation and is updated as the method changes. Last updated 2026-11-18.*

Grid gross loss (GGL) is the share of electricity that enters the grid and is lost in transmission and distribution before it reaches an end user. We calculate it directly from the utility-level data reported on EIA-861 instead of from EIA's State Electricity Profiles (which eGRID uses). That's because the State Electricity Profiles apply a single national loss rate to every state, so they carry no regional information about losses (see [Differences from eGRID](#differences-from-egrid)).

The current method is defined at the national level. A regional allocation is not defined yet.

## Formula

```
GGL = total_losses / grid_energy_input

total_losses = reported_losses + nonreporter_losses

grid_energy_input = (total_disposition - sales_for_resale - behind_the_meter_sales)
                    + nonreporter_sales + nonreporter_losses
                  = (retail_sales - behind_the_meter_sales) + furnished_without_charge + consumed_by_respondent
                    + reported_losses + nonreporter_sales + nonreporter_losses

nonreporter_sales = adjustment_sales - customer_sited_sales

nonreporter_losses = nonreporter_sales × reported_losses / (retail_sales - behind_the_meter_sales)
```

| Variable | Description | Source |
|---|---|---|
| `reported_losses` | Total energy losses reported by EIA-861 respondents | `core_eia861__yearly_operational_data_misc.total_energy_losses_mwh` |
| `total_disposition` | Total disposition | `core_eia861__yearly_operational_data_misc.total_disposition_mwh` |
| `sales_for_resale` | Sales for resale | `core_eia861__yearly_operational_data_misc.sales_for_resale_mwh` |
| `retail_sales` | Retail sales by respondents | `core_eia861__yearly_operational_data_misc.retail_sales_mwh` |
| `behind_the_meter_sales` | Retail sales by third-party behind-the-meter providers | `core_eia861__yearly_operational_data_misc.retail_sales_mwh` where `entity_type = 'Behind the Meter'` |
| `furnished_without_charge` | Furnished without charge | `core_eia861__yearly_operational_data_misc.furnished_without_charge_mwh` |
| `consumed_by_respondent` | Consumed by respondent without charge | `core_eia861__yearly_operational_data_misc.consumed_by_respondent_without_charge_mwh` |
| `adjustment_sales` | EIA's state-level adjustment to retail sales | `core_eia861__yearly_sales.sales_mwh` where `utility_id_eia = 99999` |
| `customer_sited_sales` | Sales to ultimate customers by customer-sited nonutility plants | EIA-923 nonutility source and disposition schedule, `Retail Sales` column (also published as `Sales_Ult_Cust_CS_{year}.xlsx` in the EIA-861 files) |
| `nonreporter_sales` | Grid-delivered retail sales by utilities that do not report operational data | Calculated |
| `nonreporter_losses` | Estimated losses on `nonreporter_sales` | Calculated |

All values are annual totals in MWh, summed over all EIA-861 respondents in the 50 states and DC.

## Elements of the formula

### Reported losses (`reported_losses`)

Losses are taken directly from the total energy losses that each utility reports in Schedule 2, Part B (energy sources and disposition) of EIA-861. The form instructions define this as electricity "lost from transmission, distribution, and/or unaccounted for", calculated as the difference between total sources and the sum of the other dispositions. In other words, reported losses are a residual of each utility's energy accounts, so they also absorb any metering, billing-period, or accounting mismatches. This is probably why reported losses vary more from year to year than physical losses should (see [Validation](#validation)).

Every respondent that reports losses is included, not only utilities with retail customers. That includes wholesale power marketers, transmission companies, and federal power marketing administrations, which together report about 15 TWh of losses (about 8% of the total) on transmission-level flows. Not every respondent reports losses, though (see [Coverage of the operational data](#coverage-of-the-operational-data)).

### Total disposition minus sales for resale (`total_disposition - sales_for_resale`)

Each utility reports how all of the electricity it obtained was disposed of:

```
total_disposition = retail_sales + sales_for_resale + furnished_without_charge
                    + consumed_by_respondent + reported_losses
```

Summed nationally, total disposition is about 7.1 PWh (2023), well above the roughly 4.1 PWh that actually flows through the grid. That's because sales for resale double count: power sold from a generating utility to a distribution utility appears once as the seller's sale for resale, and again in the buyer's disposition when it is sold to retail customers. Power that passes through several intermediaries can be counted several times. This double counting is why using total disposition directly as the denominator (as in the original draft of PR #161) gives a national GGL of about 2.4% in 2023 instead of 4.36%.

Subtracting sales for resale removes these intermediate transactions and leaves only the final uses of electricity plus losses. Nationally, every sale for resale is (in principle) a purchase by another respondent, so the subtraction doesn't lose any delivered energy, with two exceptions: power sold to utilities that don't file operational data (handled by `nonreporter_sales` below), and power exported to Canada and Mexico (see [Open questions](#open-questions)).

The remaining components of the denominator are:

- **Retail sales (`retail_sales`)**: electricity sold to end-use customers, including both bundled sales and energy-only sales by retail power marketers. Delivery-only service (where a wires utility delivers energy that a customer bought from a separate retailer) is not part of disposition, because the wires utility never owns that energy. This avoids counting retail-choice customers twice.
- **Furnished without charge (`furnished_without_charge`)**: electricity delivered to customers without charge, such as municipal street lighting. It travels over the grid like a retail sale, so it is part of the energy that can be lost.
- **Consumed by respondent without charge (`consumed_by_respondent`)**: electricity a utility uses "in its electric and other departments" (for example, offices and service buildings). Station use at power plants is not part of this, because utilities report net generation (gross generation minus station use) as their source, so station use never enters their energy accounts. Note that this is *not* the same as "direct use" in the State Electricity Profiles (see below).
- **Reported losses (`reported_losses`)**: included so that the denominator represents all energy that entered the grid, not only what was delivered. This makes GGL a share of gross energy supplied (losses / (delivered + losses)) rather than a ratio to delivered energy (losses / delivered).

### Behind-the-meter sales (`behind_the_meter_sales`)

About a dozen EIA-861 respondents (for example, Sunrun, Tesla and Sunnova) own rooftop solar on their customers' premises and sell its output to those customers. They report these as retail sales (about 12 TWh in 2023), but the electricity never enters the grid, so it is subtracted from the denominator.

### Sales by utilities that don't report operational data (`nonreporter_sales`)

The EIA-861 sales file includes a record for each state named "Adjustment {year}" with utility ID 99999. EIA describes this as an adjustment "for non-respondents and for customer-sited respondents", and these records are flagged as imputed rather than observed. Comparing these records to the other EIA-861 files shows that they are made up of three parts (2024 shown):

| Part of the 99999 adjustment | TWh (2024) | Grid-delivered? |
|---|---|---|
| Bundled sales by short-form (EIA-861S) utilities | 90.1 | Yes |
| Customer-sited sales (nonutility plants selling directly to an end user on the same site) | 17.4 | No |
| Energy-only sales adjustment in retail-choice states | 20.0 | Yes |
| Other | 0.1 | |
| **Total** | **127.7** | |

The bundled part of the adjustment equals short-form sales plus customer-sited sales in every state, to within 0.1 TWh. The energy-only part appears to be an adjustment for retail power marketers. Besides a net addition of sales, it also moves large volumes between customer sectors (for example, in Ohio it subtracts about 20 TWh of commercial sales and adds about 23 TWh of industrial sales), which suggests EIA reclassifies some marketers' customers.

**Short-form utilities.** Utilities may file the short form (EIA-861S) instead of the full form if they sell less than 200,000 MWh per year, provide only bundled service, are not needed for the quality of EIA's statistical estimates, and don't report on the monthly form (EIA-861M). The short form only collects retail sales, revenue and customer counts by state and balancing authority, so these utilities report no energy sources, dispositions or losses. The short form is offered in seven out of every eight years, and in the eighth year (most recently 2019) all utilities file the full form. Its coverage has changed over time:

| Years | Short-form utilities | Short-form sales (TWh) |
|---|---|---|
| 2012–2018 | about 1,100 | about 27 |
| 2019 | none (all utilities filed the full form) | 0 |
| 2020–2025 | about 1,690 | 87–93 |

These utilities buy most of their power from full respondents, so that power is subtracted as a sale for resale in `total_disposition - sales_for_resale`. But because they don't report operational data, their retail sales never come back into the denominator. Adding `nonreporter_sales` back is required to close the national energy balance (see [Validation](#validation)).

**Customer-sited sales** are excluded from `nonreporter_sales`. They are equal to the retail sales reported by nonutility plants on EIA-923 (17.4 TWh in 2024, matching exactly by state), and EIA describes them as customer-sited, meaning the plant sells directly to an end user at the same site. Like direct use and behind-the-meter sales, this electricity doesn't travel over the grid.

### Estimated losses on those sales (`nonreporter_losses`)

Utilities that don't report operational data also don't report losses, but their customers' electricity is still lost in transmission and distribution. We estimate these losses by assuming that these utilities have the same ratio of losses to grid-delivered retail sales as full respondents. The estimated losses are added to both the numerator and the denominator.

This is a simple assumption. Short-form utilities are small municipal and cooperative distribution utilities, whose reported loss rates (about 3.6% to 5.0% of retail sales in 2023) are similar to the national average, so the choice of rate has a small effect on the national result (using 3.6% or 5.0% instead of the national ratio changes 2023 GGL by about 0.03 percentage points).

## What is not included, and why

- **Direct use**: Electricity that commercial and industrial facilities generate and consume on site (about 125–137 TWh per year) never enters the grid, so it cannot be lost. Because no utility buys or sells it, it never appears in EIA-861 disposition, so `grid_energy_input` excludes it automatically and the GGL formula needs no direct use adjustment. Direct use only has to be subtracted when the grid's energy is calculated from the supply side (from generation), as in eGRID and in our [validation](#validation).
- **Customer-sited and behind-the-meter sales**: electricity sold to an end user on the same site as the generator never enters the grid (see above).
- **Station use**: Electricity consumed within power plants is already excluded from net generation, which is the source utilities report.
- **Delivery-only sales**: excluded to avoid double counting retail-choice customers, as described under retail sales.
- **U.S. territories**: PUDL's EIA-861 tables cover the 50 states and DC only.

## Coverage of the operational data

The EIA-861 frame file lists every entity that reports to EIA-861 and which schedules each one files. In 2024:

| Entities (2024) | Count | Reports retail sales? | Reports operational data (including losses)? |
|---|---|---|---|
| Full-form utilities with retail customers | 1,521 | Yes | Yes |
| Full-form entities without retail customers (wholesale marketers, transmission companies, federal agencies, etc.) | 192 | No | Yes |
| Short-form utilities | 1,687 | Yes (short form) | No |
| Other (for example, demand-side management administrators) | 13 | No | No |

All operational data records are flagged as observed (none are imputed). Full-form respondents account for about 97% of national retail sales (3,755 of 3,874 TWh in 2023). The remaining 3% is the 99999 adjustment described above.

However, not every full-form respondent reports losses. In 2023, respondents reporting zero or missing losses accounted for about 707 TWh of retail sales (19%):

| Entity type (2023) | Respondents | Reporting zero or missing losses | Their retail sales (TWh) |
|---|---|---|---|
| Retail power marketers | 258 | 173 | 571 |
| Investor-owned utilities | 166 | 29 | 55 |
| Community choice aggregators | 26 | 18 | 34 |
| Cooperatives | 604 | 37 | 14 |
| Municipal utilities | 463 | 22 | 13 |
| Behind-the-meter providers | 13 | 13 | 12 |
| Other (political subdivisions, state and federal agencies, wholesale marketers, transmission companies) | 192 | 96 | 7 |

Most of these are not missing losses, because the losses on that electricity are reported by someone else. For example, in most retail-choice states, the wires utility reports losses on all the energy it delivers, including energy sold by retail marketers and community choice aggregators. But there are some real gaps:

- **Texas wires utilities report no operational data.** Oncor, CenterPoint, AEP Texas and Texas-New Mexico Power own most of the distribution system in ERCOT, but their operational data records are empty. In ERCOT, losses are reported (if at all) by retail electric providers and wholesale entities, and 66 Texas retail marketers with 267 TWh of retail sales report zero losses.
- **Some utilities balance purchases against sales with no losses.** For example, Central Maine Power reports purchases exactly equal to its retail sales plus sales for resale.
- **Some losses may be reported in other categories.** For example, Consolidated Edison reports no losses but reports 1.4 TWh of energy consumed by respondent, much more than other utilities of its size. Losses reported in other categories end up in the denominator instead of the numerator.

## Validation

### National energy balance

Physically, the electricity that enters the grid must equal the electricity delivered from it plus what is lost along the way:

```
net_generation + net_international_imports - direct_use - customer_sited_sales
    = delivered_energy + total_losses
```

Each side of this balance can be calculated from different EIA sources, so it is an independent check on the losses we use. Losses are a small difference between two large numbers (about 5% of supply), so the balance cannot check losses precisely, but it can show whether reported losses are roughly right and whether any terms are missing.

| Term | Source |
|---|---|
| `net_generation` | EIA-923 net generation for all plants in the 50 states and DC, including the state-level adjustment rows (plant ID 99999) |
| `net_international_imports` | EIA-930 hourly interchange between U.S. balancing authorities and Canadian and Mexican balancing authorities, summed over the year |
| `direct_use` | Self-supplied direct use from the EIA-923 nonutility source and disposition schedule (see below) |
| `customer_sited_sales` | EIA-923 nonutility source and disposition schedule, `Retail Sales` column |
| `delivered_energy` | EIA-861 retail sales (including utility 99999 and excluding delivery-only service), minus `customer_sited_sales` and `behind_the_meter_sales`, plus `furnished_without_charge` and `consumed_by_respondent` |
| `total_losses` | Reported losses plus estimated losses on nonreporter sales |

#### Direct use

The raw `Direct_Use` column in the EIA-923 nonutility source and disposition schedule (about 330 TWh in 2023) is not the right value for `direct_use`. It is all electricity the facility used, including electricity it bought from the grid (reported as incoming electricity, about 220 TWh in 2023). Purchased electricity was delivered over the grid and appears in a utility's retail sales, so it must stay in the supply. Only the self-supplied part should be subtracted. Candidate definitions, calculated per plant and summed nationally, are:

| `direct_use` definition (TWh) | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|---|
| `max(Direct_Use - Incoming, 0)` | 133 | 133 | 128 | 127 | 128 | 125 | 122 |
| `min(Direct_Use, Gross - Station_Use)` | 141 | 140 | 136 | 136 | 137 | 135 | 133 |
| eGRID / State Electricity Profiles | 144 | 143 | 139 | 139 | 140 | 137 | 135 |

We use the first definition, which is the physically consistent one for this balance (whatever a facility used and did not buy from the grid, it supplied itself). EIA's value is closer to the second.

#### Results

| TWh | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|---|
| `net_generation` | 4,181 | 4,131 | 4,010 | 4,110 | 4,231 | 4,183 | 4,309 |
| `+ net_international_imports` | 47 | 45 | 58 | 46 | 49 | 21 | 10 |
| `- direct_use` | 133 | 133 | 128 | 127 | 128 | 125 | 122 |
| `- customer_sited_sales` | 14 | 15 | 14 | 15 | 15 | 17 | 17 |
| **Supply into the grid** | **4,081** | **4,029** | **3,925** | **4,014** | **4,137** | **4,063** | **4,179** |
| `delivered_energy` | 3,860 | 3,803 | 3,716 | 3,799 | 3,916 | 3,864 | 3,963 |
| `total_losses` | 202 | 208 | 213 | 184 | 219 | 176 | 191 |
| **Residual** (supply − delivered − losses) | **19** | **18** | **−4** | **30** | **2** | **23** | **25** |
| Residual as % of supply | 0.5% | 0.4% | −0.1% | 0.8% | 0.0% | 0.6% | 0.6% |
| Residual with EIA's direct use definition | 11 | 10 | −12 | 21 | −7 | 13 | 14 |

The balance closes within 1% of supply in every year. However, the residual is positive in most years (about 15 TWh on average), meaning that more energy enters the grid than is accounted for by deliveries plus reported losses.

#### What the balance says about losses

The balance gives a second, top-down estimate of losses (`supply - delivered_energy`), which doesn't depend on what utilities report:

| | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|---|
| Implied by the balance | 5.4% | 5.6% | 5.3% | 5.3% | 5.3% | 4.9% | 5.2% |
| Implied, with EIA's direct use definition | 5.2% | 5.4% | 5.1% | 5.1% | 5.1% | 4.7% | 4.9% |
| Reported (`total_losses` / supply) | 4.9% | 5.2% | 5.4% | 4.6% | 5.3% | 4.3% | 4.6% |

Three things stand out:

- **Reported losses may be slightly too low.** In most years, implied losses are 0.2–0.7 percentage points higher than reported losses. This is consistent with the gaps in loss reporting described in [Coverage of the operational data](#coverage-of-the-operational-data), especially the Texas wires utilities. However, the gap is about the same size as the uncertainty from the direct use definition, so it isn't conclusive.
- **The balance is too uncertain to replace reported losses.** The residual (up to about 30 TWh) is up to 15% of losses, and the choice of direct use definition alone moves implied losses by 8–11 TWh. The uncertainty would be larger at the regional level, where interchange between balancing authorities adds its own errors.
- **Reported losses vary more from year to year than implied losses.** For example, the reported rate drops from 5.4% in 2020 to 4.6% in 2021 and rises back to 5.3% in 2022, while the implied rate stays near 5.3%. Because reported losses are a residual of each utility's energy accounts, some of this variation is probably accounting noise rather than physical change.

### National GGL

| | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|---|---|
| Current method | 4.97% | 5.18% | 5.42% | 4.63% | 5.29% | 4.36% | 4.59% |
| eGRID (U.S.) | 4.87% | 5.11% | 5.27% | 4.45% | 5.12% | 4.20% | 4.40% |

eGRID values are calculated from its reported losses, total disposition and direct use without rounding.

## Differences from eGRID

eGRID calculates GGL from Table 10 of EIA's State Electricity Profiles. We checked how Table 10 is built and found that each state's "estimated losses" is that state's retail sales multiplied by a single national ratio (all EIA-861 losses divided by all EIA-861 retail sales). This holds exactly for every state we checked (Alaska, California, Maine and Texas) in every year from 2001 to 2024. As a result, eGRID's regional GGL values do not reflect regional differences in losses; they differ across interconnects only because of differences in the denominator.

The current method gives a national GGL about 0.1 to 0.2 percentage points higher than eGRID's. Most of this difference (about 0.11 percentage points in 2023) comes from eGRID including the 99999 adjustment in its denominator without any matching losses. The remainder comes from other items in eGRID's denominator (customer-sited sales, unaccounted electricity and international exports).

## Open questions

- **Missing losses in Texas**: Texas wires utilities report no operational data, and most Texas retail marketers report zero losses. An energy balance for ERCOT (which has very limited interchange with the rest of the country) could show how large this gap is.
- **International exports**: power exported to Canada and Mexico is reported as sales for resale, so it is removed from the denominator, while any losses incurred moving it to the border remain in the numerator. EIA-930 shows gross exports of 4–17 TWh per year (14 TWh in 2023), so adding them back to the denominator would lower national GGL by at most about 0.02 percentage points.
- **Year-to-year noise in reported losses**: reported national loss rates swing by up to 0.8 percentage points between years, more than the top-down estimate does. A multi-year average of the loss ratio might give a more stable estimate.
- **Regional allocation**: not yet defined. The main challenges are allocating losses reported by entities without retail sales (wholesale marketers, transmission companies, federal power marketing administrations), and the different ways losses are reported in restructured markets (for example, Texas retail providers report losses, while California wires utilities report losses on energy they deliver for community choice aggregators).
