---
stoplight-id: fossil_only_residual_mix
---

## Background on residual mixes
The GHG Protocol Scope 2 Guidance describes two methods for calculating emissions from electricity consumption: location-based and market-based accounting. Location-based accounting uses the emission rate of the electricity physically delivered by the grid, while market-based accounting uses the emission rate of the specific generation that an entity has a contractual claim to, for example through energy attribute certificates (EACs) such as renewable energy certificates (RECs).

In market-based accounting, a residual mix emission rate is needed for electricity that has no specific claim attached to it, for example:
- Electricity purchased from an unknown source, such as unspecified wholesale market purchases or unspecified imports from another region.
- "Null power", where the environmental attributes of the generation have been stripped away and sold separately as unbundled certificates.

The residual mix represents the "leftover" mix of generation that no other entity has claimed. Using a residual mix instead of the grid average for this electricity prevents clean energy from being claimed twice: once by the owner of the certificate, and again by every consumer using a grid average that includes that clean energy.

The Center for Resource Solutions (CRS) [Guidance for Calculating Residual Mix (March 2024)](https://resource-solutions.org/document/030624/) describes two types of residual mix:
- **Type A**: the mix of all unclaimed generation attributes in a region and time period, calculated as total emissions minus the emissions from all specified transactions (EAC retirements, wholesale market attributions, bilateral power purchase agreements, and retained distributed generation).
- **Type B**: a Type A residual mix plus Renewable Portfolio Standard (RPS) generation and other potentially shared specified attributes. It is closer to a utility's "standard delivery" mix.

## Fossil-only residual mix as a proxy
Calculating a Type A residual mix requires data on all contractual claims and transactions in a market. This data is generally not publicly available, especially at the hourly resolution needed for hourly market-based accounting, so it is impractical for most users to calculate an actual residual mix.

When the data needed to calculate a residual mix is not available, the CRS guidance recommends a data hierarchy of alternative sources (Tables 5 and 6 of the guidance), which includes "regional fossil-only resource mixes and emissions factors". The CRS guidance describes this approach as the "most conservative" alternative, since all renewables and other clean energy are "assumed to be transacted and removed", leaving only fossil generation as unclaimed. It is less accurate than an actual residual mix, but avoids double counting clean energy claims.

Because the fossil-only residual mix can be calculated from public generation and emissions data, OGE can provide it at the same regional and temporal resolution as the rest of the power sector data. Singularity's [blog post on proxy residual emission rates](https://medium.com/singularity-energy/making-proxy-residual-emission-rates-data-more-accessible-6e245d7fea0b) describes this approach in more detail.

## Calculating the fossil-only residual mix
The power sector data files for each BA and for the US include a `residual` fuel category, in addition to the rows for each individual fuel category and the `total` row. The `residual` row sums the generation, fuel consumption, and emissions from the following fuel categories:
- `coal`
- `natural_gas`
- `petroleum`

The emission rates in the `residual` row are calculated in the same way as for every other fuel category, by dividing the total emissions from these fuel categories by their total net generation. For example, for CO2:

fossil-only residual mix CO2 rate (lb/MWh) = (coal CO2 + natural gas CO2 + petroleum CO2) &frasl; (coal net generation + natural gas net generation + petroleum net generation)

The `residual` row is provided at hourly, monthly, and annual resolution, for all pollutants, and for each type of emissions adjustment (`_for_electricity`, `_adjusted`, and `_for_electricity_adjusted`). Since the residual mix does not include any biomass fuel categories, the biomass-adjusted (`_adjusted`) values of the residual mix are generally equal to the unadjusted values, except for small amounts of biomass co-fired at subplants whose primary fuel is coal, natural gas, or petroleum.

The `residual` row is included for every BA and every time period, so that all BA files have the same rows. If a BA has no generation from these fuel categories in a time period (for example, a BA with only hydro or wind generation), the values in the `residual` row are zero, including the emission rates.

### Fuel categories excluded from the residual mix
Following the CRS guidance, the fossil-only residual mix only includes coal, oil, and gas. In particular, the `waste` fuel category is not included, even though it includes some fossil-derived fuels such as non-biogenic municipal solid waste and tire-derived fuels, because:
- The `waste` category also includes biogenic municipal solid waste, and the biogenic and non-biogenic components cannot be separated at the fuel category level.
- Waste-to-energy generation can generate RECs in some states. Since the fossil-only residual mix is meant to be a conservative estimate that excludes any generation that could be claimed through RECs, it should not be included.

The `other` fuel category is also not included. It consists of purchased steam and waste heat, which OGE treats as clean fuels, and generation reported with an unknown or "other" fuel type, which cannot be identified as fossil.

In the power sector data, generation and emissions are assigned to a fuel category based on the primary fuel of each subplant (see [Assigning Plant Primary Fuel](../Data%20Aggregation/Plant%20Primary%20Fuel.md#assigning-fuel-categories)). The residual mix may therefore include small amounts of secondary fuels consumed at subplants whose primary fuel is coal, natural gas, or petroleum, and may exclude small amounts of these fuels consumed at subplants with a different primary fuel.

## Using the fossil-only residual mix
- **Regional boundaries**: The CRS guidance recommends using the residual mix of the region where the electricity was generated or purchased from, and narrowing geographic boundaries where this improves the quality of the data. The BA-level residual mix is likely to be the most appropriate in most cases. The national (`US.csv`) residual mix is provided for completeness, but covers a much larger area than the markets in which electricity is generally transacted.
- **Temporal resolution**: The CRS guidance focuses on annual residual mixes. Hourly fossil-only rates can vary substantially from hour to hour, so they are only appropriate for hourly market-based accounting. Monthly or annual rates should be used for annual accounting.
- **Generated rates**: The residual mix emission rates are generated emission rates, based on the fossil generation located within each region. They do not account for fossil electricity imported into or exported from the region.
- **Small BAs**: As described in the [emissions overview](../Emissions%20Calculations/Overview.md), emission rates in small BAs can be very high when a fossil plant is starting up and consuming fuel while generating little electricity. This can have a larger effect on the residual mix than on the total, since the residual mix only includes fossil generation.

## Future Work, Known Issues, and Open Questions
- Consider providing a consumption-based fossil-only residual mix, which would account for fossil generation imported into each region.
