---
stoplight-id: biomass_adjustments
---

## Background on adjusting biomass emissions
The combustion of biomass releases greenhouse gases and other air pollutants into the atmosphere. However, the EPA's eGRID database includes a legacy calculation of biomass-adjusted emissions values. As they explain in the eGRID technical support document:
> Prior editions of eGRID applied a biomass adjustment to the annual emission values based on an assumption of zero emissions from biomass combustion. This assumes that the amount of carbon sequestered during biomass growth equals the amount released during combustion, without consideration of other factors. For reasons of consistency, the same approach is applied by eGRID on an ongoing basis, even though this approach is not consistent with the most recent science on the topic.

This approach of assuming zero emissions from biomass combustion is problematic for several reasons:
1. There is much debate in the academic literature about the assumption of zero net emissions from biomass. This is far from a comprehensive literature review on the topic, but for example see: [Johnson 2009](https://www.sciencedirect.com/science/article/pii/S0195925508001637), [Cherubini et al 2011](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1757-1707.2011.01102.x), [Haberl et al 2012](https://www.sciencedirect.com/science/article/pii/S0301421512001681), [Downie et al 2014](https://www.sciencedirect.com/science/article/pii/S0961953413004820)
2. This approach selectively applies a partial life-cycle accounting approach to biomass fuels (as it consideres the upstream emissions impacts of the fuel), which is inconsistent with the treatment of other fuels in this dataset

Based on our current understanding of this topic, it may not be appropriate to use biomass-adjusted emissions data for carbon accounting or other general uses, unless they are being used in a specific policy or regulatory context that treats biomass emissions as carbon neutral. Thus, biomass-adjusted emissions are only included in the OGE dataset for consistency with eGRID and for use in these niche cases. All emissions data that has been adjusted for biomass emissions will include `_adjusted` in the name of the column.

## Calculating biomass-adjusted emissions
Adjusted CO2 emissions are set to zero for all biomass fuel consumption, including agricultural byproducts (AB), black liquor (BLQ), landfill gas (LFG), biogenic municipal solid waste (MSB), other biomass gas (OBG), other biomass liquids (OBL), other biomass solids (OBS), sludge waste (SLW), wood and wood waste solids (WDS), and wood waste liquids (WDL).

No adjustment is applied to CH4, N2O, NOx, or SO2 emissions, so the adjusted values of these pollutants are equal to their unadjusted values. Because CO2e is calculated from CO2, CH4, and N2O, adjusted CO2e only differs from unadjusted CO2e because of the CO2 adjustment.

Prior versions of OGE also adjusted CH4, N2O, NOx, and SO2 emissions from landfill gas (LFG), following the eGRID methodology at the time. That methodology assumed landfills would otherwise flare the gas, so it set adjusted CH4, N2O, and SO2 from LFG to zero, and subtracted a flaring baseline from NOx. Starting in 2023, eGRID removed these landfill gas adjustments for all pollutants except CO2, and OGE now does the same. This approach was likely updated due to problems related to applying a consequential LCA concept to an otherwise direct emissions attributional inventory of power sector emissions.

## Future Work, Known Issues, and Open Questions
- Determine whether we should continue publishing biomass-adjusted emissions ([details](https://github.com/singularity-energy/open-grid-emissions/issues/130))
- Consider removing the `_adjusted` and `_for_electricity_adjusted` columns for CH4, N2O, NOx, and SO2, which are now identical to the unadjusted and `_for_electricity` columns, respectively. 
- Look into updated emission factors for "other biomass" fuels ([details](https://github.com/singularity-energy/open-grid-emissions/issues/69))
- Consider the biogenic and nonbiogenic components of MSW fuel when adjusting emissions from CEMS ([details](https://github.com/singularity-energy/open-grid-emissions/issues/51))