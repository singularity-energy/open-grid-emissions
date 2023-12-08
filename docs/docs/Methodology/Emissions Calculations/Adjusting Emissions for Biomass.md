---
stoplight-id: biomass_adjustments
---

## Background on adjusting biomass emissions
The combustion of biomass releases greenhouse gases and other air pollutants into the atmosphere. However, the EPA's eGRID database includes a legacy calculation of biomass-adjusted emissions values. As they explain in the eGRID technical support document:
> Prior editions of eGRID applied a biomass adjustment to the annual emission values based on an assumption of zero emissions from biomass combustion. This assumes that the amount of carbon sequestered during biomass growth equals the amount released during combustion, without consideration of other factors. For reasons of consistency, the same approach is applied in eGRID2020.

However, this approach of assuming zero emissions from biomass combustion is problematic for several reasons:
1. There is much debate in the academic literature about the assumption of zero net emissions from biomass. This is far from a comprehensive literature review on the topic, but for example see: [Johnson 2009](https://www.sciencedirect.com/science/article/pii/S0195925508001637), [Cherubini et al 2011](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1757-1707.2011.01102.x), [Haberl et al 2012](https://www.sciencedirect.com/science/article/pii/S0301421512001681), [Downie et al 2014](https://www.sciencedirect.com/science/article/pii/S0961953413004820)
2. This approach selectively applies a partial life-cycle accounting approach to biomass fuels (as it consideres the upstream emissions impacts of the fuel), which is inconsistent with the treatment of other fuels in this dataset

Based on our current understanding of this topic, it may not be appropriate to use biomass-adjusted emissions data for carbon accounting or other general uses, unless they are being used in a specific policy or regulatory context that treats biomass emissions as carbon neutral. Thus, biomass-adjusted emissions are only included in the OGEI dataset for consistency with eGRID and for use in these niche cases. All emissions data that has been adjusted for biomass emissions will include `_adjusted` in the name of the column.

## Calculating biomass-adjusted emissions
Adjusted CO2 emissions are set to zero for all biomass fuel consumption, including agricultural byproducts (AB), black liquor (BLQ), landfill gas (LFG), biogenic municipal solid waste (MSB), other biomass gas (OBG), other biomass liquids (OBL), other biomass solids (OBS), sludge waste (SLW), wood and wood waste solids (WDS), and wood waste liquids (WDL).

All other emissions (CH4, N2O, NOx, and SO2) are only adjusted for landfill gas (LFG), based on the assumption that landfills would otherwise flare (burn) the gas anyway even if they did not use it for electricity generation (see [this issue](https://github.com/singularity-energy/open-grid-emissions/issues/73)). Adjusted CH4, N2O, and SO2 for LFG consumption is set to zero, based on the assumption (per the eGRID methodology) that "the gas would have been combusted in a flare and would have produced some emissions of NOx, SO2, CH4, and N2O anyway."

Adjusting NOx emissions from landfill gas is done by subtracting baseline flared NOx emissions (calculated using an emission factor for flaring of landfill gas of 0.078 lb NOx/mmbtu, derived from table 2.4-4 of the EPA's AP-42 document) from the combustion NOx emissions (see [this issue](https://github.com/singularity-energy/open-grid-emissions/issues/73) for further discussion of this process).

This methodology for adjusting landfill gas emissions is potentially problematic for multiple reasons. As opposed to other biomass fuels, where the arguement for these adjustments is based in attributional life-cycle accounting and involves considering upstream carbon sequestration by the biomass, the argument for adjusting landfill gas emissions is based in a consequential emissions framework that assumes that no additional emissions are occuring compared to a counterfactual baseline scenario. Selectively applying a consequential emissions approach to landfill gas emissions is not necessarily appropriate for inventorying attributional emissions from the power sector.

## Future Work, Known Issues, and Open Questions
- Determine whether we should continue publishing biomass-adjusted emissions ([details](https://github.com/singularity-energy/open-grid-emissions/issues/130))
- Update methodology for adjusting landfill gas emissions ([details](https://github.com/singularity-energy/open-grid-emissions/issues/73))
- Look into updated emission factors for "other biomass" fuels ([details](https://github.com/singularity-energy/open-grid-emissions/issues/69))
- Consider the biogenic and nonbiogenic components of MSW fuel when adjusting emissions from CEMS ([details](https://github.com/singularity-energy/open-grid-emissions/issues/51))