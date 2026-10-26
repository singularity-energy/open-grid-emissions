---
stoplight-id: partial_cems_plant
---

### Shaping EIA subplants with partial CEMS plant data

In certain cases, there is hourly CEMS data reported for some, but not all subplants at a single plant. In this case, we use the combined hourly profile of the subplants that do report to CEMS to shape the monthly data for the subplants that do not report any data to CEMS. This approach assumes that the hourly operational of all subplants at a single plant will be similar, although further research is needed to establish the accuracy of this assumption.

We exclude any non-combustion plants (hydro, wind, solar, nuclear, geothermal) from this method. This prevents the following cases:

1) Cases where a backup fossil generator reports to CEMS at a plant whose other subplants are renewable
2) Cases where a fossil-fuel plant has a single rooftop solar array installed resulting in one solar subplant

Battery storage at these plants is also excluded, since it isn't operated like the plant's combustion units (for example, the batteries at Moss Landing in CISO and DeCordova in ERCO). It is shaped with storage charging and discharging profiles instead (see [Energy storage profiles](Shaping%20Using%20Fleet-Specific%20Profiles.md#energy-storage-profiles)). Compressed air storage that burns natural gas when discharging (such as McIntosh) reports to CEMS as a combustion unit, so it uses its own CEMS data.

To calculate the profiles used to shape the data, the pipeline calculates the total net generation and fuel consumption for each plant-hour, and converts each hourly value to a percentage of the monthly total value. This percentage is then multiplied by the monthly total value for the subplants that we are shaping to get the hourly profile.

Net generation is shaped using the net generation profile, and all fuel and emissions-related totals are shaped using the fuel consumption profile. If one of these profiles is missing, the other is used in its place (for example, if there is no net generation profile, the fuel consumption profile is used to shape the net generation data). Additionally, if the monthly net generation total is negative, it is assigned a flat hourly profile.