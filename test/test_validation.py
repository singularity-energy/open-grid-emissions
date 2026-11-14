import pandas as pd

from oge.validation import hourly_profile_source_metric

DATA_METRICS = [
    "net_generation_mwh",
    "co2_mass_lb",
    "co2_mass_lb_for_electricity",
    "nox_mass_lb",
    "nox_mass_lb_for_electricity",
    "so2_mass_lb",
    "so2_mass_lb_for_electricity",
]


def test_hourly_profile_source_metric_keeps_method_names():
    """Profile methods should keep their full names, in the order of the methods,
    including those that aren't in the order (which go last)."""
    cems = pd.DataFrame({"plant_id_eia": [1], **{c: [1.0] for c in DATA_METRICS}})
    no_partial_cems = cems.iloc[0:0]
    shaped_eia_data = pd.DataFrame(
        {
            "ba_code": "CISO",
            "profile_method": [
                "flat_negative_generation",
                "net_load_proxy",
                "a_new_method",
                "residual_profile",
            ],
            **{c: 1.0 for c in DATA_METRICS},
        }
    )
    plant_attributes = pd.DataFrame({"plant_id_eia": [1], "ba_code": ["CISO"]})

    result = hourly_profile_source_metric(
        cems, no_partial_cems, no_partial_cems, shaped_eia_data, plant_attributes
    )

    expected = [
        "cems_reported",
        "eia_shaped_residual_profile",
        "eia_shaped_net_load_proxy",
        "eia_shaped_flat_negative_generation",
        "eia_shaped_a_new_method",
    ]
    for ba_code in ["CISO", "US Total"]:
        methods = result.loc[result["ba_code"] == ba_code, "profile_method"]
        assert methods.to_list() == expected
    assert result["net_generation_mwh"].to_list() == [0.2] * 10


def test_summarize_cems_measurement_quality_maps_mats_codes():
    """MATS startup and shutdown values count as measured, and manually calculated
    values count as imputed."""
    from oge.validation import summarize_cems_measurement_quality

    cems = pd.DataFrame(
        {
            "co2_mass_lb": [1.0, 1.0, 1.0],
            "co2_mass_measurement_code": ["Measured", "Measured", "Substitute"],
            "nox_mass_lb": [1.0, 1.0, 1.0],
            "nox_mass_measurement_code": ["Measured", "LME", "Measured"],
            "so2_mass_lb": [1.0, 1.0, 1.0],
            "so2_mass_measurement_code": ["Measured", "Measured", "Measured"],
            "hg_mass_lb": [0.6, 0.2, 0.2],
            "hg_mass_measurement_code": [
                "Measured",
                "Startup or Shutdown",
                "Manually Calculated",
            ],
        }
    )

    result = summarize_cems_measurement_quality(cems).set_index("index")

    assert result.loc["Measured", "hg_mass_lb"] == 0.8
    assert result.loc["Imputed", "hg_mass_lb"] == 0.2
