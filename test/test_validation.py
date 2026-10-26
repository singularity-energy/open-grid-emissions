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
