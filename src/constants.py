CLEAN_FUELS = [
    "SUN",  # solar
    "MWH",  # electricity/storage
    "WND",  # wind
    "WAT",  # water/hydro
    "WH",  # waste heat
    "PUR",  # purchased steam
    "NUC",  # nuclear
]

UNIT_CONVERSIONS = {"lb": ("kg", 0.453592), "mmbtu": ("GJ", 1.055056)}

TIME_RESOLUTIONS = {"hourly": "H", "monthly": "M", "annual": "A"}
