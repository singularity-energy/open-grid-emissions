import os
import argparse
import sys

sys.path.append("src")


import consumed as consumed


def main(args):
    year = int(snakemake.wildcards.year)
    path_prefix = f"{snakemake.wildcards.run_type}/{year}/"

    print("18. Calculating and exporting consumption-based results")
    hourly_consumed_calc = consumed.HourlyConsumed(
        snakemake.input[0],
        path_prefix,
        year=year,
        small=args.small,
        skip_outputs=snakemake.config["data_pipeline"]["skip_outputs"],
    )
    hourly_consumed_calc.run()
    hourly_consumed_calc.output_results()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--small",
        type=bool,
        default=False,
    )

    args = parser.parse_args()

    main(args)
