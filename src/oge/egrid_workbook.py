"""Formats and writes the eGRID-format workbooks.

Each workbook has a contents sheet followed by one sheet for each sheet of the eGRID
data file, in eGRID's layout: column headers in the first row and column codes in the
second, colored by category, with eGRID's column widths and number formats. The headers
of columns that OGE does not calculate yet are grayed out.
"""

import re
from datetime import date
from importlib.metadata import version

import pandas as pd
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell

from oge.egrid_schema import (
    EGRID_CATEGORIES,
    EGRID_SHEET_DESCRIPTIONS,
    EGRID_SHEETS,
    GJ_UNITS,
    METRIC_UNITS,
    REGION_SHEETS,
    get_egrid_columns,
    get_unit,
)

# header fill color of each category, from the eGRID data file. Columns that are not
# in a category use the color of the None key
CATEGORY_COLORS = {
    None: "#F2F2F2",
    1: "#F2DCDB",
    2: "#E6B8B7",
    3: "#D58785",
    4: "#EBF1DE",
    5: "#C4D79B",
    6: "#76933C",
    7: "#DAEEF3",
    8: "#92CDDC",
    9: "#31869B",
    10: "#FDE9D9",
    11: "#FABF8F",
    12: "#E26B0A",
    13: "#FFFFCC",
    14: "#FFFF99",
    15: "#E4DFEC",
    16: "#B1A0C7",
    17: "#60497A",
}
# categories whose header fill is dark enough to need white text
DARK_CATEGORIES = [6, 9, 17]
# Layout of each sheet, from the eGRID data file: the number of columns at the left
# that stay visible when scrolling (None for no frozen panes) and the height of the
# header row
FROZEN_COLUMNS = {
    "UNT": 6,
    "GEN": 6,
    "PLNT": 5,
    "ST": 3,
    "BA": 3,
    "SRL": 3,
    "NRL": 3,
    "US": None,
    "GGL": None,
}
HEADER_HEIGHTS = {"UNT": 60.75, "GEN": 60.75, "GGL": 60.75}
DEFAULT_HEADER_HEIGHT = 67.5
# Column widths, from the eGRID data file. Columns not listed here have the default
# width. Widths of the regional sheets' data columns are keyed by the code without the
# sheet's prefix
DEFAULT_COLUMN_WIDTH = 13
COLUMN_WIDTHS = {
    "UNT": {"PNAME": 35.45, "PRGCODE": 17.73},
    "GEN": {"PNAME": 35.45, "CFACT": 13.82, "GENERSRC": 30.18, "GENYRRET": 15.82},
    "PLNT": {
        "PNAME": 35.45,
        "OPRNAME": 32.82,
        "UTLSRVNM": 35.73,
        "SECTOR": 15.54,
        "BANAME": 32.73,
        "SRNAME": 18,
        "CAPFAC": 13.82,
        "NAMEPCAP": 14,
        "PWRTOHT": 14,
        "UNCO2E": 15.73,
        "CHPCHTIOZ": 15.18,
        "BION2O": 14.73,
        "PLGENATO": 17.27,
        "PLGENACO": 17.27,
        "PLTOPR": 17,
        "PLCOPR": 17,
    },
    "ST": {"PSTATABB": 14.73, "FIPSST": 19.18},
    "BA": {"BANAME": 73.82, "BACODE": 19.18},
    "SRL": {"SUBRGN": 14.73, "SRNAME": 19.18},
    "NRL": {"NERC": 14.73, "NERCNAME": 30.18},
    "US": {"YEAR": 14.82},
    "GGL": {"YEAR": 10.73, "REGION": 22, "ESTLOSS": 11.82, "GGRSLOSS": 23.73},
}
REGION_DATA_WIDTHS = {
    "NAMEPCAP": 14.73,
    "GENATO": 17.27,
    "GENATH": 14.73,
    "GENACO": 17.27,
    "CLPR": 14.73,
    "TOPR": 17.27,
    "THPR": 14.73,
    "COPR": 17.27,
    "NBGNOL": 14.73,
    "NBGNBM": 15.73,
    "NBGNWI": 14.73,
    "NBGNOP": 16.54,
    "NBCLPR": 14.73,
    "NBBMPR": 15.73,
    "NBWIPR": 14.73,
    "NBOPPR": 16.27,
}
# Number formats, from the eGRID data file. Quantities are whole numbers with negative
# values in parentheses, emission rates have three decimals, resource mixes are
# percentages, and the exceptions below are keyed by the code on the unit, generator
# and plant sheets
QUANTITY_FORMAT = "#,##0;\\ \\(#,##0\\)"
RATE_FORMAT = "#,##0.000"
PERCENT_FORMAT = "0.0%"
NUMBER_FORMATS = {
    "NAMEPCAP": "#,##0.0",
    "HRSOP": "#,##0.0",
    "NUMUNT": "#,##0",
    "NUMGEN": "#,##0",
    "PWRTOHT": "#,##0",
    "CAPFAC": "#,##0.0000",
    "NBFACTOR": "#,##0.0000",
    "ELCALLOC": "#,##0.0000",
    "CFACT": "#,##0.000;\\ \\(#,##0.000\\)",
}
# the unit sheet reports emissions with one decimal, and the unit sheet's heat input,
# the plant sheet's unadjusted values and the grid gross loss sheet are whole numbers
# without parentheses
UNIT_MASS_FORMAT = "#,##0.0"
# plant Hg masses and the CH4 and N2O emission rates by fuel type are small, so they
# have more decimals
HG_MASS_FORMAT = "#,##0.000"
SMALL_RATE_FORMAT = "#,##0.0000"
# Hg columns of the plant and regional sheets, without the sheet's prefix. OGE calculates
# Hg, but leaves it blank where its Hg data are incomplete (e.g. every region in 2024),
# so their headers are not grayed out even when they have no data
HG_CODES = [
    "HGAN",  # annual Hg
    "HGRTA",  # output rate
    "HGRA",  # input rate
    "HGCRT",  # combustion output rate
    "CHGRT",  # coal output rate
    "FSHGRT",  # fossil output rate
    "CHGR",  # coal input rate
    "FSHGR",  # fossil input rate
    "UNHG",  # plant unadjusted Hg
    "UNHGSRC",  # plant unadjusted Hg source
]
# links on the contents sheet
OGE_REPOSITORY_URL = "https://github.com/singularity-energy/open-grid-emissions"
OGE_DOCUMENTATION_URL = "https://docs.singularity.energy/docs/open-grid-emissions"
OGE_ISSUES_URL = "https://github.com/singularity-energy/open-grid-emissions/issues"


def write_egrid_workbook(
    tables: dict[str, pd.DataFrame], year: int, path: str, metric: bool
):
    """Writes the tables of each sheet to an Excel workbook in the layout of the eGRID
    data file, after a contents sheet.

    The first row of each sheet holds the column headers and the second row the column
    codes, colored by category. The headers of columns that OGE does not calculate,
    i.e. columns with no data other than Hg, are grayed out.

    Args:
        tables (dict[str, pd.DataFrame]): table of each sheet, keyed by its name in
            `EGRID_SHEETS`, with the columns of `get_egrid_columns()`.
        year (int): data year.
        path (str): path of the workbook.
        metric (bool): whether the tables are in metric units.
    """
    with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
        formats = WorkbookFormats(writer.book)
        write_contents_sheet(writer.book, formats, year, metric)
        for sheet in EGRID_SHEETS:
            sheet_name = f"{sheet}{year % 100:02d}"
            table = tables[sheet]
            table.to_excel(
                writer, sheet_name=sheet_name, startrow=2, header=False, index=False
            )
            worksheet = writer.sheets[sheet_name]
            number_formats = get_number_formats(sheet, year, metric)
            for i, (code, header, category) in enumerate(
                get_egrid_columns(sheet, year, metric)
            ):
                calculated = table[code].notna().any() or is_hg_column(sheet, code)
                worksheet.write(0, i, header, formats.header(category, calculated))
                worksheet.write(
                    1, i, code, formats.header(category, calculated, code_row=True)
                )
                worksheet.set_column(
                    i,
                    i,
                    get_column_width(sheet, code),
                    formats.get({"num_format": number_formats[code]}),
                )
            # filter on the column codes
            worksheet.autofilter(1, 0, len(table) + 1, len(table.columns) - 1)
            worksheet.set_row(0, HEADER_HEIGHTS.get(sheet, DEFAULT_HEADER_HEIGHT))
            if FROZEN_COLUMNS[sheet] is not None:
                worksheet.freeze_panes(2, FROZEN_COLUMNS[sheet])


class WorkbookFormats:
    """Creates the formats of an eGRID workbook, reusing each distinct format.

    Every format uses the Arial font at eGRID's size unless the properties set another.
    """

    FONT = {"font_name": "Arial", "font_size": 8.5}

    def __init__(self, workbook: xlsxwriter.Workbook):
        self.workbook = workbook
        self.formats = {}

    def get(self, properties: dict) -> xlsxwriter.format.Format:
        """Returns the format with these properties, creating it if needed."""
        properties = {**self.FONT, **properties}
        key = tuple(sorted(properties.items()))
        if key not in self.formats:
            self.formats[key] = self.workbook.add_format(properties)
        return self.formats[key]

    def header(
        self, category: int | None, calculated: bool = True, code_row: bool = False
    ) -> xlsxwriter.format.Format:
        """Returns the format of a column header or code, colored by the column's
        category. The headers of columns that OGE does not calculate are grayed out.
        """
        properties = {
            "bold": True,
            "text_wrap": True,
            "border": 1,
            "border_color": "#000000",
        }
        if not code_row:
            properties |= {"align": "center", "valign": "vcenter"}
        if not calculated:
            properties |= {"bg_color": "#FFFFFF", "font_color": "#A6A6A6"}
        else:
            properties |= {
                "bg_color": CATEGORY_COLORS[category],
                "font_color": "#FFFFFF" if category in DARK_CATEGORIES else "#000000",
            }
        return self.get(properties)


def is_hg_column(sheet: str, code: str) -> bool:
    """Whether a column of a plant or regional sheet is one of `HG_CODES`, including
    its metric GJ twin.

    Args:
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        code (str): column code.

    Returns:
        bool: whether the column holds Hg data.
    """
    if sheet == "PLNT":
        prefix = "PL"
    elif sheet in REGION_SHEETS:
        prefix = REGION_SHEETS[sheet][0]
    else:
        return False
    base = code.removeprefix(prefix)
    return base in HG_CODES or base.removesuffix("2") in HG_CODES


def get_number_formats(sheet: str, year: int, metric: bool = False) -> dict[str, str]:
    """Returns the Excel number format of each column of a sheet of the eGRID data
    file, following the formats of the eGRID data file.

    The GJ columns of the metric file have the format of the column they follow.

    Args:
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        year (int): data year.
        metric (bool, optional): whether to return the formats of the metric units
            file. Defaults to False.

    Returns:
        dict[str, str]: number format of each column code. Columns that are not
            numbers have the "General" format.
    """
    formats = {}
    for code, header, category in get_egrid_columns(sheet, year):
        units = re.findall(r"\(([^()]+)\)", header)
        unit = units[-1] if units else None
        if code in NUMBER_FORMATS and sheet not in REGION_SHEETS:
            number_format = NUMBER_FORMATS[code]
        elif unit == "resource mix" or code == "GGRSLOSS":
            number_format = PERCENT_FORMAT
        elif unit not in [*METRIC_UNITS, *GJ_UNITS, "MW", "Btu/kWh"] or unit == "feet":
            number_format = "General"
        elif "/" in unit and unit != "Btu/kWh":
            small = category in [13, 14] and ("CH4" in code or "N2O" in code)
            number_format = SMALL_RATE_FORMAT if small else RATE_FORMAT
        elif sheet == "UNT":
            number_format = UNIT_MASS_FORMAT if unit in ["tons", "lbs"] else "#,##0"
        elif sheet == "PLNT" and "Hg" in header:
            number_format = HG_MASS_FORMAT
        elif sheet == "GGL" or (sheet == "PLNT" and category == 2):
            number_format = "#,##0"
        else:
            number_format = QUANTITY_FORMAT
        formats[code] = number_format
        if metric and get_unit(header) in GJ_UNITS:
            formats[f"{code}2"] = number_format
    return formats


def get_column_width(sheet: str, code: str) -> float:
    """Returns the width of a column of a sheet of the eGRID data file.

    Args:
        sheet (str): sheet name without the year, one of `EGRID_SHEETS`.
        code (str): column code.

    Returns:
        float: column width, in Excel's character units.
    """
    if code in COLUMN_WIDTHS[sheet]:
        return COLUMN_WIDTHS[sheet][code]
    if sheet in REGION_SHEETS:
        return REGION_DATA_WIDTHS.get(
            code.removeprefix(REGION_SHEETS[sheet][0]).removesuffix("2"),
            DEFAULT_COLUMN_WIDTH,
        )
    return DEFAULT_COLUMN_WIDTH


def write_contents_sheet(
    workbook: xlsxwriter.Workbook,
    formats: WorkbookFormats,
    year: int,
    metric: bool,
):
    """Adds a contents sheet to an eGRID workbook, laid out like the contents sheet of
    the eGRID data file.

    The sheet has a table of contents, links to OGE's documentation, a legend of the
    header colors with a link to each category on each sheet, notes and conversion
    factors. Each section is a gray box with a thick outline. The table of contents
    uses the eGRID data file's text, except that it says the file was produced with
    OGE.

    Args:
        workbook (xlsxwriter.Workbook): workbook to add the sheet to.
        formats (WorkbookFormats): formats of the workbook.
        year (int): data year.
        metric (bool): whether the workbook is in metric units.
    """
    worksheet = workbook.add_worksheet("Contents")
    yy = f"{year % 100:02d}"
    link = {"underline": 1, "font_color": "#0563C1"}

    title = {"font_size": 14, "bold": True, "align": "center", "text_wrap": True}
    heading = {"font_size": 12, "bold": True}
    large = {"font_size": 12}

    # table of contents
    sheets = [
        [
            (1, 1, i, {**large, "align": "center"}),
            (
                2,
                3,
                f"internal:'{sheet}{yy}'!A1",
                {**large, **link, "string": f"{sheet}{yy}"},
            ),
            (4, 10, f"{EGRID_SHEET_DESCRIPTIONS[sheet]} year {year} data", large),
        ]
        for i, sheet in enumerate(EGRID_SHEETS, start=1)
    ]
    write_contents_box(
        worksheet,
        formats,
        1,
        [
            [(1, 10, "Table of Contents", title)],
            [(1, 10, "for", title)],
            [(1, 10, f"eGRID{year}_Data{'_Metric' if metric else ''}.xlsx", title)],
            [
                (
                    1,
                    10,
                    f"eGRID{year} {describe_sheets()} Data Files",
                    title,
                )
            ],
            [
                (
                    1,
                    10,
                    date.today().strftime("%B %d, %Y"),
                    {**heading, "align": "center"},
                )
            ],
            [
                (1, 4, "Produced with", {**heading, "align": "right", "valign": "top"}),
                (
                    5,
                    10,
                    OGE_REPOSITORY_URL,
                    {
                        **large,
                        **link,
                        "valign": "top",
                        "string": f"Open Grid Emissions {version('oge')}.",
                    },
                ),
            ],
            [
                (1, 1, "Sheet", {**large, "underline": 1, "align": "center"}),
                (2, 3, "Name", {**large, "underline": 1}),
                (4, 10, "Description", {**large, "underline": 1}),
            ],
            *sheets,
            [(2, 10, "", {})],
        ],
        last_col=10,
    )
    worksheet.set_row(4, 43.9)
    for row in [5, 6]:
        worksheet.set_row(row, 25.9)

    # documentation and feedback
    feedback_row = 1 + len(sheets) + 9
    write_contents_box(
        worksheet,
        formats,
        feedback_row,
        [
            [(1, 10, "Documentation and Feedback", heading)],
            [(1, 10, OGE_DOCUMENTATION_URL, {**link, "string": "OGE documentation"})],
            [
                (
                    1,
                    10,
                    OGE_ISSUES_URL,
                    {**link, "string": "Report an issue or ask a question"},
                )
            ],
        ],
        last_col=10,
    )

    # color coding legend, with a link to the first column of each category on each
    # sheet
    legend_row = feedback_row + 4
    first_columns = {}
    for col, sheet in enumerate(EGRID_SHEETS[:-1], start=9):
        for i, (_, _, category) in enumerate(get_egrid_columns(sheet, year, metric)):
            if category is not None and (category, col) not in first_columns:
                first_columns[(category, col)] = (sheet, i)
    legend = [
        [
            (1, 7, "Category", {**heading, "bottom": 2}),
            (8, 8, "Color", {**heading, "align": "center", "bottom": 2}),
            (
                9,
                16,
                "Link to sheet and category",
                {**heading, "align": "center", "bottom": 2},
            ),
        ]
    ]
    thin = {"left": 1, "right": 1, "top": 1, "bottom": 1}
    for category, name in EGRID_CATEGORIES.items():
        row = [
            (1, 7, f"{category}) {name}", thin),
            (8, 8, "", {**thin, "bg_color": CATEGORY_COLORS[category]}),
        ]
        for col in range(9, 17):
            if (category, col) in first_columns:
                sheet, i = first_columns[(category, col)]
                cell = xl_rowcol_to_cell(0, i)
                row.append(
                    (
                        col,
                        col,
                        f"internal:'{sheet}{yy}'!{cell}",
                        {**thin, **link, "align": "center", "string": sheet},
                    )
                )
            else:
                row.append((col, col, "", thin))
        legend.append(row)
    legend.append(
        [
            (1, 7, "Not calculated by OGE yet", thin),
            (8, 8, "", {**thin, "bg_color": "#FFFFFF"}),
            *[(col, col, "", thin) for col in range(9, 17)],
        ]
    )
    write_contents_box(
        worksheet,
        formats,
        legend_row,
        [[(1, 16, "Color Coding Legend", {**heading, "align": "center"})], *legend],
        last_col=16,
    )

    # notes
    notes_row = legend_row + len(legend) + 2
    notes = [
        "Values in parentheses are negative numbers.",
        "Rows are the plants and balancing authorities in OGE's results. Plants are "
        "identified by their EIA plant ID (ORISPL).",
        "Columns that OGE does not calculate yet are blank, and their headers are "
        "grayed out.",
    ]
    write_contents_box(
        worksheet,
        formats,
        notes_row,
        [[(1, 10, "Notes", heading)], *[[(1, 10, note, {})] for note in notes]],
        last_col=10,
    )

    # conversion factors
    factors_row = notes_row + len(notes) + 2
    factors = [
        ("1 megawatt-hour (MWh)", "1,000 kilowatt-hour (kWh)"),
        ("1 short ton", "2,000 pounds (lb)"),
        ("1.10231 short ton", "1 metric ton"),
        ("2.2046 lb", "1 kilogram (kg)"),
        ("0.9478 MMBtu", "1 Gigajoule (GJ)"),
        ("0.27778 MWh", "1 GJ"),
    ]
    write_contents_box(
        worksheet,
        formats,
        factors_row,
        [
            [(1, 5, "Conversion Factors", heading)],
            *[[(1, 2, left, {}), (3, 4, right, {})] for left, right in factors],
        ],
        last_col=5,
    )

    for col, width in {
        0: 1.73,
        1: 15,
        2: 11,
        4: 12.54,
        5: 18.73,
        6: 11,
        9: 5.73,
    }.items():
        worksheet.set_column(col, col, width)
    worksheet.set_column(10, 16, 6)


def describe_sheets() -> str:
    """Returns a list of the sheets of the eGRID data file for the contents sheet's
    title, e.g. "Unit, Generator, ..., and Grid Gross Loss (%)", in title case.
    """
    descriptions = [
        " ".join(
            word.capitalize() if word.islower() else word
            for word in description.split()
        )
        for description in EGRID_SHEET_DESCRIPTIONS.values()
    ]
    return f"{', '.join(descriptions[:-1])}, and {descriptions[-1]}"


def write_contents_box(
    worksheet: xlsxwriter.worksheet.Worksheet,
    formats: WorkbookFormats,
    first_row: int,
    rows: list[list[tuple]],
    last_col: int,
):
    """Writes a section of the contents sheet as a gray box with a thick outline,
    starting at column B.

    Args:
        worksheet (xlsxwriter.worksheet.Worksheet): contents sheet.
        formats (WorkbookFormats): formats of the workbook.
        first_row (int): first row of the box, from 0.
        rows (list[list[tuple]]): cells of each row of the box, as (first column, last
            column, value, format properties), merged across their columns. Values
            that start with "internal:" or "https:" are written as links, whose text
            is the "string" property.
        last_col (int): last column of the box, from 0.
    """
    last_row = first_row + len(rows) - 1

    def cell_format(row: int, first_col: int, cell_last_col: int, properties: dict):
        properties = {
            "font_size": 11,
            "bg_color": CATEGORY_COLORS[None],
            "valign": "vcenter",
            "border_color": "#000000",
            **properties,
        }
        for side, on_edge in {
            "left": first_col == 1,
            "right": cell_last_col == last_col,
            "top": row == first_row,
            "bottom": row == last_row,
        }.items():
            if on_edge:
                properties[side] = 2
        return formats.get(properties)

    for row, cells in enumerate(rows, start=first_row):
        for col in range(1, last_col + 1):
            worksheet.write_blank(row, col, None, cell_format(row, col, col, {}))
        for first_col, cell_last_col, value, properties in cells:
            properties = dict(properties)
            string = properties.pop("string", None)
            merged_format = cell_format(row, first_col, cell_last_col, properties)
            if cell_last_col > first_col:
                worksheet.merge_range(
                    row, first_col, row, cell_last_col, "", merged_format
                )
            if isinstance(value, str) and value.startswith(("internal:", "https:")):
                worksheet.write_url(row, first_col, value, merged_format, string=string)
            else:
                worksheet.write(row, first_col, value, merged_format)
