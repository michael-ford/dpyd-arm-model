#!/usr/bin/env python3
"""Convert all CSV results to a single Excel workbook with multiple sheets."""

from pathlib import Path

import pandas as pd

# Files to exclude from export
EXCLUDE_FILES = set()

# Mapping of CSV filenames to clean sheet names
SHEET_NAME_MAP = {
    # WT Unified model
    "wt_unified_absolute_risks.csv": "Absolute Risks",
    "wt_unified_odds_ratios.csv": "Odds Ratios",
    "wt_unified_model_summary.csv": "Model Summary",
    "wt_unified_publication_summary.csv": "Publication Summary",
    # Pairwise analysis
    "pairwise_probability_results.csv": "Pairwise Probability",
}


# Columns to drop from specific sheet types
DROP_COLUMNS = {
    "publication_summary": ["P_Best"],
}


def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Replace underscores with spaces in column names and apply title case."""
    df.columns = [col.replace("_", " ").title() for col in df.columns]
    return df


def get_sheet_name(csv_path: Path) -> str:
    """Get clean sheet name from mapping or generate one."""
    filename = csv_path.name

    if filename in SHEET_NAME_MAP:
        return SHEET_NAME_MAP[filename]

    # Fallback: clean up the filename
    name = filename.replace(".csv", "").replace("_", " ").title()
    return name[:31]  # Excel 31-char limit


def export_to_xlsx(output_dir: str = "output", xlsx_name: str = "dpyd_nma_results.xlsx"):
    """Export all CSV files in output directory to a single Excel workbook."""
    output_path = Path(output_dir)
    xlsx_path = output_path / xlsx_name

    # Find all CSV files
    csv_files = sorted(output_path.rglob("*.csv"))

    if not csv_files:
        print(f"No CSV files found in {output_dir}")
        return

    print(f"Found {len(csv_files)} CSV files")

    # Filter out excluded files
    csv_files = [f for f in csv_files if f.name not in EXCLUDE_FILES]

    # Create Excel writer
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for csv_file in csv_files:
            sheet_name = get_sheet_name(csv_file)
            print(f"  Adding: {csv_file.name} -> '{sheet_name}'")

            # Read CSV and drop unwanted columns
            df = pd.read_csv(csv_file)
            for key, cols in DROP_COLUMNS.items():
                if key in csv_file.name:
                    df = df.drop(columns=[c for c in cols if c in df.columns], errors="ignore")

            # Clean column names and write to sheet
            df = clean_column_names(df)
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"\nExported to: {xlsx_path}")
    return xlsx_path


if __name__ == "__main__":
    export_to_xlsx()
