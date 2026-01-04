#!/usr/bin/env python3
"""Convert all CSV results to a single Excel workbook with multiple sheets."""

from pathlib import Path

import pandas as pd

# Files to exclude from export
EXCLUDE_FILES = {"model_comparison_13hetho.csv"}

# Mapping of CSV filenames to clean sheet names
SHEET_NAME_MAP = {
    # Model comparisons
    "model_comparison_dic.csv": "Model Comparison DIC",
    # WT Binary model
    "wt_binary_absolute_risks.csv": "Binary Absolute Risks",
    "wt_binary_odds_ratios.csv": "Binary Odds Ratios",
    "wt_binary_model_summary.csv": "Binary Model Summary",
    "wt_binary_publication_summary.csv": "Binary Publication Summary",
    # WT Unified model
    "wt_unified_absolute_risks.csv": "Unified Absolute Risks",
    "wt_unified_odds_ratios.csv": "Unified Odds Ratios",
    "wt_unified_model_summary.csv": "Unified Model Summary",
    "wt_unified_publication_summary.csv": "Unified Publication Summary",
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

            # Read CSV, clean column names, and write to sheet
            df = pd.read_csv(csv_file)
            df = clean_column_names(df)
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"\nExported to: {xlsx_path}")
    return xlsx_path


if __name__ == "__main__":
    export_to_xlsx()
