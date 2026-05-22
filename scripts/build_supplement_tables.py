"""
Build publication-ready supplement tables S8, S9, S10a, S10b
from raw analysis CSVs.

Outputs to /tmp/supplement-tables/:
  Table_S8_genotyping.xlsx / .html
  Table_S9_qba.xlsx        / .html
  Table_S10a_ld_breakdown.xlsx / .html
  Table_S10b_ld_summary.xlsx   / .html
"""
# /// script
# requires-python = ">=3.10"
# dependencies = ["openpyxl", "pandas"]
# ///
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

SRC = Path("/Users/mikeford/dpyd-arm-model/output/sensitivity")
OUT = Path("/tmp/supplement-tables")
OUT.mkdir(exist_ok=True, parents=True)

VARIANT_NAME = {
    "WT_Clean": "WT (clean)",
    "WT_Biased": "WT (biased)",
    "HapB3": "HapB3",
    "2846hetho": "c.2846A>T",
    "2Ahetho": "*2A",
    "13hetho": "*13",
}

HEADER_FILL = PatternFill("solid", fgColor="1F3864")
HEADER_FONT = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
TITLE_FONT = Font(name="Calibri", size=12, bold=True)
CAPTION_FONT = Font(name="Calibri", size=9, italic=True, color="404040")
FOOTNOTE_FONT = Font(name="Calibri", size=9, color="404040")
BODY_FONT = Font(name="Calibri", size=10)
TOTAL_FILL = PatternFill("solid", fgColor="EDEDED")
TOTAL_FONT = Font(name="Calibri", size=10, bold=True)

THIN = Side(border_style="thin", color="808080")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
LEFT = Alignment(horizontal="left", vertical="center", wrap_text=True)
TITLE_ALIGN = Alignment(horizontal="left", vertical="center", wrap_text=True)


def style_header_row(ws, row, n_cols):
    for c in range(1, n_cols + 1):
        cell = ws.cell(row=row, column=c)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = CENTER
        cell.border = BORDER


def style_body_cell(cell, align="center", bold=False):
    cell.font = TOTAL_FONT if bold else BODY_FONT
    cell.alignment = CENTER if align == "center" else LEFT
    cell.border = BORDER


def autosize(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


def write_title_block(ws, title, caption, n_cols):
    ws.cell(row=1, column=1, value=title).font = TITLE_FONT
    ws.cell(row=1, column=1).alignment = TITLE_ALIGN
    ws.merge_cells(start_row=1, end_row=1, start_column=1, end_column=n_cols)
    ws.cell(row=2, column=1, value=caption).font = CAPTION_FONT
    ws.cell(row=2, column=1).alignment = TITLE_ALIGN
    ws.merge_cells(start_row=2, end_row=2, start_column=1, end_column=n_cols)
    ws.row_dimensions[2].height = 32


def write_footnotes(ws, start_row, n_cols, notes):
    for i, note in enumerate(notes):
        r = start_row + i
        ws.cell(row=r, column=1, value=note).font = FOOTNOTE_FONT
        ws.cell(row=r, column=1).alignment = TITLE_ALIGN
        ws.merge_cells(start_row=r, end_row=r, start_column=1, end_column=n_cols)


HTML_CSS = """
<style>
body { font-family: 'Calibri', 'Helvetica', sans-serif; max-width: 1100px; margin: 24px auto; color: #222; }
h2 { font-size: 14pt; margin-bottom: 4px; }
p.caption { font-size: 10pt; color: #555; margin-top: 0; margin-bottom: 12px; font-style: italic; }
table { border-collapse: collapse; font-size: 10pt; margin-bottom: 8px; }
th { background: #1F3864; color: white; padding: 6px 10px; border: 1px solid #555; text-align: center; vertical-align: middle; }
td { padding: 5px 10px; border: 1px solid #aaa; text-align: center; vertical-align: middle; }
td.left { text-align: left; }
tr.total td { background: #EDEDED; font-weight: bold; }
p.footnote { font-size: 9pt; color: #555; margin: 4px 0; }
</style>
"""


def html_page(title, body):
    return (
        f"<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>{title}</title>{HTML_CSS}</head><body>{body}</body></html>"
    )


def yn(v):
    s = str(v).strip().upper()
    if s == "TRUE":
        return "Yes"
    if s == "FALSE":
        return "No"
    return "—"


def build_s8():
    df = pd.read_csv(SRC / "genotyping" / "genotyping_combined_summary.csv")
    diag_start = df.index[df["comparison"] == "---DIAGNOSTICS---"][0]
    diagnostics = df.iloc[diag_start + 1:].copy()
    body = df.iloc[:diag_start].copy()

    def label_pair(s):
        suffix = ""
        s_clean = s
        if "(" in s:
            i = s.index("(")
            s_clean, suffix = s[:i].strip(), s[i:].strip()
        a, b = [t.strip() for t in s_clean.split(" vs ")]
        out = f"{VARIANT_NAME.get(a, a)} vs {VARIANT_NAME.get(b, b)}"
        if suffix:
            out += f" {suffix}"
        return out

    body["Comparison"] = body["comparison"].map(label_pair)

    def fmt_sucra(v):
        try:
            f = float(v)
        except (TypeError, ValueError):
            return "—"
        if f != f:  # NaN check
            return "—"
        return f"{f:.1f}"

    display_headers = [
        "Comparison",
        "OR", "95% CrI", "Excludes 1?", "SUCRA (%)",
        "OR", "95% CrI", "Excludes 1?", "SUCRA (%)",
    ]

    def render_row(r):
        return [
            r["Comparison"],
            f"{float(r['OR_loose']):.2f}",
            str(r["CrI_loose"]),
            yn(r["CrI_excludes_1_loose"]),
            fmt_sucra(r["SUCRA_loose"]),
            f"{float(r['OR_strict']):.2f}",
            str(r["CrI_strict"]),
            yn(r["CrI_excludes_1_strict"]),
            fmt_sucra(r["SUCRA_strict"]),
        ]

    rows = [render_row(r) for _, r in body.iterrows()]
    diag = {row["comparison"]: (row["OR_loose"], row["OR_strict"])
            for _, row in diagnostics.iterrows()}

    title = "Table S8. Genotyping Misclassification Sensitivity Analysis (WT-Split Model)"
    caption = (
        "Network meta-analysis after splitting WT into WT (clean) and WT (biased) nodes to test whether "
        "genotyping completeness drives the primary finding. Two WT (clean) definitions are reported: "
        "loose (n = 9 cohorts; includes studies with pan-negative-by-design WT arms plus three studies "
        "with functionally pan-negative WT controls — Amstutz_2009, Jennings_2013, Froehlich_2015) and "
        "strict (n = 6 cohorts; restricted to studies with pre-screening exclusion of variant carriers). "
        "WT (biased) comprises studies that did not test for all four clinically actionable variants. "
        "OR = odds ratio; CrI = credible interval; SUCRA = surface under the cumulative ranking curve."
    )

    wb = Workbook()
    ws = wb.active
    ws.title = "Table S8"
    n_cols = len(display_headers)
    write_title_block(ws, title, caption, n_cols)

    ws.merge_cells(start_row=3, end_row=4, start_column=1, end_column=1)
    cell = ws.cell(row=3, column=1, value="Comparison")
    cell.fill = HEADER_FILL
    cell.font = HEADER_FONT
    cell.alignment = CENTER
    cell.border = BORDER

    ws.cell(row=3, column=2, value="Loose definition (n = 9 WT-clean cohorts)")
    ws.merge_cells(start_row=3, end_row=3, start_column=2, end_column=5)
    ws.cell(row=3, column=6, value="Strict definition (n = 6 WT-clean cohorts)")
    ws.merge_cells(start_row=3, end_row=3, start_column=6, end_column=9)
    for c in (2, 6):
        cell = ws.cell(row=3, column=c)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = CENTER
        cell.border = BORDER

    for i, h in enumerate(display_headers, start=1):
        if i == 1:
            continue
        cell = ws.cell(row=4, column=i, value=h)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = CENTER
        cell.border = BORDER

    r = 5
    for row in rows:
        for c, v in enumerate(row, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            style_body_cell(cell, align="left" if c == 1 else "center")
        r += 1

    notes = [
        (
            f"Convergence (loose): DIC = {diag['DIC'][0]}, max PSRF = {diag['max_PSRF'][0]}, "
            f"converged = {yn(diag['converged'][0])}.  "
            f"Convergence (strict): DIC = {diag['DIC'][1]}, max PSRF = {diag['max_PSRF'][1]}, "
            f"converged = {yn(diag['converged'][1])}."
        ),
        ("WT (clean) vs WT (biased) is reported as a diagnostic comparison to assess whether genotyping "
         "completeness drives outcome heterogeneity; SUCRA is not meaningful for this two-WT contrast."),
        ("PSRF = potential scale reduction factor (Gelman-Rubin diagnostic, threshold ≤ 1.10). "
         "DIC = deviance information criterion."),
    ]
    write_footnotes(ws, r + 1, n_cols, notes)
    autosize(ws, [38, 9, 14, 12, 11, 9, 14, 12, 11])
    ws.row_dimensions[3].height = 22
    ws.row_dimensions[4].height = 22
    wb.save(OUT / "Table_S8_genotyping.xlsx")

    html_rows = [
        "<tr>"
        "<th rowspan='2'>Comparison</th>"
        "<th colspan='4'>Loose definition (n = 9 WT-clean cohorts)</th>"
        "<th colspan='4'>Strict definition (n = 6 WT-clean cohorts)</th>"
        "</tr>",
        "<tr>" + "".join(f"<th>{h}</th>" for h in display_headers[1:]) + "</tr>",
    ]
    for row in rows:
        cells = "".join(
            f"<td class='left'>{v}</td>" if i == 0 else f"<td>{v}</td>"
            for i, v in enumerate(row)
        )
        html_rows.append(f"<tr>{cells}</tr>")
    table_html = "<table>" + "".join(html_rows) + "</table>"
    body_html = f"<h2>{title}</h2><p class='caption'>{caption}</p>{table_html}"
    for n in notes:
        body_html += f"<p class='footnote'>{n}</p>"
    (OUT / "Table_S8_genotyping.html").write_text(html_page(title, body_html))


def build_s9():
    df = pd.read_csv(SRC / "genotyping" / "qba_table.csv")

    def check(v):
        return "✓" if int(v) == 1 else "—"

    df_disp = pd.DataFrame({
        "Study": df["Study"],
        "WT arm N": df["WT_N"].astype(int),
        "WT events": df["WT_R"].astype(int),
        "WT event rate": df["WT_Event_Rate"].map(lambda x: f"{x:.3f}"),
        "HapB3": df["Tested_HapB3"].map(check),
        "*2A": df["Tested_2A"].map(check),
        "c.2846A>T": df["Tested_2846"].map(check),
        "*13": df["Tested_13"].map(check),
        "Untested variants": df["Untested_Variants"].map(
            lambda x: "—" if str(x).strip() in ("None", "nan", "") else str(x)),
        "Expected misclassified (n)": df["Expected_Misclassified"].map(lambda x: f"{float(x):.1f}"),
        "% of WT arm": df["Pct_of_WT_Arm"].map(lambda x: f"{float(x):.1f}"),
        "Max event-rate shift": df["Max_Event_Rate_Shift"].map(lambda x: f"{float(x):.4f}"),
    })

    total_n = int(df["WT_N"].sum())
    total_r = int(df["WT_R"].sum())
    total_misclassified = float(df["Expected_Misclassified"].sum())
    total_row = {
        "Study": "Total",
        "WT arm N": str(total_n),
        "WT events": str(total_r),
        "WT event rate": f"{total_r / total_n:.3f}",
        "HapB3": "", "*2A": "", "c.2846A>T": "", "*13": "",
        "Untested variants": "",
        "Expected misclassified (n)": f"{total_misclassified:.1f}",
        "% of WT arm": f"{100 * total_misclassified / total_n:.1f}",
        "Max event-rate shift": "",
    }

    title = "Table S9. Quantitative Bias Analysis for Wild-Type Genotyping Completeness"
    caption = (
        "Per-study estimate of patients potentially misclassified as wild-type because the study did not "
        "test for all four clinically actionable variants (HapB3, *2A, c.2846A>T, *13). Expected misclassified count "
        "= WT arm N × Σ carrier frequency of each untested variant. The maximum event-rate shift quantifies the "
        "expected upward movement in the study's observed WT event rate if all expected misclassified patients "
        "experienced events at the variant-specific rates implied by the primary network meta-analysis."
    )
    headers = list(df_disp.columns)
    n_cols = len(headers)

    wb = Workbook()
    ws = wb.active
    ws.title = "Table S9"
    write_title_block(ws, title, caption, n_cols)
    for i, h in enumerate(headers, start=1):
        ws.cell(row=3, column=i, value=h)
    style_header_row(ws, 3, n_cols)

    r = 4
    for _, row in df_disp.iterrows():
        for c, h in enumerate(headers, start=1):
            cell = ws.cell(row=r, column=c, value=row[h])
            style_body_cell(cell, align="left" if h == "Study" else "center")
        r += 1
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=r, column=c, value=total_row[h])
        cell.fill = TOTAL_FILL
        cell.font = TOTAL_FONT
        cell.alignment = LEFT if h == "Study" else CENTER
        cell.border = BORDER
    r += 1

    notes = [
        "✓ indicates the variant was genotyped in the WT arm; — indicates not tested. "
        "European carrier frequencies used to compute expected misclassified counts (CPIC 2017): "
        "HapB3 = 4.7%, *2A = 1.6%, c.2846A>T = 0.7%, *13 = 0.1%.",
        "Expected misclassified = WT arm N × Σᵥ (carrier frequency of untested variant v). "
        "% of WT arm = expected misclassified / WT arm N × 100. "
        "Max event-rate shift = Σᵥ [carrier frequency of v × (expected event rate of v − pooled "
        "WT event rate, 0.323)], where expected event rate of v is derived by applying the primary "
        "NMA odds ratio for variant v to the pooled WT baseline.",
    ]
    write_footnotes(ws, r + 1, n_cols, notes)
    autosize(ws, [26, 10, 10, 12, 9, 7, 12, 7, 28, 14, 11, 16])
    wb.save(OUT / "Table_S9_qba.xlsx")

    html_rows = ["<tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr>"]
    for _, row in df_disp.iterrows():
        cells = "".join(
            f"<td class='left'>{row[h]}</td>" if h == "Study" else f"<td>{row[h]}</td>"
            for h in headers
        )
        html_rows.append(f"<tr>{cells}</tr>")
    tot_cells = "".join(
        f"<td class='left'>{total_row[h]}</td>" if h == "Study" else f"<td>{total_row[h]}</td>"
        for h in headers
    )
    html_rows.append(f"<tr class='total'>{tot_cells}</tr>")
    table_html = "<table>" + "".join(html_rows) + "</table>"
    body_html = f"<h2>{title}</h2><p class='caption'>{caption}</p>{table_html}"
    for n in notes:
        body_html += f"<p class='footnote'>{n}</p>"
    (OUT / "Table_S9_qba.html").write_text(html_page(title, body_html))


def build_s10a():
    df = pd.read_csv(SRC / "ld" / "ld_study_breakdown.csv")

    def i(x):
        try:
            return str(int(round(float(x))))
        except (TypeError, ValueError):
            return str(x)

    df_disp = pd.DataFrame({
        "Study": df["Study"],
        "SNP genotyped": df["SNP_Used"],
        "Original N": df["HapB3_N"].map(i),
        "Original events": df["HapB3_R"].map(i),
        "Event rate": df["Event_Rate"].map(lambda x: f"{x:.3f}"),
        "Expected discordant": df["Expected_Discordant"].map(lambda x: f"{float(x):.2f}"),
        "Best-case N": df["Best_N"].map(i),
        "Best-case events": df["Best_R"].map(i),
        "Worst-case N": df["Worst_N"].map(i),
        "Worst-case events": df["Worst_R"].map(i),
        "LD adjustment applied": df["LD_Adjustment_Applied"],
    })

    total_n = int(df["HapB3_N"].sum())
    total_r = int(df["HapB3_R"].sum())
    total_disc = float(df["Expected_Discordant"].sum())
    total_best_n = int(df["Best_N"].sum())
    total_best_r = int(df["Best_R"].sum())
    total_worst_n = int(df["Worst_N"].sum())
    total_worst_r = int(df["Worst_R"].sum())
    total_row = {
        "Study": "Total",
        "SNP genotyped": "",
        "Original N": str(total_n),
        "Original events": str(total_r),
        "Event rate": f"{total_r / total_n:.3f}",
        "Expected discordant": f"{total_disc:.2f}",
        "Best-case N": str(total_best_n),
        "Best-case events": str(total_best_r),
        "Worst-case N": str(total_worst_n),
        "Worst-case events": str(total_worst_r),
        "LD adjustment applied": "",
    }

    title = "Table S10a. Linkage Disequilibrium Sensitivity — Per-Study HapB3 Breakdown"
    caption = (
        "HapB3 arm characteristics for the 16 studies contributing HapB3 data. "
        "Studies genotyped the causal SNP c.1129−5923C>G, the tag SNP c.1236G>A, or both. "
        "Expected discordant = HapB3 N × 1/300 (tag-SNP discordance rate). "
        "Best-case adjustment removes the discordant patient and treats them as a non-event; "
        "worst-case treats them as an event. Adjustment was applied only when expected discordant ≥ 0.5."
    )
    headers = list(df_disp.columns)
    n_cols = len(headers)
    left_cols = {"Study", "SNP genotyped", "LD adjustment applied"}

    wb = Workbook()
    ws = wb.active
    ws.title = "Table S10a"
    write_title_block(ws, title, caption, n_cols)
    for i_, h in enumerate(headers, start=1):
        ws.cell(row=3, column=i_, value=h)
    style_header_row(ws, 3, n_cols)

    r = 4
    for _, row in df_disp.iterrows():
        for c, h in enumerate(headers, start=1):
            cell = ws.cell(row=r, column=c, value=row[h])
            style_body_cell(cell, align="left" if h in left_cols else "center")
        r += 1
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=r, column=c, value=total_row[h])
        cell.fill = TOTAL_FILL
        cell.font = TOTAL_FONT
        cell.alignment = LEFT if h in left_cols else CENTER
        cell.border = BORDER
    r += 1

    notes = [
        "Discordant patients allocated to Medwid_2023 (scenario A) and Wigle_2021 (scenario B) — the two largest "
        "tag-SNP studies — for the LD-adjusted re-analyses (see Table S10b).",
        "Tag SNP c.1236G>A is in tight LD with causal HapB3 (D' = 1, r² ≈ 1 in European populations); "
        "approximately 1 in 300 carriers is discordant.",
    ]
    write_footnotes(ws, r + 1, n_cols, notes)
    autosize(ws, [26, 30, 11, 13, 11, 14, 12, 13, 13, 14, 22])
    wb.save(OUT / "Table_S10a_ld_breakdown.xlsx")

    html_rows = ["<tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr>"]
    for _, row in df_disp.iterrows():
        cells = "".join(
            f"<td class='left'>{row[h]}</td>" if h in left_cols else f"<td>{row[h]}</td>"
            for h in headers
        )
        html_rows.append(f"<tr>{cells}</tr>")
    tot_cells = "".join(
        f"<td class='left'>{total_row[h]}</td>" if h in left_cols else f"<td>{total_row[h]}</td>"
        for h in headers
    )
    html_rows.append(f"<tr class='total'>{tot_cells}</tr>")
    table_html = "<table>" + "".join(html_rows) + "</table>"
    body_html = f"<h2>{title}</h2><p class='caption'>{caption}</p>{table_html}"
    for n in notes:
        body_html += f"<p class='footnote'>{n}</p>"
    (OUT / "Table_S10a_ld_breakdown.html").write_text(html_page(title, body_html))


def build_s10b():
    df = pd.read_csv(SRC / "ld" / "ld_summary.csv")
    primary_or, primary_cri_low, primary_cri_high = 2.005, 1.294, 3.186
    primary_psrf = 1.02854
    primary_dic = 151.492

    scen_name = {
        "best_medwid": "Best case (Medwid_2023, no event)",
        "best_wigle": "Best case (Wigle_2021, no event)",
        "worst_medwid": "Worst case (Medwid_2023, event)",
        "worst_wigle": "Worst case (Wigle_2021, event)",
    }

    rows = [[
        "Primary (no LD adjustment)",
        f"{primary_or:.3f}",
        f"{primary_cri_low:.2f}–{primary_cri_high:.2f}",
        "—",
        f"{primary_psrf:.3f}",
        f"{primary_dic:.2f}",
        "Yes",
        "No",
    ]]
    for _, r in df.iterrows():
        rows.append([
            scen_name.get(r["scenario"], r["scenario"]),
            f"{float(r['HapB3_OR']):.3f}",
            f"{float(r['CrI_low']):.2f}–{float(r['CrI_high']):.2f}",
            f"{float(r['delta_from_primary']):+.3f}",
            f"{float(r['PSRF']):.3f}",
            f"{float(r['DIC']):.2f}",
            yn(r["converged"]),
            yn(r["escalated"]),
        ])

    headers = [
        "Scenario",
        "HapB3 OR",
        "95% CrI",
        "Δ log-OR from primary",
        "Max PSRF",
        "DIC",
        "Converged",
        "Escalated",
    ]

    title = "Table S10b. Linkage Disequilibrium Sensitivity — HapB3 Odds Ratio Across LD Scenarios"
    caption = (
        "HapB3 vs wild-type odds ratio under best-case and worst-case allocation of the single expected discordant "
        "patient identified via tag-SNP genotyping (see Table S10a). Allocation tested against the two largest "
        "tag-SNP studies (Medwid_2023 and Wigle_2021). All four LD-adjusted scenarios converged (PSRF ≤ 1.10) "
        "and remained within ±0.04 log-OR of the primary estimate."
    )
    n_cols = len(headers)

    wb = Workbook()
    ws = wb.active
    ws.title = "Table S10b"
    write_title_block(ws, title, caption, n_cols)
    for i_, h in enumerate(headers, start=1):
        ws.cell(row=3, column=i_, value=h)
    style_header_row(ws, 3, n_cols)
    r = 4
    for row in rows:
        is_primary = row[0].startswith("Primary")
        for c, v in enumerate(row, start=1):
            cell = ws.cell(row=r, column=c, value=v)
            if is_primary:
                cell.fill = TOTAL_FILL
                cell.font = TOTAL_FONT
            else:
                cell.font = BODY_FONT
            cell.alignment = LEFT if c == 1 else CENTER
            cell.border = BORDER
        r += 1
    notes = [
        "Δ log-OR from primary = log(HapB3 OR scenario) − log(HapB3 OR primary). Values within ±0.10 indicate "
        "consistency with the primary analysis.",
        "PSRF = potential scale reduction factor (Gelman-Rubin diagnostic, convergence threshold ≤ 1.10). "
        "DIC = deviance information criterion.",
        "Escalated = scenario was re-run with the primary model's full MCMC sampling configuration after the "
        "reduced sensitivity-run sampling failed convergence; none of the LD scenarios required escalation.",
    ]
    write_footnotes(ws, r + 1, n_cols, notes)
    autosize(ws, [38, 11, 14, 18, 11, 11, 12, 12])
    wb.save(OUT / "Table_S10b_ld_summary.xlsx")

    html_rows = ["<tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr>"]
    for row in rows:
        is_primary = row[0].startswith("Primary")
        cls = " class='total'" if is_primary else ""
        cells = "".join(
            f"<td class='left'>{v}</td>" if i == 0 else f"<td>{v}</td>"
            for i, v in enumerate(row)
        )
        html_rows.append(f"<tr{cls}>{cells}</tr>")
    table_html = "<table>" + "".join(html_rows) + "</table>"
    body_html = f"<h2>{title}</h2><p class='caption'>{caption}</p>{table_html}"
    for n in notes:
        body_html += f"<p class='footnote'>{n}</p>"
    (OUT / "Table_S10b_ld_summary.html").write_text(html_page(title, body_html))


def main():
    build_s8()
    build_s9()
    build_s10a()
    build_s10b()
    for f in sorted(OUT.iterdir()):
        print(f.name, f.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
