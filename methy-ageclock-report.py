#!/usr/bin/env python3

import base64
import csv
import heapq
import html
import math
import re
from collections import Counter
from pathlib import Path

CHIP_ID = "209648820170"
SAMPLE_TYPE = "Buffy coat"
PLATFORM = "Infinium MethylationEPIC v2"

PANDU_LOGO = "pandu-logo.png"
PANDU_SLOGAN = "pandu-slogan.png"
METADATA_FILE = "metadata.csv"
REPORT_NAME = "Methylation_Analysis_Report.html"

CSS = r"""
:root {
  color-scheme: light;
  font: 15px/1.5 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  --red: #d24645;
  --ink: #252321;
  --muted: #586067;
  --line: #e3e6e8;
  --soft: #f7f8f7;
}
* { box-sizing: border-box; }
body { margin: 0; color: var(--ink); background: #f4f5f4; }
.sheet {
  max-width: 900px; margin: 32px auto; background: #fff; border: 1px solid var(--line);
  border-radius: 12px; box-shadow: 0 12px 35px #25313a0b; overflow: hidden;
}
.brandline {
  padding: 24px 38px; display: flex; justify-content: space-between; align-items: center;
  gap: 20px; border-bottom: 1px solid var(--line);
}
.partner-logo { display: block; width: 96px; height: auto; }
.brand-logo { display: block; width: 164px; height: auto; }
.hero { padding: 31px 38px 27px; border-top: 3px solid var(--red); }
h1 {
  font-size: clamp(1.85rem, 4.2vw, 2.6rem); letter-spacing: -.035em;
  line-height: 1.15; margin: 0 0 10px; max-width: 720px;
}
.sample {
  display: flex; width: fit-content; flex-wrap: wrap; align-items: baseline; gap: 8px 14px;
  margin-top: 12px; border: 1px solid var(--line); border-radius: 7px;
  padding: 9px 12px; font-size: .83rem; overflow-wrap: anywhere;
}
.sample span { color: var(--muted); }
.sample strong { font-weight: 680; }
.hero-cards {
  display: grid; grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 12px; margin-top: 24px;
}
.hero-cards .sample {
  width: 100%; min-width: 0; margin-top: 0;
  flex-direction: column; align-items: flex-start; gap: 4px;
}
.content { padding: 0 38px 34px; }
.section {
  border: 1px solid var(--line); border-top: 3px solid var(--red); border-radius: 9px;
  margin-top: 16px; overflow: hidden; break-inside: avoid;
}
.section-head { padding: 18px 20px 14px; }
.section h2 { margin: 0; font-size: 1.15rem; letter-spacing: -.02em; line-height: 1.25; }
.cards {
  display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px;
  padding: 0 20px 20px;
}
.card { min-width: 0; border: 1px solid var(--line); border-radius: 8px; padding: 14px; background: var(--soft); }
.card .label { color: var(--muted); font-size: .72rem; font-weight: 720; letter-spacing: .04em; }
.card .value { margin-top: 4px; font-size: .98rem; font-weight: 690; overflow-wrap: anywhere; }
.alert {
  margin: 0 20px 20px; padding: 12px 14px; border: 1px solid var(--line);
  border-left-width: 3px; border-radius: 7px; background: var(--soft); font-size: .88rem;
}
.table-wrap { width: 100%; overflow-x: auto; padding: 0 20px 20px; }
table { width: 100%; border-collapse: collapse; font-size: .82rem; }
th {
  text-align: left; color: var(--muted); font-size: .7rem;
  letter-spacing: .035em; padding: 9px 8px; border-bottom: 1px solid var(--line);
}
td { padding: 10px 8px; border-bottom: 1px solid var(--line); vertical-align: top; }
tbody tr:last-child td { border-bottom: 0; }
.num { text-align: right; white-space: nowrap; font-variant-numeric: tabular-nums; }
.center { text-align: center; }
.word { overflow-wrap: anywhere; }
.analysis-summary-table th:first-child, .analysis-summary-table td:first-child {
  width: 230px; min-width: 230px; overflow-wrap: normal; word-break: normal;
}
.analysis-summary-table td:nth-child(2) {
  text-align: justify;
  text-align-last: left;
  text-justify: inter-word;
  hyphens: auto;
  overflow-wrap: anywhere;
}
details.report {
  margin: 0 20px 20px; border: 1px solid var(--line); border-radius: 8px; overflow: hidden;
}
details.report summary { cursor: pointer; padding: 11px 13px; font-size: .84rem; font-weight: 650; }
details.report summary::-webkit-details-marker { display: none; }
details.report summary::after { content: "+"; float: right; color: var(--muted); }
details.report[open] summary::after { content: "−"; }
details.report .table-wrap { padding: 0 13px 13px; }
.software-versions-table th { text-transform: uppercase; }
details.report .software-versions-wrap { padding: 0 20px 20px; }
.figure { margin: 0 20px 20px; border: 1px solid var(--line); border-radius: 8px; padding: 12px; }
.figure img { display: block; width: 100%; height: auto; }
.figure-guide { margin-top: 12px; padding-top: 10px; border-top: 1px solid var(--line); color: var(--muted); font-size: .85rem; }
.figure-guide strong { color: var(--ink); }
@media (max-width: 700px) {
  .sheet { margin: 0; border: 0; border-radius: 0; box-shadow: none; }
  .brandline, .hero, .content { padding-left: 20px; padding-right: 20px; }
  .partner-logo { width: 84px; }
  .brand-logo { width: 138px; }
  .cards { grid-template-columns: 1fr; }
  .hero-cards { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  /* Stack long method descriptions on small screens without breaking words. */
  .analysis-summary-table, .analysis-summary-table tbody,
  .analysis-summary-table tr, .analysis-summary-table td {
    display: block; width: 100%;
  }
  .analysis-summary-table thead { display: none; }
  .analysis-summary-table td:first-child {
    min-width: 0; font-weight: 650; border-bottom: 0; padding-bottom: 2px;
  }
  .analysis-summary-table td:nth-child(2) { padding-top: 2px; }
  .analysis-summary-table tbody tr { border-bottom: 1px solid var(--line); }
  .analysis-summary-table tbody tr:last-child { border-bottom: 0; }
}
@media (max-width: 390px) {
  .brandline { align-items: flex-start; flex-direction: column; gap: 8px; }
  .hero-cards { grid-template-columns: 1fr; }
}
@media print {
  @page { size: A4; margin: 10mm; }
  body { background: white; }
  .sheet { margin: 0; max-width: none; border: 0; box-shadow: none; }
  .brandline { padding-top: 4px; padding-bottom: 14px; }
  .hero { padding-top: 21px; padding-bottom: 20px; }
  .content { padding-bottom: 17px; }
  .section, .card, details.report, .figure {
    print-color-adjust: exact; -webkit-print-color-adjust: exact;
  }
  .section { break-inside: auto; }
  .figure { break-inside: avoid; }
  details.report > .table-wrap { display: block !important; }
}
"""

# Report labels
QC_TITLE = "Quality Control & Preprocessing"
SAMPLE_QC_TITLE = "Sample QC"
BETA_DISTRIBUTION_TITLE = "Beta-value Distribution"
PCA_TITLE = "Principal Component Analysis"
CLOCK_TITLE = "DNA Methylation Clock"
AGE_DIFF_TITLE = "Epigenetic Age Difference"

FDR_CUTOFF = 0.05
TOP_ROWS = 10

GUIDE_SAMPLE_QC = (
    "Call Rate ≥95% indicates good probe detection. "
    "Valid CpGs ≥95% indicates high data completeness after preprocessing."
)
GUIDE_BETA_DISTRIBUTION = (
    "DNA methylation typically shows a bimodal distribution, with two peaks "
    "near β = 0 (unmethylated) and β = 1 (methylated)."
)
GUIDE_PCA = "Samples with similar DNA methylation patterns tend to cluster together in Principal Component Analysis."

FILES = {
    "sample_qc": "1_qc_preprocessing/01_sample_qc.png",
    "beta_distribution": "1_qc_preprocessing/02_beta_distribution.png",
    "pca": "1_qc_preprocessing/03_pca.png",
    "dmp": "2_diff_methylation/04_dmp_results.csv",
    "volcano": "2_diff_methylation/05_dmp_volcano.png",
    "dmr": "2_diff_methylation/06_dmr_results.csv",
    "comparison_info": "2_diff_methylation/comparison_info.csv",
    "clock": "3_epigenetic_age/08_methylation_clock.png",
    "age_difference": "3_epigenetic_age/09_age_difference.png",
    "versions": "software_versions/software_versions.csv",
}

def esc(value):
    return html.escape(str(value), quote=True)

def image_data_url(path):
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"

def info_cards(items):
    return (
        '<div class="cards">'
        + "".join(
            f'<div class="card"><div class="label">{esc(label)}</div>'
            f'<div class="value">{esc(value)}</div></div>'
            for label, value in items
        )
        + '</div>'
    )

def notice(message):
    return f'<div class="alert">{esc(message)}</div>'

def section(title, body):
    return (
        '<section class="section">'
        f'<div class="section-head"><h2>{esc(title)}</h2></div>'
        + body + '</section>'
    )

def collapsible_html(title, body):
    return (
        '<details class="report">'
        f'<summary>{esc(title)}</summary>'
        f'{body}</details>'
    )

def display_cell_value(value):
    if value is None or str(value).strip().casefold() in {"", "na", "n/a", "nan", "none"}:
        return ""
    return value

def html_table(headers, rows):
    head = ''.join(f'<th class="{rows[0][i][1]}">{esc(h)}</th>' for i, h in enumerate(headers))
    body = ''.join(
        '<tr>' + ''.join(
            f'<td class="{cls}">{esc(display_cell_value(value))}</td>'
            for value, cls in row
        ) + '</tr>'
        for row in rows
    )
    return (
        '<div class="table-wrap"><table><thead><tr>' + head + '</tr></thead>'
        '<tbody>' + body + '</tbody></table></div>'
    )

def figure(path, alt_text, guide=None):
    caption = f'<figcaption class="figure-guide"><strong>Guide:</strong> {esc(guide)}</figcaption>' if guide else ''
    return (
        '<figure class="figure">'
        f'<img src="{image_data_url(path)}" alt="{esc(alt_text)}">'
        f'{caption}</figure>'
    )

def rows_from_csv(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        yield from reader

def number(value):
    try:
        val = float(value)
        return val if math.isfinite(val) else math.nan
    except (TypeError, ValueError):
        return math.nan

def fmt(value, digits=3):
    val = number(value)
    return f'{val:.{digits}f}' if math.isfinite(val) else 'N/A'

def fmt_p(value):
    val = number(value)
    if not math.isfinite(val):
        return 'N/A'
    return f'{val:.2e}' if 0 < abs(val) < 0.0001 else f'{val:.5f}'

def load_comparison_info(path):
    rows = list(rows_from_csv(path, ("reference_group", "comparison_group")))
    reference_group = (rows[0]['reference_group'] or '').strip()
    comparison_group = (rows[0]['comparison_group'] or '').strip()
    return reference_group, comparison_group

def load_metadata(path):
    metadata = list(rows_from_csv(path, ("sample", "group", "age")))
    return metadata

def load_dmp(path):
    required = ("CpG", "seqnames", "position", "associated.genes",
                "gene.regions", "mean_beta_reference", "mean_beta_comparison",
                "delta_beta", "logFC_M", "p_value", "FDR")
    total, significant, lowest, highest = 0, 0, [], []
    for row in rows_from_csv(path, required):
        total += 1
        fdr = number(row["FDR"])
        if fdr < FDR_CUTOFF:
            significant += 1
            logfc = number(row["logFC_M"])
            if logfc < 0:
                rank = (-logfc, -total)
                if len(lowest) < TOP_ROWS:
                    heapq.heappush(lowest, (rank, row))
                elif rank > lowest[0][0]:
                    heapq.heapreplace(lowest, (rank, row))
            elif logfc > 0:
                rank = (logfc, -total)
                if len(highest) < TOP_ROWS:
                    heapq.heappush(highest, (rank, row))
                elif rank > highest[0][0]:
                    heapq.heapreplace(highest, (rank, row))
    bottom = [row for _, row in sorted(lowest, key=lambda item: item[0], reverse=True)]
    top = [row for _, row in sorted(highest, key=lambda item: item[0], reverse=True)]
    return total, significant, bottom, top

def load_dmr(path):
    required = ("DMR", "seqnames", "start", "end", "no.cpgs",
                "meandiff", "min_smoothed_fdr", "overlapping.genes")
    total, lowest, highest = 0, [], []
    for row in rows_from_csv(path, required):
        total += 1
        mean_diff = number(row["meandiff"])
        if mean_diff < 0:
            rank = (-mean_diff, -total)
            if len(lowest) < TOP_ROWS:
                heapq.heappush(lowest, (rank, row))
            elif rank > lowest[0][0]:
                heapq.heapreplace(lowest, (rank, row))
        elif mean_diff > 0:
            rank = (mean_diff, -total)
            if len(highest) < TOP_ROWS:
                heapq.heappush(highest, (rank, row))
            elif rank > highest[0][0]:
                heapq.heapreplace(highest, (rank, row))
    bottom = [row for _, row in sorted(lowest, key=lambda item: item[0], reverse=True)]
    top = [row for _, row in sorted(highest, key=lambda item: item[0], reverse=True)]
    return total, bottom, top

def metadata_section(metadata):
    counts = Counter((row['group'] or '(empty)').strip() for row in metadata)
    cards = [(f'{group} Samples', count) for group, count in sorted(counts.items())]

    columns = [('sample', 'Sample ID'), ('group', 'Group'), ('age', 'Age (years)'), ('sex', 'Sex')]
    if any(row.get('sample_name(opsional)', '').strip() for row in metadata):
        columns.insert(1, ('sample_name(opsional)', 'Sample Name'))
    if any(row.get('batch(opsional)', '').strip() for row in metadata):
        columns.append(('batch(opsional)', 'Batch'))
    cells = [[(row.get(key, '') or '',
               'center' if key in ('age', 'sex', 'batch(opsional)') else 'word')
              for key, _ in columns]
             for row in metadata]
    return section('Sample Information',
                   info_cards(cards)
                   + collapsible_html('Sample Metadata',
                                      html_table([title for _, title in columns], cells)))

def qc_section(folder):
    body = (
        figure(folder / FILES['sample_qc'], f'{SAMPLE_QC_TITLE}: Call Rate and Valid CpGs', GUIDE_SAMPLE_QC)
        + figure(folder / FILES['beta_distribution'], BETA_DISTRIBUTION_TITLE, GUIDE_BETA_DISTRIBUTION)
        + figure(folder / FILES['pca'], PCA_TITLE, GUIDE_PCA)
    )
    return section(QC_TITLE, body)

def dmp_section(folder, result, reference_group, comparison_group):
    tested, significant, lowest, highest = result
    body = info_cards([
        ('CpGs Tested for DMPs', f'{tested:,}'),
        (f'Significant DMPs (FDR < {FDR_CUTOFF:g})', f'{significant:,}'),
    ])
    for title, selected, direction in (
        (f'Top {TOP_ROWS} Lowest logFC M (Lower in {comparison_group})', lowest, 'negative'),
        (f'Top {TOP_ROWS} Highest logFC M (Higher in {comparison_group})', highest, 'positive'),
    ):
        if selected:
            rows = [[
                (r['CpG'], 'word'),
                (r['seqnames'], 'word'),
                (r['position'], 'num'),
                (fmt(r['delta_beta']), 'num'),
                (fmt(r['logFC_M']), 'num'),
                (fmt_p(r['FDR']), 'num'),
                (r['associated.genes'] or '', 'word'),
            ] for r in selected]
            body += collapsible_html(
                title,
                html_table(('CpG', 'Chr', 'Position (hg38)', 'Δβ', 'logFC M',
                            'FDR', 'Associated Genes'), rows)
            )
        else:
            comparison_direction = 'lower' if direction == 'negative' else 'higher'
            body += notice(
                f'No significant DMPs with {direction} logFC M '
                f'({comparison_direction} in {comparison_group}).'
            )
    body += figure(folder / FILES['volcano'], f'DMP Volcano Plot — {comparison_group} vs {reference_group}')
    return section('Differentially Methylated Positions (DMP)', body)

def dmr_section(result, comparison_group):
    total, lowest, highest = result
    body = info_cards([('DMRs Identified', f'{total:,}')])
    if total == 0:
        body += notice('No DMRs were reported in the results file.')
    else:
        for title, selected, direction in (
            (f'Top {TOP_ROWS} Lowest Mean Diff (Lower in {comparison_group})', lowest, 'negative'),
            (f'Top {TOP_ROWS} Highest Mean Diff (Higher in {comparison_group})', highest, 'positive'),
        ):
            if selected:
                rows = [[
                    (r['DMR'], 'word'),
                    (r['seqnames'], ''),
                    (f"{r['start']}–{r['end']}", 'word'),
                    (r['no.cpgs'], 'num'),
                    (fmt(r['meandiff']), 'num'),
                    (fmt_p(r['min_smoothed_fdr']), 'num'),
                    (r['overlapping.genes'] or '', 'word'),
                ] for r in selected]
                body += collapsible_html(
                    title,
                    html_table(('DMR', 'Chr', 'Coordinates (hg38)', 'CpG Count',
                                'Mean Diff', 'Min Smoothed FDR', 'Overlapping Genes'), rows)
                )
            else:
                comparison_direction = 'lower' if direction == 'negative' else 'higher'
                body += notice(
                    f'No reported DMRs with {direction} Mean Diff '
                    f'({comparison_direction} in {comparison_group}).'
                )
    return section('Differentially Methylated Regions (DMR)', body)

def clock_section(folder):
    body = figure(folder / FILES['clock'], f'{CLOCK_TITLE}: DNAm Age vs Actual Age')
    return section(CLOCK_TITLE, body)

def age_diff_section(folder):
    body = figure(folder / FILES['age_difference'], f'{AGE_DIFF_TITLE}: DNAm Age − Actual Age (years)')
    return section(AGE_DIFF_TITLE, body)

def normalize_version(value):
    raw = str(value).strip()
    if not raw or raw.upper() == "N/A":
        return "N/A"
    match = re.search(
        r"(?<![0-9A-Za-z])v?(\d+(?:\.\d+)+(?:[-+._]?[0-9A-Za-z]+)*)",
        raw,
        re.IGNORECASE,
    )
    if not match:
        return "N/A"
    return f"v{match.group(1)}"

def software_versions_section(path):
    software_functions = {
        'SeSAMe': QC_TITLE,
        'ggplot2': 'Visualization',
        'ggrepel': 'Non-overlapping Labels',
        'sesameData': 'SeSAMe Reference Data',
        'limma': 'DMP Analysis',
        'DMRcate': 'DMR Analysis',
        'IlluminaHumanMethylationEPICv2anno.20a1.hg38': 'CpG Annotation',
        'methylclock': 'Horvath Age Estimation',
        'methylclockData': 'Horvath Coefficients & Reference Methylation Values',
    }
    versions = list(rows_from_csv(path, ('Software', 'Version')))

    entries = []
    for row in versions:
        name = row['Software']
        role = software_functions.get(name)
        label = f"{name} ({role})" if role else name
        entries.append(
            "<tr>"
            f"<td><strong>{esc(label)}</strong></td>"
            f"<td>{esc(normalize_version(row['Version']))}</td>"
            "</tr>"
        )
    body = "".join(entries)
    table = (
        '<div class="table-wrap software-versions-wrap"><table class="software-versions-table">'
        '<thead><tr><th>Software</th><th>Version</th></tr></thead>'
        f'<tbody>{body}</tbody></table></div>'
    )
    return collapsible_html('Software Versions', table)

def analysis_summary_section(versions_path):
    rows = [
        [
            (QC_TITLE, 'word'),
            ('SeSAMe preprocesses IDAT data and performs quality control using '
             'array reference data from sesameData. It also consolidates probe-level '
             'beta values into CpG-level beta values. ggplot2 generates the '
             f'{SAMPLE_QC_TITLE}, {BETA_DISTRIBUTION_TITLE}, and {PCA_TITLE} plots, '
             f'while ggrepel improves sample label placement in the {PCA_TITLE} plot.', 'word'),
        ],
        [
            ('Differential Methylation', 'word'),
            ('limma identifies differentially methylated positions (DMPs) from '
             'CpG-level beta values using statistical modeling. DMRcate identifies '
             'differentially methylated regions (DMRs) from probe-level beta values '
             'with EPIC v2 remapping and probe consolidation. '
             'IlluminaHumanMethylationEPICv2anno.20a1.hg38 provides genomic '
             'coordinates and gene annotations (hg38) for DMPs. '
             'ggplot2 generates the DMP volcano plot.', 'word'),
        ],
        [
            ('Epigenetic Age', 'word'),
            ('methylclock estimates DNA methylation age from CpG-level beta values '
             'using the Horvath Pan-Tissue Clock, supported by Horvath coefficients '
             'and reference methylation values from methylclockData. '
             f'ggplot2 generates the {CLOCK_TITLE} and {AGE_DIFF_TITLE} plots, '
             f'while ggrepel improves sample label placement in the {CLOCK_TITLE} plot.', 'word'),
        ],
    ]
    analysis_table = html_table(('Analysis', 'Description'), rows).replace(
        '<table>', '<table class="analysis-summary-table">', 1
    )
    analysis_table = analysis_table.replace(
        'IlluminaHumanMethylationEPICv2anno.20a1.hg38',
        'IlluminaHuman<wbr>Methylation<wbr>EPICv2<wbr>anno.20a1.hg38',
        1
    )
    return section('Analysis Summary',
                   collapsible_html('Analysis Methods', analysis_table)
                   + software_versions_section(versions_path))

def build_report(folder, slogan, logo, metadata, dmp, dmr, reference_group, comparison_group):
    brand = (
        f'<img class="partner-logo" src="{image_data_url(slogan)}" alt="Your Omics Partner">'
        f'<img class="brand-logo" src="{image_data_url(logo)}" alt="Pandu Biosains">'
    )
    parts = [
        metadata_section(metadata),
        qc_section(folder),
        dmp_section(folder, dmp, reference_group, comparison_group),
        dmr_section(dmr, comparison_group),
        clock_section(folder),
        age_diff_section(folder),
        analysis_summary_section(folder / FILES['versions']),
    ]
    return (
        '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>DNA Methylation Analysis Report</title>'
        f'<style>{CSS}</style></head>'
        '<body><main class="sheet">'
        f'<header class="brandline">{brand}</header>'
        '<section class="hero"><h1>DNA Methylation Analysis Report</h1>'
        '<div class="hero-cards">'
        f'<div class="sample"><span>Platform</span><strong>{esc(PLATFORM)}</strong></div>'
        f'<div class="sample"><span>Sample Type</span><strong>{esc(SAMPLE_TYPE)}</strong></div>'
        f'<div class="sample"><span>Chip ID</span><strong>{esc(CHIP_ID)}</strong></div>'
        f'<div class="sample"><span>Total Samples</span><strong>{len(metadata)}</strong></div>'
        '</div></section><div class="content">'
        + ''.join(parts)
        + '</div></main></body></html>\n'
    )

def generate(folder):
    slogan = folder / PANDU_SLOGAN
    logo = folder / PANDU_LOGO
    metadata = load_metadata(folder / METADATA_FILE)
    dmp = load_dmp(folder / FILES['dmp'])
    dmr = load_dmr(folder / FILES['dmr'])
    reference_group, comparison_group = load_comparison_info(folder / FILES['comparison_info'])
    destination = folder / REPORT_NAME
    destination.write_text(build_report(folder, slogan, logo, metadata, dmp, dmr,
                                        reference_group, comparison_group),
                           encoding='utf-8')
    print(f'Report generated: {destination}')

if __name__ == '__main__':
    generate(Path(__file__).resolve().parent)