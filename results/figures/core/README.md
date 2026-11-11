# Fixed-study core figures

Rendered from the reference fixed-study core reports and primary network geometry only.
No pcnetmeta refit, ranking, or absolute-risk output is used.

- `primary_forest.{png,pdf}` — all five modeled nodes versus wild-type.
- `primary_network.{png,pdf}` — all 10 observed study-co-occurrence edges.
- `hapb3_sensitivity_forest.{png,pdf}` — 12 prespecified HapB3 scenarios.
- `core_figure_data.csv` — machine-readable plotted data, including derived node coordinates and direct edges.
- `input_checksums.csv` — SHA-256 checksums of the four committed input artifacts.

Renderer: `R/render_core_figures.R`; execution environment: pinned Docker image `runtime-recorded`.
