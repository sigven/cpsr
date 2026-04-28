// =============================================================================
// CPSR Typst preamble
// Defines brand colors, page chrome, show rules, and reusable components.
// Included via `include-before-body` in cpsr_report_pdf.qmd.
// =============================================================================

// ---------------------------------------------------------------------------
// Brand colors
// ---------------------------------------------------------------------------
#let cpsr-teal    = rgb("#007a74")  // primary brand
#let cpsr-dark    = rgb("#2c313c")  // headings, body chrome
#let cpsr-neutral = rgb("#8B8989")  // empty / no-finding state

// Pathogenicity scale (mirrors color_palette$pathogenicity in R)
#let cpsr-col-benign  = rgb("#077009")
#let cpsr-col-lbenign = rgb("#6FB572")
#let cpsr-col-vus     = rgb("#2c313c")
#let cpsr-col-lpath   = rgb("#9C3948")
#let cpsr-col-path    = rgb("#9E0142")

// ---------------------------------------------------------------------------
// Page chrome — header and footer only; margins/paper come from Quarto YAML
// ---------------------------------------------------------------------------
#set page(
  header: context {
    set text(size: 7.5pt, fill: cpsr-teal)
    grid(
      columns: (1fr, auto),
      align: (left, right),
      [*CPSR* — Cancer Predisposition Sequencing Reporter],
      text(fill: luma(140))[#datetime.today().display("[day] [month repr:long] [year]")],
    )
    v(-6pt)
    line(length: 100%, stroke: 0.4pt + cpsr-teal)
  },
  footer: context {
    line(length: 100%, stroke: 0.4pt + cpsr-teal)
    v(-6pt)
    set text(size: 7.5pt, fill: luma(140))
    grid(
      columns: (1fr, auto),
      align: (left, right),
      [For research use only],
      counter(page).display("1 of 1", both: true),
    )
  },
  // Suppress header/footer on the cover page
  header-ascent: 30%,
  footer-descent: 30%,
)

// First page: no header (cover occupies the full top)
#set page(header: context {
  if counter(page).get().first() == 1 { none }
  else {
    set text(size: 7.5pt, fill: cpsr-teal)
    grid(
      columns: (1fr, auto),
      align: (left, right),
      [*CPSR* — Cancer Predisposition Sequencing Reporter],
      text(fill: luma(140))[#datetime.today().display("[day] [month repr:long] [year]")],
    )
    v(-6pt)
    line(length: 100%, stroke: 0.4pt + cpsr-teal)
  }
})

// ---------------------------------------------------------------------------
// Typography
// ---------------------------------------------------------------------------
#set text(size: 10pt, fill: cpsr-dark, font: "Source Sans Pro", fallback: true)
#set par(leading: 0.65em)

// ---------------------------------------------------------------------------
// Heading show rules
// ---------------------------------------------------------------------------
#show heading.where(level: 1): it => block(above: 1.6em, below: 0.6em)[
  #set text(size: 13pt, weight: "bold", fill: cpsr-teal)
  #it.body
  #v(-4pt)
  #line(length: 100%, stroke: 1.2pt + cpsr-teal)
]

#show heading.where(level: 2): it => block(above: 1.2em, below: 0.4em)[
  #set text(size: 11pt, weight: "bold", fill: cpsr-dark)
  #it.body
]

#show heading.where(level: 3): it => block(above: 1em, below: 0.3em)[
  #set text(size: 10pt, weight: "bold", fill: cpsr-dark)
  #it.body
]

// ---------------------------------------------------------------------------
// Reusable components
// ---------------------------------------------------------------------------

// Pathogenicity / evidence badge (inline pill)
#let cpsr-badge(label, color: cpsr-neutral) = box(
  fill: color,
  inset: (x: 5pt, y: 2pt),
  radius: 3pt,
  text(fill: white, size: 8pt, weight: "bold", label),
)

// KPI summary box used on the cover page
#let cpsr-kpi-box(title, value, color: cpsr-neutral) = rect(
  fill: color,
  radius: 4pt,
  inset: (x: 10pt, y: 9pt),
  width: 100%,
)[
  #set text(fill: white)
  #text(size: 7.5pt)[#title]
  #linebreak()
  #text(size: 15pt, weight: "bold")[#value]
]

// Metadata label+value pair (used inside the cover page info grid)
#let cpsr-meta(label, value) = (
  text(weight: "bold", size: 8.5pt, fill: cpsr-dark)[#label],
  text(size: 8.5pt, fill: cpsr-dark)[#value],
)

// ---------------------------------------------------------------------------
// Cover + executive summary page
// Call this function at the very top of cpsr_report_pdf.qmd via a raw Typst
// block emitted from R so that all values are resolved at render time.
// ---------------------------------------------------------------------------
#let cpsr-cover(
  sample_id:        "—",
  genome_assembly:  "—",
  gene_panel:       "—",
  cpsr_version:     "—",
  // KPI values (strings, pre-formatted by R)
  genes_path:       "None",
  genes_bm_pgx:     "None",
  genes_sf:         "Not determined",
  n_path:           "N = 0",
  n_vus:            "N = 0",
  n_benign:         "N = 0",
  // KPI box colors (hex strings, resolved by R from color_palette)
  col_genes_path:   "#8B8989",
  col_bm_pgx:       "#8B8989",
  col_sf:           "#8B8989",
  col_path:         "#8B8989",
  col_vus:          "#8B8989",
  col_benign:       "#8B8989",
) = {

  // --- Banner ---
  rect(fill: cpsr-teal, width: 100%, inset: (x: 18pt, y: 14pt))[
    #set text(fill: white)
    #grid(
      columns: (1fr, auto),
      align: (left + horizon, right + horizon),
      column-gutter: 1em,
      [
        #text(size: 17pt, weight: "bold")[Cancer Predisposition Sequencing Reporter]
        #linebreak()
        #text(size: 8.5pt)[cpsr v#cpsr_version]
      ],
      text(size: 8.5pt)[#datetime.today().display("[day] [month repr:long] [year]")],
    )
  ]

  v(1.4em)

  // --- Sample metadata ---
  grid(
    columns: (auto, 1fr, auto, 1fr),
    column-gutter: 1.2em,
    row-gutter: 0.55em,
    ..cpsr-meta("Sample ID",         sample_id),
    ..cpsr-meta("Genome assembly",   genome_assembly),
    ..cpsr-meta("Gene panel",        gene_panel),
    grid.cell(colspan: 3)[],   // blank right cells on second row
  )

  v(1.8em)

  // --- KPI heading ---
  text(size: 9pt, weight: "bold", fill: cpsr-dark)[Summary of findings]
  v(0.6em)

  // --- KPI grid: 3 columns × 2 rows ---
  grid(
    columns: (1fr, 1fr, 1fr),
    column-gutter: 7pt,
    row-gutter: 7pt,
    cpsr-kpi-box("Genes with pathogenic variants",      genes_path,  color: rgb(col_genes_path)),
    cpsr-kpi-box("Genes with biomarkers / PGx",        genes_bm_pgx, color: rgb(col_bm_pgx)),
    cpsr-kpi-box("Genes with secondary findings",      genes_sf,    color: rgb(col_sf)),
    cpsr-kpi-box("Pathogenic / Likely pathogenic",     n_path,      color: rgb(col_path)),
    cpsr-kpi-box("Variants of uncertain significance", n_vus,       color: rgb(col_vus)),
    cpsr-kpi-box("Benign / Likely benign",             n_benign,    color: rgb(col_benign)),
  )

  v(1em)
  line(length: 100%, stroke: 0.4pt + luma(200))
  pagebreak()
}
