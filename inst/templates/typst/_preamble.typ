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
#let cpsr-teal-dark = rgb("#00504c") // table headers, section headings on cover
#let cpsr-rule    = rgb("#cbd5d4")  // hairlines
#let cpsr-stripe  = rgb("#f3f7f7")  // zebra rows
#let cpsr-muted   = rgb("#566262")  // secondary text (AA contrast on white)

// Page margins - must match the `margin` block in cpsr_report_pdf.qmd; used
// to let the cover band run full-bleed
#let cpsr-mx = 1.5cm
#let cpsr-bar = 0.32cm  // height of the full-bleed top/bottom bars

// Logo (copied together with the templates to the render directory)
#let cpsr-logo(height) = image("typst/cpsr_logo.png", height: height)

// Report-level values set by the cover, read by header/footer on later pages
#let cpsr-sample = state("cpsr-sample", "")
#let cpsr-ver    = state("cpsr-ver", "")
#let cpsr-today  = datetime.today().display("[day] [month repr:long] [year]")

// Pathogenicity scale (mirrors color_palette$pathogenicity in R)
#let cpsr-col-benign  = rgb("#077009")
#let cpsr-col-lbenign = rgb("#6FB572")
#let cpsr-col-vus     = rgb("#2c313c")
#let cpsr-col-lpath   = rgb("#9C3948")
#let cpsr-col-path    = rgb("#9E0142")

// ---------------------------------------------------------------------------
// Page chrome - full-bleed teal bands frame every page. Top band (page 2+)
// names the section(s) on the page and carries sample ID + date; bottom band
// closes the page; footer (white) holds logo, disclaimer/links and page count.
// Margins/paper come from Quarto YAML
// ---------------------------------------------------------------------------
#let cpsr-band-top = 0.95cm   // top band, pages 2+
#let cpsr-band-bot = 0.55cm   // bottom band, all pages

// Section label for the current page: the level-1 heading(s) on the page; if
// the first heading starts below the top part of the page, the section
// continued from the previous page is listed first. Call inside `context`.
#let cpsr-section-label() = {
  let pg = here().page()
  let hs = query(heading.where(level: 1))
  let on = hs.filter(h => h.location().page() == pg)
  let before = hs.filter(h => h.location().page() < pg)
  let names = ()
  if before.len() > 0 and (
      on.len() == 0 or on.first().location().position().y > 5cm) {
    names.push(before.last().body)
  }
  for h in on { names.push(h.body) }
  if names.len() == 0 { none } else { names.join([ · ]) }
}

#set page(
  background: context {
    let pg = here().page()
    if pg == 1 {
      // cover: thin bar, the cover has its own large logo and band
      place(top, rect(width: 100%, height: cpsr-bar, fill: cpsr-teal, stroke: none))
    } else {
      place(top, rect(
        width: 100%, height: cpsr-band-top, fill: cpsr-teal, stroke: none,
        inset: (x: cpsr-mx),
      )[
        #set text(fill: white, size: 8.5pt)
        #align(horizon, grid(
          columns: (1fr, auto),
          align: (left + horizon, right + horizon),
          text(weight: "bold", size: 9pt, tracking: 0.1em)[#upper(cpsr-section-label())],
          [#text(weight: "bold")[#cpsr-sample.final()] · #cpsr-today],
        ))
      ])
    }
    place(bottom, rect(
      width: 100%, height: cpsr-band-bot, fill: cpsr-teal, stroke: none))
  },
  header: none,
  footer: context {
    line(length: 100%, stroke: 0.5pt + cpsr-rule)
    v(3pt)
    set text(size: 7.5pt, fill: cpsr-muted)
    grid(
      columns: (auto, 1fr, auto),
      column-gutter: 8pt,
      align: (left + horizon, left + horizon, right + horizon),
      // cover page already carries the large logo
      if counter(page).get().first() == 1 { none } else { cpsr-logo(0.8cm) },
      [For research use only · cpsr v#cpsr-ver.final() ·
        #link("https://sigven.github.io/cpsr")[sigven.github.io/cpsr]],
      counter(page).display("1 of 1", both: true),
    )
  },
  // footer kept clear of the bottom band
  footer-descent: 8%,
)

// ---------------------------------------------------------------------------
// Typography
// ---------------------------------------------------------------------------
// Sans-serif font in line with the HTML report - Source Sans Pro is shipped
// with the rmarkdown R package, DejaVu Sans with the conda environment
// (fallback for symbols not covered by Source Sans Pro); both are made
// available to Typst through TYPST_FONT_PATHS (see cpsr::write_cpsr_output)
#set text(
  font: ("Source Sans Pro", "DejaVu Sans"),
  size: 10pt, fill: cpsr-dark)
#set par(leading: 0.65em)
// tables (wrapped in figures by Quarto) may break across pages
#show figure: set block(breakable: true)
// table cells: ragged-right text without hyphenation (narrow columns)
#show table.cell: set par(justify: false)
#show table.cell: set text(hyphenate: false)
// tables set slightly smaller than body text (10pt)
#show table: set text(size: 8.5pt)
// links in teal (clickable in the PDF)
#show link: set text(fill: cpsr-teal)

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

// Metadata label+value pair (used inside the cover page info grid)
#let cpsr-meta(label, value) = (
  text(weight: "bold", size: 8.5pt, fill: cpsr-muted)[#upper(label)],
  text(size: 10pt, weight: "bold", fill: cpsr-dark)[#value],
)

#let cpsr-cover-heading(body) = text(
  size: 12pt, weight: "bold", fill: cpsr-teal-dark)[#body]

// ---------------------------------------------------------------------------
// Cover page (sample metadata), followed by key findings and contents
// Called at the very top of cpsr_report_pdf.qmd via a raw Typst block emitted
// from R (cpsr:::pdf_cover_typst) so that all values are resolved at render
// time.
// ---------------------------------------------------------------------------
#let cpsr-cover(
  sample_id:        "—",
  genome_assembly:  "—",
  gene_panel:       "—",
  sex_age:          "—",
  cpsr_version:     "—",
) = {

  cpsr-sample.update(sample_id)
  cpsr-ver.update(cpsr_version)

  // --- Logo (white) and report date ---
  grid(
    columns: (auto, 1fr),
    align: (left + horizon, right + horizon),
    cpsr-logo(2.4cm),
    align(right)[
      #set text(size: 8.5pt, fill: cpsr-muted)
      Report date
      #linebreak()
      #text(size: 10pt, weight: "bold", fill: cpsr-dark)[#cpsr-today]
    ],
  )

  v(0.6cm)

  // --- Full-bleed band: sample name as the headline ---
  pad(x: -cpsr-mx, rect(
    fill: cpsr-teal, width: 100%, height: 5.2cm, stroke: none,
    inset: (x: cpsr-mx, bottom: 0.9cm, top: 0.9cm),
  )[
    #set text(fill: white)
    #align(bottom + left)[
      #text(size: 9pt, weight: "bold", tracking: 0.12em)[#upper[Cancer predisposition sequencing report]]
      #v(3pt)
      // sample ID on a single line: 26pt, scaled down to fit the band for
      // long IDs (max. 40 characters; min. 14pt)
      #layout(size => {
        let fs = 26pt
        let w = measure(text(size: fs, weight: "bold", sample_id)).width
        if w > size.width { fs = calc.max(14pt, fs * (size.width / w)) }
        text(size: fs, weight: "bold", sample_id)
      })
      #v(3pt)
      #text(size: 9pt)[Cancer Predisposition Sequencing Reporter · cpsr v#cpsr_version]
    ]
  ])

  v(0.8cm)

  // --- Sample metadata ---
  grid(
    columns: (auto, 1fr),
    column-gutter: 1.4em,
    row-gutter: 0.6em,
    ..cpsr-meta("Genome assembly", genome_assembly),
    ..cpsr-meta("Sex / age", sex_age),
    ..cpsr-meta("Virtual gene panel", gene_panel),
  )

  v(0.5cm)
  line(length: 100%, stroke: 0.5pt + cpsr-rule)
  v(0.4cm)
}

// Key findings on the cover: a framed panel with a colored title band, a
// tinted lead section (most important findings) and a white list of
// supporting findings. Content from cpsr::report_highlights(); the R code
// emits #cpsr-highlights[ #cpsr-hl-lead[..markdown..] #cpsr-hl-support[..markdown..] ]
// Pathogenicity label (e.g. Pathogenic) as a colored pill - white bold text
// on the pathogenicity color, as in the Clinical significance column
#let cpsr-hl-pill(label, fill) = box(
  fill: rgb(fill), radius: 3pt, inset: (x: 6pt, y: 2.5pt), baseline: 1.5pt,
  text(fill: white, weight: "bold", size: 0.92em)[#label],
)
#let cpsr-hl-marker = (
  text(fill: cpsr-teal, size: 8pt, weight: "bold")[#sym.bullet],
  text(fill: cpsr-teal, weight: "bold")[–],
)
#let cpsr-highlights(body) = block(
  width: 100%, breakable: false, clip: true, radius: 5pt,
  stroke: 1pt + cpsr-teal, fill: white,
)[
  #set block(spacing: 0pt)
  #block(width: 100%, fill: cpsr-teal, inset: (x: 14pt, y: 8pt))[
    #set text(fill: white)
    #grid(
      columns: (auto, 1fr), column-gutter: 8pt, align: horizon,
      polygon.regular(vertices: 6, size: 11pt, fill: white),
      text(size: 10.5pt, weight: "bold", tracking: 0.14em)[#upper[Key findings]],
    )
  ]
  #body
]
#let cpsr-hl-lead(body) = block(
  width: 100%, inset: (x: 14pt, y: 12pt), fill: cpsr-teal.lighten(90%),
)[
  #set text(size: 10.5pt)
  #set par(justify: false, leading: 0.62em)
  #set block(spacing: 0.8em)
  #set list(spacing: 0.65em, marker: cpsr-hl-marker)
  #show strong: set text(fill: cpsr-teal-dark)
  #body
]
#let cpsr-hl-support(body) = block(
  width: 100%, inset: (x: 14pt, y: 11pt),
  stroke: (top: 0.75pt + cpsr-rule),
)[
  #set text(size: 9.5pt)
  #set par(justify: false, leading: 0.6em)
  #set block(spacing: 0.8em)
  #set list(spacing: 0.7em, marker: cpsr-hl-marker)
  #show strong: set text(fill: cpsr-teal-dark)
  #body
]

// Contents (level-1 headings, page numbers resolved by Typst) and end of
// the cover page
#let cpsr-cover-end() = {
  v(0.6cm)
  cpsr-cover-heading[Contents]
  v(0.2cm)
  {
    set text(size: 10pt)
    outline(title: none, depth: 1, indent: 0pt)
  }
  pagebreak()
}
