# omakase 0.1.0

First release.

## Reading

* `sashimi_from_bam()` reads coverage and spliced junctions from BAM/CRAM files
  in a single pass, with strand-specific protocols, mapping-quality and
  junction-count filters, and per-sample or per-group tracks.
* `sashimi_from_rmats()` turns an rMATS `SE`/`A5SS`/`A3SS`/`MXE`/`RI` event
  table into one figure per event, with the event's two isoforms drawn as
  transcript models.
* `sashimi_from_junctions()` reads STAR `SJ.out.tab` and regtools junction BED
  files, for projects that kept junctions but not alignments.
* `sashimi_from_tags()` reads 5′-end tag data (CAGE, STRT, CamoTSS) and turns
  single-base start sites into readable tracks.
* `sashimi_from_tables()` and `read_sashimi_dir()` build the object from tidy
  tables, so counts from anywhere can be drawn.
* Sample manifests and palette files use the same formats `ggsashimi` accepts.

## The data contract

* `sashimi_data()` is one documented object — six tidy tables — that every
  reader produces and every plotting function consumes. Writing a new reader
  means producing those tables and nothing else.

## Drawing

* `plot_sashimi()` returns a `patchwork` of `ggplot` panels rather than writing
  a file, so figures compose and can be modified after the fact.
* Six junction arc geometries (`sine`, `parabola`, `bezier`, `xspline`,
  `elbow`, `arch`), with height rules that stagger a handful of arcs and
  switch to span-scaling when there are many.
* Intron compression through an explicit, exactly invertible piecewise-linear
  map, exposed as `intron_map()` and as a `scales` transform, so axis labels
  report true genomic coordinates.
* Coverage normalisation — CPM, RPM, RPKM, DESeq2 median-of-ratios size
  factors, or your own — which no comparable tool offers.

## Transcript consequences

* `classify_consequence()` says what an alternative start site does to the
  transcript: 5′ UTR change, promoter swap, or N-terminal/CDS impact, from the
  open reading frames, with a guard against distal mispairings.
* `find_orf()`, `orf_table()` and `splice_mrna()` call reading frames on
  spliced mRNA.
* `plot_consequence()` draws the composition as a radial donut, bar, lollipop
  or stacked column.

## Elsewhere

* A command-line interface at `inst/scripts/omakase` for shell pipelines.
* 230 unit tests covering the arc geometry, the compression map's round trip,
  normalisation, PSI, and the consequence decision tree.
