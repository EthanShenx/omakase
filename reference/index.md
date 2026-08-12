# Package index

## Reading data

Every reader returns the same object, so the input format and the figure
are independent choices.

- [`sashimi_from_bam()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_bam.md)
  : Read coverage and junctions from alignment files
- [`sashimi_from_rmats()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_rmats.md)
  : Build a sashimi data object from rMATS events and alignments
- [`sashimi_from_junctions()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_junctions.md)
  : Build a sashimi data object from junction files
- [`sashimi_from_tags()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_tags.md)
  : Build a sashimi data object from 5'-tag data
- [`sashimi_from_tables()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_tables.md)
  : Build a sashimi data object from tidy tables
- [`read_sashimi_dir()`](https://EthanShenx.github.io/omakase/reference/read_sashimi_dir.md)
  : Read a directory of sashimi tables
- [`as_sashimi_data()`](https://EthanShenx.github.io/omakase/reference/as_sashimi_data.md)
  : Coerce an object to sashimi data
- [`read_manifest()`](https://EthanShenx.github.io/omakase/reference/read_manifest.md)
  : Read a sample manifest
- [`read_annotation()`](https://EthanShenx.github.io/omakase/reference/read_annotation.md)
  : Read a GTF or GFF3 annotation
- [`read_junctions()`](https://EthanShenx.github.io/omakase/reference/read_junctions.md)
  : Read a splice junction file
- [`write_junctions()`](https://EthanShenx.github.io/omakase/reference/write_junctions.md)
  : Write junctions to a BED file
- [`read_bed()`](https://EthanShenx.github.io/omakase/reference/read_bed.md)
  [`read_bed12()`](https://EthanShenx.github.io/omakase/reference/read_bed.md)
  : Read a BED file
- [`read_tag_bed()`](https://EthanShenx.github.io/omakase/reference/read_tag_bed.md)
  : Read a 5'-tag BED file
- [`read_rmats()`](https://EthanShenx.github.io/omakase/reference/read_rmats.md)
  : Read an rMATS event table
- [`read_palette()`](https://EthanShenx.github.io/omakase/reference/read_palette.md)
  : Read a palette file

## The data contract

The tidy tables every reader produces and every plot consumes.

- [`sashimi_data()`](https://EthanShenx.github.io/omakase/reference/sashimi_data.md)
  : Construct a sashimi data object
- [`validate_sashimi_data()`](https://EthanShenx.github.io/omakase/reference/validate_sashimi_data.md)
  : Validate a sashimi data object
- [`loci()`](https://EthanShenx.github.io/omakase/reference/loci.md) :
  Locus identifiers held by a sashimi data object
- [`combine_sashimi()`](https://EthanShenx.github.io/omakase/reference/combine_sashimi.md)
  : Combine sashimi data objects
- [`` `[`( ``*`<sashimi_data>`*`)`](https://EthanShenx.github.io/omakase/reference/sub-.sashimi_data.md)
  : Subset a sashimi data object by locus
- [`write_sashimi_data()`](https://EthanShenx.github.io/omakase/reference/write_sashimi_data.md)
  : Write a sashimi data object to disk
- [`add_models()`](https://EthanShenx.github.io/omakase/reference/add_models.md)
  : Add transcript models to a sashimi data object

## Drawing sashimi plots

- [`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md)
  : Draw a sashimi figure
- [`plot_sashimi_all()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi_all.md)
  : Draw a sashimi figure for every locus
- [`sashimi_track()`](https://EthanShenx.github.io/omakase/reference/sashimi_track.md)
  : Draw a single coverage panel with its junction arcs
- [`sashimi_annotation()`](https://EthanShenx.github.io/omakase/reference/sashimi_annotation.md)
  : Draw the annotation panel for one locus
- [`save_sashimi()`](https://EthanShenx.github.io/omakase/reference/save_sashimi.md)
  : Save a sashimi figure
- [`presets()`](https://EthanShenx.github.io/omakase/reference/presets.md)
  : Plotting presets

## Drawing genome tracks

The browser-style view of a locus: stacked tracks over a shared
coordinate axis, built from composable pieces.

- [`plot_tracks()`](https://EthanShenx.github.io/omakase/reference/plot_tracks.md)
  : Draw a genome-track figure
- [`track_models()`](https://EthanShenx.github.io/omakase/reference/track_models.md)
  : Build a transcript model track
- [`track_coverage()`](https://EthanShenx.github.io/omakase/reference/track_coverage.md)
  : Build a coverage track
- [`track_features()`](https://EthanShenx.github.io/omakase/reference/track_features.md)
  : Build a feature track
- [`track_axis()`](https://EthanShenx.github.io/omakase/reference/track_axis.md)
  : Build a coordinate axis track
- [`track_spacer()`](https://EthanShenx.github.io/omakase/reference/track_spacer.md)
  : Build a spacer track
- [`save_tracks()`](https://EthanShenx.github.io/omakase/reference/save_tracks.md)
  : Save a genome-track figure

## Transforming

Normalisation, aggregation, PSI and intron compression.

- [`normalize_tracks()`](https://EthanShenx.github.io/omakase/reference/normalize_tracks.md)
  : Normalise coverage tracks
- [`normalize_methods()`](https://EthanShenx.github.io/omakase/reference/normalize_methods.md)
  : Coverage normalisation methods
- [`aggregate_tracks()`](https://EthanShenx.github.io/omakase/reference/aggregate_tracks.md)
  : Aggregate replicate tracks into one track per group
- [`compute_psi()`](https://EthanShenx.github.io/omakase/reference/compute_psi.md)
  : Compute PSI (percent spliced in) per locus and group
- [`delta_psi()`](https://EthanShenx.github.io/omakase/reference/delta_psi.md)
  : Difference in PSI between two groups
- [`intron_map()`](https://EthanShenx.github.io/omakase/reference/intron_map.md)
  : Build a genome-to-plot coordinate map that compresses introns
- [`compress_coords()`](https://EthanShenx.github.io/omakase/reference/compress_coords.md)
  [`expand_coords()`](https://EthanShenx.github.io/omakase/reference/compress_coords.md)
  : Apply or invert an intron compression map
- [`intron_trans()`](https://EthanShenx.github.io/omakase/reference/intron_trans.md)
  : A ggplot2 axis transform for compressed genomic coordinates
- [`shrink_methods()`](https://EthanShenx.github.io/omakase/reference/shrink_methods.md)
  : Intron compression methods

## Arc geometry

- [`arc_path()`](https://EthanShenx.github.io/omakase/reference/arc_path.md)
  : Build the path of a junction arc
- [`arc_shapes()`](https://EthanShenx.github.io/omakase/reference/arc_shapes.md)
  : Arc shapes available for junction curves
- [`arc_heights()`](https://EthanShenx.github.io/omakase/reference/arc_heights.md)
  : Compute arc apex heights from junction counts
- [`arc_height_rules()`](https://EthanShenx.github.io/omakase/reference/arc_height_rules.md)
  : Arc height rules
- [`arc_widths()`](https://EthanShenx.github.io/omakase/reference/arc_widths.md)
  : Compute arc line widths from junction counts

## Transcript consequences

What an alternative start site does to the transcript, from the open
reading frames.

- [`classify_consequence()`](https://EthanShenx.github.io/omakase/reference/classify_consequence.md)
  : Classify the consequence of an alternative start site
- [`consequence_summary()`](https://EthanShenx.github.io/omakase/reference/consequence_summary.md)
  : Summarise a classified switch table
- [`consequence_levels()`](https://EthanShenx.github.io/omakase/reference/consequence_levels.md)
  : Consequence categories and subtypes
- [`consequence_labels()`](https://EthanShenx.github.io/omakase/reference/consequence_labels.md)
  : Human-readable labels for consequence subtypes
- [`plot_consequence()`](https://EthanShenx.github.io/omakase/reference/plot_consequence.md)
  : Draw the consequence composition
- [`find_orf()`](https://EthanShenx.github.io/omakase/reference/find_orf.md)
  : Call the open reading frame of an mRNA
- [`orf_table()`](https://EthanShenx.github.io/omakase/reference/orf_table.md)
  : Call ORFs for a set of transcript models
- [`splice_mrna()`](https://EthanShenx.github.io/omakase/reference/splice_mrna.md)
  : Splice exons into an mRNA sequence
- [`count_uatg()`](https://EthanShenx.github.io/omakase/reference/count_uatg.md)
  : Count upstream AUGs in a 5' UTR

## Appearance

- [`theme_omakase()`](https://EthanShenx.github.io/omakase/reference/theme_omakase.md)
  : The omakase theme
- [`theme_omakase_axes()`](https://EthanShenx.github.io/omakase/reference/theme_omakase_axes.md)
  : A themed variant that keeps the axes
- [`omakase_palette()`](https://EthanShenx.github.io/omakase/reference/omakase_palette.md)
  : Colour palettes
- [`omakase_palettes()`](https://EthanShenx.github.io/omakase/reference/omakase_palettes.md)
  : Palette names available
- [`scale_colour_omakase()`](https://EthanShenx.github.io/omakase/reference/scale_colour_omakase.md)
  [`scale_color_omakase()`](https://EthanShenx.github.io/omakase/reference/scale_colour_omakase.md)
  [`scale_fill_omakase()`](https://EthanShenx.github.io/omakase/reference/scale_colour_omakase.md)
  : Discrete colour and fill scales using the omakase palettes

## Utilities

- [`parse_region()`](https://EthanShenx.github.io/omakase/reference/parse_region.md)
  : Parse a genomic region string
- [`format_activity()`](https://EthanShenx.github.io/omakase/reference/format_activity.md)
  : Format junction or activity values with adaptive precision
- [`format_count()`](https://EthanShenx.github.io/omakase/reference/format_count.md)
  : Format an integer count
- [`format_coord()`](https://EthanShenx.github.io/omakase/reference/format_coord.md)
  : Format a genomic coordinate
- [`format_percent()`](https://EthanShenx.github.io/omakase/reference/format_percent.md)
  : Format a percentage
- [`event_types()`](https://EthanShenx.github.io/omakase/reference/event_types.md)
  : Alternative splicing event types
- [`strand_modes()`](https://EthanShenx.github.io/omakase/reference/strand_modes.md)
  : Strand specificity settings
