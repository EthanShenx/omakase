# Read a sample manifest

Reads the tab-separated manifest used to describe a set of alignment
files. The format matches the one `ggsashimi` accepts, so an existing
manifest works unchanged: column 1 is a sample identifier, column 2 the
path to the file, and any further columns are metadata that can be named
in `group_col` or `color_col`.

## Usage

``` r
read_manifest(path, group_col = NULL, label_col = NULL, base_dir = NULL)
```

## Arguments

- path:

  Path to the manifest, or a character vector of file paths, or a single
  alignment file.

- group_col:

  Column (name or 1-based index) that assigns samples to groups.
  Defaults to column 3 when present, otherwise every sample is its own
  group.

- label_col:

  Column carrying the label to print on each panel. Defaults to the
  sample identifier.

- base_dir:

  Directory that relative paths in the manifest are resolved against;
  defaults to the manifest's own directory.

## Value

A data frame with columns `sample`, `path`, `group`, `label`, plus any
extra manifest columns.

## Details

A header row is optional. When absent, columns are named `sample`,
`path`, `V3`, `V4` and so on.

## Examples

``` r
m <- tempfile(fileext = ".tsv")
writeLines(c("s1\ta.bam\tEndothelial", "s2\tb.bam\tEpithelial"), m)
read_manifest(m)
#>   sample                  path          V3       group label
#> 1     s1 /tmp/Rtmp8fD4O2/a.bam Endothelial Endothelial    s1
#> 2     s2 /tmp/Rtmp8fD4O2/b.bam  Epithelial  Epithelial    s2
```
