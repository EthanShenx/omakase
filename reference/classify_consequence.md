# Classify the consequence of an alternative start site

Compares a reference ("main") isoform with an alternative one and
reports what the switch between them does to the transcript.

## Usage

``` r
classify_consequence(x, min_overlap = 0)
```

## Arguments

- x:

  A data frame with one row per switch. Required columns:
  `main_protein`, `alt_protein`, `main_utr5_len`, `alt_utr5_len`,
  `main_first_exon_start`, `main_first_exon_end`,
  `alt_first_exon_start`, `alt_first_exon_end`, `main_tss`, `alt_tss`,
  `strand`. Optional: `gene_id`, `gene_name`, `main_tx_start`,
  `main_tx_end`, `alt_tx_start`, `alt_tx_end` (needed for the overlap
  guard), `main_n_uATG`, `alt_n_uATG`, `main_aa_len`, `alt_aa_len`.

- min_overlap:

  Pairs whose reciprocal transcript-span overlap is at or below this are
  excluded as distal mispairings.

## Value

The input with added columns: `category`, `subtype`, `body_overlap`,
`d_utr5_bp`, `aa_delta`, `n_uATG_gained` and `uorf_gained`.

## Details

Writing \\P_m\\ and \\P_a\\ for the peptides the two isoforms encode,
and \\\ell_m\\, \\\ell_a\\ for their 5' UTR lengths, the decision tree
is:

1.  \\P_m \ne \emptyset,\\ P_a = \emptyset\\ gives **ORF loss**; \\P_m =
    \emptyset,\\ P_a \ne \emptyset\\ gives **ORF gain**; both empty is
    unclassified. Calling a pair with no ORF in either isoform an "ORF
    loss" would blame the start site for a failure that was already
    there.

2.  \\P_a = P_m\\, so the protein is untouched and the change is purely
    5': if the two first exons are disjoint **and** the alternative site
    lies downstream, it is a **promoter swap**; otherwise a **5' UTR
    change**, `longer`, `shorter` or `equal` by the sign of \\\ell_a -
    \ell_m\\.

3.  \\P_a \ne P_m\\: if \\P_m\\ ends with \\P_a\\ the alternative
    protein is a truncation; if \\P_a\\ ends with \\P_m\\ it is an
    extension; otherwise the N-terminus is simply different.

A guard runs first. The reciprocal overlap of the two transcript spans,
\$\$o = \frac{\max(0, \min(e_1, e_2) - \max(s_1, s_2))}{\min(e_1 -
s_1,\\ e_2 - s_2)},\$\$ is zero when the two "isoforms" do not overlap
at all. That is not an alternative start site of one transcription unit,
it is two separate transcripts sharing a gene label, and such pairs are
marked `"Distal (excluded)"` rather than classified.

## Examples

``` r
sw <- data.frame(
  gene_name = c("A", "B", "C"),
  strand = "+",
  main_protein = c("MKVLA", "MKVLA", "MKVLA"),
  alt_protein  = c("MKVLA", "VLA",   ""),
  main_utr5_len = c(100, 100, 100), alt_utr5_len = c(200, 100, 100),
  main_first_exon_start = 1000, main_first_exon_end = 1200,
  alt_first_exon_start = c(1050, 1050, 1050),
  alt_first_exon_end = c(1250, 1250, 1250),
  main_tss = 1000, alt_tss = 1050,
  main_tx_start = 1000, main_tx_end = 5000,
  alt_tx_start = 1050, alt_tx_end = 5000
)
classify_consequence(sw)[, c("gene_name", "category", "subtype")]
#>   gene_name       category           subtype
#> 1         A   5'UTR change            longer
#> 2         B N-terminal/CDS N-term truncation
#> 3         C N-terminal/CDS          ORF loss
```
