# Example data

Everything in this directory is public. No unpublished data is distributed with
omakase.

## `bams/` — RNA-seq alignments

Six BAM files (with indexes) covering `chr10:27,035,000-27,050,000`, the locus
of the human gene *ABI1*. Two libraries from each of three cell types:

| Sample | File | Cell type |
|---|---|---|
| ENCLB024ZZZ | `ENCFF088HTJ.bam` | Endothelial |
| ENCLB025ZZZ | `ENCFF450AIU.bam` | Endothelial |
| ENCLB271TJH | `ENCFF841RJE.bam` | Epithelial |
| ENCLB459IUG | `ENCFF756PUW.bam` | Epithelial |
| ENCLB008ZZZ | `ENCFF871MPV.bam` | Mesenchymal |
| ENCLB009ZZZ | `ENCFF603IAA.bam` | Mesenchymal |

The reads are from the **ENCODE Project**, whose data are released without
restriction for research use. The region-restricted subsets were prepared by the
[ggsashimi](https://github.com/guigolab/ggsashimi) authors and are redistributed
here under that project's MIT licence. Using the same files means an omakase
figure can be compared directly against a ggsashimi one.

* ENCODE Project Consortium (2012). An integrated encyclopedia of DNA elements
  in the human genome. *Nature* **489**, 57–74.
* Garrido-Martín D., Palumbo E., Guigó R., Breschi A. (2018). ggsashimi:
  Sashimi plot revised for browser- and annotation-independent splicing
  visualization. *PLOS Computational Biology* **14**(8), e1006360.

## `annotation.gtf`

GENCODE annotation restricted to the same window. Redistributed from the
ggsashimi examples. GENCODE data are released under the terms described at
<https://www.gencodegenes.org/>.

## `samples.tsv`

A sample manifest in the format [`read_manifest()`][omakase::read_manifest]
reads: sample identifier, path, then any number of metadata columns. The layout
is deliberately the one ggsashimi accepts, so an existing manifest works
unchanged.

## `palette.txt`

A one-colour-per-line palette file, again in ggsashimi's format, read by
[`read_palette()`][omakase::read_palette].

## `demo_consequence.tsv`

**Synthetic.** 180 invented start-site switches, each a made-up pair of
peptides, 5′ UTR lengths, first exons and transcript spans, put through
[`classify_consequence()`][omakase::classify_consequence]. The table is
therefore the classifier's genuine output, and the numbers describe nothing
real. Regenerate it with `data-raw/make_demo_consequence.R`.
