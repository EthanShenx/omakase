# omakase: Elegant Sashimi Plots and Transcript-Consequence Figures

A chef's-choice toolkit for RNA-seq splicing and
transcription-start-site figures. 'omakase' reads coverage and junctions
straight from BAM/SAM files, rMATS event tables, STAR junction files,
5'-tag BEDs or tidy data frames, and renders them as publication-ready
sashimi plots built on 'ggplot2' and 'patchwork'. Every panel is
returned as a plot object, so figures compose and can be modified after
the fact rather than being written straight to disk. Beyond plotting,
the package classifies the consequence of an alternative transcription
start site on the transcript it produces - 5' UTR change, promoter swap
or N-terminal/CDS impact - from open reading frames called on the
spliced mRNA, and visualises that composition as a radial donut. Intron
compression, coverage normalisation, junction arc geometry and PSI are
all exposed as documented, swappable pieces.

## See also

Useful links:

- <https://github.com/EthanShenx/omakase>

- <https://EthanShenx.github.io/omakase/>

- Report bugs at <https://github.com/EthanShenx/omakase/issues>

## Author

**Maintainer**: Ethan Shen <ethanshen111@gmail.com>
