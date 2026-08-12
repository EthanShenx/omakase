# Consequence categories and subtypes

Consequence categories and subtypes

## Usage

``` r
consequence_levels()
```

## Value

A named list with `categories` and `subtypes`, in the order they are
drawn.

## Examples

``` r
consequence_levels()
#> $categories
#> [1] "5'UTR change"      "Promoter swap"     "N-terminal/CDS"   
#> [4] "Unclassified"      "Distal (excluded)"
#> 
#> $subtypes
#>  [1] "longer"                   "shorter"                 
#>  [3] "equal"                    "alt first exon"          
#>  [5] "alt N-terminus"           "N-term extension"        
#>  [7] "N-term truncation"        "ORF loss"                
#>  [9] "ORF gain"                 "no ORF in either"        
#> [11] "non-overlapping isoforms"
#> 
```
