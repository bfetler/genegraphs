# Gene Graphs

## Synteny Blocks in Genetics

Genomic sequencing and mapping have enabled comparison of the
general structures of genomes of many different species. The 
general finding is that organisms of relatively recent divergence 
show similar blocks of genes in the same relative positions in
the genome. This situation is called synteny, translated roughly 
as possessing common chromosome sequences. For example, many of 
the genes of humans are syntenic with those of other mammals — 
not only apes but also cows, mice, and so on. Study of synteny 
can show how the genome is cut and pasted in the course of 
evolution.

--from https://en.wikipedia.org/wiki/Synteny

## Exercise

The exercise here is to make an algorithm to create
a breakpoint graph comparing two genes.  I used a DiGraph 
with a DFS variation to identify Connected Components.
