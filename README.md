# GCPSR (Genealogical Concordance Phylogenetic Species Recognition)

This package contains perl scripts for the implementation of the
genealogical concordance phylogenetic species recognition method.

## Prerequisites

You need to have perl 5 on your computer (on Mac and Linux this is by
default installed).

The scripts in this repository use the `Bio::Tree::Tree` package.

    # To install using cpan
    cpan Bio::Tree::Tree

## Perl scripts

To print the usage information for the scripts, run

    perl concordance_non-discordance.pl -h

and

    perl exhaustive_subdivision.pl -h

To run both steps with a minimal support of two

    ./concordance_non-discordance.pl -count=2 *.nwk | ./exhaustive_subdivision.pl - >gcpsr.nwk

### Step 1: Concordance and non-discordance analysis

The script `concordance_non-discordance.pl` takes **newick** or
**nexus** tree files (_one tree per file!_) as input and outputs a
tree by keeping clades that are highly supported by single gene trees
and are not conflicting with other clades with the same level of
support.

> **Important:** the names of the individuals/strains in the tree
> files have to agree across files and all the strains have to be
> present in all of the tree files.

The script has two parameters: the minimum support value
(`-min=<int>`) to keep a clade as a potential concordant clade and the
minimum number (`-count=<int>`) of single gene trees containing the
given clade with sufficient support. (Default values are `-min=95` and
`-count=1`: a clade has at least **95** as support value and is found
in at least **1** of the trees with that support.)

### Step 2. Exhaustive subdivision

The script `exhaustive_subdivision.pl` takes a single **newick** or
**nexus** tree as input and does the exhaustive subdivision analysis
and classifies all the individuals/strains into well (defined by
`-count=<int>`) supported phylogenetic species.

## GCPSR as implemented in this repository

### Step 1. (Concordance and non-discordance analysis)
1. Individual single locus gene trees are searched for clades with
   high support (_support_ >= `<int>`, `<int>` comes from the
   `-min=<int>` argument).
2. Highly supported clades are screened for concordance. Only those
   highly supported clades are considered as concordant that are
   present in at least `<int>` single locus gene trees (`<int>` comes
   from the `-count=<int>` argument).
3. Highly supported clades that are concordant are screened for
   discordance. All clades that are in conflict with each other are
   removed. This way all clades that are kept are non-discordant.
   
This produces a tree that has clades that are both concordant and
non-discordant across the single locus phylogenies. The support value
for each of the clades is the number of single locus phylogenies
containing the given clade with at least the minimum support values.

### Step 2. (Exhaustive subdivision)
1. The tree produced by step one is screened for clades with high the
   high support (_support_ >= `<int>`, `<int>` comes from
   `-count=<int>` argument). This is the cutoff value for recognizing
   a concordant clade as a phylogenetic species.
2. Each individual/strain is placed into a potential phylogenetic
   species that has the fewest members (smallest clade) and includes
   the given strain. All subclades of the given clade are removed
   (clades that would specify a species within this potential
   phylogenetic species). This ensures monophyly and that each species
   contains at least two strains.

After all strains are placed into the least inclusive clade, the
phylogenetic species tree is printed in newick format to the STDOUT
(standard output).


## Theory behind this method

You can read more about the underlying theory in the [md/theory.md](md/theory.md)
file.
