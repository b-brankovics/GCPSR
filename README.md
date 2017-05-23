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

## Theory behind this method

You can read more about the underlying theory in the (md/theory.md)
file.
