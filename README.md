# GCPSR (Genealogical Concordance Phylogenetic Species Recognition)

This package contains perl scripts for the implementation of the
genealogical concordance phylogenetic species recognition method.

## Prerequisites

You need to have perl 5 on your computer (on Mac and Linux this is by
default installed).

The scripts in this repository use the `Bio::Tree::Tree` package.

    # To install using cpan
    cpan Bio::Tree::Tree

## Data preprocessing 

It is important to ensure that the right input trees are used for the
GCPSR analysis. The majority-rule consensus trees produced from single
locus phylogenies have to meet the following requirements:
- Trees are rooted using the same outgroup
- The outgroup contains at least two strains
- The clades have bootstrap (BS) or Bayesian posterior probability
  (BPP) support values 
- Support values appear as internal node names
  (e.g. _(((A,B)**80**,(C,D)**92**)**100**,(X,Y)**100**)_ or with
  branch lengths: _(((A:0.03,B:0.02)**80**:0.07,(C:0.04,D:0.01)**92**:0.05)**100**:0.10,(X:0.07,Y:0.04)**100**:0.2)_)
- Trees are in either **newick** ("`*.nwk`") or **nexus**
  ("`*.nex.*.t`" or "`*.nex.*.tre`") format and have one of the file
  extensions recognized by the scripts

## Perl scripts

To print the usage information for the scripts, run

    perl concordance_non-discordance.pl -h

and

    perl exhaustive_subdivision.pl -h

To run both steps with a minimal support of two

    ./concordance_non-discordance.pl -count=2 *.nwk | ./exhaustive_subdivision.pl - -count=2 >gcpsr.nwk



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

### Parameters

#### Minimum support

Specified by `-min=<int>` using
`concordance_non-discordance.pl`. Default value is **95**.
Recommended value for trees with bootstrap support at least is **80**,
for trees with Bayesian posterior probability **0.95** (or **95** in percentage
representation).

Ensures that only highly supported clades are recognized as concordant
or phylogenetic species.

#### Concordance cutoff

Specified by `-count=<init>` using
`concordance_non-discordance.pl`. Default value is **2**.
It is recommended to keep this value as low as possible.

Ensures that only those clades have to be compared for discordance
that are concordant through multiple single locus phylogenies. This is
needed when analyzing large numbers of loci for which we do not know
whether they are under _balancing selection_ or other influences that
would compromise phylogenetic species recognition. Before identifying
definitive phylogenetic species all trees that suggest conflicting
groupings, should be carefully analyzed why they should be _"ignored"_
(e.g. mating type loci and TRI gene cluster in _Fusarium graminearum_
species complex are under balancing selection; or if the locus is
difficult to align, which produces uncertainty in the correct homology
estimation in the alignment.)

#### Phylogenetic species cutoff

Specified by `-count=<init>` using
`exhaustive_subdivision.pl`. Default value is **1**.
It is recommended have the majority of the loci supporting the
recognition of the phylogenetic species (e.g. 3 out of 4).

Ensures that only those lineages are recognized as phylogenetic
species that are supported by a large number of the single locus
phylogenies (preferably by the majority).

## Example
Example analyses using 4 test tree files:

    ./concordance_non-discordance.pl  files-for-testing/* -min=95 -count=1 | ./exhaustive_subdivision.pl - -count=3

This analysis recognizes three phylogenetic species in the test data
set:
- outgroup: (X,Y)
- PS1: (A,B,C)
- PS2: (D,E,F,G)

Because these clades were highly supported (>= 95) in at least 3
single locus phylogenies (out of the 4). Although, both (D,E) and
(F,G) clades were highly supported in 3 trees these are removed,
because the following clades (D,F) and (E,G) were also highly
supported that were in conflict with clades (D,E) and (F,G).

    ./concordance_non-discordance.pl  files-for-testing/* -min=95 -count=2 | ./exhaustive_subdivision.pl - -count=3

This analysis recognizes three phylogenetic species in the test data
set:
- outgroup: (X,Y)
- PS1: (A,B,C)
- PS2: (D,E)
- PS3: (F,G)

Because these clades were highly supported (>= 95) in at least 3
single locus phylogenies (out of the 4). Since a clade was only
considered as concordant if it is supported by at least two trees,
clades were (D,E) and (F,G) were not contradicted by other concordant
clades.


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

## How to use this method for phylogenetic species recognition

### Phase 1: looking for dominant/general trends

The first question, whether there are multiple loci supporting the
separation of the species in your data.

For this you can use a relatively large concordance cutoff
(`-count=<int>`, first script), so a small
number of conflicting loci do not destroy your general pattern.

### Phase 2: find which loci are in conflict with the general pattern

First, save the the result of the analysis in phase 1 to a file
(e.g. `trend.nwk`). This after using a relatively large concordance
cutoff for the first script and using that result for exhaustive
subdivision. This result contains potential phylogenetic species.

Check whether you get the same result if you set concordance cutoff to
1 (`-count=1`).

If not, then you have gene trees that contradict the general
trend/pattern in your data set (gene tree forest). These should be
identified. To do this, use the following command (this assumes that
the reference tree is `trend.nwk` and all the individual single locus
trees are in the `single_locus_phylogenies` folder and have the `.nwk`
extension):

    perl find_conflicting_tree.pl -ref=trend.nwk single_locus_phylogenies/*.nwk

Examine all the conflicting loci and try to identify why they are
showing a conflict. Could it be that they are under balancing
selection, which helps maintain ancestral polymorphism even through
speciation, blocking lineage sorting. Lineage sorting is the one of
the phenomenon that GCPSR is exploiting. Another option could be that
the locus has an unreliable alignment, which produces a tree with high
support, but should not be considered as a reliable tree. (Tree
estimation is based on assuming homologous relationship between
nucleotides in the same column in the alignments. If the alignment is
incorrect, then the tree does not reflect evolutionary relation, but
the error itself.)

Some conflicts are not under these influences and demonstrate true
discordance. Discordance indicates that we are examining below the
phylogenetic species level. "_The transition from concordance among
branches to incongruity among branches can be used to diagnose
species._" ([Taylor, J. W., Jacobson, D. J., Kroken, S., Kasuga, T., Geiser, D. M., Hibett, D. S., and Fisher, M. C. 2000. Phylogenetic species recognition and species concepts in fungi. _Fungal Genetics and Biology_ 31, 21–32.)](http://taylorlab.berkeley.edu/sites/default/files/taylorlab/publications/taylor2000.pdf)

### Phase 3: final analysis on the refined locus set

If all conflicting trees could be reasoned away, then run the analysis
once again on the trees that are kept with `-count=1` for the first
step. For the second script use `-count=<int>` where `<int>` should be
a number representing the majority of single locus gene trees used
(e.g. 6 out of 10, so `-count=6`).

The final result contains the phylogenetic species identified in your
data set.


## Theory behind this method

You can read more about the underlying theory in the [md/theory.md](md/theory.md)
file.
