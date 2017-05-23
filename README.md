# GCPSR (Genealogical Concordance Phylogenetic Species Recognition)

This package contains perl scripts for the implementation of the
genealogical concordance phylogenetic species recognition method.

## GCPSR _sensu_ [Dettman _et al._ (2003)](http://dx.doi.org/10.1554/03-073)

1. Grouping of individuals

    A clade is recognized as an independent evolutionary lineage if
    it satisfied either of two criteria:

    - Genealogical concordance:

        The clade is present in the majority of the
        single-locus genealogies.

    - Genealogical nondiscordance:

        The clade is well supported in at least one single-locus
        genealogy, as judged by both MP bootstrap proportions
        and Bayesian posterior probabilities, and is not
        contradicted in any other single-locus genealogy at the same
        level of support.

2. Ranking of groups

    This steps identifies which independent evolutionary lineages
    could be considered as phylogenetic species.
	Two ranking criteria were applied:

    - Genetic differentiation

        To prevent minor tip clades from being recognized,
        phylogenetic species have to be relatively distinct and well
        differentiated from other species.

    - Exhaustive subdivision:

        All individuals have to be placed within a phylogenetic species.

## Perl scripts

To run both steps with a minimal support of two

    ./concordance_non-discordance.pl -min=2 *.nwk | ./exhaustive_subdivision.pl - >gcpsr.nwk
