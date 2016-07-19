# GCPSR (Genealogical Concordance Phylogenetic Species Recognition)

This package contains perl scripts for the implementation of the
genealogical concordance phylogenetic species recognition method.

# GCPSR _sensu_ [Dettman _et al._ (2003)](http://dx.doi.org/10.1554/03-073)

1. Grouping of individuals

    A clade was recognized as an independent evolutionary lineage if
    it satisfied either of two criteria:

    - Genealogical concordance:

        The clade was present in the majority of the
        single-locus genealogies

    - Genealogical nondiscordance:

        The clade was well supported in at least one single-locus
        genealogy, as judged by both MP bootstrap proportions
        and Bayesian posterior probabilities, and was not
        contradicted in any other single-locus genealogy at the same
        level of support

2. Ranking of groups

    When deciding which independent evolutionary lineages represented
    phylogenetic species, characteristics of lineages in combined
    data analyses were also considered. Two ranking criteria were
    applied:

    - Genetic differentiation

        to prevent minor tip clades from being recognized,
        phylogenetic species had to be relatively distinct and well
        differentiated from other species

    - Exhaustive subdivision:

        all individuals had to be placed within a phylogenetic species
