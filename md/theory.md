# Theoretical Background

The operational criteria for GCPSR for fungi were introduced in the theoretical
paper of [Taylor _et al._ (2000)]((http://dx.doi.org/10.1006/fgbi.2000.1228),
and first implemented as an operational framework using a two-step methodology by
[Dettman _et al._ (2003)](http://dx.doi.org/10.1554/03-073).
After genetic isolation, two populations (or species) undergo the following
stages: shared polymorphism (apparent polyphyly), loss of shared polymorphism
(after fixation in one of the species) and reciprocal monophyly (after fixation in
both species). According to GCPSR, the genetic isolation between populations
(phylogenetic species) can be detected by a combination of genealogical
concordance and non-discordance. Genealogical concordance is meant by the
concordant reciprocal monophyly of multiple gene genealogies. Genealogical
non-discordance means that no grouping supported by high support values for one of
the genes is contradicted by another gene with the same level of support.

"The strength of GCPSR lies in its comparison of more than one gene genealogy.
A requirement of each gene genealogy is that recombination does not occur within
the gene, and in practice, parts of genes are often used to construct the
genealogies. Where the different gene trees are concordant they have the same tree
topology due to fixation of formerly polymorphic loci following genetic isolation;
these concordant branches connect species. Conflict among the gene trees is likely
to be due to recombination among individuals within a species, and the transition
from concordance to conflict determines the limits of species (Fig. 2)." -
[Taylor _et al._ (2000)]((http://dx.doi.org/10.1006/fgbi.2000.1228)

![Fig. 2](md/figures/Transition.png)
**Fig. 2.** Simultaneous analysis of three gene genealogies shows how the transition from concordance among branches to incongruity among branches can be used to diagnose species. -
[Taylor _et al._ (2000)]((http://dx.doi.org/10.1006/fgbi.2000.1228)

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
