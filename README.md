Genomic surveillance of a resurgence of COVID-19 in Guangzhou, China


# Summary
We conducted a genomic surveillance to a COVID-19 outbreak related with imported cases between March to April 2020 in Guangzhou, China.

# Dependencies

MAFFT 7.458
IQ-Tree2 rc2
RAxML 8.2.12
ggtree 1.14.6
python 3.7
R version 3.5.1


# Phylogenetic analysis
Multiple alignments of genome sequences were performed by using MAFFT v7.458(Katoh et al., 2002) and manually inspected by using MEGA v10.1.8(Sudhir et al., 2018). Given the bias of genome coverage of public genome and sequences in this study, part of the 5’ and 3’ untranslated region was removed and the aligned genome length was 29697 nucleotides. We explored the phylogenetic structure with maximum likelihood (ML)method. ML Phylogenies of large alignment (>6000 genomes) were inferred by using IQ-Tree2(Minh et al., 2020) (rc2) with the best-fitting substitution model parameters (GTR+F+R2) estimated by Model Finder and 1000 rapid bootstrapping replicates. Phylogenetic analyses of <200 viral genomes were performed by using RAxML v8.2.12(Stamatakis, 2014)with 1000 bootstrap replicates and employing the GTRGAMMA+I model. The generated phylogenetic trees were visualized with the R package ggtree v1.14.6(Yu et al., 2017).



# Figure scripts
