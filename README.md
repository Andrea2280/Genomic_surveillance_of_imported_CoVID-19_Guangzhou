# Genomic surveillance of a resurgence of COVID-19 in Guangzhou, China



## Summary
While China experienced a decline in daily growth rate of COVID-19, importations of new cases from other countries and their related local infections caused a rapid rise. In Guangzhou, a transport hub in South China, we conducted a genomic surveillance to a COVID-19 outbreak related with imported cases in March and April 2020. The SARS-CoV-2 of 125 imported cases from 25 countries were widely distributed in the global phylogeny, with novel lineages identified in countries lacking publicly available viral genomic data. Nonetheless, local outbreaks were traced to few specific lineages imported before the enhanced surveillance and control of COVID-19 in Guangzhou. The transmission networks were revealed, including cases without epidemiological information. Additionally, in ten cases, we found in-frame deletions on viral genome, without enrichment in the S protein. Intra-host variations analysis shows a stronger purifying selection on shared variations than the singletons, but the C-U substitutions have advantage in selection.

## Highlights
- A recent COVID-19 spike induced by foreign imported cases in Guangzhou city provides window into the genomic evolution and transmission of SARS-CoV-2.
- The imported viral strains are diverse and distributed broadly in the global phylogenetics tree.
- Our viral genome sequencing provides phylogenetic data on infected individuals from African and South/West Asia that are previously less reported.
- Local related cases could be traced back to two specific small lineages from Africa.
- Deletion mutations could be identified in the viral genome.
- C-U substitution might have selective advantage during viral transmissions.


## Dependencies

- MAFFT 7.458
- IQ-Tree2 rc2
- RAxML 8.2.12
- ggtree 1.14.6
- Biopython
- python 3.7
- R version 3.5.1


## Data
1. sample_info.csv : Sequencing sample information.
2. Public_SARS-CoV-2_geneomes_info.csv: Public SARS-CoV-2 genome sequences ID used in the analyses.  


## Phylogenetic analysis
Multiple alignments of genome sequences were performed by using MAFFT v7.458(Katoh et al., 2002) and manually inspected by using MEGA v10.1.8(Sudhir et al., 2018). Given the bias of genome coverage of public genome and sequences in this study, part of the 5’ and 3’ untranslated region was removed and the aligned genome length was 29697 nucleotides. We explored the phylogenetic structure with maximum likelihood (ML)method. ML Phylogenies of large alignment (>6000 genomes) were inferred by using IQ-Tree2(Minh et al., 2020) (rc2) with the best-fitting substitution model parameters (GTR+F+R2) estimated by Model Finder and 1000 rapid bootstrapping replicates. Phylogenetic analyses of <200 viral genomes were performed by using RAxML v8.2.12(Stamatakis, 2014)with 1000 bootstrap replicates and employing the GTRGAMMA+I model. The generated phylogenetic trees were visualized with the R package ggtree v1.14.6(Yu et al., 2017).


~~~
# Multi-Alignment
mafft --thread 10 --auto $fasta_file > $align_result_file

# Trim Alignment file
python Trim_UTR_multialign_fasta.py $align_result_file $align_trim_file

# Iqtree phylogenetic analysis
iqtree -s $align_trim_file --prefix iqtree_result -m GTR+F+R2 -B 1000 -cmax 30 -redo -T 24

# RAXML phylogenetic analysis
raxmlHPC-PTHREADS -s $align_trim_file -n raxmltree_result -m GTRGAMMAI -f a -x 12345 -N 1000 -p 123456 -T 24 -k
~~~
## Genome assembly and variations calling

1. To obtain consensus sequences and deletion mutations: `snakemake -s nCov_assembly_pipeline.py  -p`

2. iSNV and SNP mutations calling uses iSNV_calling.sh[]: bash iSNV_calling.sh  3,4,5

## variations annotation

1. To convert iSNV table into vcf format(VCFv4.1) that can be recognized  by SnpEff software: python iSNVTable_2_vcf_allsample.py  -i  ./data/   -o ./allsamplesVcfPath/   -r   MN908947.3
2. Annoting the mutation by using SnpEff software : `java -jar snpEff.jar ann  MN908947   ${sample}.vcf    >./{sample}.snpeffAnno.vcf`
3. To Integrate all the samples' mutations annotion files into a single file: `python Snpeff_results_integrateTo1file -i  ./bigtable2vcf/outsnpeff/   -o  .allsampMutatAnno.txt.txt`

4. To annotate the indel that were filted:`java -jar snpEff.jar ann  MN908947   ${sample}.indel.vcf    >./{sample}.indel.snpeffAnno.vcf`

## Figure scripts
