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
- pangolin V1.1.8(https://github.com/hCoV-2019/pangolin)
- Biopython
- python 2.7
- python 3.6
- R version 3.5.1

## Repo Contents
1. data: preprocessed data used for analysis.
2. result: result files for figure visualization.
3. scripts: python script and R script for analysis and visulization.

## Phylogenetic analysis
1. Multi-Alignment of the SARS-CoV-2 genomes : `mafft --thread 10 --auto $fasta_file > $align_result_file`

2. Trim Alignment file: `python ./scripts/Trim_UTR_multialign_fasta.py $align_result_file $align_trim_file`

3. Iqtree phylogenetic analysis: `iqtree -s $align_trim_file --prefix iqtree_result -m GTR+F+R2 -B 1000 -cmax 30 -redo -T 24`

4. RAXML phylogenetic analysis: `raxmlHPC-PTHREADS -s $align_trim_file -n raxmltree_result -m GTRGAMMAI -f a -x 12345 -N 1000 -p 123456 -T 24 -k`

5. pangolin lineage analysis: `pangolin $fasta_file -t 24`


## Genome assembly and variations calling

1. To obtain consensus sequences and deletion mutations: `snakemake -s nCov_assembly_pipeline.py  -p`

2. iSNV and SNP mutations calling uses iSNV_calling.sh[]: bash iSNV_calling.sh  3,4,5

## variations annotation

1. To convert iSNV table into vcf format(VCFv4.1) that can be recognized  by SnpEff software: python iSNVTable_2_vcf_allsample.py  -i  ./data/   -o ./allsamplesVcfPath/   -r   MN908947.3
2. Annoting the mutation by using SnpEff software : `java -jar snpEff.jar ann  MN908947   ${sample}.vcf    >./{sample}.snpeffAnno.vcf`
3. To Integrate all the samples' mutations annotion files into a single file: `python Snpeff_results_integrateTo1file -i  ./bigtable2vcf/outsnpeff/   -o  .allsampMutatAnno.txt.txt`

4. To annotate the indel that were filted:`java -jar snpEff.jar ann  MN908947   ${sample}.indel.vcf    >./{sample}.indel.snpeffAnno.vcf`


## Citation
If you use data, results or conclusion from this work, please cite:

## Acknowledgement
This study was supported by National Natural Science Foundation of China (31870079, 91953122, 31871326), National Science and Technology Major Project of the Ministry of Science and Technology of China (2017ZX10103011, 2018ZX10305410, 2018ZX10201001), Guangdong Provincial Novel Coronavirus Scientific and Technological Project (2020111107001), Guangdong Basic and Applied Basic Research Foundation (2020A1515010776 and 2020B1515020057) and the Beijing Nova Program (Z181100006218114 and Z181100006218110) to M.N. and P.L..

