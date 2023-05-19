# NIPT-FF

A pipeline utilizes maternal and parental iSNPS to obtain fetal fraction in cfDNA of pregnant women

## Prerequisite

- Mother, father and pregnant women NGS sequencing fastq data
- Basic NGS analysis tools including the GATK bundle
- HPC environment and an open-source workload manager and job scheduler like slurm (optional)

## Main Steps

## Mapping

Since this is a relatively simple pipeline and BWA automatically performs adapter trimming, no adapter trimming and QC tools are implemented here.

## Marking duplicates

Hybridization capture is used for building the library, and deduplication is necessary here

## Calling SNPs

GATK HaplotypeCaller and Vardict are applied in this pipeline. Note that the results of the two tools are not merged. For pregnant women's cfDNA, it is better to run Vardict for calling genotyping and allele frequency because GATK's algorithm is complex and not suitable to detect low frequency AF with its default parameters.

## Merging GVCF and performing genotyping

For GATK HaplotypeCaller results, the parental and maternal GVCF should be merged first, and performed joint genotyping on the two samples pre-called with HaplotypeCaller.

## Calculating fetal fraction





