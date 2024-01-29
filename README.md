# **OPTIMIZATION OF WHOLE EXOME SEQUENCE PIPELINE: HUMAN REFERENCE PAN/GENOME**
Master Thesis (Masterarbeit)

Submitted in Partial Fulfillment of the Degree Requirements Master of Science (M.Sc.) in Biotechnology

at the Department of Biotechnology, University of Applied Sciences Mannheim (Hochschule Mannheim), Germany.

by Quynh Nhu Nguyen

Under the Guidance of Dr. Phuc Loi Luu of “on-site Supervisor” and Prof. Dr. Markus Gumbel, Hochschule Mannheim.

Mannheim, February 2024 

## Abstract

### Introduction
In genomics research, haploid linear reference genomes of superior quality are indispensable. However, their use can introduce biases where DNA fragments carrying the reference allele are more likely to map accurately or receive higher quality scores. The influence of bias may cause a minor distortion in the outcomes of diverse genetic analyses, including the estimation of heterozygosity, the investigation of population affinities, and the possible neglect of specific elements in clinical and research contexts. To get around this problem, new methods have come up that try to replace the linear reference with a variation graph containing known alternative variants at each genetic locus. In this context, we examine the effectiveness of applying the Human Whole Exome Sequence (WES) Pangenome in Developmental and Epileptic Encephalopathy syndrome (DEE) and Hereditary Breast and Ovarian Cancer syndrome (HBOC) in improving mapping accuracy. This approach aims to address reference or alternative bias, provide better coverage for pathogenic single-nucleotide polymorphisms (SNPs), as well as small insertions and deletions (Indels).

### Methods
We advocated employing Whole Exome Sequencing (WES) for the creation of a pangenome, which we subsequently applied to clinical sample data. Our approach involved utilizing phased, diploid short reads obtained from 5 WES epilepsy blood samples (4 for training, 1 for testing) sourced from the Children's Hospital 2 in Ho Chi Minh City, Vietnam. Additionally, we incorporated 12 WES breast carcinoma samples (9 for training, 3 for testing) from the BioProject PRJNA489865 into our study. The training sets were used for constructing the pangenome, while the test sets were for evaluating read alignment against the pangenome efficiently.
Through the nf-core/pangenome pipeline, we generated the DEE-WES Pangenome and HBOC-WES Pangenome. This process commenced with a reference assembly GRCh38, and we incorporated additional assemblies, created from VCFs obtained from nf-core sarek pipeline, into a variation graph. Subsequently, we developed the Up/Downstream Pangenome WDL workflow, involving 4 main steps: 1) Pangenome short read mapping, 2) Germline variant calling, 3) Annotation, and 4) Comprehensive reporting. 

### Results
Pangenome analysis integrated into both DEE and HBOC patients brought substantial advancements in mapping accuracy, allele representation and bias reduction, notably benefiting pathogenic heterozygous alleles. These improvements were more pronounced in cancer cases, underscoring the efficacy of pangenome analysis in managing the genetic diversity present in cancer genomes.

#### Mapping Accuracy
In DEE, mapped and paired reads increased by 6.963%, accompanied by a 6.020% rise in properly paired reads and a 1.177% decrease in reads mapped with MQ=0. HBOC showed a 3.018% rise in mapped and paired reads, with a 6.583% increase in properly paired reads and a 1.145% decrease in MQ=0 reads. Mean coverage surged by 6.638% (DEE) and 51.777% (HBOC).
Allele Representation: Coverage Enhancement
The DEE-WES Pangenome demonstrated elevated reads for pathogenic heterozygous alleles, with an average increase of 2.247 alternate allele reads and 4.419 reference allele reads. In HBOC, alternate allele reads experienced a surge of 20.828, and reference allele reads saw a rise of 19.142. Pathogenic homozygous alleles displayed a noteworthy increase, with a rise of 16.243 reads in DEE and 13.887 reads in HBOC. 

#### Allelic Balance 
Allelic balance improved by 0.094 in cancer, indicating a more balanced representation. Epilepsy maintained stability in allelic balance.

#### Bias Reduction
Pangenome analysis reduced both reference and alternate biases in cancer (DEE and HBOC). In cancer, alternate bias (β < 0.45) reduced by 0.022, and reference bias (β > 0.55) decreased by 0.054. DEE saw a decrease in alternate bias by 0.043.

#### Comparison
Differences in pangenome impact between DEE and HBOC stem from unique genetic structures. Cancer genomes, with greater heterogeneity, benefit more from pangenome analysis. Conversely, genetic stability in disorders like Epilepsy diminishes pangenome influence.

#### Conclusion
WES Pangenome analysis shows promise in refining mapping accuracy and reducing biases, particularly in cancer genomes. It proves valuable for clinical diagnostics, enabling tailored treatments and revolutionizing personalized medicine.

Keywords: allelic balance, alternate alleles, alternate bias, breast cancer, coverage, epilepsy, heterozygous alleles, pangenome, reference bias, reference alleles, target drug, whole exome sequencing (WES).

