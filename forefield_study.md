# LongMeta application

We used the longMeta pipeline to analyze three environmental whole metagenomics datasets. The samples were collected during two summer seasons (2014 and 2015) from three different proglacial systems sited in front of Midtre Lovénbreen (Svalbard), the Storglaciären (Sweden) and the Greenlandic ice sheet in proximity of point 601. This dataset has 65 soil whole shotgun metagenomics samples for a total of 3,531,237,696 trimmed and quality checked paired-end Illumina reads (2x100 bp and 2x150 bp) and 433 Gb of data.

The aim of this project was to explore taxonomy profiles of the microbial succession to study the rock weathering and nitrogen fixing biological processes and their protagonists.

In this page we will show how we integrated LongMeta scripts, alignment and assembly software to obtain gene and taxonomy profiles from this complex dataset.

#### 1. Assembly step

We first performed a co-assembly of all the 65 samples with the software MEGAHIT v1.1.3 (https://github.com/voutcn/megahit).

```bash
megahit -1 M1_R1_trimmed_paired.fastq.gz,M2_R1_trimmed_paired.fastq.gz,M3_R1_trimmed_paired.fastq.gz,M4_R1_trimmed_paired.fastq.gz,M5_R1_trimmed_paired.fastq.gz,M6_R1_trimmed_paired.fastq.gz,M7_R1_trimmed_paired.fastq.gz,M8_R1_trimmed_paired.fastq.gz,M9_R1_trimmed_paired.fastq.gz,M10_R1_trimmed_paired.fastq.gz,M11_R1_trimmed_paired.fastq.gz,M12_R1_trimmed_paired.fastq.gz,M13_R1_trimmed_paired.fastq.gz,M14_R1_trimmed_paired.fastq.gz,M15_R1_trimmed_paired.fastq.gz,M16_R1_trimmed_paired.fastq.gz,M17_R1_trimmed_paired.fastq.gz,M18_R1_trimmed_paired.fastq.gz,M19_R1_trimmed_paired.fastq.gz,M20_R1_trimmed_paired.fastq.gz,M21_R1_trimmed_paired.fastq.gz,M22_R1_trimmed_paired.fastq.gz,M23_R1_trimmed_paired.fastq.gz,S11c_R1_trimmed_paired.fastq.gz,S2c_R1_trimmed_paired.fastq.gz,S4a_R1_trimmed_paired.fastq.gz,S5c_R1_trimmed_paired.fastq.gz,S10_R1_trimmed_paired.fastq.gz,S12_R1_trimmed_paired.fastq.gz,S3a_R1_trimmed_paired.fastq.gz,S4b_R1_trimmed_paired.fastq.gz,S6a_R1_trimmed_paired.fastq.gz,S11a_R1_trimmed_paired.fastq.gz,S151617_R1_trimmed_paired.fastq.gz,S3b_R1_trimmed_paired.fastq.gz,S4c_R1_trimmed_paired.fastq.gz,S6b_R1_trimmed_paired.fastq.gz,S11b_R1_trimmed_paired.fastq.gz,S2a_R1_trimmed_paired.fastq.gz,S3c_R1_trimmed_paired.fastq.gz,S5a_R1_trimmed_paired.fastq.gz,S6c_R1_trimmed_paired.fastq.gz,G1a_R1_trimmed_paired.fastq.gz,G3a_R1_trimmed_paired.fastq.gz,G5a_R1_trimmed_paired.fastq.gz,G7a_R1_trimmed_paired.fastq.gz,G1b_R1_trimmed_paired.fastq.gz,G3b_R1_trimmed_paired.fastq.gz,G5b_R1_trimmed_paired.fastq.gz,G7b_R1_trimmed_paired.fastq.gz,G1c_R1_trimmed_paired.fastq.gz,G3c_R1_trimmed_paired.fastq.gz,G5c_R1_trimmed_paired.fastq.gz,G7c_R1_trimmed_paired.fastq.gz,G2a_R1_trimmed_paired.fastq.gz,G4a_R1_trimmed_paired.fastq.gz,G6a_R1_trimmed_paired.fastq.gz,G8a_R1_trimmed_paired.fastq.gz,G2b_R1_trimmed_paired.fastq.gz,G4b_R1_trimmed_paired.fastq.gz,G6b_R1_trimmed_paired.fastq.gz,G2c_R1_trimmed_paired.fastq.gz,G4c_R1_trimmed_paired.fastq.gz,G6c_R1_trimmed_paired.fastq.gz,G8c_R1_trimmed_paired.fastq.gz -2 M1_R2_trimmed_paired.fastq.gz,M2_R2_trimmed_paired.fastq.gz,M3_R2_trimmed_paired.fastq.gz,M4_R2_trimmed_paired.fastq.gz,M5_R2_trimmed_paired.fastq.gz,M6_R2_trimmed_paired.fastq.gz,M7_R2_trimmed_paired.fastq.gz,M8_R2_trimmed_paired.fastq.gz,M9_R2_trimmed_paired.fastq.gz,M10_R2_trimmed_paired.fastq.gz,M11_R2_trimmed_paired.fastq.gz,M12_R2_trimmed_paired.fastq.gz,M13_R2_trimmed_paired.fastq.gz,M14_R2_trimmed_paired.fastq.gz,M15_R2_trimmed_paired.fastq.gz,M16_R2_trimmed_paired.fastq.gz,M17_R2_trimmed_paired.fastq.gz,M18_R2_trimmed_paired.fastq.gz,M19_R2_trimmed_paired.fastq.gz,M20_R2_trimmed_paired.fastq.gz,M21_R2_trimmed_paired.fastq.gz,M22_R2_trimmed_paired.fastq.gz,M23_R2_trimmed_paired.fastq.gz,S11c_R2_trimmed_paired.fastq.gz,S2c_R2_trimmed_paired.fastq.gz,S4a_R2_trimmed_paired.fastq.gz,S5c_R2_trimmed_paired.fastq.gz,S10_R2_trimmed_paired.fastq.gz,S12_R2_trimmed_paired.fastq.gz,S3a_R2_trimmed_paired.fastq.gz,S4b_R2_trimmed_paired.fastq.gz,S6a_R2_trimmed_paired.fastq.gz,S11a_R2_trimmed_paired.fastq.gz,S151617_R2_trimmed_paired.fastq.gz,S3b_R2_trimmed_paired.fastq.gz,S4c_R2_trimmed_paired.fastq.gz,S6b_R2_trimmed_paired.fastq.gz,S11b_R2_trimmed_paired.fastq.gz,S2a_R2_trimmed_paired.fastq.gz,S3c_R2_trimmed_paired.fastq.gz,S5a_R2_trimmed_paired.fastq.gz,S6c_R2_trimmed_paired.fastq.gz,G1a_R2_trimmed_paired.fastq.gz,G3a_R2_trimmed_paired.fastq.gz,G5a_R2_trimmed_paired.fastq.gz,G7a_R2_trimmed_paired.fastq.gz,G1b_R2_trimmed_paired.fastq.gz,G3b_R2_trimmed_paired.fastq.gz,G5b_R2_trimmed_paired.fastq.gz,G7b_R2_trimmed_paired.fastq.gz,G1c_R2_trimmed_paired.fastq.gz,G3c_R2_trimmed_paired.fastq.gz,G5c_R2_trimmed_paired.fastq.gz,G7c_R2_trimmed_paired.fastq.gz,G2a_R2_trimmed_paired.fastq.gz,G4a_R2_trimmed_paired.fastq.gz,G6a_R2_trimmed_paired.fastq.gz,G8a_R2_trimmed_paired.fastq.gz,G2b_R2_trimmed_paired.fastq.gz,G4b_R2_trimmed_paired.fastq.gz,G6b_R2_trimmed_paired.fastq.gz,G2c_R2_trimmed_paired.fastq.gz,G4c_R2_trimmed_paired.fastq.gz,G6c_R2_trimmed_paired.fastq.gz,G8c_R2_trimmed_paired.fastq.gz --k-min 27 --k-max 147 --k-step 8 -t 64 -m 0.95 -o Forefield_Assembly

##input files:
#*_R1_trimmed_paired.fastq.gz : forward Illumina reads
#*_R2_trimmed_paired.fastq.gz : reverse Illumina reads

##output file:

#Forefield_Assembly/final.contigs.fa : assembly file
````

#### 2. Assembly Check

We checked the assembly with the script longMeta-summary. Thanks to this script we were also able to trim the assembly retrieving only sequences longer than 300 bp. This script also outputs a file reporting the length of each contig and will be used in the next pipeline steps.

```bash
longMeta-summary --assembly-input Forefield_Assembly/final.contigs.fa --minimum-length 300 --assembly-output assembly_forefield.fasta --length-output assembly_forefield_length.txt

# |  File              |  Contig number   |  Total bp       |  Mean  |  Median  |  Max length   |  Min length |  n50  |  GC%  |
# |  final.contigs.fa  |  41841465        |  30378371287    |  726   |  17208   |  561967       |  300        |  841  |  62   |

##input file:

#final.contigs.fa : assembly input file

##output files:

#assembly_forefield.fasta: assembly output file, file containing only contigs longer than -min
#assembly_forefield_length.txt : length output file, tab delimited file reporting contigs and the correspondent lengths
```

#### 3. Diamond protein assignment

We used DIAMOND 0.9.22 (http://www.diamondsearch.org) to map known proteins to the assembly. You can download the database of known protein from the LongMeta GitHUB page.

```bash

#Format the nr.gz database for Diamond software
zcat nr.gz | diamond makedb -d nr

#Run Diamond
mkdir tmp
diamond blastx -d nr -q assembly_forefield.fasta -o assembly_forefield_nr.m8 -t tmp -e 0.000001 -F 15 --range-culling --range-cover 20 --id 50 --top 10 -f 6 -p 55 -c1 -b4.0
gzip assembly_forefield_nr.m8

##input files:

#assembly_forefield.fasta : assembly file
#nr : Diamond-formatted database; it can be downloaded from the main GitHUB LongMeta page

##output file:

#assembly_forefield_nr.m8.gz : Diamond output (blast tabular format)
```

#### 4. Assign taxonomy and coding region

We assigned taxonomy and coding regions to each contig with longMeta-assignment. This script screens the Diamond output file (tabular format). It then identifies the coding regions by assigning non-overlapping Diamond proteins to contigs. Taxonomy is inferred with a modified Lowest Common Ancestor algorithm (LCA).

```bash
longMeta-assignment --assembly assembly_forefield_nr.m8.gz --threads 3 --taxonomy h --acc2taxid accession_taxid.txt --taxid2taxon taxid_taxonomy.txt --gene-output assembly_forefield_nr_gene.txt --taxonomy-output assembly_forefield_nr_taxa.txt --tmp tmp > assembly_forefield_nr_assignment.log

##input files:

#assembly_forefield_nr.m8.gz : Diamond file for assembly (blast tabular format)
#accession_taxid.txt : tab delimited file reporting the protein accession numbers and the associated taxids
#taxid_taxonomy.txt : tab delimited file reporting the taxids and the associated taxonomical paths

##output files:

#assembly_forefield_nr_gene.txt : gene output file, Diamond file reporting only best non overlapping assigned proteins
#assembly_forefield_nr_taxa.txt : taxonomy output file,  file with taxonomy path associated to each contig

head assembly_forefield_nr_assignment.log
#The number sequences in the Diamond file is 33796900
#The average number of assigned genes per sequence is 1.28
#27573258 of sequences have 1 gene
#4774988 of sequences have 2 genes
#844816 of sequences have 3 genes
#280086 of sequences have 4 genes
#123960 of sequences have 5 genes
#65042 of sequences have 6 genes
#37919 of sequences have 7 genes
#24350 of sequences have 8 genes

tail assembly_forefield_nr_assignment.log
#1 of sequences have 257 genes
#1 of sequences have 259 genes
#1 of sequences have 274 genes
#1 of sequences have 300 genes
#1 of sequences have 324 genes
#1 of sequences have 349 genes
#1 of sequences have 362 genes
#1 of sequences have 442 genes
#The number of taxonomy assigned sequences with the BEST algorithm is 288266
#The number of taxonomy assigned sequences with the LCA algorithm is 33508634
```

#### 5. Chimera detection

Chimeric contigs can be an issue in high-complexity and low-coverage assemblies. We used the script longMeta-chimera to detect and split contigs identified as chimeric.

```bash
longMeta-chimera --gene-input assembly_forefield_nr_gene.txt --gene-output assembly_forefield_nr_gene_chim.txt  --taxonomy-input assembly_forefield_nr_taxa.txt --taxonomy-output assembly_forefield_nr_taxa_chim.txt --assembly-input assembly_forefield.fasta --assembly-output assembly_forefield_chim.fasta --length-input assembly_forefield_length.txt --length-output assembly_forefield_length_chim.txt --acc2taxid accession_taxid.txt --taxid2taxon taxid_taxonomy.txt

#The number of contigs is 33796900
#The contigs to check are 446452
#Chimeric contigs were searched at the genus level
#Total number of chimeric contig was 1582 which represent the 0.005% of the imported contigs

##input files:

#assembly_forefield_nr_gene.txt : gene input file
#assembly_forefield_nr_taxa.txt : taxonomy input file
#assembly_forefield.fasta : assembly input file
#assembly_forefield_length.txt : length input file
#accession_taxid.txt : tab delimited file reporting the protein accession numbers and the associated taxids
#taxid_taxonomy.txt : tab delimited file reporting the taxids and the associated taxonomical paths

##output files:

#assembly_forefield_nr_gene_chim.txt : gene output file
#assembly_forefield_nr_taxa_chim.txt : taxonomy output file
#assembly_forefield_chim.fasta : assembly output file
#assembly_forefield_length_chim.txt : length output file
```

#### 6. Read mapping

Illumina reads were mapped back to the assembly using bowtie2 v 2.3.4.3 (https://github.com/BenLangmead/bowtie2).

```bash
#build bowtie2 index
mkdir Map
bowtie2-build --threads 30 -o 3 assembly_forefield_chim.fasta Map/assembly_forefield

#run bowtie2 with all the samples
for file in *R1_trimmed_paired.fastq.gz
do
  name=${file%%_R1_trimmed_paired.fastq.gz}
  echo ${name}.sam 2>&1 | tee -a alignment_to_assembly.log
  bowtie2 --phred33 --local -I 100 -X 800 --no-hd --no-unal -D 30 -R 3 -N 0 -L 20 -i S,1,0.25 --non-deterministic --threads 20 -x Map/assembly_forefield -1 ${name}_R1_trimmed_paired.fastq.gz -2 ${name}_R2_trimmed_paired.fastq.gz -S ${name}.sam 2>&1 | tee -a alignment_to_assembly.log
  gzip ${name}.sam
done

##input files:

# assembly_forefield_chim.fasta : assembly file
#*_R1_trimmed_paired.fastq.gz : forward Illumina reads
#*_R2_trimmed_paired.fastq.gz : reverse Illumina reads

##output files:

#*.sam.gz : mapping files (sam format)
```

#### 7. Calculate taxonomy and gene abundance

We calculated taxonomy and gene coverage with the script longMeta-coverage. This script calculates base coverage for each gene and taxon present in each sample (i.e. each sam file).

```bash
mkdir AllForefield
longMeta-coverage --taxonomy-input assembly_forefield_nr_taxa_chim.txt --sam-input-tab listSam.txt --length-input assembly_forefield_length_chim.txt --output-folder AllForefield --alignment-AS 40 --gene-input assembly_forefield_nr_gene_chim.txt --acc2gene accession_protein.txt --acc2go accession_GO.txt

##input files:

#assembly_forefield_nr_taxa_chim.txt : taxonomy input file
#assembly_forefield_length_chim.txt : length input file
#assembly_forefield_nr_gene_chim.txt : gene input file
#listSam.txt : tab delimited file reporting sam files in the first column, read length can be specified in second column
#accession_protein.txt : tab delimited file reporting the protein accession numbers and the associated protein names
#accession_GO.txt : tab delimited file reporting the protein accession numbers and the associated Gene Ontology (GO) accession codes

##output files:

#AllForefield/taxon_profile/genus.txt : tab delimited file reporting taxon base coverage at genus level
#AllForefield/taxon_profile/family.txt : tab delimited file reporting taxon base coverage at family level
#AllForefield/taxon_profile/order.txt: tab delimited file reporting taxon base coverage at order level
#AllForefield/taxon_profile/class.txt : tab delimited file reporting taxon base coverage at class level
#AllForefield/taxon_profile/phylum.txt : tab delimited file reporting taxon base coverage at phylum level
#AllForefield/taxon_profile/domain.txt : tab delimited file reporting taxon base coverage at domain level

#AllForefield/gene_profile/gene_total.txt : tab delimited file reporting gene base coverage
#AllForefield/gene_profile/gene_to_go.txt : tab delimited file reporting correspondences between Gene Ontology entries and gene names

#AllForefield/gene_profile/genus/* : tab delimited files reporting gene base coverage at genus level
#AllForefield/gene_profile/family/* : tab delimited files reporting gene base coverage at family level
#AllForefield/gene_profile/order/* : tab delimited files reporting gene base coverage at order level
#AllForefield/gene_profile/class/* : tab delimited files reporting gene base coverage at class level
#AllForefield/gene_profile/phylum/* : tab delimited files reporting gene base coverage at phylum level
#AllForefield/gene_profile/domain/* : tab delimited files reporting gene base coverage at domain level
```

#### 8. Relative abundance for taxonomy data

Relative taxonomical abundances for each samples were calculated with the script longMeta-relative. This script bases the abundance calculation on the base coverages at genus level that were calculated by longMeta-coverage.

```bash
longMeta-relative --input-folder AllForefield --taxid2taxon taxid_taxonomy.txt

##input files:

#AllForefield/taxon_profile/genus.txt : tab delimited file reporting taxon base coverage at genus level
#AllForefield/taxon_profile/family.txt : tab delimited file reporting taxon base coverage at family level
#AllForefield/taxon_profile/order.txt: tab delimited file reporting taxon base coverage at order level
#AllForefield/taxon_profile/class.txt : tab delimited file reporting taxon base coverage at class level
#AllForefield/taxon_profile/phylum.txt : tab delimited file reporting taxon base coverage at phylum level
#AllForefield/taxon_profile/domain.txt : tab delimited file reporting taxon base coverage at domain level

#taxid_taxonomy.txt : tab delimited file reporting the taxids and the associated taxonomical paths

##output files:

#AllForefield/taxon_profile/genus_relative.txt : relative tab delimited file reporting taxon base coverage at genus level
#AllForefield/taxon_profile/family_relative.txt : relative tab delimited file reporting taxon base coverage at family level
#AllForefield/taxon_profile/order_relative.txt: relative tab delimited file reporting taxon base coverage at order level
#AllForefield/taxon_profile/class_relative.txt : relative tab delimited file reporting taxon base coverage at class level
#AllForefield/taxon_profile/phylum_relative.txt : relative tab delimited file reporting taxon base coverage at phylum level
#AllForefield/taxon_profile/domain_relative.txt : relative tab delimited file reporting taxon base coverage at domain level
```

#### 9. Explore gene profiles

Explore the dataset from the gene and taxonomical point of view. In our study we wanted to identify the organisms contributing to nitrogen fixation and rock weathering processes. Therefore we investigated i) nitrogenase genes, ii) oxalate biosinthesis genes (obcA), iii) cyanide synthase genes and iv) siderophore-related genes.

```bash
#i) nitrogenase genes
#total
longMeta-explore --input-folder AllForefield --gene 'nitrogenase' --rank genus --output-file nitrogenase_genus.txt

##input files:

#AllForefield/gene_profile/genus/* : tab delimited files reporting gene base coverage at genus level

##output files:

#nitrogenase_genus.txt : tab delimited files reporting nitrogenase base coverage for each genus associated with at least one nitrogenase coding region
#nitrogenase_genus.txt.all : average coding region base coverage for each sample

head nitrogenase_genus.txt
# taxon	gene	S6c.sam.gz	S6b.sam.gz	S6a.sam.gz	S5c.sam.gz	S5a.sam.gz	S4c.sam.gz	S4b.sam.gz	S4a.sam.gz	S3c.sam.gz	S3b.sam.gz	S3a.sam.gz	S2c.sam.gz	S2a.sam.gz	S151617.sam.gz	S12.sam.gz	S11c.sam.gz	S11b.sam.gz	S11a.sam.gz	S10.sam.gz	M9.sam.gz	M8.sam.gz	M7.sam.gz	M6.sam.gz	M5.sam.gz	M4.sam.gz	M3.sam.gz	M23.sam.gz	M22.sam.gz	M21.sam.gz	M20.sam.gz	M2.sam.gz	M19.sam.gz	M18.sam.gz	M17.sam.gz	M16.sam.gz	M15.sam.gz	M14.sam.gz	M13.sam.gz	M12.sam.gz	M11.sam.gz	M10.sam.gz	M1.sam.gz	G8c.sam.gz	G8b.sam.gz	G8a.sam.gz	G7c.sam.gz	G7b.sam.gz	G7a.sam.gz	G6c.sam.gz	G6b.sam.gz	G6a.sam.gz	G5c.sam.gz	G5b.sam.gz	G5a.sam.gz	G4c.sam.gz	G4b.sam.gz	G4a.sam.gz	G3c.sam.gz	G3b.sam.gz	G3a.sam.gz	G2c.sam.gz	G2b.sam.gz	G2a.sam.gz	G1c.sam.gz	G1b.sam.gz	G1a.sam.gz
# Achromobacter	dinitrogenase reductase	0.50000	0	0.28571	0	0.14286	0	0.02381	0	4.16667	2.14286	3.11905	0	0.14286	0	0	1.11905	1.92857	1.50000	0	0.09524	0	0	0	0	0.09524	0	0	0	0	0.14286	0.14286	0.14286	0	0	0	0	0	0	0	0	0	0.28571	0.28570.14286	0.14286	0
# Acidobacteria_Unclassified	Nitrogenase (molybdenum-iron)-specific transcriptional regulator NifA	0	0.72586	0.31646	0	0	0.34965	0.31646	0	0.34965	0	0.34960.19814	0	0	0	0	0	0	0	0.23310	0.23310	0	0	0	0	0	0	0.35865	0.42194	0	0.17483	0	0.31646	0	0	0.62037	0.61336	0.34965	0	0	0	0.69930	0.52743	0	0	0.52743	0	0	0	0.17483	0	0
# Actinobacteria_Unclassified	dinitrogenase reductase activating glycohydrolase	0	0	0	0	0	0	0	0	0	0	0	0	0	0.18727	0	0.28090	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1.07678	1.68539	0	0.28090
# Actinobacteria_Unclassified	nitrogenase molybdenum-iron protein subunit beta	0	0	0	0	0	0	0	0	0	0	0	0	0	2.52294	0.45872	0.91743	0	0.30581	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.45872	0	0	0	0	0.38226	0	0	0	0	0	0	0	0	0	0	0	0	0	0
# Alphaproteobacteria_Unclassified	dinitrogenase reductase	0	0	0	0	0	0	0	0	1.77596	1.36612	1.50273	0	0	0	0	0	1.09290	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.40984	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
# Bacteria_Unclassified	Nitrogenase iron protein 2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1.93798	2.03480.58140	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
# Bacteria_Unclassified	dinitrogenase reductase	0.58716	0	0.32986	0.67157	2.43679	0	0	0	0.27903	1.21060	0.80156	0.37974	0.72338	0	0	1.08333	0.62232	0.50000.35417	0	0	0.11111	0	0.39216	0	0	0.62500	0	0.30581	0	0.20833	0	0.11111	0	0	0	0.30581	0.22222	0.06250	0	0	0	0.16667	0.49020	0.34396	0	0.92593	0.92593	0	0	0	0.23039	0.72046	0.49020	0.43860	0.53365	0.62500	0	0	0	0	0.53180	0
# Bacteria_Unclassified	dinitrogenase reductase activating glycohydrolase	0	0.15244	0.38314	0	0.07622	0.15244	0.30488	0.26779	0.42373	0.30488	0	0	0	0.15241.61023	0	0	0	0.29612	0	0	0	0	0	0	0	0.88857	0.50441	0	0	0	0	0.10163	0	0	0	0	0	0.10163	0	0.11467	0	0	0	0	0.15244	0.15244	0	0	0	0	0	0.19055	0	0.15244	0	0	0.45732	0.15244	5.27287	7.49725	3.14980.59260	0.13974	0.38314
# Bacteria_Unclassified	dinitrogenase reductase activating glycohydrolase (draG)	0	0	0.81605	0	0	0	1.61290	1.14247	0.73925	0.67204	0	0	0	0.40323	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.81845	0	0	0	3.05060	0.81845	1.33929

#genus level
longMeta-explore --input-folder AllForefield --gene 'nitrogenase' --output-file nitrogenase.txt

##input file:

#AllForefield/gene_profile/gene_total.txt : tab delimited file reporting gene base coverage

#output files:

#nitrogenase_genus.txt : tab delimited files reporting nitrogenase base coverage for each sample
#nitrogenase.txt.all : average coding region base coverage for each sample

head nitrogenase.txt
# gene	S6c.sam.gz	S6b.sam.gz	S6a.sam.gz	S5c.sam.gz	S5a.sam.gz	S4c.sam.gz	S4b.sam.gz	S4a.sam.gz	S3c.sam.gz	S3b.sam.gz	S3a.sam.gz	S2c.sam.gz	S2a.sam.gz	S151617.sam.gz	S12.sam.gz	S11c.sam.gz	S11b.sam.gz	S11a.sam.gz	S10.sam.gz	M9.sam.gz	M8.sam.gz	M7.sam.gz	M6.sam.gz	M5.sam.gz	M4.sam.gz	M3.sam.gz	M23.sam.gz	M22.sam.gz	M21.sam.gz	M20.sam.gz	M2.sam.gz	M19.sam.gz	M18.sam.gz	M17.sam.gz	M16.sam.gz	M15.sam.gz	M14.sam.gz	M13.sam.gz	M12.sam.gz	M11.sam.gz	M10.sam.gz	M1.sam.gz	G8c.sam.gz	G8b.sam.gz	G8a.sam.gz	G7c.sam.gz	G7b.sam.gz	G7a.sam.gz	G6c.sam.gz	G6b.sam.gz	G6a.sam.gz	G5c.sam.gz	G5b.sam.gz	G5a.sam.gz	G4c.sam.gz	G4b.sam.gz	G4a.sam.gz	G3c.sam.gz	G3b.sam.gz	G3a.sam.gz	G2c.sam.gz	G2b.sam.gz	G2a.sam.gz	G1c.sam.gz	G1b.sam.gz	G1a.sam.gz
# ADP-ribosyl-(nitrogenase)-activating glycohydrolase	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.84034	1.19040.77031	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
# Dinitrogenase reductase ADP-ribosyltransferase (DRAT)	0	0	0	0	0	0	0	0	2.20798	0	2.42165	0	0	0.35613	0	0	0	0
# Fe-only nitrogenase accessory protein AnfO	0.44429	0	0.52654	0	0	0	0	0	0	0	0	0.30120	0.85341	0	1.40562	0	0.71496	0.26040.48769	0	0	0	0	0	0	0	0.27401	0.20080	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
# Mo-dependent nitrogenase	0.58480	0	0.07184	0	0	0	0	0	0	0	0	0	0	1.33294	9.36594	0.07184	0	0	1.45062	0.28731.20690	0.28736	0	0.40898	0	0.19941	0.68942	2.43116	0	0	0.45045	0	0	0	0	0.86207	0.57471	0	0	0	0	0.75355	0	0	6.38017	3.61381	10.34736	0	0	0
# Mo-dependent nitrogenase C-terminal domain family protein	0	0	0	0	0	0	0	0	0	0	0	0	0	1.74419	3.87597	0	0.47954	0.70241	0	0	0	0	0	0	0	0	0	0	0	0	0	0.74074	0	0	0	0	0	0	0	0	4.99785	2.93497	8.83506	0	0	0
# Mo-dependent nitrogenase family protein	0	0.37968	0.36873	0	0.88496	1.40118	0	0	0	0	0	0	0	13.53801	29.32815	0	0	0.79745	1.17994	0	0.36232	1.03245	0.14749	0.22140	0.38217	13.06514	3.23241	0	0	0.15949	0	0.28612	0.55447	0.38760	1.12827	0.48981	0	0.29499	0.29601	0.17653.76130	0	0	0	0	0	0	0	0	0	0	0.36873	0	0	0.36873	0	0	0	0	9.74989	8.74455	14.77808	0.19530.39062	0.39062
# Mo-nitrogenase MoFe protein subunit NifD precursor	0.89400	1.25438	1.18135	1.59414	0.72386	0.68257	1.89899	0	0.58003	0	1.83309	0	2.62353	0.73455	1.02041	1.21334	1.61201.15733	0.41279	0.09183	0.66788	0.25918	1.23858	0.31702	0.23812	0	0	0.69878	0	0.12771	0.20576	0.25147	0.25543	0.12771	0	0	0.39010	0	0.74713	0.51086	0.87740.25543	0.19157	0.51020	0.41850	0	0	0	0.85034	0.15432	0.92593	0.86412	0.15432	0.67319	0.30864	0.55096	0	0	0	0	0	0	0	0	0	0
# Mo-nitrogenase MoFe protein subunit NifK	0.55975	1.17210	0.62798	1.36225	0.29183	0.62528	2.69113	0.35311	0.14441	0.19960	0.19455	0.09728	0.35311	0.12286	0	1.41750	0.47725	0.13900.61987	0.08431	0.45733	0.06653	0.06653	0.06485	0.06485	0.13054	0	0	0	0	0.19455	0.13307	0.29273	0.12970	0	0	0.06485	0	0.39266	0.28603	1.40096	0.44900.19960	0	0.09728	0.53196	0	0.72957	0.19455	0	0	0	0	0	0	0	0	0	0	0
# Mo-nitrogenase iron protein subunit NifH	0	0	0	0	0	0.17361	0	0	1.12360	0	0.69444	0	1.76505	0.34722	0	1.39098	2.46338	2.04410.55221	0.19403	1.30257	0	1.06522	0.72097	0	0	0	0.18727	0	0	0.92593	0	0	0	0	0.18727	1.47082	0.58989	0	0	1.53895	1.04931.43120	0	1.01128	0	0.28090	0	0	0	0.30120	0	0	0	0	0	0

#ii) oxalate biosynthesis genes
#total
longMeta-explore --input-folder AllForefield --gene '3-keto-5-aminohexanoate cleavage' --rank genus --output-file obc_genus.txt
#genus level
longMeta-explore --input-folder AllForefield --gene '3-keto-5-aminohexanoate cleavage' --output-file obc.txt

#iii) cyanide synthase genes
#total
longMeta-explore --input-folder AllForefield --gene 'cyanide synthase' --rank genus --output-file cyanide_genus.txt
#genus level
longMeta-explore --input-folder AllForefield --gene 'cyanide synthase' --output-file cyanide.txt

#iv) siderophore related genes
#total
longMeta-explore.pl --input-folder AllForefield --gene 'siderophore' --rank genus --output-file siderophore_genus.txt
#genus level
longMeta-explore --input-folder AllForefield --gene 'siderophore' --output-file siderophore.txt
```

