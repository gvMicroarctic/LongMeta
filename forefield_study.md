# LongMeta application

We used the longMeta pipeline to analyze three environmental whole metagenomics datasets. The samples were collected during two summer seasons (2014 and 2015) from three different proglacial systems sited in front of Midtre Lovénbreen (Svalbard), the Storglaciären (Sweden) and the Greenlandic ice sheet in proximity of point 601. This dataset has 65 soil whole shotgun metagenomics samples for a total of 3,531,237,696 trimmed and quality checked paired-end Illumina reads (2x100 bp and 2x150 bp) and 433 Gb of data.

The aim of this project was to explore taxonomy profiles of the microbial succession to study the rock weathering and nitrogen fixing biological processes and their protagonists.

In this page we will show how we integrated LongMeta scripts, alignment and assembly software to obtain gene and taxonomy profiles from this complex dataset.

#### 1. Assembly step

We first performed a co-assembly of all the 65 samples with the software MEGAHIT v1.1.3 (https://github.com/voutcn/megahit).

```bash
megahit -1 M1_R1_trimmed_paired.fastq.gz,M2_R1_trimmed_paired.fastq.gz,M3_R1_trimmed_paired.fastq.gz,M4_R1_trimmed_paired.fastq.gz,M5_R1_trimmed_paired.fastq.gz,M6_R1_trimmed_paired.fastq.gz,M7_R1_trimmed_paired.fastq.gz,M8_R1_trimmed_paired.fastq.gz,M9_R1_trimmed_paired.fastq.gz,M10_R1_trimmed_paired.fastq.gz,M11_R1_trimmed_paired.fastq.gz,M12_R1_trimmed_paired.fastq.gz,M13_R1_trimmed_paired.fastq.gz,M14_R1_trimmed_paired.fastq.gz,M15_R1_trimmed_paired.fastq.gz,M16_R1_trimmed_paired.fastq.gz,M17_R1_trimmed_paired.fastq.gz,M18_R1_trimmed_paired.fastq.gz,M19_R1_trimmed_paired.fastq.gz,M20_R1_trimmed_paired.fastq.gz,M21_R1_trimmed_paired.fastq.gz,M22_R1_trimmed_paired.fastq.gz,M23_R1_trimmed_paired.fastq.gz,S11c_R1_trimmed_paired.fastq.gz,S2c_R1_trimmed_paired.fastq.gz,S4a_R1_trimmed_paired.fastq.gz,S5c_R1_trimmed_paired.fastq.gz,S10_R1_trimmed_paired.fastq.gz,S12_R1_trimmed_paired.fastq.gz,S3a_R1_trimmed_paired.fastq.gz,S4b_R1_trimmed_paired.fastq.gz,S6a_R1_trimmed_paired.fastq.gz,S11a_R1_trimmed_paired.fastq.gz,S151617_R1_trimmed_paired.fastq.gz,S3b_R1_trimmed_paired.fastq.gz,S4c_R1_trimmed_paired.fastq.gz,S6b_R1_trimmed_paired.fastq.gz,S11b_R1_trimmed_paired.fastq.gz,S2a_R1_trimmed_paired.fastq.gz,S3c_R1_trimmed_paired.fastq.gz,S5a_R1_trimmed_paired.fastq.gz,S6c_R1_trimmed_paired.fastq.gz,G1a_R1_trimmed_paired.fastq.gz,G3a_R1_trimmed_paired.fastq.gz,G5a_R1_trimmed_paired.fastq.gz,G7a_R1_trimmed_paired.fastq.gz,G1b_R1_trimmed_paired.fastq.gz,G3b_R1_trimmed_paired.fastq.gz,G5b_R1_trimmed_paired.fastq.gz,G7b_R1_trimmed_paired.fastq.gz,G1c_R1_trimmed_paired.fastq.gz,G3c_R1_trimmed_paired.fastq.gz,G5c_R1_trimmed_paired.fastq.gz,G7c_R1_trimmed_paired.fastq.gz,G2a_R1_trimmed_paired.fastq.gz,G4a_R1_trimmed_paired.fastq.gz,G6a_R1_trimmed_paired.fastq.gz,G8a_R1_trimmed_paired.fastq.gz,G2b_R1_trimmed_paired.fastq.gz,G4b_R1_trimmed_paired.fastq.gz,G6b_R1_trimmed_paired.fastq.gz,G2c_R1_trimmed_paired.fastq.gz,G4c_R1_trimmed_paired.fastq.gz,G6c_R1_trimmed_paired.fastq.gz,G8c_R1_trimmed_paired.fastq.gz -2 M1_R2_trimmed_paired.fastq.gz,M2_R2_trimmed_paired.fastq.gz,M3_R2_trimmed_paired.fastq.gz,M4_R2_trimmed_paired.fastq.gz,M5_R2_trimmed_paired.fastq.gz,M6_R2_trimmed_paired.fastq.gz,M7_R2_trimmed_paired.fastq.gz,M8_R2_trimmed_paired.fastq.gz,M9_R2_trimmed_paired.fastq.gz,M10_R2_trimmed_paired.fastq.gz,M11_R2_trimmed_paired.fastq.gz,M12_R2_trimmed_paired.fastq.gz,M13_R2_trimmed_paired.fastq.gz,M14_R2_trimmed_paired.fastq.gz,M15_R2_trimmed_paired.fastq.gz,M16_R2_trimmed_paired.fastq.gz,M17_R2_trimmed_paired.fastq.gz,M18_R2_trimmed_paired.fastq.gz,M19_R2_trimmed_paired.fastq.gz,M20_R2_trimmed_paired.fastq.gz,M21_R2_trimmed_paired.fastq.gz,M22_R2_trimmed_paired.fastq.gz,M23_R2_trimmed_paired.fastq.gz,S11c_R2_trimmed_paired.fastq.gz,S2c_R2_trimmed_paired.fastq.gz,S4a_R2_trimmed_paired.fastq.gz,S5c_R2_trimmed_paired.fastq.gz,S10_R2_trimmed_paired.fastq.gz,S12_R2_trimmed_paired.fastq.gz,S3a_R2_trimmed_paired.fastq.gz,S4b_R2_trimmed_paired.fastq.gz,S6a_R2_trimmed_paired.fastq.gz,S11a_R2_trimmed_paired.fastq.gz,S151617_R2_trimmed_paired.fastq.gz,S3b_R2_trimmed_paired.fastq.gz,S4c_R2_trimmed_paired.fastq.gz,S6b_R2_trimmed_paired.fastq.gz,S11b_R2_trimmed_paired.fastq.gz,S2a_R2_trimmed_paired.fastq.gz,S3c_R2_trimmed_paired.fastq.gz,S5a_R2_trimmed_paired.fastq.gz,S6c_R2_trimmed_paired.fastq.gz,G1a_R2_trimmed_paired.fastq.gz,G3a_R2_trimmed_paired.fastq.gz,G5a_R2_trimmed_paired.fastq.gz,G7a_R2_trimmed_paired.fastq.gz,G1b_R2_trimmed_paired.fastq.gz,G3b_R2_trimmed_paired.fastq.gz,G5b_R2_trimmed_paired.fastq.gz,G7b_R2_trimmed_paired.fastq.gz,G1c_R2_trimmed_paired.fastq.gz,G3c_R2_trimmed_paired.fastq.gz,G5c_R2_trimmed_paired.fastq.gz,G7c_R2_trimmed_paired.fastq.gz,G2a_R2_trimmed_paired.fastq.gz,G4a_R2_trimmed_paired.fastq.gz,G6a_R2_trimmed_paired.fastq.gz,G8a_R2_trimmed_paired.fastq.gz,G2b_R2_trimmed_paired.fastq.gz,G4b_R2_trimmed_paired.fastq.gz,G6b_R2_trimmed_paired.fastq.gz,G2c_R2_trimmed_paired.fastq.gz,G4c_R2_trimmed_paired.fastq.gz,G6c_R2_trimmed_paired.fastq.gz,G8c_R2_trimmed_paired.fastq.gz --k-min 27 --k-max 147 --k-step 8 -t 64 -m 0.95 -o Forefield_Assembly

#input file:
#*_R1_trimmed_paired.fastq.gz : forward Illumina read files
#*_R2_trimmed_paired.fastq.gz : reverse Illumina read files

#output file:
#Forefield_Assembly/k147.contigs.fa : assembly file
````

#### 2. Assembly Check

We checked the assembly with the script longMeta-summary. Thanks to this script we were also able to trim the assembly retrieving only sequences longer than 300 bp. This script also outputs a file reporting the length of each contig and will be used in next pipeline steps.

```bash
longMeta-summary --assembly-input Forefield_Assembly/final.contigs.fa --minimum-length 300 --assembly-output assembly_forefield.fasta --length-output assembly_forefield_length.txt


# |  File             |  Contig number   |  Total bp       |  Mean  |  Median  |  Max length   |  Min length |  n50  |  GC%  |
# |  k147.contigs.fa  |  41841465        |  30378371287    |  726   |  17208   |  561967       |  300        |  841  |  62   |

#input file:
#k147.contigs.fa : assembly input file (fasta or fastq format)

#output files:
#assembly_forefield.fasta: assembly output file, file containing only contigs longer than -min
#assembly_forefield_length.txt : length output file, tab-separated file reporting contigs and the correspondent lengths
```

#### 3. Diamond protein assignment

We used DIAMOND 0.9.22 (http://www.diamondsearch.org) to map known proteins to the assembly. You can download the database of known protein from the longMeta GitHUB page.

```bash

#Diamond format the nr.gz database
zcat nr.gz | diamond makedb -d nr

#Run Diamond
mkdir tmp
diamond blastx -d nr -q assembly_forefield.fasta -o assembly_forefield_nr.m8 -t tmp -e 0.000001 -F 15 --range-culling --range-cover 20 --id 50 --top 10 -f 6 -p 55 -c1 -b4.0
gzip assembly_forefield_nr.m8

#input files
#assembly_forefield.fasta : assembly file
#nr : Diamond-formatted database; it can be downloaded from the main GitHUB LongMeta page

#output file:
#assembly_forefield_nr.m8.gz : tabular blast file
```

#### 4. Assign taxonomy and coding region

We assigned taxonomy and coding regions to each contig with longMeta-assignment. This script screens the Diamond output file (tabular format). It then assigns non-overlapping Diamond proteins to the sequences indentifying in this way the coding regions. Taxonomy is inferred with a modified Lowest Common Ancestor algorithm (LCA).

```bash
longMeta-assignment --assembly assembly_forefield_nr.m8.gz --threads 3 --taxonomy h --acc2taxid accession_taxid.txt --taxid2taxon taxid_taxonomy.txt --gene-output assembly_forefield_nr_gene.txt --taxonomy-output assembly_forefield_nr_taxa.txt --tmp tmp > assembly_forefield_nr_assignment.log

##input files
#assembly_forefield_nr.m8.gz : Diamond file for assembly (blast tabular format)
#accession_taxid.txt : database file reporting taxids and taxonomy information
#taxid_taxonomy.txt : database file reporting accession numbers and taxids

##output files:
#assembly_forefield_nr_gene.txt : gene output file, Diamond file reporting only best non overlapping assigned proteins
#assembly_forefield_nr_taxa.txt : taxonomy output file, tab-separated file with taxonomy path associated to each contig

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

#input files:
#assembly_forefield_nr_gene.txt : gene input file
#assembly_forefield_nr_taxa.txt : taxonomy input file
#assembly_forefield.fasta : assembly input file
#assembly_forefield_length.txt : length input file
#accession_taxid.txt : database file reporting taxids and taxonomy information
#taxid_taxonomy.txt : database file reporting accession numbers and taxids

#output files:
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

#input files:
# assembly_forefield_chim.fasta : assembly file
#*_R1_trimmed_paired.fastq.gz : forward Illumina reads
#*_R2_trimmed_paired.fastq.gz : reverse Illumina reads

#output files:
#*.sam.gz : mapping files (sam format)

```

#### 7. Calculate taxonomy and gene abundance

We calculated taxonomy and gene coverage with the script longMeta-abundance. This script calculates base coverage for each gene and taxon present in the dataset.

```bash
mkdir AllForefield
longMeta-abundance --taxonomy-input assembly_forefield_nr_taxa_chim.txt --sam-input-tab listSam.txt --length-input assembly_forefield_length_chim.txt --output-folder AllForefield --alignment-AS 40 --gene-input assembly_forefield_nr_gene_chim.txt --acc2gene accession_protein.txt --acc2go accession_GO.txt

#Input files
#assembly_forefield_nr_taxa_chim.txt : taxonomy input file
#assembly_forefield_length_chim.txt : length input file
#assembly_forefield_nr_gene_chim.txt : gene input file
#listSam.txt : tab-sepated file with sam files in the first column, read length can be specified in second column

#accession_protein.txt : XX
#accession_GO.txt : XX

#Output files

#AllForefield/taxon_profile/genus.txt : base coverage at genus level
#AllForefield/taxon_profile/family.txt : base coverage at family level
#AllForefield/taxon_profile/order.txt: base coverage at order level
#AllForefield/taxon_profile/class.txt : base coverage at class level
#AllForefield/taxon_profile/phylum.txt : base coverage at phylum level
#AllForefield/taxon_profile/domain.txt : base coverage at domain level

#AllForefield/gene_profile/gene_to_go.txt : correspondences between Gene Ontology entries and gene names
#AllForefield/gene_profile/gene_total.txt : coding gene coverages

#AllForefield/gene_profile/genus/* : gene files at genus level
#AllForefield/gene_profile/family/* : gene files at family level
#AllForefield/gene_profile/order/* : gene files at order level
#AllForefield/gene_profile/class/* : gene files at class level
#AllForefield/gene_profile/phylum/* : gene files at phylum level
#AllForefield/gene_profile/domain/* : gene files at domain level

```

#### 8. Relative abundance for taxonomy data

Relative taxonomical abundances were calculated with the script longMeta-abundance where all the taxonomy relative abundances at high taxon level are assigned as the genus average coverage multiplied by the number of different genera and divided by the number of all the genera.

I could include the figure here

```bash
longMeta-relative --input-folder AllForefield --taxid2taxon taxid_taxonomy.txt


#input files

#AllForefield/taxon_profile/genus.txt : base coverage at genus level
#AllForefield/taxon_profile/family.txt : base coverage at family level
#AllForefield/taxon_profile/order.txt: base coverage at order level
#AllForefield/taxon_profile/class.txt : base coverage at class level
#AllForefield/taxon_profile/phylum.txt : base coverage at phylum level
#AllForefield/taxon_profile/domain.txt : base coverage at domain level

#output files

#AllForefield/taxon_profile/genus_relative.txt : relative base coverage at genus level
#AllForefield/taxon_profile/family_relative.txt : relative base coverage at family level
#AllForefield/taxon_profile/order_relative.txt: relative base coverage at order level
#AllForefield/taxon_profile/class_relative.txt : relative base coverage at class level
#AllForefield/taxon_profile/phylum_relative.txt : relative base coverage at phylum level
#AllForefield/taxon_profile/domain_relative.txt : relative base coverage at domain level

```

#### 9. Explore gene profiles

Explore the dataset from the gene and taxonomical point of view. In our study we wanted to explore nitrogenase and rock weathering genes at taxonomical level and therefore we investigated i) nitrogenase, ii) oxalate biosinthesis genes (obcA), iii) cyanide synthase genes and iv) siderophore-related genes.

```bash
#i) nitrogenase genes
#total
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene 'nitrogenase' --rank genus --output-file nitrogenase_genus.txt

#input files
#AllForefield/gene_profile/genus/* : gene files at genus level

#output files
#nitrogenase_genus.txt : output file
#nitrogenase_genus.txt.all

head nitrogenase_genus.txt


#genus level
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene 'nitrogenase' --output-file nitrogenase.txt

#input files
#AllForefield/gene_profile/gene_total.txt : coding gene coverages

#output files
#nitrogenase.txt : output file
#nitrogenase.txt.all

head nitrogenase.txt


#ii) oxalate byosynthesis genes
#total
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene '3-keto-5-aminohexanoate cleavage' --rank genus --output-file obc_genus.txt
#genus level
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene '3-keto-5-aminohexanoate cleavage' --output-file obc.txt

#iii) cyanide sinthase genes
#total
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene 'cyanide synthase' --rank genus --output-file cyanide_genus.txt
#genus level
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene 'cyanide synthase' --output-file cyanide.txt

#iv) siderophore related genes
#total
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene 'siderophore' --rank genus --output-file siderophore_genus.txt
#genus level
~/Script/toGITHUB/longMeta-explore.pl --input-folder AllForefield --gene 'siderophore' --output-file siderophore.txt

```
