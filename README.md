# LongMeta

This is a pipeline for the annotation of long reads such as assembly contigs. We suggest to use this pipeline only on shallow environmental assemblies where a gene prediction approach for functional annotation is not applicable. Some of this steps can also be applied for Nanopore read annotation.

Click [here](forefield_study.md) to see how the LongMeta pipeline was applied to environmental data.

#### Database

Different databases are used in this pipeline. All the database files can be downloaded from the following links. The files are big so make sure you have enough space on your system.

These are the formatted datasets for the NCBI-nr database v5. These datasets are updated every two months (if case the nr.gz was updated on the the NCBI website).

- [nr.gz](https://www.cerealsdb.uk.net/LongMeta/nr.gz) [67 Gb] : NCBI blast protein database.
- [accession_taxid.txt](https://www.cerealsdb.uk.net/LongMeta/accession_taxid.txt) [4.5 Gb] : tab-separated file reporting the protein accession numbers and the accociated taxids.
- [taxid_taxonomy.txt](https://www.cerealsdb.uk.net/LongMeta/taxid_taxonomy.txt) [239 Mb]: tab-separated file reporting the taxids and the associated taxonomical paths.
- [accession_protein.txt](https://www.cerealsdb.uk.net/LongMeta/accession_protein.txt) [9.9 Gb] : tab-separated file reporting the protein accession numbers and the associated protein names.
- [accession_GO.txt](https://www.cerealsdb.uk.net/LongMeta/accession_GO.txt) [1.6 Gb] : tab-separated file reporting the protein accession numbers and the associated Gene Ontology (GO) accession codes.
- [GO.txt](https://www.cerealsdb.uk.net/LongMeta/GO.txt) [3.6 Mb] : tab-separated file reporting the GO codes, the GO and the associated GO names and ontologies.

Here's a snapshot of the format of these profiles

```bash
zcat nr.gz | head
# >KJX92028.1 hypothetical protein TI39_contig5958g00003 [Zymoseptoria brevis]
# MAWTRQLVPLMLLFCGAHGLQRSSTATDQLSNSALQALGSHADLAAFVNDVEAVPEIANVILAHRGITIMAPVDSAWLRV
# DAIKRRNPAFLAWHIMNANVLTSDVPLVQYEQHPGITIPTFLSGSKNWTYSGEPASLISGGQSLTAITLKTEDNVIWVSG
# ASNVSYIKQANISYDRGIIHKIDPALQFPTSAYETAFAVGLYSYCWAVFTAGLDQEIRRIPNSTFLLPINEAFHAALPFL
# LGASREEFKRIVYRHVIPGRVLWSHEFYNASHETFEGSIVQIRGGNGRRWFVDDAMILDGSDKPLYNGVGHVVNKVLLPT
# >PYI97175.1 lysine 2,3-aminomutase [Verrucomicrobia bacterium]PYJ33862.1 lysine 2,3-aminomutase [Verrucomicrobia bacterium]
# MITPVSEEGNGKRFVSHAPGFWPQTPTELWNDWKWQLKNRVTSLAHLEQHLDLSDEERSGVLLSGDKLALAVTPHFFNLV
# PRNNPEDPIRRQVIPRIEETWTSPYDMADPCGEDSHMPVPGLVHRYPDRVLFLVTDRCASYCRYCTRSRVVSGVGEQELH
# TNFEEAFRYLQQHNEVRDVLLSGGDALIFSDDKIDKLLSRLRSIKHIEFVRIGTRVPIFLPQRITPDLCALLAKHHPLWM
# SVHVNHPRELTIEVKEALERLVNAGIPLGNQSVLLAGVNDDLETMKTLVHKLLMCRVRPYYIYQCDLINGSSHLRTSVAK

head accession_taxid.txt
# A0A023GPI8.1	232300
# A0A023GS28.1	37565
# A0A023GS29.1	37565
# A0A023IWD9.2	262245
# A0A023IWE0.1	580329
# A0A023IWE1.1	67723
# A0A023IWE2.1	1324310
# A0A023IWE3.1	580329
# A0A023IWG1.1	580329
# A0A023IWG2.1	67723

head taxid_taxonomy.txt
# 161256	Eukaryota	Arthropoda	Insecta	Coleoptera	Chrysomelidae	Galerucella	Galerucella lineola
# 2554519	Eukaryota	Arthropoda	Insecta	Lepidoptera	Erebidae	Ostha	Ostha sp. INCTA270-10
# 871351	Bacteria	Firmicutes	Bacilli	Lactobacillales	Lactobacillaceae	Lactobacillus	Lactobacillus sp. MF44(2010)
# 455508	Unclassified	Unclassified	Unclassified	Unclassified	Unclassified	Unclassified	Unclassified
# 1060012	Eukaryota	Arthropoda	Insecta	Lepidoptera	Lepidoptera_Unclassified	Lepidoptera_Unclassified	Lepidoptera sp. BOLD:AAV9997
# 314295	Eukaryota	Chordata	Mammalia	Primates	Primates_Unclassified	Primates_Unclassified	Primates_Unclassified
# 43255	Eukaryota	Chordata	Actinopteri	Perciformes	Harpagiferidae	Harpagifer	Harpagifer_Unclassified
# 1288524	Bacteria	Actinobacteria	Actinobacteria	Micrococcales	Microbacteriaceae	Curtobacterium	Curtobacterium sp. I12A-00125
# 1318028	Viruses	Negarnaviricota	Insthoviricetes	Articulavirales	Orthomyxoviridae	Alphainfluenzavirus	Influenza A virus
# 2466516	Eukaryota	Arthropoda	Insecta	Hymenoptera	Diapriidae	Diapriidae_Unclassified	Diapriidae sp. BIOUG10490-E09

head accession_protein.txt
#WP_003522443.1	hypothetical protein H009_14803
#XP_022177924.1	uncharacterized protein LOC111038960 isoform X3
#WP_145841639.1	hypothetical protein B9N43_07350
#WP_129747056.1	hypothetical protein NU08_2170
#PRQ37177.1	hypothetical protein RchiOBHm_Chr4g0399661
#TMS39860.1	hypothetical protein L596_006323
# QHN97289.1	hypothetical protein Ahy_B08g091168
# WP_003241080.1	hypothetical protein CVV77_19315
# WP_005204715.1	hypothetical protein F902_03102
# ROV92053.1	hypothetical protein VPNG_09832

head accession_go.txt
# WP_046691645.1	GO:0016021,GO:0006355,GO:0016020,GO:0003677,
# XP_007388560.1	GO:0022900,GO:0016020,GO:0004129,GO:0016021,GO:1902600,
# XP_015192253.1	GO:0016021,GO:0004672,GO:0005524,GO:0016020,GO:0006468,
# WP_048329662.1	GO:0051060,GO:0005975,GO:0030246,GO:0004553,
# WP_004397554.1	GO:0000155,GO:0016740,GO:0016772,GO:0016310,GO:0016301,GO:0000160,GO:0004673,GO:0023014,GO:0007165,GO:0018106,
# WP_006036850.1	GO:0016874,GO:0003824,
# WP_056333144.1	GO:0016740,GO:0003824,GO:0016757,GO:0009435,GO:0016763,GO:0004514,
# WP_056419841.1	GO:0005886,GO:0016021,GO:0016020,GO:0055085,GO:0022857,
# XP_010741415.1	GO:0016020,GO:0016021,GO:0055114,GO:0047057,GO:0042373,
# WP_015963251.1	GO:0016020,GO:0007049,GO:0032153,GO:0005886,GO:0016021,GO:0043093,GO:0005887,GO:0051301,

head go.txt
# GO:0000001	mitochondrion inheritance	biological_process
# GO:0000002	mitochondrial genome maintenance	biological_process
# GO:0000003	reproduction	biological_process
# GO:0019952	reproduction	biological_process
# GO:0050876	reproduction	biological_process
# GO:0000005	obsolete ribosomal chaperone activity	molecular_function
# GO:0000006	high-affinity zinc transmembrane transporter activity	molecular_function
# GO:0000007	low-affinity zinc ion transmembrane transporter activity	molecular_function
# GO:0000008	obsolete thioredoxin	molecular_function
# GO:0000013	obsolete thioredoxin	molecular_function
```

#### Getting started

1. Download the folder longMeta
```bash
wget XXXlink to gitHUB ##add gitHUB link
```
2. Give the correct permissions to the scripts
```bash
cd GitHUB
chmod +x longMeta/* #check if they come with correct permissions already
```
3. Move the scripts to your bin or set it in your PATH
```bash
export PATH={path-to-longMeta}/longMeta:$PATH
```

#### Detailed option command line

LongMeta consists in 5 different scripts:

1. longMeta-summary
2. longMeta-assignment
3. longMeta-chimera
4. longMeta-abundance
5. longMeta-explore
6. longMeta-relative


1. Assembly summary : longMeta-summary

The script longMeta-summary is used to check and trim the assembly (if needed).

```bash
longMeta-summary [--help] [--sequence-input INPUT_FILE] [--type {fasta, fastq}] [--minimum-length POS_INTEGER] [--length-output OUTPUT_FILE] [--sequence-output OUTPUT_FILE]

-h, --help
	show help message

--sequence-input, -s INPUT_FILE
	assembly file

OPTIONAL

--type {fasta, fastq}
	format of input file. default: fasta
--minimum-length, -min POS_INTEGER
	minimum length to consider. default: 0
--length-output OUTPUT_FILE
	file reporting contigs and the correspondent lengths, only for contigs longer than -min
--sequence-output OUTPUT_FILE
	file reporting contigs, only for contigs longer than -min
```

2. Assignment : longMeta-assignment

The script longMeta-assignment assign taxonomy and functionality to each contig. This script works with Diamond which can be found XXXX. It outputs only Diamond matches that have been identified as the best matches and a file with contig and associated taxonomy.

```bash
longMeta-assignment [--help] [--tmp TEMPORARY_FOLDER] [--threads POS_INTEGER] [--assembly INPUT_FILE] [--paired-read INPUT_FILE] [--unpaired-read INPUT_FILE] [--gene-output FILE_OUTPUT] [--taxonomy-output FILE_OUTPUT] [--taxonomy {h, b, l}] [--taxid2taxon DATABASE_FILE] [--acc2taxid DATABASE_FILE] [--max-equal POS_INTEGER] [--cutoff-ID-best POS_INTEGER] [--min-best POS_INTEGER] [--all-lca {mixed, gene, all}] [--min-lca POS_INTEGER] [--cutoff-ID-lca POS_INTEGER]

--help, -h
	show help message
--tmp TEMPORARY_FOLDER
	temporary folder
--threads POS_INTEGER
	number of threads

--assembly, -a INPUT_FILE
	Diamond file for assembly
--paired-read, -p INPUT_FILE
	Diamond file for paired reads (comma separated files)
--unpaired-read, -u INPUT_FILE
	Diamond file for unpaired reads (comma separated files)

--gene-output FILE_OUTPUT
	output gene file

--taxonomy-output FILE_OUTPUT
	output taxonomy file
--taxonomy {h, b, l}
	pipeline for the taxonomy assignment. h stands for hybrid, l for lca and b for best. default: h
--taxid2taxon DATABASE_FILE
	database file reporting taxids and taxonomy information
--acc2taxid DATABASE_FILE
	database file reporting accession numbers and taxids

Options for taxonomy assignment:
--max-equal POS_INTEGER
	maximum number of sequences to be used for taxonomical classification when sequences with the same Diamond identity score and bit-score were assigned to the same sequence area. If not specified, LongMeta uses all the sequences
--cutoff-ID-best POS_INTEGER
	identity score value to consider a Diamond alignment for the BEST algorithm. default: 80
--min-best POS_INTEGER	minimum number of sequences assigned to a sequence to run the BEST algorithm. default: 3
--all-lca {all, gene, mixed}
	LCA algorithm can use all the sequences that Diamond assigned to a sequence (all), only the sequences that were assigned as genes (gene) or a mixed approach (mixed). In the latter, all the sequences are used only if less than a n number of sequences is present. N can be set with the option --min-lca. Default: mixed
--min-lca POS_INTEGER
	minimum number of sequences needed to run the LCA algorithm on only gene-assigned matches. default: 3. (to use with --all-lca mixed)
--cutoff-ID-lca POS_INTEGER
	minimum identity score to give more weight to a certain Diamond alignment. default: 80

Options for gene assignment:
--overlap POS_INTEGER
	maximum number of overlapping bases between gene coding regions

```

3. Chimera detection (optional) : longMeta-chimera

The script longMeta-chimera screens for chimeric contigs.

```bash
longMeta-chimera [--help] [--gene-input INPUT_FILE] [--taxonomy-input INPUT_FILE] [--sequence-input INPUT_FILE] [--length-input INPUT_FILE] [--gene-output OUTPUT_FILE] [--taxonomy-output OUTPUT_FILE] [--sequence-output OUTPUT_FILE] [--length-output OUTPUT_FILE] [--taxid2taxon DATABASE_FILE] [--acc2taxid DATABASE_FILE] [--taxon-rank {domain, phylum, class, order, family, genus, species}] [--ID-limit POS_INTEGER] [--cluster-limit POS_INTEGER]

-h, --help
	show help message

--gene-input, -g INPUT_FILE
	gene file
--taxonomy-input, -t INPUT_FILE
	taxonomy file
--taxid2taxon DATABASE_FILE
	database file reporting taxids and taxonomy information
--acc2taxid DATABASE_FILE
	database file reporting accession numbers and taxids

Optional:

--sequence-input, -s INPUT_FILE
	contig file
--length-input, -l INPUT_FILE
	length file
--gene-output, -go OUTPUT_FILE
	output gene file
--taxonomy-output, -to OUTPUT_FILE
	output taxonomy file
--sequence-output, -so OUTPUT_FILE
	output contig file
--length-output, -lo OUTPUT_FILE
	output length file
--taxon-rank {domain, phylum, class, order, family, genus, species}
	taxonomical level for chimera finding. default: genus
--ID-limit POS_INTEGER
	minimum identity score to consider a gene coding region . default: 80
--cluster-limit POS_INTEGER
	number of consecutive gene coding region assigned to the same taxon to divide the contig. default: 2
```

4. Abundance calculation : longMeta-abundance

The script longMeta-abundance works with output files from longMeta-assignment and reads aligned to the assembly. We used bwahaha or bowtie2 to align Illumina reads back to assembly.

```bash
longMeta-abundance [--help] [--taxonomy-input INPUT_FILE] [--sam-input INPUT_FILE] [--length-input INPUT_FILE] [--output-folder OUTPUT_FOLDER] [--ignore_uncl {yes,no}] [--perc-limit POS_NUMBER] [--phyla-exclusion {yes,no}] [--gene-input INPUT_FILE] [--acc2gene DATABASE_FILE] [--acc2go DATABASE_FILE] [--read-length POS_INTEGER] [--alignment-AS NUMBER]

--help, -h
	show help message

--taxonomy-input, -t INPUT_FILE
	taxonomy file
--sam-input, -s INPUT_FILE
	sam files (comma-separated)
--sam-input-tab, -st INPUT_FILE
	text file reporting sam files in first column, read length can also be specified in second column
--length-input, -l INPUT_FILE
	length file
--output-folder, -o OUTPUT_FOLDER
	output folder

Perform gene profiling:

--gene-input, -g INPUT_FILE
	gene file
--acc2gene DATABASE_FILE
	accession to gene database
--acc2go DATABASE_FILE
	accession to gene ontology database

Optional:

--ignore-uncl {yes,no}
	consider contigs that did not get any matches from Diamond. default: no
--phyla-exclusion {yes,no}
	exclude contigs assigned to Metazoa and Streptophyta. default: yes
--perc-limit POS_NUMBER
	sensitivity threshold, percentage of bases and contigs that need to be assigned to a genus to be considered. default: 0.1
--read-length POS_INTEGER
	Length of the mapped reads. default: 150
--alignment-AS NUMBER
	minimum alignment score to consider a match (AS)
```


5. Explore data : longMeta-explore

The script longMeta-explore uses the folder data with all longMeta data and allows an easy exploration of the data

```bash
longMeta-explore [--help] [--input-folder INPUT_FOLDER] [--gene STRING] [--go {go-numbers}] [--go-name STRING] [--go-ontology {molecular_function, biological_process, cellular_component}] [--taxon STRING] [--taxon-rank {domain, phylum, class, order, family, genus}] [--output-file OUTPUT_FILE] [--go2def DATABASE] [--output-type {gene,GO}]

--help, -h
	show help message

--input-folder, -in INPUT_FOLDER
	longMeta folder created in the previous steps of the pipeline
--gene, -g STRING
	partial or entire gene name to use for gene retrieval. If the name is made up from more than one word, they must be enclosed in quotation marks. Special characters (e.g. ` or /) must be escaped with ‘\’. Write ‘all’ if all the genes must be retrieved.
--go-accession, -ga GO_CODE
	gene ontology code to use for gene retrival
--go-name, -gn STRING
	gene ontology name (or part) to use for gene retrival

--go-ontology, -go {molecular_function, biological_process, cellular_component}
	gene ontology category
--go2def, DATABASE
	gene ontology database
--taxon, -t STRING
	taxa (comma-separated)
--taxon-rank, -r {domain, phylum, class, order, family, genus}
	rank of --taxon
--output-file, -out OUTPUT_FILE
	output file reporting the genes of interest
--output-type {gene,GO}
	output type
```

6. Calculate relative abundance : longMeta-relative

The script longMeta-relative gives the taxonomic relative abundance

```bash

longMeta-relative [--help] [--input-folder INPUT_FOLDER] [--taxid2taxon DATABASE_FILE]

--help, -h
	show help message

--input_folder, -in INPUT_FOLDER
	longMeta folder created in the previous steps of the pipeline
--taxid2taxon DATABASE_FILE
	database file reporting taxids and taxonomy information

```
