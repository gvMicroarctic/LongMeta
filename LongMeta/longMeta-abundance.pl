#!/usr/bin/perl

use strict;
use POSIX;

#Assign taxonomy and functionality to sam files.

#usage: longMeta-abundance [--help] [--taxonomy-input INPUT_FILE] [--sam-input INPUT_FILE] [--length-input INPUT_FILE] [--output-folder OUTPUT_FOLDER] [--ignore_uncl {yes,no}] [--perc-limit POS_NUMBER] [--phyla-exclusion {yes,no}] [--gene-input INPUT_FILE] [--acc2gene DATABASE_FILE] [--acc2go DATABASE_FILE] [--read-length POS_INTEGER] [--alignment-AS NUMBER]

#help message
if (($ARGV[0] eq '--help') or ($ARGV[0] eq '-h')) {
    print "usage: longMeta-abundance [--help] [--taxonomy-input INPUT_FILE] [--sam-input INPUT_FILE] [--length-input INPUT_FILE] [--output-folder OUTPUT_FOLDER] [--ignore_uncl {yes,no}] [--perc-limit POS_NUMBER] [--phyla-exclusion {yes,no}] [--gene-input INPUT_FILE] [--acc2gene DATABASE_FILE] [--acc2go DATABASE_FILE] [--read-length POS_INTEGER] [--alignment-AS NUMBER]\n\n";
    
    print "--help, -h\n\tshow help message\n\n";
    
    print "--taxonomy-input, -t INPUT_FILE\n\ttaxonomy file\n";
    print "--sam-input, -s INPUT_FILE\n\tsam files (comma-separated)\n";
    print "--sam-input-tab, -st INPUT_FILE\n\ttext file reporting sam files in first column, read length can also be specified in second column\n";
    print "--length-input, -l INPUT_FILE\n\tlength file\n";
    print "--output-folder, -o OUTPUT_FOLDER\n\toutput folder\n\n";
    
    print "Perform gene profiling:\n\n";
    
    print "--gene-input, -g INPUT_FILE\n\tgene file\n";
    print "--acc2gene DATABASE_FILE\n\taccession to gene database\n";
    print "--acc2go DATABASE_FILE\n\taccession to gene ontology database\n\n";
    
    print "Optional:\n\n";
    
    print "--ignore-uncl {yes,no}\n\tconsider contigs that did not get any matches from Diamond. default: no\n";
    print "--phyla-exclusion {yes,no}\n\texclude contigs assigned to Metazoa and Streptophyta. default: yes\n";
    print "--perc-limit POS_NUMBER\n\tsensitivity threshold, percentage of bases and contigs that need to be assigned to a genus to be considered. default: 0.1\n";
    print "--read-length POS_INTEGER\n\tLength of the mapped reads. default: 150\n";
    print "--alignment-AS NUMBER\n\tminimum alignment score to consider a match (AS)\n\n";
    
    exit;
}

my %args = @ARGV;

#check options
my %allOption;
$allOption{'--taxonomy-input'} = ''; $allOption{'-t'} = '';
$allOption{'--sam-input'} = ''; $allOption{'-s'} = '';
$allOption{'--sam-input-tab'} = ''; $allOption{'-st'} = '';
$allOption{'--length-input'} = ''; $allOption{'-l'} = '';
$allOption{'--output-folder'} = ''; $allOption{'-o'} = '';

$allOption{'--gene-input'} = ''; $allOption{'-g'} = '';
$allOption{'--acc2gene'} = '';
$allOption{'--acc2go'} = '';

$allOption{'--ignore-uncl'} = '';
$allOption{'--phyla-exclusion'} = '';
$allOption{'--perc-limit'} = '';
$allOption{'--read-length'} = '';
$allOption{'--alignment-AS'} = '';

$allOption{'--help'} = ''; $allOption{'-h'} = '';

foreach my $a (keys %args) {
    if (!(defined($allOption{$a}))) {
        print "Error: Invalid command option: $a. To print help message: longMeta-abundance --help\n";
        exit;
    }
}

#import options

#output folder
my $out = $args{'--output-folder'};
if ($out eq "") {
    $out = $args{'-o'};
    if ($out eq "") {
        print "Error: --output-folder must be specified. To print help message: longMeta-abundance --help\n";
        exit;
    }
}
if (!(-e $out and -d $out)) {
    print "Error: longMeta-abundance cannot open $out. To print help message: longMeta-abundance --help\n";
    exit;
}

#taxonomy input file
my $inLCA = $args{'--taxonomy-input'};
if ($inLCA eq "") {
    $inLCA = $args{'-t'};
    if ($inLCA eq "") {
        print "Error: --taxonomy-input must be specified. To print help message: longMeta-abundance --help\n";
        exit;
    }
}
if ($inLCA =~ /.gz$/) {
    if (-e $inLCA) {
        open(LCA, "gunzip -c $inLCA |");
    } else {
        print "Error: longMeta-abundance cannot open $inLCA. To print help message: longMeta-abundance --help\n";
        exit;
    }
} else {
    if (!(open(LCA, "<$inLCA"))) {
        print "Error: longMeta-abundance cannot open $inLCA. To print help message: longMeta-abundance --help\n";
        exit;
    }
}

`mkdir ./${out}/taxon_profile/ 2> /dev/null`;

#gene file and databases
my $geneClassification = 'no';
my $goClassification = 'no';
my $inGENE = $args{'--gene-input'}; #gene input file
if ($inGENE eq "") {
    $inGENE = $args{'-g'};
}
if ($inGENE ne "") {
    if ($inGENE =~ /.gz$/) {
        if (-e $inGENE) {
            open(CUT, "gunzip -c $inGENE |");
        } else {
            print "Error: longMeta-abundance cannot open $inGENE. To print help message: longMeta-abundance --help\n";
            exit;
        }
    } else {
        if (!(open(CUT, "<$inGENE"))) {
            print "Error: longMeta-abundance cannot open $inGENE. To print help message: longMeta-abundance --help\n";
            exit;
        }
    }
    #gene database
    my $acc2gene = $args{'--acc2gene'}; #gene input file
    if ($acc2gene eq "") {
        print "Error: --acc2gene must be specified. To print help message: longMeta-abundance --help\n";
        exit;
    }
    
    if ($acc2gene =~ /.gz$/) {
        if (-e $inGENE) {
            open(GENE, "gunzip -c $acc2gene |");
        } else {
            print "Error: longMeta-abundance cannot open $acc2gene. To print help message: longMeta-abundance --help\n";
            exit;
        }
    } else {
        if (!(open(GENE, "<$acc2gene"))) {
            print "Error: longMeta-abundance cannot open $acc2gene. To print help message: longMeta-abundance --help\n";
            exit;
        }
    }
    `mkdir ./${out}/gene_profile/`;
    `mkdir ./${out}/gene_profile/domain`;
    `mkdir ./${out}/gene_profile/phylum`;
    `mkdir ./${out}/gene_profile/class`;
    `mkdir ./${out}/gene_profile/order`;
    `mkdir ./${out}/gene_profile/family`;
    `mkdir ./${out}/gene_profile/genus`;
    $geneClassification = 'yes';
    
    #GO database
    my $acc2go = $args{'--acc2go'}; #gene input file
    if ($acc2go ne "") {
        if ($acc2go =~ /.gz$/) {
            if (-e $acc2go) {
                open(GO, "gunzip -c $acc2go |");
            } else {
                print "Error: longMeta-abundance cannot open $acc2go. To print help message: longMeta-abundance --help\n";
                exit;
            }
        } else {
            if (!(open(GO, "<$acc2go"))) {
                print "Error: longMeta-abundance cannot open $acc2go. To print help message: longMeta-abundance --help\n";
                exit;
            }
        }
        $goClassification = 'yes';
    }
}

#alignment score
my $alignment_check = $args{'--alignment-AS'}; #alignment check
if ($alignment_check ne "") {
    if ($alignment_check !~ /^-?\d+$/) {
        print "Error: --alignment-AS must be an integer. To print help message: longMeta-abundance --help\n";
        exit;
    }
}

#sam files
my @sam = split(/,/, $args{'--sam-input'});
my %read_length;
if ($sam[0] eq "") {
    @sam = split(/,/, $args{'-s'});
    if ($sam[0] eq "") {
        #open text file
        my $text = $args{'--sam-input-tab'}; #text file
        if ($text eq '') {
            $text = $args{'-st'}; #text file
        }
        if ($text eq '') {
            print "Error: either --sam-input or --sam-input-tab must be specified. To print help message: longMeta-abundance --help\n";
            exit;
        } else {
            if (!(open(TEXT, "<$text"))) {
                print "Error: longMeta-abundance cannot open $text. To print help message: longMeta-abundance --help\n";
                exit;
            } else {
                while(defined(my $input = <TEXT>)) {
                    chomp($input);
                    if ($input =~ /\t/) {
                        my @info = split(/\t/, $input);
                        push @sam, $info[0];
                        if ($info[1] !~ /^?\d+$/) {
                            print "Error: there was an error in $text, the second file column must specify the read length of the sam file. To print help message: longMeta-abundance --help\n";
                        } else {
                            $read_length{$info[0]} = $info[1];
                        }
                    } else {
                        push @sam, $input;
                    }
                }
                close(TEXT);
            }
        }
    }
}

foreach my $sample (@sam) {
    if ($sample =~ /.gz$/) {
        if (-e $sample) {
            open(SAM, "gunzip -c $sample |");
        } else {
            print "Error: longMeta-abundance cannot open $sample. To print help message: longMeta-abundance --help\n";
            exit;
        }
    } else {
        if (!(open(SAM, "<$sample"))) {
            print "Error: longMeta-abundance cannot open $sample. To print help message: longMeta-abundance --help\n";
            exit;
        }
    }
    close(SAM);
}

#read length
my $read_lengthTmp = $args{'--read-length'};
if ($read_lengthTmp eq '') {
    my $defNum = keys %read_length;
    if ($defNum == 0) {
        foreach my $s (@sam) {
            $read_length{$s} = 150;
        }
    }
} elsif ($read_lengthTmp !~ /^\d+$/) {
    print "Error: --read-length must be a positive integer. To print help message: longMeta-abundance --help\n";
    exit;
} else {
    foreach my $s (@sam) {
        $read_length{$s} = $read_lengthTmp;
    }
}

#length file
my $inLEN = $args{'--length-input'};
if ($inLEN eq "") {
    $inLEN = $args{'-l'};
    if ($inLEN eq "") {
        print "Error: --length-input must be specified. To print help message: longMeta-abundance --help\n";
        exit;
    }
}
if ($inLEN =~ /.gz$/) {
    if (-e $inLEN) {
        open(LEN, "gunzip -c $inLEN |");
    } else {
        print "Error: longMeta-abundance cannot open $inLEN. To print help message: longMeta-abundance --help\n";
        exit;
    }
} else {
    if (!(open(LEN, "<$inLEN"))) {
        print "Error: longMeta-abundance cannot open $inLEN. To print help message: longMeta-abundance --help\n";
        exit;
    }
}

my $phylaExclusion = $args{'--phyla-exclusion'};
if ($phylaExclusion eq "") {
    $phylaExclusion = "yes";
} elsif (($phylaExclusion ne "yes") && ($phylaExclusion ne "no") ) {
    print "Error: --phyla-exclusion must be either yes or no. To print help message: longMeta-abundance --help\n";
    exit;
}

my %trimPhyla;
if ($phylaExclusion eq "yes") {
    $trimPhyla{'Streptophyta'} = '';
    $trimPhyla{'Acanthocephala'} = '';
    $trimPhyla{'Annelida'} = '';
    $trimPhyla{'Arthropoda'} = '';
    $trimPhyla{'Brachiopoda'} = '';
    $trimPhyla{'Bryozoa'} = '';
    $trimPhyla{'Chaetognatha'} = '';
    $trimPhyla{'Chordata'} = '';
    $trimPhyla{'Cnidaria'} = '';
    $trimPhyla{'Ctenophora'} = '';
    $trimPhyla{'Cycliophora'} = '';
    $trimPhyla{'Echinodermata'} = '';
    $trimPhyla{'Entoprocta'} = '';
    $trimPhyla{'Gastrotricha'} = '';
    $trimPhyla{'Gnathostomulida'} = '';
    $trimPhyla{'Hemichordata'} = '';
    $trimPhyla{'Kinorhyncha'} = '';
    $trimPhyla{'Loricifera'} = '';
    $trimPhyla{'Mollusca'} = '';
    $trimPhyla{'Nematoda'} = '';
    $trimPhyla{'Nematomorpha'} = '';
    $trimPhyla{'Nemertea'} = '';
    $trimPhyla{'Onychophora'} = '';
    $trimPhyla{'Orthonectida'} = '';
    $trimPhyla{'Placozoa'} = '';
    $trimPhyla{'Platyhelminthes'} = '';
    $trimPhyla{'Porifera'} = '';
    $trimPhyla{'Priapulida'} = '';
    $trimPhyla{'Rhombozoa'} = '';
    $trimPhyla{'Rotifera'} = '';
    $trimPhyla{'Tardigrada'} = '';
    $trimPhyla{'Xenacoelomorpha'} = '';
} else {
    $trimPhyla{'noExclusion'} = '';
}

my @rank = ('domain', 'phylum', 'class', 'order', 'family','genus');

my $datestring;

my $unclassified = $args{'--ignore-uncl'};
if ($unclassified eq "") {
    $unclassified = "no";
} elsif (($unclassified ne "yes") && ($unclassified ne "no") ) {
    print "Error: --ignore-uncl must be either yes or no. To print help message: longMeta-abundance --help\n";
    exit;
}

my $perc_limit = $args{'--perc-limit'};
if ($perc_limit eq "") {
    $perc_limit = 0.1;
} elsif ($perc_limit !~ /^[\d\.]+$/) {
    print "Error: --perc-limit must be a positive number. To print help message: longMeta-abundance --help\n";
    exit;
}

my %assigned;
my %offs;
my %gephy;
my %countScale;
my $countScaleTot = 0;

while (defined(my $input = <LCA>)) {
    chomp($input);
    my @taxa = split(/\t/, $input);
    if (!(defined($trimPhyla{$taxa[2]}))) { #only if not Metazoa or Streptophyta (by default)
        $assigned{$taxa[0]} ='';
        $countScale{$taxa[6]}++;
        $countScaleTot++;
        $offs{$taxa[0]}{'domain'} = $taxa[1]; #hash with contig - rank = taxon
        $offs{$taxa[0]}{'phylum'} = $taxa[2];
        $offs{$taxa[0]}{'class'} = $taxa[3];
        $offs{$taxa[0]}{'order'} = $taxa[4];
        $offs{$taxa[0]}{'family'} = $taxa[5];
        $offs{$taxa[0]}{'genus'} = $taxa[6];
        $gephy{$taxa[6]}{'domain'} = $taxa[1];
        $gephy{$taxa[6]}{'phylum'} = $taxa[2];
        $gephy{$taxa[6]}{'class'} = $taxa[3];
        $gephy{$taxa[6]}{'order'} = $taxa[4];
        $gephy{$taxa[6]}{'family'} = $taxa[5];
        $gephy{$taxa[6]}{'genus'} = $taxa[6];
    } else {
        $assigned{$taxa[0]} ='';
        $countScale{$taxa[6]}++;
        $countScaleTot++;
        $offs{$taxa[0]}{'domain'} = 'Unclassified';
        $offs{$taxa[0]}{'phylum'} = 'Unclassified';
        $offs{$taxa[0]}{'class'} = 'Unclassified';
        $offs{$taxa[0]}{'order'} = 'Unclassified';
        $offs{$taxa[0]}{'family'} = 'Unclassified';
        $offs{$taxa[0]}{'genus'} = 'Unclassified';
        $gephy{'Unclassified'}{'domain'} = 'Unclassified';
        $gephy{'Unclassified'}{'phylum'} = 'Unclassified';
        $gephy{'Unclassified'}{'class'} = 'Unclassified';
        $gephy{'Unclassified'}{'order'} = 'Unclassified';
        $gephy{'Unclassified'}{'family'} = 'Unclassified';
        $gephy{'Unclassified'}{'genus'} = 'Unclassified';
    }
}
close(LCA);

my %len;
my %lenScale;
my $lenScaleTot = 0;
while(defined(my $input=<LEN>)) {
    chomp($input);
    my $contig;
    my @info = split(/\t/, $input);
    if ($info[0] =~ /\s/) { #if contig header contains spaces
        my @infoC = split(/\s/, $info[0]);
        $contig = $infoC[0];
    } else {
        $contig = $info[0];
    }
    if (defined($offs{$contig})) {
        $lenScale{$offs{$contig}{'genus'}} += $info[1];
        $lenScaleTot += $info[1];
    }
    $len{$contig} = $info[1]; #contig = len
}
close(LEN);

if ($unclassified eq "yes") { #count Unclassified component
    foreach my $contig (keys %len) {
        if (!(defined($assigned{$contig}))) {
            $countScale{'Unclassified'}++;
            $offs{$contig}{'domain'} = 'Unclassified';
            $offs{$contig}{'phylum'} = 'Unclassified';
            $offs{$contig}{'class'} = 'Unclassified';
            $offs{$contig}{'order'} = 'Unclassified';
            $offs{$contig}{'family'} = 'Unclassified';
            $offs{$contig}{'genus'} = 'Unclassified';
            $gephy{'Unclassified'}{'domain'} = 'Unclassified';
            $gephy{'Unclassified'}{'phylum'} = 'Unclassified';
            $gephy{'Unclassified'}{'class'} = 'Unclassified';
            $gephy{'Unclassified'}{'order'} = 'Unclassified';
            $gephy{'Unclassified'}{'family'} = 'Unclassified';
            $gephy{'Unclassified'}{'genus'} = 'Unclassified';
            $lenScale{'Unclassified'} += $len{$contig};
            $lenScaleTot += $len{$contig};
        }
    }
}

undef %assigned;

my %sumLenTrim;
my %genusExcluded;
foreach my $genus (sort {$lenScale{$b} <=> $lenScale{$a}} keys %lenScale) {
    my $divCount = ($countScale{$genus}/$countScaleTot)*100;
    my $divLen = ($lenScale{$genus}/$lenScaleTot)*100;
    
    if (($divLen >= $perc_limit) or ($divCount >= $perc_limit)) {
        $sumLenTrim{$genus} = $lenScale{$genus}; #trimmed version of sumLen
        
    } else {
        $genusExcluded{$genus} = '';
    }
}

my %def;
my %defLen;
my %nr_protein;
my %nr_go;

if ($geneClassification eq "yes") {
    #retrieve the accession number and start/end positions on the contig for only the pre-selected contigs. So the main key is the accession number
    my %allAcc;
    while (defined(my $input = <CUT>)) {
        chomp($input);
        my @matches = split /\t/, $input;
        if (defined $offs{$matches[0]}) {
            $allAcc{$matches[1]} = '';  #protein accession number = ''
            if ($matches[6] < $matches[7]) { #$matches[6] and $matches[7] are the initial and final bp of the contig with the gene match
                $def{$matches[0]}{$matches[1]}{'S'} = $matches[6];  #contig - protein accession number - S|E = start position | end position
                $def{$matches[0]}{$matches[1]}{'E'} = $matches[7];
            } else {
                $def{$matches[0]}{$matches[1]}{'S'} = $matches[7];
                $def{$matches[0]}{$matches[1]}{'E'} = $matches[6];
            }
        }
    }
    close(CUT);
    
    my $cutNum = keys %def;
    
    print "Uploaded $cutNum contigs from gene file.\n";
    
    #Index file with accession number and proteins.
    #create hash with contig as main hash, then name of the protein and star/end bo of the contig with gene match.
    while (defined(my $input = <GENE>)) {
        chomp($input);
        my ($acc, $prot) = split /\t/, $input;
        if (defined($allAcc{$acc})) {
            $nr_protein{$acc} = $prot; #accession numb = prot name
        }
    }
    close(GENE);
    
    #gene length on each contig for gene abundance calculation
    foreach my $contig (keys %def) {
        foreach my $gene (keys %{$def{$contig}}) {
            my $diff = $def{$contig}{$gene}{'E'} - $def{$contig}{$gene}{'S'} + 1;
            $defLen{$contig}{$nr_protein{$gene}}+= $diff;
        }
    }
    
    if ($goClassification eq 'yes') {
        while (defined(my $input = <GO>)) {
            chomp($input);
            my ($acc, $go) = split /\t/, $input;
            chop($go);
            if (defined($allAcc{$acc})) {
                my @allGO = split(/,/, $go);
                foreach my $g (@allGO) {
                    $nr_go{$acc}{$g} = '';
                }
            }
        }
    }
}

my %finalTaxa;
my %finalGene;

my %final_combined;
my %final_tot;
my %final_go;

foreach my $sample (@sam) {
    
    my $datestring = localtime();
    print "Starting on $sample: $datestring\n";
    
    if ($sample =~ /.gz$/) {
        open(SAM, "gunzip -c $sample |");
    } else {
        open(SAM, "<$sample");
    }
    
    my %taxa;
    my %sumRead;
    my %finalProv;
    
    my $len = $read_length{$sample};
    
    my $div1;
    my $div2;
    
    if ($len == 150) {
        $div1 = 25;
        $div2 = 50;
    } elsif ($len == 100) {
        $div1 = 15;
        $div2 = 35;
    } elsif ($len == 300) {
        $div1 = 50;
        $div2 = 100;
    } else {
        $div1 = ceil($len/3);
        $div2 = ceil($len/6);
    }

    my $div3 = $div2*2;
    
    while (defined(my $input = <SAM>)) { #read the sam file
        my ($read,$flag,$contig,$start_al,$rest) = (split /\t/, $input,5);
        
        $taxa{$contig} += $read_length{$sample}; #number of bases per contig
        
        if ($geneClassification eq 'yes') {
            if (defined($def{$contig})) { #if there are coding regions in the contig
                my $AS;
                my $pos1 = $start_al + $div2;
                my $pos2 = $start_al + $div3;
                my $pos3 = $start_al + $len;
                if ($alignment_check ne "") {
                    ($AS) = $rest =~ /AS:i:(.*?)\t/;
                    if ($AS >= $alignment_check) {
                        $finalGene{$contig}{$start_al} += $div1;
                        $finalGene{$contig}{$pos1} += $div2;
                        $finalGene{$contig}{$pos2} += $div2;
                        $finalGene{$contig}{$pos3} += $div1;
                    }
                } else {
                    $finalGene{$contig}{$start_al} += $div1;
                    $finalGene{$contig}{$pos1} += $div2;
                    $finalGene{$contig}{$pos2} += $div2;
                    $finalGene{$contig}{$pos3} += $div1;
                }
            }
        }
    }
    close(SAM);
    
    $datestring = localtime();
    print "SAM loaded: $datestring\n";
    
    foreach my $contig (keys %finalGene) {
        my %contigGene;
        foreach my $pos (keys %{$finalGene{$contig}}) {
            foreach my $gene (keys %{$def{$contig}}) {
                if (($pos >= $def{$contig}{$gene}{'S'}) && ($pos <= $def{$contig}{$gene}{'E'})) {
                    my $geneName = $nr_protein{$gene};
                    if (defined($nr_go{$gene})) {
                        foreach my $go (keys %{$nr_go{$gene}}) {
                            $final_go{$geneName}{$go} = ''; #protein - GO = ''
                        }
                    }
                    $contigGene{$geneName} += $finalGene{$contig}{$pos};
                    last;
                }
            }
        }
    
        foreach my $gene (keys %contigGene) {
            my $coverage = $contigGene{$gene}/$defLen{$contig}{$gene};
            $final_tot{$gene}{$sample}{'coverage'} += $coverage; #hash final_tot: gene - sample = proportion
            $final_tot{$gene}{$sample}{'count'}++; #hash final_tot: gene - sample = proportion
            foreach my $ph (@rank) { #for all the ranks
                if (defined($genusExcluded{$offs{$contig}{'genus'}})) { #contigs that were excluded from taxonomy classification
                    $final_combined{$ph}{'Unclassified'}{$gene}{$sample}{'coverage'} += $coverage; #hash final_combined: rank - taxon - gene -sample = occurance.
                    $final_combined{$ph}{'Unclassified'}{$gene}{$sample}{'count'}++; #hash final_combined: rank - taxon - gene -sample = occurance.
                } else {
                    if (defined($offs{$contig}{$ph})) {
                        $final_combined{$ph}{$offs{$contig}{$ph}}{$gene}{$sample}{'coverage'} += $coverage; #hash final_combined: rank - taxon - gene -sample = occurance.
                        $final_combined{$ph}{$offs{$contig}{$ph}}{$gene}{$sample}{'count'}++;
                    } else {
                        $final_combined{$ph}{'Unclassified'}{$gene}{$sample}{'coverage'} += $coverage; #hash final_combined: rank - taxon - gene -sample = occurance.
                        $final_combined{$ph}{'Unclassified'}{$gene}{$sample}{'count'}++;
                    }
                }
            }
        }
        undef %contigGene;
    }
    undef %finalGene;
    
    $datestring = localtime();
    print "SAM analyzed: $datestring\n";
    
    foreach my $contig (keys %taxa) { #foreach contig
        $sumRead{$offs{$contig}{'genus'}} += $taxa{$contig}; #genus = number of reads
    }
    undef %taxa;
    
    foreach my $genus (sort {$sumRead{$b} <=> $sumRead{$a}} keys %sumRead) { #foreach genus
        if (defined $sumLenTrim{$genus}) {
            if ($genus ne "") { #if classified
                my $coverage = $sumRead{$genus}/$sumLenTrim{$genus};
                $finalProv{$genus} = $coverage;
            }
        }
    }
    undef %sumRead;
    
    foreach my $ph (@rank) {
        foreach my $genus (keys %finalProv) {
            $finalTaxa{$ph}{$gephy{$genus}{$ph}}{$sample}{'coverage'} += $finalProv{$genus};
            $finalTaxa{$ph}{$gephy{$genus}{$ph}}{$sample}{'count'}++;
        }
    }
    
}

undef %gephy;

my $datestring = localtime();
print "Write up taxon: $datestring\n";

foreach my $ph (@rank) { #for each rank
    open(my $file, ">./${out}/taxon_profile/${ph}.txt") or die "Coudn't open the file ./${out}/taxon_profile/${ph}.txt";
    print $file "taxa\t";
    print $file join("\t",@sam);
    foreach my $taxon (keys %{$finalTaxa{$ph}}) {
        print $file "\n$taxon";
        foreach my $sample (@sam) { #file
            if (defined($finalTaxa{$ph}{$taxon}{$sample})) {
                my $p = sprintf "%.5f", ($finalTaxa{$ph}{$taxon}{$sample}{'coverage'}/$finalTaxa{$ph}{$taxon}{$sample}{'count'});
                print $file "\t$p";
            } else {
                print $file "\t0";
            }
        }
    }
    close($file);
}

undef %finalTaxa;

$datestring = localtime();
print "Write up genus: $datestring\n";

if ($geneClassification eq 'yes') {
    
    foreach my $ph (@rank) { #foreach taxon
        foreach my $taxon (sort keys %{$final_combined{$ph}}) { #foreach taxon
            my $taxaname = $taxon;
            $taxaname =~ s/\s/__/g;
            $taxaname =~ s/\//---/g;
            open(my $file1, ">./${out}/gene_profile/${ph}/${ph}_${taxaname}.txt") or die "Coudn't open the file ./${out}/gene_profile/${ph}/${ph}_${taxaname}.txt";
            print $file1 "gene\t", join("\t", @sam);
            foreach my $gene (sort keys %{$final_combined{$ph}{$taxon}}) { #foreach gene
                if ($gene eq '') {
                    print $file1 "\nUnknown protein";
                } else {
                    print $file1 "\n$gene";
                }
                foreach my $sample (@sam) { #foreach sample
                    if (exists $final_combined{$ph}{$taxon}{$gene}{$sample}) {
                        my $p = sprintf "%.5f", ($final_combined{$ph}{$taxon}{$gene}{$sample}{'coverage'}/$final_combined{$ph}{$taxon}{$gene}{$sample}{'count'});
                        print $file1 "\t$p";
                    } else {
                        print $file1 "\t0";
                    }
                }
            }
            close($file1);
        }
    }
    
    #Create file with total gene abundance
    open(my $file2, ">./${out}/gene_profile/gene_total.txt") or die "Coudn't open the file ./${out}/gene_profile/gene_total.txt"; #file with gene names as name of the rows and samples as name columns.
    
    print $file2 "gene\t", join("\t", @sam);
    foreach my $gene (sort keys %final_tot) { #foreach gene
        if ($gene eq '') {
            print $file2 "\nUnknown protein";
        } else {
            print $file2 "\n$gene";
        }
        foreach my $sample (@sam) { #foreach sample
            if (exists $final_tot{$gene}{$sample}) {
                my $p = sprintf "%.5f", ($final_tot{$gene}{$sample}{'coverage'}/$final_tot{$gene}{$sample}{'count'});
                print $file2 "\t$p"; #here I divide by the sum adjusted with the readtrimmed proportion
            } else {
                print $file2 "\t0";
            }
        }
    }
    close($file2);
    
    #Create file with protein - gene ontology correspondences
    if ($goClassification eq 'yes') {
        open(my $file3, ">./${out}/gene_profile/gene_to_go.txt") or die "Coudn't open the file ./${out}/gene_profile/gene_to_go.txt"; #file with gene names and go
        foreach my $gene (sort keys %final_go) { #foreach gene
            if ($gene ne '') {
                print $file3 "$gene\t";
                foreach my $go (sort keys %{$final_go{$gene}}) {
                    print $file3 "$go,";
                }
                print $file3 "\n";
            }
        }
        close($file3);
    }
}
