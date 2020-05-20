#!/usr/bin/perl

use strict;

#Explore dataset created with longMeta

#usage: longMeta-explore [--help] [--taxonomy-input INPUT_FILE] [--sam-input INPUT_FILE] [--length-input INPUT_FILE] [--output-folder OUTPUT_FOLDER] [--ignore_uncl {yes,no}] [--perc-limit POS_NUMBER] [--phyla-exclusion {yes,no}] [--gene-input INPUT_FILE] [--acc2gene DATABASE_FILE] [--acc2go DATABASE_FILE] [--read-length POS_INTEGER] [--alignment-AS NUMBER] [--output-type {gene,GO}]

#help message
if (($ARGV[0] eq '--help') or ($ARGV[0] eq '-h')) {
    print "usage: longMeta-explore [--help] [--input-folder INPUT_FOLDER] [--gene STRING] [--go {go-numbers}] [--go-name STRING] [--go-ontology {molecular_function, biological_process, cellular_component}] [--taxon STRING] [--taxon-rank {domain, phylum, class, order, family, genus}] [--output-file OUTPUT_FILE] [--go2def DATABASE] [--output-type {gene,GO}]\n\n";
    
    print "--help, -h\n\tshow help message\n\n";
    
    print "--input-folder, -in INPUT_FOLDER\n\tlongMeta folder created in the previous steps of the pipeline\n";
    
    print "--gene, -g STRING\n\tpartial or entire gene name to use for gene retrieval. If the name is made up from more than one word, they must be enclosed in quotation marks. Special characters (e.g. ` or /) must be escaped with \‘\\’. Write ‘all’ if all the genes must be retrieved.\n";
    print "--go-accession, -ga GO_CODE\n\tgene ontology code to use for gene retrival\n";
    print "--go-name, -gn STRING\n\tgene ontology name (or part) to use for gene retrival\n\n";
    print "--go-ontology, -go {molecular_function, biological_process, cellular_component}\n\tgene ontology category\n";
    
    print "--go2def, DATABASE\n\tgene ontology database\n";
    
    print "--taxon, -t STRING\n\ttaxa (comma-separated)\n";
    print "--taxon-rank, -r {domain, phylum, class, order, family, genus}\n\trank of --taxon\n";
    
    print "--output-file, -out OUTPUT_FILE\n\toutput file reporting the genes of interest\n";
    print "--output-type {gene,GO}\n\toutput type\n\n";
    
    exit;
}

my %args = @ARGV;

#check options
my %allOption;
$allOption{'--input-folder'} = ''; $allOption{'-in'} = '';
$allOption{'--gene'} = ''; $allOption{'-g'} = '';

$allOption{'--go-accession'} = ''; $allOption{'-ga'} = '';
$allOption{'--go-name'} = ''; $allOption{'-gn'} = '';
$allOption{'--go-ontology'} = ''; $allOption{'-go'} = '';

$allOption{'--go2def'} = ''; $allOption{'-gc'} = '';

$allOption{'--taxon'} = ''; $allOption{'-t'} = '';
$allOption{'--taxon-rank'} = ''; $allOption{'-r'} = '';

$allOption{'--output-file'} = ''; $allOption{'-out'} = '';
$allOption{'--output-type'} = '';

$allOption{'--help'} = ''; $allOption{'-h'} = '';

foreach my $a (keys %args) {
    if (!(defined($allOption{$a}))) {
        print "Error: Invalid command option: $a. To print help message: longMeta-explore --help\n";
        exit;
    }
}

#import options

#input folder
my $in = $args{'--input-folder'};
if ($in eq "") {
    $in = $args{'-in'};
    if ($in eq "") {
        print "Error: --input-folder must be specified. To print help message: longMeta-explore --help\n";
        exit;
    }
}
if (!(-e $in and -d $in)) {
    print "Error: longMeta-explore cannot open $in. To print help message: longMeta-explore --help\n";
    exit;
}

my $define = "no";

#gene
my $gene = $args{'--gene'};
if ($gene eq "") {
    $gene = $args{'-g'};
    if ($gene ne "") {
        $define = "yes";
    }
}

#goAcc
my $goAcc = $args{'--go-accession'};
if ($goAcc eq "") {
    $goAcc = $args{'-ga'};
}

#goName
my $goName = $args{'--go-name'};
if ($goName eq "") {
    $goName = $args{'-gn'};
}

#goOnt
my %goOntAll;
$goOntAll{'molecular_function'} = '';
$goOntAll{'biological_process'} = '';
$goOntAll{'cellular_component'} = '';
my $goOnt = $args{'--go-ontology'};

if ($goOnt eq "") {
    $goOnt = $args{'-go'};
} else {
    if (!(defined($goOntAll{$goOnt}))) {
        print "Error: --go-ontology must be either molecular_function, biological_process or cellular_component. To print help message: longMeta-explore --help\n";
        exit;
    }
}

my %goAllDB;
my @goAll;
if (($goAcc ne '') or ($goName ne '') or ($goOnt ne '')) {
    my $goDB = $args{'--go2def'};
    if ($goDB eq "") {
        print "Error: --go2def must be defined. To print help message: longMeta-explore --help\n";
        exit;
    } else {
        if ($goDB =~ /.gz$/) {
            if (-e $goDB) {
                open(GO, "gunzip -c $goDB |");
            } else {
                print "Error: longMeta-explore cannot open $goDB. To print help message: longMeta-explore --help\n";
                exit;
            }
        } else {
            if (!(open(GO, "<$goDB"))) {
                print "Error: longMeta-explore cannot open $goDB. To print help message: longMeta-explore --help\n";
                exit;
            }
        }
    }
}

if ($goAcc ne '') {
    while(defined(my $input = <GO>)) {
        chomp($input);
        my @info = split(/\t/, $input);
        $goAllDB{$info[0]} = '';
    }
    close(GO);
    #load go file
    @goAll = split(/,/, $goAcc);
    foreach my $go (@goAll) {
        if (!(defined($goAllDB{$go}))) {
            print "Error: $go is not present in the Gene Ontology database. To print help message: longMeta-explore --help\n";
            exit;
        }
    }
}

if ($goName ne '') {
    while(defined(my $input = <GO>)) {
        chomp($input);
        if ($input =~ /[^A-Za-z]$goName[^A-Za-z]/i) {
            my @info = split(/\t/, $input);
            if ($goOnt ne "") {
                if ($info[2] eq $goOnt) {
                    push @goAll, $info[0];
                }
            } else {
                push @goAll, $info[0];
            }
            
        }
    }
    close(GO);
    
    if ($goAll[0] eq "") {
        print "Error: $goName is not present in any of the the Gene Ontology names. To print help message: longMeta-explore --help\n";
        exit;
    }
}

my %nameGO;

if (($goOnt ne '') && ($goName eq '')) {
    while(defined(my $input = <GO>)) {
        chomp($input);
        my @info = split(/\t/, $input);
        if ($info[2] eq $goOnt) {
            push @goAll, $info[0];
            $nameGO{$info[0]} = $info[1];
        }
    }
    close(GO);
}


if (($gene eq '') && ($goAcc eq "") && ($goName eq "") && ($goOnt eq "")) {
    print "Error: either --gene, --go-accession, --go_name or --go-ontology must be specified. To print help message: longMeta-explore --help\n";
    exit;
}

#taxon
my $taxon = $args{'--taxon'};
if ($taxon eq "") {
    $taxon = $args{'-t'};
}
$taxon =~ s/\//---/g;
$taxon =~ s/ /--/g;
my @taxonAll;
if ($taxon ne "") {
    @taxonAll = split(/,/, $taxon);
}

#rank
my $rank = $args{'--taxon-rank'};
if ($rank eq "") {
    $rank = $args{'-r'};
}
if ($rank eq "") {
    if ($taxon ne "") {
        print "Error: --taxon-rank must be specified when --taxon is defined. To print help message: longMeta-explore --help\n";
        exit;
    }
} else {
    if (($rank ne "domain") && ($rank ne "phylum") && ($rank ne "class") && ($rank ne "order") && ($rank ne "family") && ($rank ne "genus")) {
        print "Error: --taxon-rank must be either domain, phylum, class, order, family or genus. To print help message: longMeta-explore --help\n";
        exit;
    }
}

#output file
my $out = $args{'--output-file'};
my $file;
if ($out ne "") {
    if (!(open($file, ">$out"))) {
        print "Error: longMeta-explore cannot open $out. To print help message: longMeta-explore --help\n";
        exit;
    }
} else {
    $out = $args{'-out'};
    if ($out eq "") {
        print "Error: --output-file must be defined. To print help message: longMeta-explore --help\n";
        exit;
    } else {
        if (!(open($file, ">$out"))) {
            print "Error: longMeta-explore cannot open $out. To print help message: longMeta-explore --help\n";
            exit;
        }
    }
}

#output type
my $type = $args{'--output-type'};
if ($type eq '') {
    $type = 'gene';
} elsif (($type ne 'gene') && ($type ne 'GO')) {
    print "Error: --output-type must be either gene or GO. To print help message: longMeta-explore --help\n";
    exit;
}

if ($gene ne "") {
    my $geneDef = $gene;
    if ($gene eq 'all') {
        $geneDef = '';
    }
    if (($taxon eq "") && ($rank eq "")) { #both taxon and rank are undefined
        my $count = 0;
        open(TMP, "<${in}/gene_profile/gene_total.txt");
        while(defined(my $input=<TMP>)) {
            chomp($input);
            if ($count == 0) {
                print $file "$input\n";
            }
            if ($input =~ /$geneDef/i) {
                print $file "$input\n";
            }
            $count++;
        }
        close(TMP);
    } elsif (($taxon eq "") && ($rank ne "")) { #rank is undefined
        my @fileAll = `ls ${in}/gene_profile/$rank/`;
        my $count = 0;
        foreach my $taxonTmp (@fileAll) {
            chop($taxonTmp);
            open(TMP, "<${in}/gene_profile/$rank/$taxonTmp");
            while(defined(my $input=<TMP>)) {
                chomp($input);
                if ($count == 0) {
                    print $file "taxon\t$input\n";
                }
                if ($input =~ /$geneDef/i) {
                    my ($taxon) = $taxonTmp =~ /${rank}_(.*).txt/;
                    $taxon =~ s/---/\//g;
                    $taxon =~ s/--/ /g;
                    print $file "$taxon\t$input\n";
                }
                $count++;
            }
            close(TMP);
        }
    } else {
        my $count = 0;
        foreach my $taxonTmp (@taxonAll) {
            open(TMP, "<${in}/gene_profile/$rank/${rank}_${taxonTmp}.txt");
            while(defined(my $input=<TMP>)) {
                chomp($input);
                if ($count == 0) {
                    print $file "taxon\t$input\n";
                }
                if ($input =~ /$geneDef/i) {
                    print $file "$taxonTmp\t$input\n";
                }
                $count++;
            }
            close(TMP);
        }
    }
}


if (($goAcc ne "") or ($goName ne "") or ($goOnt ne "")) {
    
    open(TMP, "<${in}/gene_profile/gene_to_go.txt");
    my %goReverse;
    while(defined(my $input=<TMP>)) {
        my @info = split(/\t/, $input);
        chop($info[1]);
        my @infoGO = split(/,/, $info[1]);
        foreach my $go (@infoGO) {
            $goReverse{$go}{$info[0]} = '';
        }
        
    }
    close(TMP);
    
    my %geneRet;
    my %geneRetGO;
    
    foreach my $go (@goAll) {
        my $count = 0;
        if (defined($goReverse{$go})) {
            foreach my $gene (keys %{$goReverse{$go}}) {
                if ($type eq 'gene') {
                    $geneRet{$gene} = '';
                    $count++;
                } else {
                    $geneRet{$gene}{$go} = '';
                }
            }
        }
        if ($goAcc ne "") {
            if ($count == 1) {
                print "There was 1 gene associated with $go\n";
            } elsif ($count > 1) {
                print "There were $count gene associated with $go\n";
            } else {
                print "There were not genes associated with $go\n";
            }
        }
    }
    
    
    #retrieve only proteins associated with GO
    if (($taxon eq "") && ($rank eq "")) { #both taxon and rank are undefined
        my $count = 0;
        open(TMP, "<${in}/gene_profile/gene_total.txt");
        
        if ($type eq 'gene') {
            while(defined(my $input=<TMP>)) {
                chomp($input);
                my ($gene, $rest) = split(/\t/, $input, 2);
                if ($count == 0) {
                    print $file "$input\n";
                }
                if (defined($geneRet{$gene})) {
                    print $file "$input\n";
                }
                $count++;
            }
        } else {
            my @sample;
            my %sumSample;
            my %countSample;
            
            while(defined(my $input=<TMP>)) {
                chomp($input);
                if ($count == 0) {
                    my ($rest,$sample) = split(/\t/, $input,2);
                    @sample = split(/\t/, $sample);
                    print $file "GO\tname\t$sample\n";
                } else {
                    my ($gene,$ab) = split(/\t/, $input,2);
                    if (defined($geneRet{$gene})) {
                        my @abundance = split(/\t/, $ab);
                        my $s = 0;
                        foreach my $ab (@abundance) {
                            if ($ab > 0) {
                                foreach my $go (keys %{$geneRet{$gene}}) {
                                    $sumSample{$go}{$sample[$s]}+=$ab;
                                    $countSample{$go}{$sample[$s]}++;
                                }
                            }
                            $s++;
                        }
                    }
                }
                $count++;
            }
            close(TMP);
            
            foreach my $go (sort keys %sumSample) {
                print $file "$go\t$nameGO{$go}";
                foreach my $s (@sample) {
                    if (defined($sumSample{$go}{$s})) {
                        my $div = $sumSample{$go}{$s}/$countSample{$go}{$s};
                        $div = sprintf "%.5f", $div;
                        print $file "\t$div";
                    } else {
                        print $file "\t0";
                    }
                }
                print $file "\n";
            }
        
        }
    } elsif (($taxon eq "") && ($rank ne "")) { #rank is undefined
        my @fileAll = `ls ${in}/gene_profile/$rank/`;
        my $count = 0;
        foreach my $taxonTmp (@fileAll) {
            chop($taxonTmp);
            open(TMP, "<${in}/gene_profile/$rank/$taxonTmp");
            while(defined(my $input=<TMP>)) {
                chomp($input);
                my ($gene, $rest) = split(/\t/, $input, 2);
                if ($count == 0) {
                    print $file "taxon\t$input\n";
                }
                if (defined($geneRet{$gene})) {
                    my ($taxon) = $taxonTmp =~ /${rank}_(.*).txt/;
                    $taxon =~ s/---/\//g;
                    $taxon =~ s/--/ /g;
                    print $file "$taxon\t$input\n";
                }
                $count++;
            }
            close(TMP);
        }
    } else {
        my $count = 0;
        foreach my $taxonTmp (@taxonAll) {
            open(TMP, "<${in}/gene_profile/$rank/${rank}_${taxonTmp}.txt");
            while(defined(my $input=<TMP>)) {
                chomp($input);
                my ($gene, $rest) = split(/\t/, $input, 2);
                if ($count == 0) {
                    print $file "taxon\t$input\n";
                }
                if (defined($geneRet{$gene})) {
                    print $file "$taxonTmp\t$input\n";
                }
                $count++;
            }
            close(TMP);
        }
    }
}
close($file);



#check if longMeta found any match with the user's criteria
my $count = `wc -l $out`;
my @count = split(/\s/, $count, 2);
if ($count < 2) {
    print "longMeta-explore did not find any match with these criteria\n";
    `rm $out`;
    exit;
}

my @sample;
my %sumSampleAll;
my %countSampleAll;
my $count = 0;

open(TMP, "<${in}/gene_profile/gene_total.txt");
open(my $file, ">${out}.all");

while(defined(my $input=<TMP>)) {
    chomp($input);
    if ($count == 0) {
        my ($rest,$sample) = split(/\t/, $input,2);
        @sample = split(/\t/, $sample);
        print $file "all\t$sample\n";
    } else {
        my ($gene,$ab) = split(/\t/, $input,2);
        my @abundance = split(/\t/, $ab);
        my $s = 0;
        foreach my $ab (@abundance) {
            if ($ab > 0) {
                $sumSampleAll{$sample[$s]}+=$ab;
                $countSampleAll{$sample[$s]}++;
            }
            $s++;
        }
    }
    $count++;
}
close(TMP);

print $file "all";
foreach my $s (@sample) {
    my $div = $sumSampleAll{$s}/$countSampleAll{$s};
    $div = sprintf "%.5f", $div;
    print $file "\t$div";
}
close($file)
