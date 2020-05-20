#!/usr/bin/perl

use strict;

#Calculate relative abundance for taxonomy profiles

#usage: longMeta-relative [--help] [--input-folder INPUT_FOLDER] [--taxid2taxon DATABASE_FILE]

##help message
if (($ARGV[0] eq '--help') or ($ARGV[0] eq '-h')) {
    print "usage: longMeta-relative [--help] [--input-folder INPUT_FOLDER] [--taxid2taxon DATABASE_FILE]\n\n";
    
    print "--help, -h\n\tshow help message\n\n";
    
    print "--input_folder, -in INPUT_FOLDER\n\tlongMeta folder created in the previous steps of the pipeline\n";
    print "--taxid2taxon DATABASE_FILE\n\tdatabase file reporting taxids and taxonomy information\n\n";
    
    exit;
}

my %args = @ARGV;

my %allOption;
$allOption{'--input-folder'} = ''; $allOption{'-in'} = '';
$allOption{'--taxid2taxon'} = '';
$allOption{'--help'} = ''; $allOption{'-h'} = '';

foreach my $a (keys %args) {
    if (!(defined($allOption{$a}))) {
        print "Error: Invalid command option: $a. To print help message: longMeta-relative --help\n";
        exit;
    }
}

#input folder
my $in = $args{'--input-folder'};
if ($in eq "") {
    $in = $args{'-in'};
    if ($in eq "") {
        print "Error: --input-folder must be specified. To print help message: longMeta-relative --help\n";
        exit;
    }
}
if (!(-e $in and -d $in)) {
    print "Error: longMeta-relative cannot open $in. To print help message: longMeta-relative --help\n";
    exit;
}

my $taxid2taxon = $args{'--taxid2taxon'};
if ($taxid2taxon ne "") {
    if ($taxid2taxon =~ /.gz$/) {
        if (-e $taxid2taxon) {
            open(TAXADB, "gunzip -c $taxid2taxon |");
        } else {
            print "Error: longMeta-relative cannot open $taxid2taxon. To print help message: longMeta-relative --help\n";
            exit;
        }
    } else {
        if (!(open(TAXADB, "<$taxid2taxon"))) {
            print "Error: longMeta-relative cannot open $taxid2taxon. To print help message: longMeta-relative --help\n";
            exit;
        }
    }
} else {
    print "Error: --taxid2taxon must be specified. To print help message: longMeta-relative --help\n";
    exit;
}

my %gephy;
while (defined(my $input = <TAXADB>)) {
    chomp($input);
    my @taxa = split(/\t/, $input);
    if ($taxa[6] ne 'Unclassified') {
        $gephy{$taxa[6]}{'domain'} = $taxa[1];
        $gephy{$taxa[6]}{'phylum'} = $taxa[2];
        $gephy{$taxa[6]}{'class'} = $taxa[3];
        $gephy{$taxa[6]}{'order'} = $taxa[4];
        $gephy{$taxa[6]}{'family'} = $taxa[5];
        $gephy{$taxa[6]}{'genus'} = $taxa[6];
    }
}
close(TAXADB);

$gephy{'Unclassified'}{'domain'} = 'Unclassified';
$gephy{'Unclassified'}{'phylum'} = 'Unclassified';
$gephy{'Unclassified'}{'class'} = 'Unclassified';
$gephy{'Unclassified'}{'order'} = 'Unclassified';
$gephy{'Unclassified'}{'family'} = 'Unclassified';
$gephy{'Unclassified'}{'genus'} = 'Unclassified';

my @rank = ('domain', 'phylum', 'class', 'order', 'family', 'genus');

open(GENUS, "<${in}/taxon_profile/genus.txt") or die;
my $count = 1;
my @sample;
my %sumSample;
my %genusSample;

while (my $input = <GENUS>) {
    chomp($input);
    if ($count == 1) {
        my ($rest,$sample) = split(/\t/, $input,2);
        @sample = split(/\t/, $sample);
    } else {
        my ($genus,$ab) = split(/\t/, $input,2);
        my @abundance = split(/\t/, $ab);
        my $s = 0;
        foreach my $ab (@abundance) {
            $sumSample{$sample[$s]}+=$ab;
            if ($ab > 0) {
                foreach my $ph (@rank) {
                    $genusSample{$sample[$s]}{$ph}{$gephy{$genus}{$ph}}++;
                }
            }
            $s++;
        }
    }
    $count++;
}
close(GENUS);

foreach my $ph (@rank) {
    open(TAXON, "<${in}/taxon_profile/${ph}.txt") or die;
    open(my $file, ">${in}/taxon_profile/${ph}_relative.txt") or die;
    
    print "Relative abundances at $ph level were saved in ${in}/taxon_profile/${ph}_relative.txt\n";
    
    $count = 1;
    
    while(my $input=<TAXON>) {
        chomp($input);
        if ($count == 1) {
            my ($rest,$sample) = split(/\t/, $input,2);
            @sample = split(/\t/, $sample);
            print $file "$input";
        } else {
            my ($taxon,$ab) = split(/\t/, $input,2);
            my @abundance = split(/\t/, $ab);
            my $s = 0;
            print $file "$taxon";
            foreach my $ab (@abundance) {
                my $div = (($ab*$genusSample{$sample[$s]}{$ph}{$taxon})/$sumSample{$sample[$s]})*100;
                $div = sprintf "%.5f", $div;
                $s++;
                print $file "\t$div";
            }
        }
        $count++;
        print $file "\n";
    }
    close(TAXON);
    close($file);
}
