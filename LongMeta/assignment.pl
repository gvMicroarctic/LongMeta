#!/usr/bin/perl

#this script is called from longMeta-assignment script

use strict;
use POSIX;

my %args = @ARGV;

#import options

my @diamondpair = split(/,/, $args{'-p'}); #PE files - paired reads (e.g. PE Illumina reads)
my @diamondunpair = split(/,/, $args{'-u'}); #SE files - unpaired reads or contigs (e.g. SE Illumina reads, contigs, MinION reads)

my $geneoutput = $args{'--gene-output'};
my $taxaoutput = $args{'--taxonomy-output'};

my $logName;
if ($geneoutput ne '') { #gene output file
    $logName = $geneoutput . ".log";
} else {
    $logName = $taxaoutput . ".log";
}

my $alsotax = $args{'--taxonomy'}; #taxonomy assignment
my $tmp = $args{'--tmp'}; #temporary folder

#taxonomy assignment parameters
my $max_match = $args{'--max-equal'}; #NA
my $cutoff_best = $args{'--cutoff-ID-best'}; #90

my $min_best = $args{'--min-best'}; #3

my $all_lca = $args{'--all-lca'}; #mixed (or all, gene)
my $min_lca = $args{'--min-lca'}; #3

my $cutoff_lca = $args{'--cutoff-ID-lca'}; #90


#gene assignment parameter
my $overlap = $args{'--overlap'}; #3

my $acc2taxid;
my $taxid2taxon;

if ($alsotax ne "") {
    $taxid2taxon = $args{'--taxid2taxon'};
    $acc2taxid = $args{'--acc2taxid'};
    
    if ($taxid2taxon =~ /.gz$/) {
        open(TAXA1, "gunzip -c $taxid2taxon |") or die "Coudn't open the file with taxid and taxon correspondences: $!";
    } else {
        open(TAXA1, "<$taxid2taxon") or die "Coudn't open the file with taxid and taxon correspondences: $!";
    }
    
    if ($acc2taxid =~ /.gz$/) {
        open(TAXA2, "gunzip -c $acc2taxid |") or die "Coudn't open the file with accession and taxid correspondences: $!";
    } else {
        open(TAXA2, "<$acc2taxid") or die "Coudn't open the file with accession and taxid correspondences: $!";
    }
}

#input files
my %all_file;
my %all_fileType;
my $c = 1;
if ($diamondpair[1] ne '') { #paired reads
    foreach my $diamond (@diamondpair) {
        $all_file{$diamond} = "2";
        if ($c == 1) {
            $all_fileType{$diamond} = "F";
        } else {
            $all_fileType{$diamond} = "R";
        }
        $c++;
    }
}
if ($diamondunpair[0] ne '') {
    foreach my $diamond (@diamondunpair) { #unpaired reads
        $all_file{$diamond} = "1";
    }
}

#gene assignment
my %all_acc; #all the accession numbers
my $tot = 0; #total of sequences
my $old = ''; #last sequence
my %perc; #
my %inf; #
my %taxa_finalAll; #all assigned genes
my %taxa_final; #all assigned genes with an id > cutoff
my %paired0; #paired sequences occurance
my %count; #how many genes per sequence
my $c = 0;
my $diamond = '';
my $old_diamond ='';
my $file;

open(my $fileLog, ">./$tmp/$logName") or die "Coudn't open the file $logName: $!";

my $gene;
if ($geneoutput ne "") {
    open($file, ">./$tmp/$geneoutput") or die "Coudn't open the file $geneoutput: $!";
    $gene = "yes";
}

foreach $diamond (sort keys %all_file) {
    if ($diamond =~ /.gz$/) {
        open(FILE, "gunzip -c $diamond |") or die "Coudn't open the file $diamond: $!";
    } else {
        open(FILE, "<$diamond") or die "Coudn't open the file $diamond: $!";
    }
    print $fileLog "FILE\t$diamond\n";
    while (defined(my $input = <FILE>)) {
        chomp($input);
        my @info = split(/\t/, $input);
        $all_acc{$info[1]} = '';
        if (($old ne $info[0]) && ($old ne '')) {
            $tot++;
            &annotation($old_diamond,\%perc,\%inf);
            undef %perc;
            undef %inf;
        }
        $perc{$info[2]}{$info[1]}=$info[11]; #sequence-identity_score-acc_num = bit score
        $inf{$info[1]} = $input; #acc_num=all_line
        $old = $info[0];
        $old_diamond = $diamond;
    }
    close(FILE);
}
$tot++;
&annotation($old_diamond,\%perc,\%inf); #last sequence entry
undef %perc;
undef %inf;

if ($gene eq "yes") {
    close($file);
    #Gene annotation info
    print $fileLog "AVERAGE:\t$c\t$tot\n";
    foreach my $key (sort {$b <=> $a} keys %count) {
        print $fileLog "DETAIL:\t$count{$key}\t$key\n";
    }
}

#taxonomy assignment (only if specified in options)
my %paired;
my %taxa;
my %accTotaxid;

my $done0 = keys %taxa_finalAll;
print $fileLog "ALL\t$done0\n";

if ($alsotax ne "") {
    my %taxid;
    my %taxid_sp;
    my %genus_new;
    my %genus_info;
    my @ranks = ('domain', 'phylum', 'class', 'order', 'family','genus','species');
    my @ranks_cut = ('domain', 'phylum', 'class', 'order', 'family');
    
    foreach my $sequence (keys %paired0) {
        my $count = 0;
        foreach my $presence (keys %{$paired0{$sequence}}) {
            $count++;
        }
        $paired{$sequence} = $count; #sequence = if in one or both paired files
    }
    undef %paired0;
    
    open($file, ">./$tmp/$taxaoutput") or die "Coudn't open the file $taxaoutput: $!";
    
    while (defined(my $input = <TAXA1>)) {
        chomp($input);
        my @info =split(/\t/,$input);
        
        $taxid_sp{$info[0]}{'domain'} = $info[1]; #taxid - domain
        $taxid_sp{$info[0]}{'phylum'} = $info[2]; #taxid - phylum
        $taxid_sp{$info[0]}{'class'} = $info[3]; #taxid - class
        $taxid_sp{$info[0]}{'order'} = $info[4]; #taxid - order
        $taxid_sp{$info[0]}{'family'} = $info[5]; #taxid - family
        $taxid_sp{$info[0]}{'genus'} = $info[6]; #taxid - genus
        $taxid_sp{$info[0]}{'species'} = $info[7]; #taxid - species
        
        push @{$taxa{$info[0]}}, ($info[1],$info[2],$info[3],$info[4],$info[5],$info[6],$info[7]); #taxid = @(taxon ranks); only taxonomy path for interesting accession numbers
        
        $genus_new{$info[6]}{'domain'} = $info[1]; #genus - domain
        $genus_new{$info[6]}{'phylum'} = $info[2]; #genus - phylum
        $genus_new{$info[6]}{'class'} = $info[3]; #genus - class
        $genus_new{$info[6]}{'order'} = $info[4]; #genus - order
        $genus_new{$info[6]}{'family'} = $info[5]; #genus - family
        
    }
    close(TAXA1);
    
    while (defined(my $input = <TAXA2>)) {
        chomp($input);
        my @info =split(/\t/,$input);
        if (defined($all_acc{$info[0]})) { #if the acc number is the files
            $accTotaxid{$info[0]} = $info[1]; #acc = taxid;
            if (defined($taxid_sp{$info[1]})) { #if there is a taxonomy associated with taxid
                foreach my $rank (@ranks) {
                    if ($taxid_sp{$info[1]}{$rank} !~ /Unclassified/) {
                        my $unc = $taxid_sp{$info[1]}{$rank} . "_Unclassified";
                        $genus_info{$taxid_sp{$info[1]}{'genus'}}{$unc} = ''; #genus - all_higher ranks = ''
                    } else {
                        $genus_info{$taxid_sp{$info[1]}{'genus'}}{$taxid_sp{$info[1]}{$rank}} = ''; #genus - all_higher ranks = ''
                    }
                    $genus_info{$taxid_sp{$info[1]}{'genus'}}{'Unclassified'} = ''; #genus - all_higher ranks = ''
                }
            } elsif ((!(defined($taxid_sp{$info[1]}{'genus'}))) or ($info[1] eq '')) { #if taxid is '' or no taxonomy associated with taxid: Unclassified
                foreach my $rank (@ranks) {
                    push @{$taxa{$info[1]}}, 'Unclassified'; #taxid = @(taxon ranks); only taxonomy path for interesting accession numbers
                }
            }
        }
    }
    undef %taxid_sp;
    undef %all_acc;
    close(TAXA2);
    
    my $done;
    
    #taxonomy assigned with BEST approach (first cycle)
    if ($alsotax ne "l") {
        foreach my $sequence (keys %taxa_final) {
            my %cla; #undef
            my %un; #undef
            my @det = (); #undef
            my $in = 0;
            my $in1 = 0;
            my $out = 0;
            my $s = 0;
            my $def = '';
            my %check; #undef
            my %check_sp; #undef
            my $totAcc = 0;
            $tot++;
            
            foreach my $acc (keys %{$taxa_final{$sequence}}) {
                $totAcc++;
                if ($taxa{$accTotaxid{$acc}}[5] eq '') { #if this acc is not present in TAXA2 file
                    $check{'Unclassified'}++; #genus = count; #do not have contig
                    $check_sp{'Unclassified'}{'Unclassified'}++; #genus - species = count
                } else {
                    $check{$taxa{$accTotaxid{$acc}}[5]}++; #genus = count; #do not have contig
                    $check_sp{$taxa{$accTotaxid{$acc}}[5]}{$taxa{$accTotaxid{$acc}}[6]}++; #genus - species = count; #do not have contig
                }
            }
            
            foreach my $genus (keys %check) {
                if ($genus =~ /Unclassified/) {
                    $un{$genus} = ''; #unclassified at genus level
                } else {
                    $cla{$genus} = ''; #classified at genus level
                }
            }
            foreach my $ge (keys %cla) {
                $in++;
                $def = $ge;
            }
            foreach my $ge (keys %un) {
                $in1++;
            }
            if (($in == 1) && ($totAcc >= $min_best)) { #if there is only a genus and more than --min-best matches
                foreach my $ge (keys %un) {
                    if (!(defined($genus_info{$def}{$ge}))) { #unclassified genus do not belong to same taxonomy as the genus
                        $out++;
                    }
                }
                if ($out == 0) { #if all unclassified belong to same taxonomy as the genus
                    push @det, $sequence; #sequence
                    foreach my $rank (@ranks_cut) {
                        push @det, $genus_new{$def}{$rank}; #intermediate taxa
                    }
                    push @det, $def; #genus
                    
                    foreach my $sp (sort {$check_sp{$def}{$b}<=> $check_sp{$def}{$a}} keys %{$check_sp{$def}}) {
                        if ($s == 0) {
                            push @det, $sp; #most abundant species (belonging to that genus
                            $s++;
                        } else {
                            last;
                        }
                    }
                    if (defined($paired{$sequence})) {
                        print $file join "\t", @det, "\n";
                        print $file join "\t", @det, "\n";
                        $done++;
                        $done++;
                    } else {
                        print $file join "\t", @det, "\n";
                        $done++;
                    }
                    delete $taxa_finalAll{$sequence};
                }
            }
        }
        print $fileLog "BEST\t$done\n";
    }
    
    undef %taxa_final;
    
    #taxonomy assigned with WEIGHTED LCA approach (second cycle)
    if ($alsotax ne "b") {
        
        my @retrieved;
        my $file_count = 0;
        my %keep;
        my @info;
        foreach my $diamond (sort {$all_file{$a} <=> $all_file{$b}} keys %all_file) {
            if ($diamond =~ /.gz$/) {
                open(FILE, "gunzip -c $diamond |") or die "Coudn't open the file $diamond: $!";
            } else {
                open(FILE, "<$diamond") or die "Coudn't open the file $diamond: $!";
            }
            $old = '';
            if ($all_file{$diamond} == 1) {
                while (defined(my $input = <FILE>)) {
                    chomp($input);
                    my @info = split(/\t/, $input);
                    if (defined($taxa_finalAll{$info[0]})) {
                        if (($info[0] ne $old) && ($old ne '')) {
                            &assignment(\@retrieved);
                            undef @retrieved;
                        }
                        my $input1 = $info[1] . "\t" . $info[2] . "\t" . $info[8] . "\t" . $info[9] . "\t" . $info[11];
                        
                        if ($all_lca = "mixed") { #mixed approach
                            my $geneCount = keys %{$taxa_finalAll{$info[0]}};
                            if ((defined($taxa_finalAll{$info[0]}{$info[1]})) && ($geneCount >= $min_lca)) { #if taxonomy is still to be assigned and there are more than --min-lca matches
                                push @retrieved, $input1;
                            } else {
                                push @retrieved, $input1;
                            }
                            
                        } elsif ($all_lca = "gene") { #use only assigned genes for lca
                            if (defined($taxa_finalAll{$info[0]}{$info[1]})) { #if taxonomy is still to be assigned
                                push @retrieved, $input1;
                            }
                            
                        } else {
                            push @retrieved, $input1;
                        }
                        $old = $info[0];
                    }
                }
                &assignment(\@retrieved); #last loop
                undef @retrieved;
                close(FILE);
            } else {
                $file_count++;
                if ($file_count == 1) {
                    while (defined(my $input = <FILE>)) {
                        chomp($input);
                        @info = split(/\t/, $input);
                        if (defined($taxa_finalAll{$info[0]})) {
                            if (($info[0] ne $old) && ($old ne '')) {
                                if ($paired{$old} == 2) {
                                    foreach my $input1 (@retrieved) {
                                        $keep{$old}{$input1} = '';
                                    }
                                } else {
                                    &assignment(\@retrieved);
                                }
                                undef @retrieved;
                            }
                            
                            
                            my $input1 = $info[1] . "\t" . $info[2] . "\t" . $info[8] . "\t" . $info[9] . "\t" . $info[11];
                            if ($all_lca = "mixed") { #mixed approach
                                my $geneCount = keys %{$taxa_finalAll{$info[0]}};
                                if ((defined($taxa_finalAll{$info[0]}{$info[1]})) && ($geneCount >= $min_lca)) { #if taxonomy is still to be assigned and there are more than --min-lca matches
                                    push @retrieved, $input1;
                                } else {
                                    push @retrieved, $input1;
                                }
                                
                            } elsif ($all_lca = "gene") { #use only assigned genes for lca
                                if (defined($taxa_finalAll{$info[0]}{$info[1]})) { #if taxonomy is still to be assigned
                                    push @retrieved, $input1;
                                }
                                
                            } else {
                                push @retrieved, $input1;
                            }
                            $old = $info[0];
                        }
                    }
                    if ($paired{$old} == 2) {
                        foreach my $input1 (@retrieved) {
                            $keep{$old}{$input1} = '';
                        }
                    } else {
                        &assignment(\@retrieved);
                    }
                    undef @retrieved;
                    close(FILE);
                } else {
                    while (defined(my $input = <FILE>)) {
                        chomp($input);
                        my @info = split(/\t/, $input);
                        if (defined($taxa_finalAll{$info[0]})) {
                            if (($info[0] ne $old) && ($old ne '')) {
                                if ($paired{$old} == 2) {
                                    foreach my $input1 (keys %{$keep{$old}}) {
                                        push @retrieved, $input1;
                                    }
                                    &assignment(\@retrieved);
                                } else {
                                    &assignment(\@retrieved);
                                }
                                undef @retrieved;
                                delete $keep{$old};
                            }
                            my $input1 = $info[1] . "\t" . $info[2] . "\t" . $info[8] . "\t" . $info[9] . "\t" . $info[11];
                            if ($all_lca = "mixed") { #mixed approach
                                my $geneCount = keys %{$taxa_finalAll{$info[0]}};
                                if ((defined($taxa_finalAll{$info[0]}{$info[1]})) && ($geneCount >= $min_lca)) { #if taxonomy is still to be assigned and there are more than --min-lca matches
                                    push @retrieved, $input1;
                                } else {
                                    push @retrieved, $input1;
                                }
                                
                            } elsif ($all_lca = "gene") { #use only assigned genes for lca
                                if (defined($taxa_finalAll{$info[0]}{$info[1]})) { #if taxonomy is still to be assigned
                                    push @retrieved, $input1;
                                }
                                
                            } else {
                                push @retrieved, $input1;
                            }
                            $old = $info[0];
                        }
                    }
                    if ($paired{$old} == 2) {
                        foreach my $input1 (keys %{$keep{$old}}) {
                            if ($all_lca = "mixed") { #mixed approach
                                my $geneCount = keys %{$taxa_finalAll{$info[0]}};
                                if ((defined($taxa_finalAll{$info[0]}{$info[1]})) && ($geneCount >= $min_lca)) { #if taxonomy is still to be assigned and there are more than --min-lca matches
                                    push @retrieved, $input1;
                                } else {
                                    push @retrieved, $input1;
                                }
                                
                            } elsif ($all_lca = "gene") { #use only assigned genes for lca
                                if (defined($taxa_finalAll{$info[0]}{$info[1]})) { #if taxonomy is still to be assigned
                                    print "$info[0] and $info[1]\n";
                                    push @retrieved, $input1;
                                }
                                
                            } else {
                                push @retrieved, $input1;
                            }
                        }
                        &assignment(\@retrieved);
                    } else {
                        &assignment(\@retrieved);
                    }
                    undef @retrieved;
                    close(FILE);
                }
            }
        }
    }
    close($file);
    my $done1 = keys %taxa_finalAll;
    print $fileLog "LCA\t$done1\n";
}

#gene annotation and taxonomy assignment
sub annotation {
    my $diamond = $_[0];
    my %perc = %{$_[1]};
    my %inf = %{$_[2]};
    my $aa = 0;
    
    my %finalprint;
    my %final;
    
    my $highest_id;
    my $highest_bit;
    
    if ($all_file{$diamond} == 2) {
        $paired0{$old}{$diamond} = ''; #save if the read is part of a pair
    }
    
    if (($gene eq "yes") or ($alsotax ne "l")) {
        foreach my $id (sort {$b <=> $a} keys %perc) { #foreach identity_scor (higher->lower)
            my $max = 0;
            foreach my $acc (sort {$perc{$id}{$b} <=> $perc{$id}{$a}} keys %{$perc{$id}}) { #foreach acc_num
                my @sec = split(/\t/, $inf{$acc});
                my $inside = 0;
                my $start;
                my $end;
                if ($sec[6] < $sec[7]) {
                    $start=$sec[6];
                    $end=$sec[7];
                } else {
                    $start=$sec[7];
                    $end=$sec[6];
                }
                if ($aa == 0) { #first entry
                    $highest_id = $id;
                    $highest_bit = $perc{$id}{$acc};
                    $aa++;
                    $final{$aa}{'S'}=$start;
                    $final{$aa}{'E'}=$end;
                    $finalprint{$start} = $inf{$acc};
                    $taxa_finalAll{$old}{$acc} = ''; #contig - acc_num (no ID score limit)
                    if (($id >= $cutoff_best) && ($alsotax ne "")) {
                        $taxa_final{$old}{$acc} = ''; #contig - acc_num
                    }
                } else {
                    if ($max_match eq 'NA') { #if a sequence area was assigned to a Diamond match but other Diamond matches had the same identity score and bit-score, all these matches are saved to assign taxonomy.
                        if (($highest_id == $id) && ($highest_bit == $perc{$id}{$acc}) && ($alsotax ne "")) { #if top score to more than one match
                            $taxa_finalAll{$old}{$acc} = ''; #contig - acc_num (no ID score limit)
                            if (($id >= $cutoff_best)) {
                                $taxa_final{$old}{$acc} = ''; #contig - acc_num
                            }
                        }
                    } else { #if a sequence area was assigned to a Diamond match but other Diamond matches had the same identity score and bit-score, only $max_match are used for taxonomy assignment.
                        $max++;
                        if ($max <= $max_match) {
                            if (($highest_id == $id) && ($highest_bit == $perc{$id}{$acc}) && ($alsotax ne "")) { #if top score to more than one match
                                $taxa_finalAll{$old}{$acc} = ''; #contig - acc_num (no ID score limit)
                                if (($id >= $cutoff_best)) {
                                    $taxa_final{$old}{$acc} = ''; #contig - acc_num
                                }
                            }
                        }
                    }
                    foreach my $key (sort keys %final) {
                        my $startCoding = $final{$key}{'S'} + $overlap;
                        my $endCoding = $final{$key}{'E'} - $overlap;
                        if ((($start < $startCoding) && ($end > $startCoding)) or (($start < $endCoding) && ($end > $endCoding)) or (($start >= $startCoding) && ($end <= $endCoding)) or (($start <= $startCoding) && ($end >= $endCoding)) && ($inside == 0)) {
                            $inside = 1;
                        }
                    }
                    if ($inside == 0) {
                        $aa++;
                        $final{$aa}{'S'}=$start;
                        $final{$aa}{'E'}=$end;
                        $finalprint{$start} = $inf{$acc};
                        $taxa_finalAll{$old}{$acc} = '';
                        if (($id >= $cutoff_best) && ($alsotax ne "")) {
                            $taxa_final{$old}{$acc} = ''; #contig - acc_num (no ID score limit)
                        }
                    }
                }
            }
        }
    }
    if ($gene eq "yes") {
        
        my $d = 0; #how many genes in the contig
        foreach my $key_start (sort {$a <=> $b} keys %finalprint) { #foreach match
            $c++;
            print $file "$finalprint{$key_start}\t$all_fileType{$diamond}\n";
            $d++;
        }
        $count{$d}++; #number of genes in the contig = how many contigs
    }
}

close($fileLog);

sub assignment {
    my @retrieved = @{$_[0]};
    my %final;
    my %taxa1;
    my %taxa2;
    my @lca;
    my $number = 0;
    my $c = 0; #number of genes total
    my $highest;
    my $value;
    my $prev;
    my $inside;
    my %final0;
    my @keep;
    foreach my $input (@retrieved) {
        my @info = split(/\t/, $input);
        $number++;
        if ($info[1] >= $cutoff_lca) { #give more weight when high identity score-cutoff
            my $new = $info[4]/($info[3]-$info[2]);
            $taxa1{$info[0]} += 2*($new); #taxa1: accession number = modified_bit_score;
            $taxa2{$info[0]} += 2*($info[1]); #taxa2: accession number = identity_score;
        } else {
            my $new = $info[4]/($info[3]-$info[2]);
            $taxa1{$info[0]} += $new; #taxa1: accession number = modified_bit_score;
            $taxa2{$info[0]} += $info[1]; #taxa2: accession number = identity_score;
        }
    }
    foreach my $key (sort {$taxa1{$b} <=> $taxa1{$a}} keys %taxa1) { #for each accession number (starts from highest identity score)
        if (defined $taxa{$accTotaxid{$key}}[0]) { #if accession number has an associated taxonomic path
            ${$final{$taxa{$accTotaxid{$key}}[0]}}[0] = ${$final{$taxa{$accTotaxid{$key}}[0]}}[0] . ' ' . $key; #final[1]: phylum = all accession numbers
            $final0{$taxa{$accTotaxid{$key}}[0]}{'bitscore'} += $taxa1{$key}; #sum of bit score for phylum; phylum = sum
            $final0{$taxa{$accTotaxid{$key}}[0]}{'identity'} += $taxa2{$key}; #sum of identity for phylum; phylum = sum
        }
    }
    foreach my $i (keys %final0) {
        ${$final{$i}}[1] = $final0{$i}{'bitscore'}*$final0{$i}{'identity'}; #fiNAL{PHYLUM} = BITSCORE*IDENTITY
    }
    my $start = 0;
    my $sum = 0;
    foreach my $high (sort {$final{$b}[1] <=> $final{$a}[1]} keys %final) {
        if ($start == 0) {
            $highest = $high;
            $value = $final{$highest}[1];
            $start++;
        } else {
            $sum += $final{$high}[1];
        }
    }
    if (($number == 1) or (($number == 2) && ($sum == 0)) or (($number > 2) && ($value > $sum))) { #in case one is higher
        push @lca, $highest;
        $prev = $highest;
        undef %final0;
        $inside = 0;
        foreach my $z (1..6) {
            @keep = split / /, ${$final{$highest}}[0]; #this is the array with the accession numbers that determined the assignement of the previous rank!
            shift(@keep);
            undef %final;
            undef %final0;
            foreach my $ke (@keep) { #for each acc num
                $final0{$taxa{$accTotaxid{$ke}}[$z]}{'bitscore'} += $taxa1{$ke}; #taxon{bitscore} = sum
                $final0{$taxa{$accTotaxid{$ke}}[$z]}{'identity'} += $taxa2{$ke}; #taxon{identity} = sum
                ${$final{$taxa{$accTotaxid{$ke}}[$z]}}[0] = ${$final{$taxa{$accTotaxid{$ke}}[$z]}}[0] . ' ' . $ke;
            }
            foreach my $i (keys %final0) {
                ${$final{$i}}[1] = $final0{$i}{'bitscore'}*$final0{$i}{'identity'};
            }
            $start = 0;
            $sum = 0;
            foreach my $high (sort {$final{$b}[1] <=> $final{$a}[1]} keys %final) {
                if ($start == 0) {
                    $highest = $high;
                    $value = $final{$highest}[1];
                    $start++;
                } else {
                    $sum += $final{$high}[1];;
                }
            }
            if (($inside == 0) && (($number == 1) or (($number == 2) && ($sum == 0)) or (($number > 2) && ($value > $sum)))) {
                push @lca, $highest;
                $prev = $highest;
            } else {
                if ($prev =~ /Unclassified/) {
                    push @lca, $prev;
                } else {
                    $prev = $prev . "_Unclassified";
                    push @lca, $prev;
                }
                $inside = 1; #if Lowest Common Ancestor
            }
        }
    } else { #more than two have same value
        foreach my $z (1..7) {
            push @lca, 'Unclassified';
        }
        undef %final0;
    }
    print $file "$old\t";
    if (defined($paired{$old})) {
        if ($lca[0] eq '') {
            @lca = ('Unclassified','Unclassified','Unclassified','Unclassified','Unclassified','Unclassified');
        }
        print $file join "\t", @lca, "\n";
        print $file "$old\t";
        print $file join "\t", @lca, "\n";
    } else {
        if ($lca[0] eq '') {
            @lca = ('Unclassified','Unclassified','Unclassified','Unclassified','Unclassified','Unclassified');
        }
        print $file join "\t", @lca, "\n";
    }
}
