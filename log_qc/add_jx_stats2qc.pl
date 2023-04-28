#!/usr/bin/env perl
#joins log QC stats together with jx stats and formats sets up proper study/sample column header names
use strict;
use warnings;

my $samples_file=shift;
my $skip_jxns=shift;

my %jx_stats;
my %rail_ids;
open(FIN,"<$samples_file");
while(my $line = <FIN>) 
{
    chomp($line);
    #rail_id	run	study...junction_count	junction_coverage	junction_avg_coverage
    #837730	ERR188431	ERP001942...26821	141957	5.29275567652213
    my @fields=split(/\t/,$line);
    #run:study
    $fields[1] = "external_id" if($fields[1]=~/^((sample_id)|(run)|(external_id))$/);
    $fields[2] = "study" if($fields[2]=~/study_id/);
    my $key=$fields[1]."\t".$fields[2];
    $rail_ids{$key} = $fields[0];
    next if(defined($skip_jxns) && $skip_jxns == 1);
    my $javg=pop(@fields);
    my $jcov=pop(@fields);
    my $jcnt=pop(@fields);
    $jx_stats{$key} = "$jcnt\t$jcov\t$javg";
}
close(FIN);

my $num_keys = scalar(keys %rail_ids);
die "no runs/samples in samples.tsv, terminating\n" if($num_keys < 2);

while(my $line = <STDIN>)
{
    chomp($line);
    #ERP001942   ERR188431\t<log qc stats_pre_adding_both_stats>
    my @fields = split(/\t/,$line);
    my $study = shift(@fields);
    my $run = shift(@fields);
    my $newline = join("\t",@fields);
    $run = "external_id" if($run eq "sample");
    my $key = $run."\t".$study;
    die "couldnt find entry in $samples_file for $key, terminating\n" if(!defined($rail_ids{$key}));
    my $rail_id = $rail_ids{$key};
    if(defined($skip_jxns) && $skip_jxns == 1) {
        print "$rail_id\t$run\t$study\t$newline\n";
        next;
    }
    my $jxs = $jx_stats{$key};
    print "$rail_id\t$run\t$study\t$newline\t$jxs\n";
}
