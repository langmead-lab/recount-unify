#!/usr/bin/env perl
use strict;
use warnings;

my $counts_per_study_file = shift;
#genes, exons, or junctions
my $ftype = $shift;
#e.g ${GENE_EXON_ANNOTATION_ROW_COUNTS}
my $row_counts_file = shift;

my %counts_per_study;
open(FIN,"<$counts_per_study_file");
while(my $line = <FIN>) 
{
    chomp($line);
    my ($count,$study)=split(/\t/,$line);
    $counts_per_study{uc($study)}=$count;
}
close(FIN);

my %counts_per_annotation;
open(FIN,"<$row_counts_file");
while(my $line = <FIN>) 
{
    chomp($line);
    my ($annot,$row_count)=split(/\t/,$line);
    $counts_per_annotation{uc($annot)}=$row_count;
}
close(FIN);

my $ERROR=undef;
while(my $line = <STDIN>)
{
    chomp($line);
    #e.g. exon_sums_per_study/42/SRP044042/sra.exon_sums.SRP044042.G026.gz
    #$line=~/^[^\/]+\/..\/([^\/]+)\/.*\.gz$/;
    #$line=~/^[^\/]+\/..\/([^\/]+)\/([^\.]+)\.([^\.]+).([^\.]+).([^\.]+).gz
    my @path_segment = split(/\//,$line);
    my $dir = shift(@path_segment);
    #low order of study
    shift(@path_segment);
    my $study = shift(@path_segment);
    $study = uc($study);
    my $base = pop(@path_segment);
    if(!defined($counts_per_study{$study}))
    {
        print STDERR "ERROR\t$study not in $counts_per_study_file\n";
        $ERROR=1;
    }
    my $study_count = $counts_per_study{$study};
    if($ftype=~/(exons)|(genes)/i)
    {
        my ($source,$type,$study,$annot)=split(/\./,$base);
        my $num_cols=`zcat $line | tail -n1 | tr \\t \\n | wc -l`
        #subtract one for first column (not a sample for genes/exons)
        $num_cols--;
        if($study_count != $num_cols)
        {
            print STDERR "ERROR\texpected column count:$study_count != column count:$num_cols in $line\n";
            $ERROR=1;
        }
        my $num_rows=`zcat $line | cut -f 1 | wc -l`;
        #subtract one for first row (header)
        $num_rows--;
        if($num_rows_expected != $num_rows)
        {
            print STDERR "ERROR\texpected row count:$num_rows_expected != row count:$num_rows in $line\n";
            $ERROR=1;
        }
    }
}
print "CHECKS_DONE\n";
