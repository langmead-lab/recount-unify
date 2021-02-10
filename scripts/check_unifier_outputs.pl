#!/usr/bin/env perl
use strict;
use warnings;

#e.g. exon_sums_per_study/42/SRP044042/sra.exon_sums.SRP044042.G026.gz
my $input_file_name = shift;
my $counts_per_study_file = shift;
#gene, exon, or junction(?)
my $ftype = shift;
#e.g ${GENE_EXON_ANNOTATION_ROW_COUNTS}
my $row_counts_file = shift;

#load expected column (sample) counts
my %counts_per_study;
open(FIN,"<$counts_per_study_file");
while(my $line = <FIN>) 
{
    chomp($line);
    my ($count,$study)=split(/\t/,$line);
    $counts_per_study{uc($study)}=$count;
}
close(FIN);

#load expected row (gene or exon) counts
my %counts_per_annotation;
open(FIN,"<$row_counts_file");
while(my $line = <FIN>) 
{
    chomp($line);
    my ($atype,$annot,$row_count)=split(/\t/,$line);
    #gene, exon, or junction
    if($atype eq $ftype)
    {
        $counts_per_annotation{uc($annot)}=$row_count;
    }
}
close(FIN);

#$line=~/^[^\/]+\/..\/([^\/]+)\/.*\.gz$/;
#$line=~/^[^\/]+\/..\/([^\/]+)\/([^\.]+)\.([^\.]+).([^\.]+).([^\.]+).gz
my @path_segment = split(/\//,$input_file_name);
my $dir = shift(@path_segment);
#low order of study
shift(@path_segment);
my $study = shift(@path_segment);
$study = uc($study);
my $base = pop(@path_segment);
my $ERROR=undef;
if(!defined($counts_per_study{$study}))
{
    print STDERR "ERROR\t$study not in $counts_per_study_file\n";
    $ERROR=1;
}

my $study_count = $counts_per_study{$study};
my $num_rows_expected;
if($ftype=~/(exon)|(gene)/i)
{
    my ($source,$type,$study,$annot)=split(/\./,$base);
    $num_rows_expected = $counts_per_annotation{uc($annot)};
}

my $last_line;
my $num_rows=0;
while($last_line = <STDIN>)
{
    $num_rows++;
    chomp($last_line);
    #if there are empty strings (blanks) in one or more columns in the last line, this will pick it up
    my @cols = split(/\t/,$last_line);
    my $num_cols = scalar(@cols);
    ##check number of columns (samples) in the *last line* of the file
    #my $num_cols=`zcat $line | tail -n1 | tr \\t \\n | wc -l`
    #subtract one for first column (not a sample for genes/exons)
    $num_cols--;
    if($study_count != $num_cols)
    {
        print STDERR "ERROR\texpected column count:$study_count != column count:$num_cols\tline#$num_rows\t$last_line\n";
        $ERROR=1;
    }
}

##check number of rows (lines) in the file
#my $num_rows=`zcat $line | cut -f 1 | wc -l`;
#subtract one for first row (header)
$num_rows--;
if($num_rows_expected != $num_rows)
{
    print STDERR "ERROR\texpected row count:$num_rows_expected != row count:$num_rows\n";
    $ERROR=1;
}

print $input_file_name."\t";
if(!defined($ERROR))
{
    print "NUM_ROWS_COLUMNS_AND_NO_BLANKS_CHECKS_PASSED\n";
}
else
{
    print "CHECKS_FAILED\n";
    exit -1;
}
