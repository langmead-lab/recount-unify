#!/usr/bin/env bash
set -exo pipefail

gene2snaptron_samplesF=$1

#srav3h (0), gtexv2 (2), tcgav2 (3)
dataset_id=$2

cat $gene2snaptron_samplesF | perl -ne 'chomp; $f=$_; ($gid,$samps_)=split(/\t/,$f,-1); next if(length($samps_) == 0); @f0=split(/,/,$samps_,-1); shift(@f0); @coverages=(); @samples=(); for $samp_ (@f0) { ($sid,$cnt)=split(/:/,$samp_,-1); push(@coverages,$cnt); push(@samples,$sid); } $sample_count=scalar @samples; $cov_count=scalar @coverages; $cov_sum=0; map { $cov_sum+=$_; } @coverages; $avg=$cov_sum/$cov_count; $median=-1; $m_=$cov_count/2; $m_=int($m_); @f2a=sort {$a<=>$b} @coverages; if(($cov_count % 2)==0) { $m2_=$cov_count/2; $m2_-=1; $median=($f2a[$m_]+$f2a[$m2_])/2; } else { $median=$f2a[$m_]; } printf("$f\t$sample_count\t$cov_sum\t%.3f\t%.3f\t'$dataset_id'\n",$avg,$median);'
