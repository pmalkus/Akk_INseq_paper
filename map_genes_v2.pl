#!/usr/bin/perl -w

###This script takes the .ptt and the filter_cpm file for mapping insertions to the gene
###The format of the .ptt file is the standard annotation file format from NCBI, or you can refer to the .ptt file in the demo section

if (!@ARGV){
  print "$0 <.ptt> <.cpm file> <gene_disrupt_percent (max=1) \n";
  exit;
}

my $ptt=shift @ARGV;
my $file=shift @ARGV;
my $percent=shift @ARGV;
my $out=$file."_mapped";    


my $logfile;
if ($file=~/(\w+\.scarf)/){
   $logfile=$1.".log";
}

open LOG,">>$logfile";

my $gene;
my $start;
my $end;
my @temp;
my @gene;
my $name;
my $strand;
my $annotation;

open PTT,$ptt;
while (<PTT>){
    chomp;
    if (m/^(\d+)\.\.(\d+)/){
      @temp=split /\t/,$_;
      $start->{$temp[5]}=$1;
      $end->{$temp[5]}=$2;
      $strand->{$temp[5]}=$temp[1];
      $length->{$temp[5]}=abs($end->{$temp[5]}-$start->{$temp[5]});
      $annotation->{$temp[5]}=$temp[-1];
    }
}

foreach my $gene (sort keys %$start){
     for (my $i=$start->{$gene};$i<=$end->{$gene};$i++){
          $target->{$i}=$gene;
     }
}

my $count;
my $hit;
my $inside_ORF=0;
my $total_in=0;
open CPM,$file;
while (<CPM>){
     @temp=split /\s+/,$_;
     $total_in++;
     if ($target->{$temp[1]}){
         if ($strand->{$target->{$temp[1]}} eq "+"){
             $position=($temp[1]-$start->{$target->{$temp[1]}})/$length->{$target->{$temp[1]}};
         }elsif ($strand->{$target->{$temp[1]}} eq "-"){
             $position=($temp[1]-$end->{$target->{$temp[1]}})/$length->{$target->{$temp[1]}};
         }
         if ($position<=$percent){
           $inside_ORF++;
           $hit->{$target->{$temp[1]}}++;
           $count->{$target->{$temp[1]}}+=$temp[2];
         }
     }
}

open OUT,">$out";
print OUT"Gene ID\tUnique insertions\tCount\tAnnotation\n";

foreach $gene (sort keys %$start){
     if ($hit->{$gene}){
        print OUT "$gene\t$hit->{$gene}\t$count->{$gene}\t$annotation->{$gene}\n";
     }else{
        $hit->{$gene}=0;
        $count->{$gene}=0;
        print OUT "$gene\t$hit->{$gene}\t$count->{$gene}\t$annotation->{$gene}\n";
    }
}


$intergenic=$total_in-$inside_ORF;
print LOG "Sample\tTotal insertions\t Insertions inside of ORF\tInsertions located intergenic region\n";
print LOG "$file\t$total_in\t$inside_ORF\t$intergenic\n";
close LOG;
