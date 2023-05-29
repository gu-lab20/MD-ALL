#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts=();
getopts( 'f:b:h:o:', \%opts );
die(
	qq/
Usage:   *.pl [options]
Description:
	Calculate reads count according to bed file regions.
Options: 
	-f STRING	input bed file name [required]
	-b STRING	input bam file name [required]
	-h STRING	input HTSeq file name [required]
	-o STRING	output new HTSeq file name [required]
\n/
)if ( !$opts{f} || !$opts{b} || !$opts{h} || !$opts{o});

my ( $fbed, $fbam, $fhtseq, $fout ) = ($opts{f}, $opts{b}, $opts{h}, $opts{o});

open(IN,"<", $fbed) || die ("Could not open file $fbed: $!\n");
my %HTSeq=();
my $DUX4_count=0;
while(<IN>){
	my ($chr, $start, $end, $name, $geneType, $ENSG)=split;
	my @reads=qx(samtools view $fbam $chr:$start-$end);
	my %reads=();
	for(@reads){
		my ($ID)=split;
		$reads{$ID}=1;
	}
	my $readsCount=scalar(keys %reads);
	$HTSeq{$ENSG}=$readsCount if ($readsCount > 0);
    $DUX4_count++;
}
close IN;

open(IN,"<", $fhtseq) || die ("Could not open file $fhtseq: $!\n");
open(OUT,">", $fout) || die ("Could not create file $fout: $!\n");
my @tailInfor=();
while(<IN>){
	if(/^__/){
		push(@tailInfor,$_);
		next;
	}
	chomp;
	my ($ENSG, $count)=split;
	if(exists $HTSeq{$ENSG}){
		my $newCount=$HTSeq{$ENSG};
		$count = $newCount if ($count < $newCount);
		delete($HTSeq{$ENSG});
	}
	print OUT "$ENSG\t$count\n";
}
close IN;

for my $k (sort keys %HTSeq){
	print OUT "$k\t$HTSeq{$k}\n";
}

map{print OUT $_} @tailInfor;
print OUT "__DUX4\t$DUX4_count\n";
close OUT;

