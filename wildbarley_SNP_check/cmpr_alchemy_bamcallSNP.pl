#!/usr/bin/perl
#20151111 by Li Lei
#compare the genotypes I called from the bam files and the one from the Alchemy files;
#usage:  
use strict;
use warnings;
use Data::Dumper;
my ($ALCHEMY, $VCF) = @ARGV;
open(IN1, $ALCHEMY) or die "Could not open this $ALCHEMY";
my $header =<IN1>;

my %hash;
foreach my $row (<IN1>){
chomp $row;
my @rtemp = split(/\t/,$row);
my $id = $rtemp[0]."_".$rtemp[1];
   $hash{$id}=$rtemp[5];
}

close IN1;
#print Dumper(\%hash);
open(IN2, $VCF) or die "Could not open this $VCF";
my $idt_SNPs=0;
my $rc_SNPs=0;
my $diff=0;
while(<IN2>) {
	next if /^\#/;
	chomp $_;
	my @rtemp = split(/\t/,$_);
	my $key = $rtemp[0]."_".$rtemp[1];
	#print "$key\n";
	my @temp1 = split(/\:/,$rtemp[9]);
	my $geno_code = $temp1[0];
	my $genotype;
	if($geno_code eq "1/1"){
       $genotype = $rtemp[4].$rtemp[4];
	}
	elsif($geno_code eq "0/1"){
       $genotype = $rtemp[3].$rtemp[4];
	}
	else{
       $genotype = $rtemp[3].$rtemp[3];
	}
	#print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\n";
   	if (exists $hash{$key}){
   		print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\n";
   		my $rc_geno = &reverse_complement($genotype);
   	    if ($hash{$key} eq $genotype){
   	    	$idt_SNPs++;
   	    	print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\tIDT\n";
   	    }
   	    elsif($hash{$key} eq $rc_geno){
    	    $rc_SNPs++;
    	    print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\tRC\n";
   	    }
   	    else{
   	    	$diff++;
 	        print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\tNO\n";
   	    }
   }
}

print "$idt_SNPs\t$rc_SNPs\t$diff\n";
close IN2;

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
