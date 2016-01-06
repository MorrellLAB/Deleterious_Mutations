#!/usr/bin/perl
#2015/11/11 by Li Lei
#compare the genotypes I called from the bam files and the one from the Alchemy files;
#usage:  ./cmpr_alchemy_bamcallSNP.pl ALCHEMY.vcf snp_called_from_BAM.vcf >compared_genotype_Alchemy_vcf.txt
use strict;
use warnings;
#use Data::Dumper;
my ($ALCHEMY, $VCF) = @ARGV; #$ALCHEMY is ALCHEMY data; $VCF is the SNP data I called from the BAM file;
open(IN1, $ALCHEMY) or die "Could not open this $ALCHEMY";
my $header =<IN1>; #skip the head line
my %hash; #define a hash to store the id key and genotype obtained from Alchemy data;
foreach my $row (<IN1>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $id = $rtemp[0]."_".$rtemp[1]; #make the unique id key for hash, which is made of sample id and SNP id;
        $hash{$id}=$rtemp[5];
}

close IN1;
#print Dumper(\%hash);
open(IN2, $VCF) or die "Could not open this $VCF";
my $idt_SNPs=0; #define identical SNPs
my $rc_SNPs=0; #define reverse complementary SNPs
my $diff=0; #define different SNPs (neither identical nor reverse complementary)
while(<IN2>) {
	next if /^\#/; #skip the header with "#" in the beginning
	chomp $_;
	my @rtemp = split(/\t/,$_);
	my $key = $rtemp[0]."_".$rtemp[1];#make the unique id key, which is made of sample id and SNP id;
	#print "$key\n";
	my @temp1 = split(/\:/,$rtemp[9]);
	my $geno_code = $temp1[0];
	my $genotype;
  #if...elsif..else structure to define the genotype based on the SNP calling and the genotype code in vcf file;
	if($geno_code eq "1/1"){#both alleles are same with alternative allele
       $genotype = $rtemp[4].$rtemp[4];
	}
	elsif($geno_code eq "0/1"){#one is same with reference allele, and the other is same with alternative alleles
       $genotype = $rtemp[3].$rtemp[4];
	}
	else{#both alleles are the same with reference allele
       $genotype = $rtemp[3].$rtemp[3];
	}
	#print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\n";
   	if (exists $hash{$key}){
   		print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\n";
   		my $rc_geno = &reverse_complement($genotype); #get the reverse compementary genotype
   	    if ($hash{$key} eq $genotype){##get the identical genotypes with Alchemy data;
   	    	$idt_SNPs++;
   	    	print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\tIDT\n";
   	    }
   	    elsif($hash{$key} eq $rc_geno){##get the reverse complementary genotypes with Alchemy data;
    	    $rc_SNPs++;
    	    print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\tRC\n";
   	    }
   	    else{##get the different genotypes from Alchemy data;
   	    	$diff++;
 	        print "$rtemp[0]\t$rtemp[1]\t$rtemp[3]\t$rtemp[4]\t$rtemp[9]\t$genotype\t$hash{$key}\tNO\n";
   	    }
   }
}

print "$idt_SNPs\t$rc_SNPs\t$diff\n";
close IN2;

#define a subrutine "reverse_complement" to obtain the reverse complementary sequence;
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
