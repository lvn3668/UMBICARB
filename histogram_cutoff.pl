#! /usr/local/bin/perl 
# Author: Lalitha Viswanathan
# Affiliation: UMBI CARB
# Finds horizontal gene tranfer in each cluster of orthologous genes 
# (one representative gene per genome ; with analyzes being perrformed on upto 316 microbial genomes)
# Universal tree refers to tree of 316 genomes ; 
# select cluster, number of unique genomes in cluster, distance (phylogenetic distance in universal tree) and loss (caluclated using clockrates / evolutionary distance of representative genomes in that cluster)
# Plot histogram of cluster, number of LGTs, loss rate and size of cluster 
# Wrapper for hgt perl module 
use strict;
use DB;
use Data::Dumper;

# program is to measure the number of LGTs in terms of different loss and size cutoff
# Perl Modules misisng 
my $dbh = &hgt();

my $TABLE = "complete_genome_cluster_loss_universal_tree";
my $sth1 = $dbh->prepare("select cluster_ac, family_size, sum_distance, loss from $TABLE");

my %count;
my %calc;
$sth1->execute();
my $MAX_size = 0;
my $MAX_loss = 0;
while (my ($cluster,$size, $distance, $loss) = $sth1->fetchrow_array){
	$count{$loss}{$distance}++;
	$calc{$loss}{$size}++;
	$MAX_size = $size if $MAX_size < $size;
	$MAX_loss = $loss if $MAX_loss < $loss;
}

my %total;
for(my $i = $MAX_loss; $i >= 0; $i--){
	for(my $j = $MAX_size - 1; $j> 1; $j--){
		$calc{$i}{$j} = 0 if not defined $calc{$i}{$j};
		for(my $k = $i; $k <= $MAX_loss;$k++){
			for(my $l=$j; $l<=$MAX_size; $l++){
				$total{$i}{$j} += $calc{$k}{$l};
			}
		}
		print "$i $j $total{$i}{$j}\n";
	}
}



