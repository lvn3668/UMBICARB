#! /usr/local/bin/perl -w
# Author Lalitha Viswanathan
# Affiliation UMBI CARB
# Musical Genomes with Noise in residuals calculation 
use lib "$ENV{HOME}/perllib"; 
use DB;
use Family;
use Tree;
use Node;
use Daily;
use Data::Dumper;
use pairwiseblastparser;
use treedist;
use math;
use databaseaccessors;
use utilities;
use readerwriter;

my $A = 10;
my $dbHandle = &scop_treework();

my %all_genome_ids = &get_all_genome_abbrs($dbHandle);
my %clusters = &get_cluster($dbHandle);
my %clockrates = &get_clockrates($dbHandle);
my $internodedistref = &get_internode_dist($dbHandle);
&nonuniquefy($internodedistref);

#open(FH, ">musicalgenomeswithnoiseoutput.txt");

my @tobeignored = ("yst", "worm", "fly", "pNGR234", "not_defined", "");
	my (%new_cluster_distances, $largest_distance);
	$largest_distance = 0;
	foreach my $clusterId(sort {$a <=> $b} keys %clusters){

	print "****ClusterId $clusterId \n";
	# need to get hash of ss per site genome-genome values 
	my $sspersitestring = &extract_column_from_row($clusterId, \%clusters, 1);
	# convert it to hash 
	my $sspersitehashref = &stringtohash($sspersitestring);
	# do substitutions specific to the ssper site protdist result
	# protdist returns 10 characters long labels
	# so from 999_1_ccre extract ccre
	# call function which does that
	$sspersitehashref = &trimhash($sspersitehashref);
	# Find phylogenetic tree per cluster with lowest residual with residual calculated with noise level
	# Lowest residual calculated by moving every genome in phylogenetic tree to every other location and finding residual (using clock rates / rate of evolution)
	&musicalgenomes($clusterId, $sspersitehashref);
	print "****Created music with ClusterId $clusterId \n";
}

sub musicalgenomes() {
	my $clusterid = shift || die "This function needs clusterid \n";
	my $sspersitehashref = shift || die "This function needs sspersite hasref \n";

	my $stmt;
	my $clockrate = $clockrates{$clusterid};
	foreach my $genomeId (keys %$sspersitehashref) {
			if($genomeId ne "yst") {
		foreach my $allIds (keys %$internodedistref) {
			if($allIds ne "yst") {
				for(my $noise = 0; $noise < 1; $noise +=0.1) {
					print "Computing residual for $genomeId at position $allIds with Noise level $noise \n";	
			$stmt = join("\t",$clusterid, $genomeId, $allIds);
			$residual = 0;
			foreach my $innergenomeId (keys %$sspersitehashref) {

				if($innergenomeId ne"yst"){
				if($genomeId ne $innergenomeId) {
		print "Current value of residual  $residual \n";
		print "Value1: Distance->{$innergenomeId}->{$allIds} = $internodedistref->{$innergenomeId}->{$allIds} \n";
	       print  "Value2: Distance->{$genomeId}->{$innergenomeId} = $internodedistref->{$genomeId}->{$innergenomeId} \n";
		print "After adding square of (Value1-(1-noise)*Value2) to residual \n";  
			$residual += ($internodedistref->{$innergenomeId}->{$allIds} - ((1-$noise)*$internodedistref->{$genomeId}->{$innergenomeId}))**2;
			} 
		}
	}	
	$stmt = join("\t",$stmt, $residual,"\n");
	print "Residual for $genomeId at position $allIds with noise level $noise is $residual \n";
	print $stmt;
	}	
				}
			}
		}
	}
}

#close(FH);
