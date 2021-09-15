#! /usr/local/bin/perl -w
# Author: Lalitha Viswanathan
# Affiliation: UMBI CARB
# Wrapper for newickModifier script; Connects to scop_treework db and pulls genome ids corresponding to set of orthologous cluster ids and in newickModifier, gets newick tree corresponding to
# each cluster and rewrites them 
# input is orthologous cluster ids 
# TBD : Where r the ids coming from 

use lib "/home/viswanat/perllib";
use strict;
use DB;
use Family;
use Tree;
use Node;
use Daily;
use Data::Dumper;
use Getopt::Long;

# connect to scop_treework database
my $dbh = &scop_treework();
my $clusterAc;
my (%genomeAbbrs, %options);
my $FH;
my $InfileName;
GetOptions(\%options,"clusterac=s");


$clusterAc = $options{"clusterac"};

# query to get ortholog genome Ids
#my $query = $dbh->prepare("select distinct genome_abbr from all_ortho_members_daniel B where B.cluster_ac  = $clusterAc and B.subfamily=1");

# query to get non-ortholog genomeIds
#my $query = $dbh->prepare("select distinct genome_abbr from all_ortho_members_daniel B where B.cluster_ac = $clusterAc and B.subfamily = (select max(E.subfamily) from domains_rpsblast_hits A, scop_1_69_detailed D,all_ortho_members_daniel E where A.cluster_ac = $clusterAc and A.cluster_ac = E.cluster_ac and A.identity > 95  and A.seq_ac2 = D.sid and A.domain_gi not in( select distinct seq_ac from all_ortho_members_daniel where cluster_ac = $clusterAc and subfamily = 1) and E.seq_ac = A.domain_gi and E.subfamily !=1)");

my $query = $dbh->prepare("select distinct genome_abbr from all_ortho_members_daniel B where B.cluster_ac = $clusterAc and B.subfamily = 24");
$query->execute();

$InfileName = $clusterAc."_genomeAbbrs_paralogs24.txt";
open(FH, ">$InfileName") ||  die "Error opening $InfileName $!\n";
while(my ($abbr) = $query->fetchrow_array) {
	print "$clusterAc $abbr\n";
	print FH $abbr,"\n";
}
close(FH);

system("perl newickModifier.pl -file=/home/viswanat/DSSP/completeNames_newtree.dnd -genomeIds=$InfileName");
