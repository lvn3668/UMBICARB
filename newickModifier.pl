# AAuthor: Lalitha Viswanathan
# Affiliation: UMBI CARB
# Utility script to annotate neewick format phylogenetic tree with genome abbreviations where the Newick tree does not contain branch lengths
#!/usr/bin/perl
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use lib "/home/viswanat/perllib";
use Data::Dumper;
use Parallel::ForkManager;
use Getopt::Long;
use DB;

my ($treeFile, $ClusterIds, $ModifiedtreeFile, $dbh, $query, @fullName, %options);


# connect to scop_treework database
$dbh = &scop_treework();

# get all the commands from the command line
# Input is orthologous cluster id file and corresponding newick tree dendogram file
# All genomes in that cluster file are read and the corresponding abbreviations pulled from db 
# the dendogram is written out with appropriate substitutions 
GetOptions(\%options, "file=s", "genomeIds=s");

my ($FH, $ORTHOLOGOUSCLUSTERfilehandle, @genome_abbrs, $tmp);

$ClusterIds = $options{"genomeIds"};
$treeFile = $options{"file"};

open(ORTHOLOGOUSCLUSTERfilehandle, $ClusterIds) || die "Error opening $ClusterIds $! !!! \n";

while(<ORTHOLOGOUSCLUSTERfilehandle>) {
	push(@genome_abbrs, $_);
}

$ModifiedtreeFile =  $treeFile.$ClusterIds.".modified.dnd";

# open file for rewriting the tree
open(FWRITE, ">$ModifiedtreeFile") || die "Cannot open $ModifiedtreeFile for writing! \n";

# open the newick file format n substitute it with the genome abbr appended with
# ***

open(FH, $treeFile) || die "Error opening $treeFile! $! \n";
while (<FH>) {
	foreach my $abbr (@genome_abbrs) {
		chomp($abbr);
		$abbr =~s/\s+//g;
		$query = $dbh->prepare("select genome_name from genome_details_from_clusterdb where genome_abbr='$abbr'");
		$query->execute();
		@fullName = $query->fetchrow_array;
		chomp($fullName[0]);
		
		if(!$fullName[0]) {
		print "Error for $abbr \n";
		exit(1);
		}

		# Appease tree rendering software
		# Input is a Newick tree with NO branch lengths and newick-format-node-names (input tree contains whole gneome names)
		# tbd: how is this tree generated 
		# remove braces and underscores; resulting Newick tree contains genome abbreviations as substitutions 
		# write out to file 
		$fullName[0] =~s/\(//g;
		$fullName[0] =~s/\)//g;
		$fullName[0] =~s/:/_/g;

		$tmp = "_".$fullName[0]."_";
		print "$abbr $fullName[0]  \n";	
		$_ =~s/$fullName[0]:/$tmp:/g;
	}
	print $_;
	print FWRITE $_;
	}
close(FH);
close(FWRITE);

=pod 
reads in a newick dendogram file and writes out filename.modified.dendogram file with genome abbrs as used in alignment files 
=cut 
