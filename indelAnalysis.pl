#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use warnings;
use diagnostics;
use lib "/home/viswanat/perllib";
use FileRead;
use Unmammoth;
use lib "/home/viswanat/bioperl";
use Bio::AlignIO;
use DB;
use Family;
use Tree;
use Node;
use utilites;

# Author: Lalitha Viswanathan
# Affiliation: UMBI CARB 
# Script to classify gaps in sequence alignments of orthologous microbial genomes using location of  corresponding structure hits 

my %options;
my $structFileName;
my $seqFileName;
my $DIR = "/home/viswanat/DSSP";
my $MYIN;
my @pdbIds;
my $pdbIdStrings;
my @sequenceIdsWithPdbHits;
my %hashAlign;
my @seqAnalysis;
my @gapdistribution;
my $counter;
my @eventsdetails;
my @genomeIds;
my $genomeIdcounter=0;
#######################################################
# takes in a mammoth output structure alignment file
#  and a clustalw sequence alignment file
# The goal of this code presently is to do the following
#  take a pairwise structure alignment
#  get pdbIds from that structure alignment
#  get corresponding domain entries from the domains_rpsblast_hits database and extract only those things
#   from sequence alignment 
#  Kick out the rest of the sequence alignment
#  for each cluster of orthologous seuqencess, a matrix of pairwise seuqence alignments and pairwise structure alignemnt is stored as clobs and
#  domain_rps_blast stores results of (scop domans used to build query database) and queried against structure hits (corresponding to sequences in orthologous cluster)
#  Cluster of 50 seuqences will have 2500 aln files (times 2; sequence and structure) as i/p to this script
#  After eliminating seuqences not present in results of rps blast (doman rps blast), the residuals are pairwise structure alignment and a part of the sequence alignment

#  Starting with structure hits , go back to sequence alignments to correlate gaps and classify the same 
#  Sequence alignment is represented as 0s and 1s 
#  Concatenate each column of the alignment as an identifier 
#  Then concat each such COLUMN of 0s and 1s as an identifier for that  
#  Every digit in column-wise-concatenated-string corresponds to a genome in the orthologous cluster
#  At this stage,
#  arrayelement1 -> id for col1 of alignment -> Id1 (all genomes in id1 with 0s or 1s as appropriate)
#  array element2 -> id for col2 of alignment ->Id2 
#  col3->id for col3 of alignment -> Id3
#
#
#  Take transpose of this array
#  Thus the genomes become the column information and the alignment itself (the Ids for each row of the alignment)
#  become the rows of the array
# 
#  Genome by Alignment matrix on transposition becomes Alignment by Genome matrix 
#  Eliminate duplicate rows (with the number of deleted rows corresponding to the length of the deletion) ; (merge the gaps step)
#  The relative number of 0s and 1s classifies the gap as an insertion or a deletion 
#  For equivalent number of deletions and insertions (or within a certain delta), classification is undetermined  
#  Then we are left with a set of Ids all unique and another array that stores how long each Id is (each indel)
#  Moving along x axis, determine genomes where a consensus deletion is otherwise, and viceversae (a consensus deletion is a segment where 0s outnumber 1s) ; repeat for consensus insertion     
#  say for an indel classified as deletion (more 0s than 1s in the ID), find which columns (genomes) have 0s in them i.e. are involved/affected by deletion
#  For the list of THESE genomes, get an age from the tree ; the age and location on tree determines outliers and hence lateral gene transfer
#
# TCOFFEE and LALIGN seuqence alignment o/p files [???]
#  Not sure where the clustalw alignment file is coming in from (determine delta within which misclassification is to be addressed and use age tree to identify LGTs in a given cluster)

my $TREE_FILE = "home/viswanat/remote_homolog_detection_project/tree_files/14family_tree.dnd";
my $tree = new Tree({file => $TREE_FILE});
GetOptions(\%options,"struct=s","seq=s");
$structFileName = $options{"struct"};
$seqFileName = $options{"seq"};

chomp($structFileName);
chomp($seqFileName);

print "Structure File : $structFileName \n";
print "Sequence File : $seqFileName \n";


##########################################################

$structFileName = $DIR . "/".$structFileName;
my $structureAlignment = new Unmammoth($structFileName);
$structureAlignment->read_alignment();
$structureAlignment->parse();

@pdbIds = @{$structureAlignment->getPdbIds()};

foreach my $pdbId (@pdbIds) {

# convert into lower case
$pdbId = lc($pdbId);

# remove the underscore before the chain Id 
$pdbId =~s/_//g;

# add a 'd' to it to make it d.....
$pdbId = 'd'.$pdbId;

# append an underscore after the chain Id
$pdbId .= '_';

}

$pdbIdStrings = '(\''. join ('\',\'', @pdbIds) . '\')';
my @secondarystructurefasta =  $structureAlignment->getssFasta(); 
my @aminoacidfasta =  $structureAlignment->getaaFasta(); 


# print secondary structure alignment
print "AA Align: \n @aminoacidfasta";
print "\nSS Align: \n @secondarystructurefasta";


my $seqObject = new Bio::SimpleAlign(join('',@secondarystructurefasta));


###################################################################
# after getting the pdbIds, the structure alignment 
# now get the seqIds corresponding to this alignment
# so that we can get only that fragment of the alignment from
# the SEQUENCE alignment

# get the struture homologs where structure hits to scop databsae (using rpsblast) > 95% 
#connect to scop_treework database
my $dbh = &scop_treework();

# get the right query and execute it
my $query = $dbh->prepare("select distinct concat(cluster_ac,'_',domain_gi,'_', genome_abbr) 
from domains_rpsblast_hits A  where seq_ac2 in $pdbIdStrings and identity >=95");
$query->execute();

# then print each result returned
while (my $id = $query->fetchrow_array) {
push(@sequenceIdsWithPdbHits, $id);
}

############################################################

$seqFileName = $DIR . "/".$seqFileName;

# read in the sequence alignment and print it
open MYIN, $seqFileName;
my $clustalwReader  = Bio::AlignIO->newFh(-format => 'clustalw', 
          -fh => \*MYIN);
my $out = Bio::AlignIO->newFh(-fh =>\*STDOUT, 
       -format => 'fasta');

# then get those seqIds, 
# then look for them in the sequence alignment
# then print only those sequences from the sequence alignment
# whose id matches the structure based id

# Correlate sequence ids from orthologous clusters (sequence clustering done using Clustalw) withh pdb ids from scop chains database
$counter=0;
my $aln;
while($aln = <$clustalwReader>) {
my $tmpId;
foreach my $seq ($aln->each_seq()) {
 print "genome id is ",$seq->id() ,"\t";

 $tmpId = $seq->id();
 $tmpId =~s/_/#/;
 $tmpId =~s/_/#/;
 #$tmpId =~s/_/#/;

 if($tmpId =~/(\S+)#(\S+)#(\S+)/) {
 $genomeIds[$genomeIdcounter] = $+; 
 print "genome id: $genomeIds[$genomeIdcounter] \n";
 $genomeIdcounter++;
 }
 $hashAlign{$seq->id()} = $seq->seq();
 my $s = $seq->seq();
 chomp($s);
 $s =~s/\s+//g;
 my @tmp = split("", $s);
 my $gapcounter=0;
 # At this stage, the sequence alignment is reduced to those sequences for which a structural homolog is present in the mammoth alignment file
 # In the resulting sequence alignment, convert gaps to 0s and bases to 1s 
 foreach $_ (@tmp) { 
  if($_ eq "-") {
  $gapdistribution[$counter][$gapcounter++] = 0;}  else 
  {
  $gapdistribution[$counter][$gapcounter++] = 1; 
  }
 }
 $seqAnalysis[$counter++] =  \@tmp;
} 
}
# get a matrix of 0s and 1s based on where the gap is (0) or isn't(1)
&IdEachColumnOfGapMatrix(\@gapdistribution);
@gapdistribution = &getTranspose(\@gapdistribution);
#
# thus after this we get a line of numbers all transposed
# from which we can merge duplicates
# i.e. if there is a deletion 10 residues long  in the alignment, delete 9 rows  from the TRANSPOSED ARRAY of IDS and keep track of number of rows deleted 
# i.e. uniquefy the indel events in the array  along with length of indel
# storing how many "duplicates" deleted (where duplicates are defined as repeating sets of gaps; 00011100110101000111 uniquefies as 0101010101 with 2 repeats ("000" repeats 2 times) and "00" repeating once ; number of duplicates = 2; number of unique deletions = 3 ("000", "00", "0")) ; process repeats for base pairs (marked as 1) 
# Then split each unique indel into individual digits and get genomes involved in indel and find their age
my $newarrref = &uniquelyIdentifyIndels(\@gapdistribution, \@eventsdetails);
&printMatrix($newarrref);

############################################################################
# analyze for gaps
# assign an Id to each COLUMN of gap matrix as concatenated version of 0s and 1s in all the rows of a given colum 
# in the gap matrix
sub IdEachColumnOfGapMatrix() {

my $matrix = shift || die "Gap analysis requires matrix as input! \n";
my $rows = scalar(@$matrix);
my $cols = scalar(@{$matrix->[0]});
my $idstring;
for(my $colInAlignment = 0 ; $colInAlignment  < $cols; $colInAlignment++) {
 $idstring = "";
 for(my $rowInAlignment = 0; $rowInAlignment < $rows; $rowInAlignment++) {
   $idstring .=$matrix->[$rowInAlignment]->[$colInAlignment]; 
 }
 $matrix->[0]->[$colInAlignment] = $idstring;
}
return $matrix;
}

############################################################################
sub uniquelyIdentifyIndels() {
my $TransposedArrayOfRowIdsForAlignment = shift;
my $IndelEventsLength = shift;
my $rows = scalar(@$TransposedArrayOfRowIdsForAlignment);
my $ColNumber = scalar(@{$TransposedArrayOfRowIdsForAlignment->[0]});
my $loop = 0;
my $eventscounter = 0;
my $lengthofindel = 0;
my $numberofdeletedrows = 0;
my @UniqueIndelEvents = ();
my $UniqueIndelEventsCounter=0;
my ($ones, $zeroes);
print "For merging duplicates: rows $rows columns $ColNumber \n";
while ( $loop < $rows-1) {

# keep going through id uniquely identifying each row of the alignment, if the event persists 
# then keep increasing lengthofindel and numberofdeletedrows
# i.e. do something different when you encounter something different
# i.e. when ID changes from one row to the other, it is not the same event persisting
# Thus we are essentially comparing every row of the array with the row above and the row below 
# An event is defined as a change  from a gap to an alignment (either match or mismatch) and vice versae 
if( ($TransposedArrayOfRowIdsForAlignment->[$loop]->[$ColNumber-1]) &&  
 ($TransposedArrayOfRowIdsForAlignment->[$loop+1]->[$ColNumber-1]) && 
 ($TransposedArrayOfRowIdsForAlignment->[$loop]->[$ColNumber-1] eq $TransposedArrayOfRowIdsForAlignment->[$loop+1]->[$ColNumber-1]) ) {
$lengthofindel++;
$numberofdeletedrows++;
}else{

 # if something is not persisting
 # then expand that ID into separate columns 0011110000 becomes 0 0 1 1 1 1 0 0 0 0  
for(my $ctr = 0; $ctr <scalar(@{$TransposedArrayOfRowIdsForAlignment->[$loop]})-1; $ctr++) {
$UniqueIndelEvents[$UniqueIndelEventsCounter][$ctr] = $TransposedArrayOfRowIdsForAlignment->[$loop]->[$ctr];
}
$UniqueIndelEvents[$UniqueIndelEventsCounter][$ColNumber] = '-';
$UniqueIndelEvents[$UniqueIndelEventsCounter][$ColNumber+1] = '-';
$IndelEventsLength->[$eventscounter] = $lengthofindel+1;

# reset length of indel 
$lengthofindel = 0;
# now decide whether it is a deletion, insertion or neither
# find no. of 1s, no. of 0s in the ID string

# Count number of zeros and ones
# if deletions (0s) greater than 1s then add "D" to beginning of indel length
# if deletions(1s) greater than 0s, add "I" to beginning of indel length
# if both are equal, call it undecided
# Thus IndelEventsLength now contains entries like I15 or D14 etc etc
$ones = $TransposedArrayOfRowIdsForAlignment->[$loop]->[$ColNumber-1] =~tr/1/1/;
$zeroes = $TransposedArrayOfRowIdsForAlignment->[$loop]->[$ColNumber-1] =~tr/0/0/;


if($zeroes < $ones) { $IndelEventsLength->[$eventscounter] = "D".$IndelEventsLength->[$eventscounter];}elsif
($zeroes > $ones) { $IndelEventsLength->[$eventscounter] = "I".$IndelEventsLength->[$eventscounter];}else {
$IndelEventsLength->[$eventscounter] = "U".$IndelEventsLength->[$eventscounter]; }
# if classified as an insertion or as a deletion 
# find the genomes invovled in that line of the alignment-transpose i.e. which are involved in insertion
# which in deletion
#
if($IndelEventsLength->[$eventscounter] !~/^U/) {
my @genomeidsarray = &getGenomeIdsFromEventsList($IndelEventsLength->[$eventscounter], $UniqueIndelEvents[$UniqueIndelEventsCounter]);
$UniqueIndelEvents[$UniqueIndelEventsCounter][$ColNumber] = scalar(@genomeidsarray);
if(scalar(@genomeidsarray) > 0) {
  $UniqueIndelEvents[$UniqueIndelEventsCounter][$ColNumber+1] = &findages(\@genomeidsarray);
 }
}
$eventscounter++;
print "No. of cols in newarray now is ",scalar(@{$UniqueIndelEvents[$UniqueIndelEventsCounter]}) ,"\n";
$UniqueIndelEventsCounter++;
}
$loop++;
}
print "No. of cols in newarray now is ",scalar(@{$UniqueIndelEvents[$UniqueIndelEventsCounter-1]}) ,"\n";
print "After merging duplicates rows :  number of deleted rows $numberofdeletedrows  size of matrix now ", scalar(@UniqueIndelEvents), " no. of events: ", $eventscounter, " no of cols: ", scalar(@{$UniqueIndelEvents[10]})," \n";
return \@UniqueIndelEvents;
}

############################################################################
sub findInStructureHit() {
my $id = shift ||  die " Can't search in structure hits with a blank Id $! \n";
my $pdbIds = shift || die "Need an array of PdbIds to search In: $! \n"; 

foreach my $pdbId (@$pdbIds) {
 if($pdbId eq $id) {
  return 1;
 } 
}
return 0;

}
############################################################################
sub getGenomeIdsFromEventsList() {
        my $id = shift || die "finding genome ids requires id \n";
        my $arr = shift || die "finding genome ids requires array in which to search for ids \n";
        my @genomearr = ();
        my $lookfor;
            if($id =~/^D/) {
                      $lookfor=0;
                } elsif ($id =~/^I/) {
                          $lookfor=1;
                    }
                     print "size of arr is ", scalar(@$arr), "\n";
                      for(my $counter=0;$counter < @$arr; $counter++) {
                                if(($arr->[$counter] ne "-")&& ($arr->[$counter] == $lookfor)) {
                                   print $genomeIds[$counter]," \n";
                                     push(@genomearr, $genomeIds[$counter]);
                                       }
                                 }
                                  return (@genomearr);
 }

