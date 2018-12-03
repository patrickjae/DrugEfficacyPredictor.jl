################################################################
#
# Usage: 
# > perl weighted_average_concordance_index.pl filename
#
# filename is the NCI-DREAM subhcallenge 1 submission file which
# is a comma separated file (csv) with the drugs listed as the 
# columns and the cell lines listed as the rows.
#
# Output:
# the ouput of this script is a tab-delimited text file that 
# lists the individual drug and overall scores for the tested file 
#
# This script was created by Jim Costello to score the NCI-DREAM
# drug sensitivity challenge.  Please contact me at jccostel@bu.edu
# for further questions.
#
# Modifications by Patrick Jähnichen to accept crucial variables as 
# cmd line arguments. Updated usage:
# > perl weighted_average_concordance_index.pl /path/to/pooled/std /path/to/test/data /path/to/zscore file_to_score /path/to/outfile /path/to/summary
#################################################################
use Math::Erf::Approx qw(erf);


######
# define some global variables
my $sdfile = $ARGV[0];
my $testdatafile = $ARGV[1];
my $zscorefile = $ARGV[2];
my $infile = lc($ARGV[3]);
my $target_file = $ARGV[4];
my $summary_file = $ARGV[5];
print STDERR "file \"$infile\" will be scored\n";


######
# instead of testing all of the randomly generated predictions to
# calculate the mean, standard deviation, and z-scores of drugs
# in the test data, simply read the precomputed values 
my %zscores = ();
my $overall_mean = 0;
my $overall_var = 0;
open(FILE, "$zscorefile") || die;
while(<FILE>) {
    chomp;
    next if $_ =~ /^\#/;
    my ($drug, $z, $normz) = split(/\t/, $_);
    if ($drug eq 'Overall') {
	$overall_mean = $z;
	$overall_var = $normz;
    } else {
	$zscores{$drug} = $normz;
    }
}
close FILE;

######
# The standard deviation for the anonymized drugs are contained
# within this file.  These information will be used later on.
my %pooledSD = ();
open(FILE, "$sdfile") || die;
while(<FILE>) {
    chomp;
    next if $_ =~ /^\#/;
    my ($anonId, $sd) = split(/\t/, $_);
    next if $anonId eq 'AnonID';
    $pooledSD{$anonId} = $sd;
}
close FILE;
print STDERR "The standard deviations for " . keys(%pooledSD) . " drugs have been read\n";

######
# The test data is parsed and used for scoring
my %gi50 = ();
my @header = ();
open(FILE, "$testdatafile") || die;
while(<FILE>) {
    chomp;
    my ($cellLine, @vals) = split(/\t/, $_);
    if ($cellLine eq 'CellLine') {
	@header = @vals;
    } else {
# read in all the values and set the NA values to 0
	for(my $i=0; $i<=$#vals; $i++) {
	    die "$header[$i]" if !defined($header[$i]);
	    if ($vals[$i] eq 'NA') {
		$gi50{$cellLine}{$header[$i]} = 0;
	    } else {
		$gi50{$cellLine}{$header[$i]} = $vals[$i];
	    }
	}
    }
}
close FILE;
print STDERR keys(%gi50) . " cell lines and $#header drugs have been read in from the test data\n";

######
# Parse and hash the data for the submission file to test
my @drugs = ();
my %preds = ();
open(FILE, "$infile") || die "Please enter a file to score\n";
while(<FILE>) {
    chomp;
    my ($cellLine, @vals) = split(/\,/, $_);
    if ($cellLine eq 'DrugAnonID') {
	@drugs = @vals;
    } else {
	for(my $i=0; $i<=$#vals; $i++) {
	    $preds{$cellLine}{$drugs[$i]} = $vals[$i];
	}
    }    
}
close FILE;

print STDERR "There are " . keys(%preds) . " cell lines and $#drugs drugs in the predictions file entered\n";

# some sanity checks
die "There is an error in your submission file.  Please make sure the first id in the first row (header row) is \"DrugAnonID\"\n" if $#drugs < 30;
# die "Please make sure the submission file is titled \"dream7_drugsensitivity1_predictions_<team name>.csv\"\n" if $infile !~ /dream7\_drugsensitivity1\_predictions/;

open(my $fh, "$target_file");

print $fh "DrugID\tprobabalistic c-index\tweighted probabalistic c-index\n";


#######
# loop through all the drugs, score each indvidual drug, and keep a running sum of these scores
# that will be used to calculate the overall score
my $overallscore = 0; # running sum 
my $weightsum = 0; # running sum 
my $weightctr = 0; # running counter
foreach my $drug (@drugs) {    

    # select the cell lines to score, for a given drug, there will be a variable number of cell lines
    my @cell_lines = ();
    foreach my $line(keys %preds) {
        next if $gi50{$line}{$drug} ==0;
	push(@cell_lines, $line);
    }

    my $score = calculate_prob_cindex(\@cell_lines,$drug,\%gi50,\%preds,\%pooledSD);
    $overallscore += $score;
    my $weight = $score * $zscores{$drug};
    $weightsum += $weight;
    $weightctr += $zscores{$drug};
    print $fh "$drug\t$score\t$weight\n";

} 
my $oas = sprintf("%.5f", ($overallscore / 31));
my $ws = sprintf("%.5f", ($weightsum/$weightctr));
my $pv = 1 -  (.5 * (erf(($ws - $overall_mean)/(sqrt(2*$overall_var))) + 1));
# print "##############################################\n";
# print "Submission: $infile\n";
print $fh "average probablistic c-index: $oas\n";
print $fh "weighted average probablistic c-index: $ws\n";
print $fh "weighted average probablistic c-index p-value: $pv\n";
print $fh "The final team rankings are based on the weighted average probablistic c-index $ws\n";
close($fh);

open(my $fh, ">>$summary_file");
print $fh "$infile\t$ws\n";
close($fh);
exit;

#######
# subroutine that calculates the probabilistic c-index
sub calculate_prob_cindex {
    my $cell_line_ref = shift;
    my $drugname = shift;
    my $gi50_ref = shift;
    my $pred_ref = shift;
    my $pooledSD_ref = shift;


    my $sum = 0; # keeps the running sum 
    my %done = (); # keeps track of if we have seen the cell line pair.  This is needed to account for the cell lines
                   # that have the same values.  They would be counted twice if they are not accounted for


    for(my $i=0; $i<=$#{$cell_line_ref}; $i++) {
        my $test_val1 = $gi50_ref->{$cell_line_ref->[$i]}{$drugname}; # concentration of cell line i from the test data
        my $pred_val1 = $pred_ref->{$cell_line_ref->[$i]}{$drugname}; # rank of cell line i from the predictions
        for(my $j=0; $j<=$#{$cell_line_ref}; $j++) {
	    next if $j == $i;
	    my $test_val2 = $gi50_ref->{$cell_line_ref->[$j]}{$drugname}; # concentration of cell line j from the test data
            my $pred_val2 = $pred_ref->{$cell_line_ref->[$j]}{$drugname}; # rank of cell line j from the predictions
	    if ($test_val1 >= $test_val2 && $pred_val1 < $pred_val2) { # the concordance case
		next if defined($done{$j}{$i});
		my $cdf = .5 * (erf(($test_val1 - $test_val2)/(sqrt(2*($pooledSD_ref->{$drugname} * $pooledSD_ref->{$drugname})))) + 1);
		$sum += $cdf; # the running sum will add a value between .5 and 1
		$done{$i}{$j}++;
	    } elsif ($test_val1 >= $test_val2 && $pred_val1 > $pred_val2) { # the discordance case
		next if defined($done{$j}{$i});
		my $cdf = 1 - (.5 * (erf(($test_val1 - $test_val2)/(sqrt(2*($pooledSD_ref->{$drugname} * $pooledSD_ref->{$drugname})))) + 1));
		$sum += $cdf; # the running sum will add a value between 0 and 1
		$done{$i}{$j}++;
	    } 
        }
    }
    
    # in the end, all cell line pairs will be tested once or, (n * n-1)/2  
    my $n = (($#{$cell_line_ref} + 1) * $#{$cell_line_ref}) / 2;
    return ($sum/$n);
}


