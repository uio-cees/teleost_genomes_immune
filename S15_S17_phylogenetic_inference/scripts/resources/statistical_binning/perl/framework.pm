#!/usr/bin/env perl

package framework;

use strict;
#use setenv;

use Cwd;

# for temp files
use File::Copy;
use File::Temp qw/ tempfile tempdir /;



# contains filenames for well-known simulation study framework files
use constant ROSE_OUTPUT_TXT => "rose.output.txt";
use constant ROSE_INTERNAL_OUTPUT_TXT => "rose.internal.output.txt";
use constant RANDOM_TREE_FILENAME => "random.tree";
use constant ROSE_MODEL_FILE => "rose.model";
use constant ROSE_INTERNAL_MODEL_FILE => "rose.internal.model";
use constant ROSE_RAW_SEQUENCES_FILE => "rose.fas";
use constant ROSE_RAW_SEQUENCES_INTERNAL_FILE => "rose.fas.internal";
use constant ROSE_TRUE_ALIGNMENT_FILE => "rose.aln.true";
use constant ROSE_TRUE_ALIGNMENT_INTERNAL_FILE => "rose.aln.true.internal";
use constant ROSE_MODEL_TREE_FILE => "rose.mt";
use constant ROSE_MODEL_TREE_INTERNAL_FILE => "rose.mt.internal";
use constant ROSE_TRUE_TREE_FILE => "rose.tt";




# contains fixed constants for simulation study framework
use constant STATS_FILENAME => "STATS";
use constant STATS_ETC_FILENAME => "STATS_ETC";
use constant RESULTS_FILENAME => "RESULTS";

# temp stats file - all one-time, throwaway stats
use constant STATS_TEMP_FILENAME => "STATS_TEMP";

# for pretty printing fasta style sequence data
use constant FASTA_COLUMN_WIDTH => 60;


# statistic marker names
use constant TREELENGTH => "TREELENGTH";
use constant LENGTH => "LENGTH";
use constant ALIGNMENT_ERROR => "ALIGNMENT_ERROR";
use constant TREE_FP_ERROR => "TREE_FP_ERROR";
use constant TREE_FN_ERROR => "TREE_FN_ERROR";
use constant RUNTIME => "RUNTIME";

# to support dynamic condor job script generation
use constant DAGMAN_JOB_ID_PREFIX => "JOB";

use constant NEWICK_OPEN_PAREN => "(";
use constant NEWICK_CLOSE_PAREN => ")";
use constant NEWICK_COMMA => ",";
use constant NEWICK_SEMICOLON => ";";
use constant NEWICK_COLON => ":";

use constant MAX_RANDOM_SEED => 100000;

# for G-C content tests
# fine since force uppercase on input alignment
use constant G_SYMBOL => "G";
use constant C_SYMBOL => "C";

# contains subroutines used for simulation study framework

# simstudy file format is:
# <model condition, with all params separated by _>/<R? run number>/<any method, alignest or treeest or align+treeest, with all params separated by _>
# model condition run dir has STATS file and so does method dir
# true alignment file stored at rundir

# collator will just have associative array for each stat - append
# to running string for each stat as process each run dir
# have one associative array for each modelcondition/method?

# trim starting and ending whitespace
sub trim {
    my $line = shift;
    chomp($line);
    $line =~ s/^\s+//; 
    $line =~ s/\s+$//;
    return ($line);
}

sub cleanExistingDir {
    my $finaldir = shift;

    if (-e $finaldir) {
	# print warning
	print "WARNING: removing existing dir $finaldir\n";
	framework::systemSafe("rm -rf $finaldir");
    }

    framework::systemSafe("mkdir $finaldir");
}

# output array entries to each line in a filexs
sub outputEntriesFile {
    my $file = shift;
    my $arrRef = shift;

    my $fh;
    open ($fh, ">$file");
    foreach my $entry (@$arrRef) {
	print $fh "$entry\n";
    }
    close ($fh);
}

# read lines in file as an array
sub readEntriesFile {
    my $file = shift;

    my @result = ();

    my $fh;
    open ($fh, $file);
    while (<$fh>) {
	my $line = $_;
	chomp($line);
	push(@result, $line);
    }
    close ($fh);

    return (\@result);
}

# extract file contents as a single string
# without newlines!
sub getStringFromFile {
    my $file = shift;

    if (-e "$file") {
	my $result = "";
	open (FILE, $file);
	while (<FILE>) {
	    my $line = $_;
	    chomp($line);
	    $line = framework::trim($line);
	    $result .= $line;
	}
	close (FILE);
	return ($result);
    }
    else {
	return ("");
    }
}

# strips internal node names from newick format treestring
sub stripInternalNodeNames {
    my $tree = shift;

    # strip off internal node names
    $tree =~ s/\)(\w+):/\):/g;
    # also strip off root name
    $tree =~ s/\)(\w+);/\);/g;

    return ($tree);
}

# if arrays A and B are treated as sets
# returns A - B
#
# requires array references to be passed in
#
# nondestructive - returns new array copy with difference
sub getSetArrayDifference {
    my $ARef = shift;
    my $BRef = shift;

    my @diff = ();

    my %BHash = ();
    my $b;
    foreach $b (@$BRef) {
	$BHash{"$b"} = 1;
    }

    my $a;
    foreach $a (@$ARef) {
	if (!(defined($BHash{"$a"}) && ($BHash{"$a"} == 1))) {
	    push (@diff, $a);
	}
    }

    return (\@diff);
}

# prints taxa list to a file
sub printTaxaList {
    my $listRef = shift;
    my $file = shift;

    open (FILE, ">$file");
    my $a;
    foreach $a (@$listRef) {
	print FILE "$a\n";
    }
    close (FILE);
}

sub getNumTaxaFromAlignmentFile {
    my $alignmentFile = shift;

    my $taxaRef = getTaxonListFromAlignmentFile($alignmentFile);

    return (scalar(@$taxaRef));
}

# get taxon list from FASTA alignment file
# returns as a reference to an array
sub getTaxonListFromAlignmentFile {
    my $alignmentFile = shift;

    my @taxa = ();

    open (FILE, $alignmentFile);
    while (<FILE>) {
	my $line = $_;
	chomp($line);
	if ($line =~ m/^>/) {
	    my $taxon = $line;
	    $taxon =~ s/>//g;
	    $taxon = trim($taxon);
	    push (@taxa, $taxon);
	}
    }
    close (FILE);

    return (\@taxa);
}

# creates space delimited list of taxon names from newick format treestring
# now just extract taxon names
sub getTaxonListFromTreestring {
    my $tree = shift;

    $tree =~ s/\(/ /g;
    $tree =~ s/\)/ /g;
    $tree =~ s/,/ /g;
    $tree =~ s/:\d+.?\d+/ /g;
    $tree =~ s/;/ /g;
    $tree =~ s/^\s+//;
    $tree =~ s/\s+$//g;
    $tree =~ s/\s+/ /g;
    
    return ($tree);
}

# get number of taxa from space delimited list of taxon names
sub getNumTaxaFromTaxonList {
    my $taxonList = shift;

    my @taxonArray = split(/\s/, $taxonList);
    my $numtaxa = $#taxonArray + 1;

    return ($numtaxa);
}

# helper function
sub extractEdgesCountFromDistanceFile {
    my $distanceFile = shift;
    my $line;
    my $count;
    
    open (DISTANCEFILE, $distanceFile);
    while (<DISTANCEFILE>) {
	$line = $_;
	chomp($line);
	if ($line =~ m/Tree pair 1:/) {
	    $line =~ s/Tree pair 1://g;
	    $line =~ s/\s+//g;
	    $count = $line;
	}
    }
    close (DISTANCEFILE);

    return ($count);
}

# count number off of root
# do this by counting number of children at each level in the tree
# faster than parsing tree topology
sub getNumChildrenOffRoot {
    my $treestring = shift;

    my @counts = ();

    my $depth = 0;
    my $letter;

    $counts[0] = 0;

    for (my $i = 0; $i < length($treestring); $i++) {
	$letter = substr($treestring, $i, 1);
	if ($letter eq NEWICK_OPEN_PAREN) {
	    $depth++;

	    # open account at this depth
	    if (!defined($counts[$depth])) {
		$counts[$depth] = 0;
	    }
	}
	elsif ($letter eq NEWICK_CLOSE_PAREN) {
	    $counts[$depth]++;

	    $depth--;
	}
	elsif ($letter eq NEWICK_COMMA) {
	    $counts[$depth]++;
	}
    }

    # testing
#    print "@counts\n";

    if (!defined($counts[1])) {
	return (0);
    }
    else {
	return ($counts[1]);
    }
}

# warning - pass in pernicious cases like (a,b,c,d), and you'll get 0 back!
# since no internal edges
# similarly for (a,b)
# edge off the root are technically leaf edges, not internal edges
#
# also doesn't check for correct newick format - ill-formed strings get undefined results!
sub getNumInternalEdgesFromTreeString {
    my $truetree = shift;
    # count number of open parens
    # == number of internal nodes
    # then -2
    # since don't count root's parent edge
    # and root's 2 child edges are really 1 edge since unrooted topologies
    # -> number of internal edges

    # good catch by sindhu - last applies ONLY if the root has 3+ children, i.e. root
    # actually has degree > 2
    # since last -1 is only to collapse all degree-2 nodes, which should only be root
    # root is special

    my $count = 0;

    for (my $i = 0; $i < length($truetree); $i++) {
	if (substr($truetree, $i, 1) eq NEWICK_OPEN_PAREN) {
	    $count++;
	}
    }

    # kliu testing
#    print "baz $count\n";

    # subtract 1 once since don't count parent edge on root
    $count -= 1;

    # subtract 1 again if root is a degree-2 node, i.e. has 2 children
    # otherwise DON'T! root is special
    if (getNumChildrenOffRoot($truetree) == 2) {
	$count -= 1;
    }

    # don't worry about whether nodes other than root might be
    # degree-2
    # crazy stuff like ...((A), (B))
    # won't encounter this with our current suite of methods
    # and this is arguably a malformed newick tree string

    # for pernicious cases above
    if ($count < 0) {
	$count = 0;
    }

    return ($count);
}


# kliu - weird - sometimes the phylip call fails???

# returns "" empty string signalling error

# warning - don't pass in trees with 0 internal edges
# otherwise you get a warning and an empty string result signalling error

# warning - we always scale by the true tree numtaxa - 3
# returns result as "<fp> <fn>"
# arguments are:
# 1. working directory pathname
# 2. true tree filename with path
# 3. estimated tree filename with path
sub getTreeError {
    # kliu - per serita, ordering is wrong here
    # fndist actually computes fn from first tree to second tree
    # should be estimated tree first, then true tree second
    # switch the order of arguments
    my $workdir = shift;
    my $tree2File = shift;
    my $tree1File = shift;

    if ((!(defined($tree1File))) || (!(defined($tree2File)))) {
	return "";
    }
	

    my $fnerror = 1;
    my $fperror = 1;

    # returns result as "<fp> <fn>"
    return ("$fperror" . " " . "$fnerror");
}

# calculate alignment sp error
# arguments are:
# 1. working directory pathname
# 2. true alignment filename with path
# 3. estimated alignment filename with path
sub getAlignmentError {
    my $workdir = shift;
    my $trueAlignmentFile = shift;
    my $estimatedAlignmentFile = shift;
    my $scriptdir = getcwd;

    if (!defined($trueAlignmentFile) || !defined($estimatedAlignmentFile)) {
	return;
    }

    # estimated alignment should come in FASTA format
    
    # true alignment file needs to be converted from PHYLIP to FASTA and then cleaned
    chdir setenv::BIOPERL_DIR;
    my $phylipToFastaAlignCmd = setenv::PERL_COMMAND . " phylip_to_fasta_align_converter.pl $trueAlignmentFile | sed -e 's/\\/.*//g' > $workdir/rose.aln.true.fasta";
    framework::systemSafe($phylipToFastaAlignCmd);
    chdir $scriptdir;

    # now use DataMatrix to compare
    
    chdir setenv::GSP_DIR;
    my $distancedircmd = " gsp.score.DataMatrix -v $workdir/rose.aln.true.fasta -f $estimatedAlignmentFile -sp > $workdir/sp.error";
    framework::systemSafe($distancedircmd);
    chdir $scriptdir;

    open (ALIGNMENTERRORFILE, "$workdir/sp.error");
    my $alignmentError = <ALIGNMENTERRORFILE>; chomp ($alignmentError);
    close (ALIGNMENTERRORFILE);

    # clean intermediate files
    # keep sp.error file
    `rm $workdir/rose.aln.true.fasta`;

    return ($alignmentError);
}




# appends specified statistic name value pair to STATS file
# in specified directory
sub outputSTATSFile {
    my $finaldir = shift;
    my $statisticName = shift;
    my $statisticValue = shift;

    my $statsFileName = $finaldir . "/" . STATS_FILENAME;

    open (STATSFILE, ">>$statsFileName");
    print STATSFILE "$statisticName $statisticValue\n";
    close (STATSFILE);
}

# appends specified statistic name value pair to STATS file
# in specified directory
sub outputSTATS_ETCFile {
    my $finaldir = shift;
    my $statisticName = shift;
    my $statisticValue = shift;

    my $statsetcFileName = $finaldir . "/" . STATS_ETC_FILENAME;

    open (STATS_ETCFILE, ">>$statsetcFileName");
    print STATS_ETCFILE "$statisticName $statisticValue\n";
    close (STATS_ETCFILE);
}

# don't do timestamping
# testing - output in some sort of standard way
# also does timestamping
# use typeglobbing to pass in filehandle
sub testingOutput {
    my $fh = shift;
    my $statisticName = shift;
    my $statisticValue = shift;
    
    print $fh "$statisticName $statisticValue \n";
}

# get directory prefix from filename with directory prefix
sub getDirectoryPrefix {
    my $str = shift;
    $str =~ s/\/[^\/]*$/\//;
    return ($str);
}

# get directory suffix - all text after last slash
sub getDirectoryLastDirectory {
    my $str = shift;
    my @toks = split(/\//, $str);
    return ($toks[$#toks]);
}

sub prettyPrint {
    my $sequence = shift;

    prettyFilePrint ($sequence, *STDOUT);
}

# pretty print a sequence fasta-style
sub prettyFilePrint {
    my $sequence = shift;
    my $fh = shift;

    for (my $i = 0; $i < length($sequence); $i += FASTA_COLUMN_WIDTH) {
        if ($i + FASTA_COLUMN_WIDTH >= length($sequence)) {
            # at last row
            print $fh substr($sequence, $i, length($sequence) - $i) . "\n";
        }
        else {
            # not at last row
            print $fh substr($sequence, $i, FASTA_COLUMN_WIDTH) . "\n";

        }
    }
}

sub filePrintPair {
    my $taxon = shift;
    my $sequence = shift;
    my $fh = shift;

    print $fh ">$taxon\n";
    prettyFilePrint($sequence, $fh);
    print $fh "\n";

}

# print a taxonname sequence pair fasta-style
sub printPair {
    my $taxon = shift;
    my $sequence = shift;

    filePrintPair ($taxon, $sequence, *STDOUT);
}


# extract likelihood score from RAXML log output from run
sub extractRAXMLLikelihoodScore {
    my $raxmlout = shift;
    my $line = "";
    my $llh;

    open (RAXMLOUT, $raxmlout);
    while (<RAXMLOUT>) {
	$line = $_;
	chomp($line);
    }
    close (RAXMLOUT);

    # last line should have likelihood score in second spot
    my @toks = split(/\s+/, $line);
    $llh = $toks[1];

    return ($llh);
}

# extract likelihood score from RAXML info file
sub extractRAXMLLikelihoodScoreFromInfoFile {
    my $raxmlout = shift;
    my $line = "";
    my $llh;

    open (RAXMLOUT, $raxmlout);
    while (<RAXMLOUT>) {
	$line = $_;
	chomp($line);
	# only way to match for both GTRGAMMA and GTRCAT runs
	# argh
	if ($line =~ m/Final.*GAMMA.*-/) {
	    $llh = $line;
	    $llh =~ s/^.* -/-/;
	}
    }
    close (RAXMLOUT);

    return ($llh);
}


# from lisan's script - used for random seed to r8s
# hmm... well, why not reuse it for rose?
# doesn't seem to be a limit anyways for rose?? wraparound ok? does
# so deterministically?
sub getRandomSeed {
    # kliu - what is rseed here?
    #rseed();
    # kliu - this is bad practice - called automatically by rand()
    # anyways
    # according to perldocs, just use rand()'s implicit srand() call
    #srand();
    my $seed = int(rand() * MAX_RANDOM_SEED);
    return ($seed);
}

# kliu - with recent ubuntu hardy heron upgrades,
# /bin/sh is no longer a symlink to bash, it instead is a symlink to dash
# and dash doesn't understand &>, it parses that into & >
# -> breaks &> I/O redirection in all system(3) calls 
# such as that done by Perl's system function
# since system(3) invokes /bin/sh -c
#
# workaround is to autoenforce bash shell usage around the command
#
#
#
# taken from perl alarm man pages
# run a system command for a fixed amount of time only
#
# $maxTime is in seconds
#
# kills systemCommand if not finished by maxTime seconds
# caller needs to do cleanup/recovery actions separately
#
# returns (<timedOut flag = 1 if timedOut, 0 else>, <result only if not timedOut, empty string otherwise>)
sub executeWithTimeout {
    my $systemCommand = shift;
    my $maxTime = shift;

    my $result;
    my $timedOutFlag = 0;

    # hm - eval statement needs semicolon at the end
    eval {
	local $SIG{ALRM} = sub { die "alarm\n"; };
	alarm $maxTime;
	my $wrappedCommand = "/bin/bash -c \"$systemCommand\"";
	$result = `$wrappedCommand`;
	alarm 0;
    };

    if ($@) {
	die unless $@ eq "alarm\n";
	$timedOutFlag = 1;
	return ($timedOutFlag, "");
    } 
    else {
	$timedOutFlag = 0;
	return ($timedOutFlag, $result);
    }
}

# extract from raxml run dir the last timestamp/tree and
# copy it to result tree
sub extractLastTimestampedTreeToResultTreeRAxML {
    my $workdir = shift;
    my $suffix = shift;
    
    my $lastCheckpointNumber;

    open (RAXMLLOG, "RAxML_log.$suffix");
    while (<RAXMLLOG>) {
	my $line = $_;
	chomp ($line);
	my @fields = split ($line);
	# third field should be checkpoint number
	$lastCheckpointNumber = $fields[2];
    }
    close (RAXMLLOG);

    `cp $workdir/RAxML_checkpoint.$suffix.$lastCheckpointNumber RAxML_result.$suffix`;
}

# dumb load state
# doesn't make any assumption about order, type, etc.
sub loadState {
    my $stateFile = shift;
    my $stateReferencesArrayReference = shift;
    my @stateReferencesArray = @$stateReferencesArrayReference;
    my $i = 0;

    if (!(-e $stateFile)) {
	# NOOP
	print STDERR "ERROR: state file $stateFile doesn't exist. loadState() failed.\n";
	return;
    }

    open (STATEFILE, $stateFile);
    while (<STATEFILE>) {
	my $line = $_;
	chomp($line);
	my $stateReference = $stateReferencesArray[$i];
	$$stateReference = $line;
	$i++;
    }
    close (STATEFILE);
}

# dumb save state
# doesn't make any assumption about order, type, etc.
sub saveState {
    my $stateFile = shift;
    my $stateReferencesArrayReference = shift;
    my $i = 0;
    my $stateReference;

    # caller should have cleaned up $stateFile
    # paranoid
    if (-e $stateFile) {
	print "WARNING: state file $stateFile wasn't cleaned up by caller prior to saveState() call. Cleaning up.\n";
	`rm $stateFile`;
    }

    open (STATEFILE, ">$stateFile");
    foreach $stateReference (@$stateReferencesArrayReference) {
	print STATEFILE "$$stateReference\n";
    }
    close (STATEFILE);
    
    
}

# extract newick tree string from nexus tree
sub extractNewickTreeStringFromNexusTreeFile {
    my $treefile = shift;
    my $line;
    my $reading = 0;
    my $treestring = "";

    open (TREEFILE, "$treefile");
    while (<TREEFILE>) {
	$line = $_;
	$line =~ tr/a-z/A-Z/;
	chomp($line);
	if ($line =~ m/TREE.*\=\s+\[.*\]/) {
	    # start reading until end;
	    $treestring .= $line;
	    $reading = 1;
	}
	elsif ($reading) {
	    $treestring .= $line;
	}
	
	# end if end;
	if ($line =~ m/end\;/) {
	    $reading = 0;
	}
    }
    close (TREEFILE);

    $treestring =~ s/TREE.*\=\s+\[.*\]//;
    $treestring =~ s/END\;//g;
    $treestring =~ s/\s+//g;

    return ($treestring);
}

# sometimes has :0?
# kliu - also strip negative distances
# also strip e's and E's and negatives - any combo ok
sub stripDistancesFromNewickTreeString {
    my $treestring = shift;
    
    # key is semicolon delimiter
    $treestring =~ s/\:(0|1|2|3|4|5|6|7|8|9|e|E|\-|\.)+//g;
    # also last distance as necessary
    $treestring =~ s/\:0\;/\;/g;

    return $treestring;
}

# custom midpoint rooting
# branchlens in, branchlens out
# doesn't use scientific notation e, but java does have finite precision
# tends to shave off last few decimal points (18-20th) - it's generally ok
sub midpointRootCustom {
    my $inputTree = shift;
    my $outputTree = shift;
    my $stripInternalNodeNamesFlag = shift;
    my $workdir = shift;

    midpointRoot ($inputTree, $outputTree, $stripInternalNodeNamesFlag, $workdir, 0);
}

sub midpointRootWithPAUP {
    my $inputTree = shift;
    my $outputTree = shift;
    my $stripInternalNodeNamesFlag = shift;
    my $workdir = shift;

    midpointRoot ($inputTree, $outputTree, $stripInternalNodeNamesFlag, $workdir, 1);
}


# kliu - with recent ubuntu hardy heron upgrades,
# /bin/sh is no longer a symlink to bash, it instead is a symlink to dash
# and dash doesn't understand &>, it parses that into & >
# -> breaks &> I/O redirection in all system(3) calls 
# such as that done by Perl's system function
# since system(3) invokes /bin/sh -c
#
# workaround is to autoenforce bash shell usage around the command
#
# check error codes of system calls
#
# returns exit value of call
sub systemSafe {
    my $args = shift;
    
    system("/bin/bash -c \"$args\"");

    if ($? != 0) {
	print STDERR "Warning: system invocation $args has nonzero exit code!\n";
    }

    return ($?);
}

sub createScript {
    my $methoddir = shift;
    my $script = shift;
    my $commandsArrayRef = shift;

    my $command;

    open (SCRIPT, ">$script");
    foreach $command (@$commandsArrayRef) {
	print SCRIPT "$command\n";
    }
    close (SCRIPT);
}

# it appears that omitting condor output/error files (equivalent
# to piping them to /dev/null) reduces filesystem impact
# may help to keep large condor_dagman jobs from overwhelming
# the filesystem
#
# to support dynamically created condor scripts
# requires nfs usage
# avoidWorkstationsFlag indicates if we need to avoid interactive workstations - InMastodon || InScout
# universeFlag == 0 for standard, else vanilla
# omitCondorOutputFileFlag == 1 (optional) to omit condor output file
# omitCondorErrorFileFlag == 1 (optional) to omit condor error file
sub createCondorFile {
    my $condorFilename = shift;
    my $executable = shift;
    my $arguments = shift;
    my $universeFlag = shift;
    my $avoidWorkstationsFlag = shift;
    my $additionalRequirementsString = shift;
    my $omitCondorOutputFileFlag = shift;
    my $omitCondorErrorFileFlag = shift;
    
    open (CONDORFILE, ">$condorFilename");
    
    print CONDORFILE "+Group = \"GRAD\"\n";
    print CONDORFILE "+Project = \"COMPUTATIONAL_BIOLOGY\"\n";
    print CONDORFILE "+ProjectDescription = \"creating initial vanilla\"\n";
    print CONDORFILE "+WantPreempt = False\n";
    
    if ($universeFlag == 0) {
	print CONDORFILE "Universe = standard\n";
    }
    else {
	print CONDORFILE "Universe = vanilla\n";
    }

    print CONDORFILE "Executable = $executable\n";
    print CONDORFILE "Arguments = $arguments\n";
    if (($universeFlag != 0) && $avoidWorkstationsFlag) {
	print CONDORFILE "Requirements = InMastodon || InScout\n";
    }
    if (defined($additionalRequirementsString) && ($additionalRequirementsString ne "")) {
	print CONDORFILE "Requirements = $additionalRequirementsString\n";
    }

    # for safety - force initialdir 
    # to be prefix of condorFilename
    # if it exists
    my $initialdir = getDirectoryPrefix($condorFilename);
    $initialdir = trim($initialdir);
    if ($initialdir ne "") {
	print CONDORFILE "InitialDir = $initialdir\n";
    }

    print CONDORFILE "getenv = True\n";
    
    # weird - InitialDir above seems to have problem if absolute paths provided
    # here
    my $condorFilenameNoDir = getDirectoryLastDirectory($condorFilename);
    if (!defined($omitCondorErrorFileFlag) || ($omitCondorErrorFileFlag != 1)) {
	print CONDORFILE "Error = $condorFilenameNoDir.err\n";
    }
    if (!defined($omitCondorOutputFileFlag) || ($omitCondorOutputFileFlag != 1)) {
	print CONDORFILE "Output = $condorFilenameNoDir.out\n";
    }
    print CONDORFILE "Log = $condorFilenameNoDir.log\n";
    
    print CONDORFILE "Notification = Error\n";
    
    print CONDORFILE "Queue\n";
    
    close (CONDORFILE);
}

# create a simple linear condor dagman script
sub createSequentialCondorDagmanFile {
    my $condorDagmanFilename = shift;
    my $condorScriptsArrayRef = shift;
    
    open (DAGMANFILE, ">$condorDagmanFilename");
    for (my $i = 0; $i < (scalar @$condorScriptsArrayRef); $i++) {
	print DAGMANFILE "JOB " . DAGMAN_JOB_ID_PREFIX . "$i $$condorScriptsArrayRef[$i]\n";
    }
    print DAGMANFILE "\n";
    for (my $i = 0; $i < (scalar @$condorScriptsArrayRef) - 1; $i++) {
	print DAGMANFILE "PARENT " . DAGMAN_JOB_ID_PREFIX . $i . " CHILD " . DAGMAN_JOB_ID_PREFIX . ($i+1) . "\n";
    }
    print DAGMANFILE "\n";
    close (DAGMANFILE);
}

sub sum {
    my $arrayRef = shift;

    my $sum = 0.0;
    my $entry;
    foreach $entry (@$arrayRef) {
	$sum += $entry;
    }
    
    return ($sum);
}

# find average of an array
# E[X]
sub average {
    my $arrayRef = shift;

    my $sum = 0.0;
    my $count = 0;
    my $entry;
    foreach $entry (@$arrayRef) {
	$sum += $entry;
	$count++;
    }
    
    return ($sum / $count);
}

# find average of an array's squares
# E[X^2]
sub squaresAverage {
    my $arrayRef = shift;

    my $sum = 0.0;
    my $count = 0;
    my $entry;
    foreach $entry (@$arrayRef) {
	$sum += $entry * $entry;
	$count++;
    }
    
    return ($sum / $count);
}

# find standard deviation of an array's entries
sub stddev {
    my $arrayRef = shift;

    my $avg = average($arrayRef);
    my $sqAvg = squaresAverage($arrayRef);
    
    return (sqrt($sqAvg - $avg * $avg));
}

# find extrema 
# set maxFlag to 1 to get max, else gets min
sub findExtrema {
    my $arrayRef = shift;
    my $maxFlag = shift;

    if (scalar (@$arrayRef) == 1) {
	return ($$arrayRef[0]);
    }
    elsif (scalar (@$arrayRef) < 1) {
	return;
    }

    my $extrema = $$arrayRef[0];

    my $entry;
    foreach $entry (@$arrayRef) {
	if (defined($maxFlag) && ($maxFlag == 1)) {
	    if ($entry > $extrema) {
		$extrema = $entry;
	    }
	}
	else {
	    if ($entry < $extrema) {
		$extrema = $entry;
	    }
	}
    }
    
    return ($extrema);
}

sub min {
    my $arrayRef = shift;
    return (findExtrema($arrayRef, 0));
}

sub max {
    my $arrayRef = shift;
    return (findExtrema($arrayRef, 1));
}

# return number of lines in a file
sub fileLineCount {
    my $file = shift;

    if (!(-e $file)) {
	return (0);
    }
    else {
	my $count = 0;
	open (FILE, $file);
	while (<FILE>) {
	    $count++;
	}
	close (FILE);

	return ($count);
    }
}

sub autodetectModelConditionStrings {
    my $infile = shift;
    my $numModelConditionFields = shift;
    my $sortFunctionalRef = shift;

    my $mcFieldsWidth = $numModelConditionFields - 1;

    # record with hash
    my %mcs = ();

    # use fixed comma-delimited fields width
    open (INFILE, "$infile");
    while (<INFILE>) {
	my $line = $_;
	chomp ($line);
	my @fields = split(/,/, $line);
	# cool - slices
	my @mcarray = @fields[0 .. $mcFieldsWidth];
	my $mcstring = join(',', @mcarray);
	$mcs{"$mcstring"} = 1;
    }
    close (INFILE);

    my @mcskeys = keys(%mcs);

    # sort em
    if (defined($sortFunctionalRef)) {
	@mcskeys = sort $sortFunctionalRef @mcskeys;
    }
    else {
	@mcskeys = sort @mcskeys;
    }

    # shouldn't ever need to worry about memory in perl
    return (\@mcskeys);
}

# checks if an alignment has sequences all of the same length
# if true returns alignment length,
# otherwise returns undef
sub checkAlignmentLength {
    my $alignmentRef = shift;

    my @taxa = keys(%$alignmentRef);

    if (scalar(@taxa) <= 0) {
	print STDERR "ERROR: alignment in checkAlignmentLength is empty. Returning undef.\n";
	return (undef);
    }

    my $firstTaxon = $taxa[0];
    my $alignmentLength = length($$alignmentRef{$firstTaxon});

    foreach my $v (values(%$alignmentRef)) {
	if (length($v) != $alignmentLength) {
	    return (undef);
	}
    }

    return ($alignmentLength);
}

# reads FASTA alignment file and returns
# as an associative array keyed by taxon name
# WARNING: potentially expensive op for
# really large alignments!
#
# added an additional optional flag to force uppercase or not
#
# returns reference to associative arrayx
sub readAlignment {
    my $alignmentFile = shift;
    my $forceUppercaseFlag = shift;

    my %alignment = ();
    my $line;
    my $name;
    my $sequence = "";

    open (IN, $alignmentFile);
    while (<IN>) {
	$line = $_;
	$line = framework::trim($line);
	if (defined($forceUppercaseFlag) && ($forceUppercaseFlag eq "1")) {
	    $line = uc($line);
	}
	chomp($line);
	if ($line =~ m/>/) {
	    if (defined($name)) {
		if (defined($alignment{"$name"})) {
		    print STDERR "WARNING: taxon $name is repeated in $alignmentFile!\n";
		}
		$alignment{"$name"} = $sequence;
		$sequence = "";
	    }
	    
	    $line =~ s/>//g;
	    $line = framework::trim($line);
	    $name = $line;
	}
	else {
	    $sequence .= $line;
	}
    }

    if (defined($name)) {
	if (defined($alignment{"$name"})) {
	    print STDERR "WARNING: taxon $name is repeated in $alignmentFile!\n";
	}
	$alignment{"$name"} = $sequence;
	$sequence = "";
    }

    close (IN);

    return (\%alignment);
}

sub checkRawNucleotideSequencesStrict {
    my $infile = shift;

    my $previousSeqFlag = 1;

    my $line;
    open (IN, $infile);
    while (<IN>) {
	$line = $_;
	chomp($line);
	if ($line =~ m/^>/) {
	    $line =~ s/^>//;
	    if (!$previousSeqFlag) {
		print STDERR "ERROR: empty sequences not allowed.\n";
		return (0);
	    }

	    if ($line =~ m/^[^A-Z]/) {
		print STDERR "ERROR: taxon names must start with uppercase A-Z letter. Make sure your taxon name lines have no whitespace.\n";
		return (0);
	    }

	    if ($line =~ m/[^A-Z0-9]/) {
		print STDERR "ERROR: taxon names can contain only uppercase alphanumeric symbols. Make sure your taxon name lines have no whitespace.\n";
		return (0);
	    }
	    
	    $previousSeqFlag = 0;
	}
	else {
	    if ($line =~ m/[^ACGTN]/) {
		print STDERR "ERROR: only ACGTN letters allowed. Make sure to delete any and all whitespace in your file, especially at the ends of lines.\n";
		return (0);
	    }

	    if ($line =~ m/^$/) {
		print STDERR "ERROR: empty lines not allowed in file.\n";
		return (0);
	    }

	    $previousSeqFlag = 1;
	}
    }
    close (IN);

    return (1);
}

# warning - expensive - linear in number of branches with branch lengths
# ensure that all branch lengths are decimal notation
# not exponential notation
# since MAFFT newick2mafft.rb script and Opal can't handle
# exponential notation in guide trees
sub decimalFormatBranchLengths {
    my $treefile = shift;
    my $outfile = shift;

    my $treestring = "";

    open (IN, $treefile);
    while (<IN>) {
	my $line = $_;
	chomp ($line);
	$treestring .= $line;
    }
    close (IN);


    my @splits = split(/[\(\)\,\;]+/, $treestring);

    my $tok;
    my $origtok;
    foreach $tok (@splits) {
	# match at start -> token
	if (($tok ne "") && ($tok =~ m/:/)) {
	    $tok =~ s/^.*://;
	    $origtok = $tok;
	    $tok = sprintf("%.20f", $tok);
	    $treestring =~ s/:\Q$origtok\E/:$tok/g;
	}
    }
    
    open (OUT, ">$outfile");
    print OUT "$treestring\n";
    close (OUT);

}

# returns an array of pairwise p-distances between the sequences in an alignment
sub computePdistances {
    my $alignmentRef = shift;

    # paranoid
    if (scalar(keys(%$alignmentRef)) <= 0) {
	print STDERR "ERROR: no sequences to computePdistances() on.\n";
	return ();
    }
    
    my @pdistanceArr = ();

    my @taxa = keys(%$alignmentRef);

    my ($ki, $kj);
    for (my $i = 0; $i < scalar(@taxa); $i++) {
	for (my $j = $i + 1; $j < scalar(@taxa); $j++) {
	    my $ki = $taxa[$i];
	    my $kj = $taxa[$j];
	    push(@pdistanceArr, computePdistance($$alignmentRef{$ki}, $$alignmentRef{$kj}));
	}
    }

    return (\@pdistanceArr);
}

sub computePdistance {
    my $x = shift;
    my $y = shift;

    if (length($x) != length($y)) {
	die "ERROR: two sequences $x and $y have differing length in computePdistance. Aborting.\n";
    }

    my $count = 0;
    for (my $i = 0; $i < length($x); $i++) {
	if (substr($x, $i, 1) ne substr($y, $i, 1)) {
	    $count++;
	}
    }
    
    my $result = $count / length($x);

    # testing
    # print "pdistance $result $x $y\n";

    return ($result);
}

sub computeGCFraction {
    my $alignmentRef = shift;

    my @sequences = values(%$alignmentRef);

    # paranoid
    if (scalar(@sequences <= 0)) {
	print STDERR "ERROR: no sequences to compute G-C fraction on.\n";
	return (0);
    }

    my $alignmentLength = length($sequences[0]);

    my $count = 0;
    my $sequence;
    for $sequence (@sequences) {
	for (my $i = 0; $i < length($sequence); $i++) {
	    if ((substr($sequence, $i, 1) eq uc(G_SYMBOL)) ||
		(substr($sequence, $i, 1) eq lc(G_SYMBOL)) ||
		(substr($sequence, $i, 1) eq uc(C_SYMBOL)) ||
		(substr($sequence, $i, 1) eq lc(C_SYMBOL))) {
		$count++;
	    }
	}
    }

    return ($count / ($alignmentLength * scalar(@sequences)));
}


1;
