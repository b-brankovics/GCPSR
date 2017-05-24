#!/usr/bin/perl -w
use strict;
use Bio::Tree::Tree;
use Bio::TreeIO;

# Verbosity
my ($short) = grep{/^-?-s(hort)?$/} @ARGV;
my ($verbose) = grep{/^-?-v(erbose)?$/} @ARGV;
$verbose = undef if $short;

# Number the clades
my $clade_count = 0;
# Store the strain ids
my $allids = {};

# Deafult minimum support for clades to be kept
my $min_sup = 95;

# Get parameters for min count and support
my ($sup_par) = grep{/^-min=(\d+(:?\.\d+)?)$/} @ARGV;
if ($sup_par) {
    $sup_par =~ /^-min=(\d+(:?\.\d+)?)$/;
    $min_sup = $1;
}

# Get reference tree
my ($ref) = grep{/^-?-ref=(\S+)$/} @ARGV;
$ref =~ s/-?-ref=// if $ref;
my $ref_input;
if (!$ref) {
    # Reference tree is not specified
    unshift(@ARGV, "-h"); # Prints help
} elsif ($ref =~ /\.nwk$/) {
    $ref_input = new Bio::TreeIO(-file   => $ref,
			     -format => "newick");
} elsif ($ref =~ /\.nex\.\S+\.t(re)$/) {
    $ref_input = new Bio::TreeIO(-file   => $ref,
			     -format => "nexus");
} elsif ("-" eq $ref || "STDIN" eq $ref) {
    $ref_input = new Bio::TreeIO(-fh   => \*STDIN,
			     -format => "newick");
}
# Read the tree in the file and save reference clades
my $ref_clades = {};
if ($ref_input) {
    my $ref_tree = $ref_input->next_tree;

    my $rootn = $ref_tree->get_root_node;
    # Collect all the strain ids
    for ($rootn->get_all_Descendents()) {
	$allids->{&get_id($_)}++ if $_->is_Leaf;
    }
    # Collect Clade data
    &collect_clades($ref_tree, $ref_clades, \$clade_count);
} else {
    # Reference tree was not read properly
    unshift(@ARGV, "-h");
    if ($ref) {
	print {*STDERR} "Error: Something went wrong when processing '$ref' file as reference tree file\n";
    } else {
	print {*STDERR} "Error: No reference files is defined!\n" unless grep{/^(-)?-h(elp)?$/} @ARGV;
    }
}


# Test for arguments that won't be processed
my @not_used = grep{$_ !~ /^-min=(\d+(:?\.\d+)?)$/} grep{$_ !~ /^-?-ref=(\S+)$/} grep{$_ !~ /\.nwk$/} grep{$_ !~ /\.nex\.\S+\.t(re)$/} grep{$_ !~ /^-?-s(hort)?$/} grep {$_ !~ /^-?-v(erbose)?$/} @ARGV;
if (@not_used) {
    for (@not_used) {
	if (/^(-)?-h(elp)?$/) {
	    last;
	} elsif ("-" eq $_) {
	    print{*STDERR} "$0 can only read reference tree from STDIN (standard input).\n";
	}
	print{*STDERR} "'$_' is an incorrect argument\n";
    }
    die "Usage:\n\t$0 [-h | --help] [-s | --short] [-v | --verbose] [-min=<num>] -ref=reference-tree.nwk tree1 tree2 ...\n" .
	"Description:\n\tA tool to screen a set of trees and find the ones that are in conflict (discordant) with the reference tree.\n" .
	"Input:\n" .
	"\treference-tree.nwk\n" .
	"\t\tA reference tree containing clades against which discordance will be tested.\n" .
	"\t\tIf you want to read the reference tree from STDIN, then use '-ref=-' or '-ref=STDIN'\n" .
	"\t(tree1, tree2, ...)\n" .
	"\t\tMajority-rule concensus tree files in either newick or nexus format with bootstrap or Bayesian posterior probabilty support.\n" .
	"\t\tEach locus should have a separate tree and each files should contain only one tree.\n" .
	"Options:\n" .
	"\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
	"\t-min=<num>\n\t\tAll clades that have support values (BS or BPP) lower then <num> will be skipped. (Default: 95)\n" .
	"\t-s | --short\n\t\tPrint only the names of the files that contain conflicting trees.\n\t\t(This setting will overwrite 'verbose' when both are set)\n" .
	"\t-v | --verbose\n\t\tPrint the conflicting branch information as well.\n" .
	"\n";
 
}


# Process the input files
my $conflict_count;
for my $file (grep{$_ !~ /^-?-ref=(\S+)$/} grep{/\.nwk$/ || /\.nex\.\S+\.t(re)$/} @ARGV) {
    #print "Here is file '$file'\n";
    my $input;
    # Open trees or store minsupport value
    if ($file =~ /\.nwk$/) {
	$input = new Bio::TreeIO(-file   => $file,
				 -format => "newick");
    } elsif ($file =~ /\.nex\.\S+\.t(re)$/) {
	$input = new Bio::TreeIO(-file   => $file,
				 -format => "nexus");
    }
    # Read the tree in the file
    next unless $input;
    my $tree = $input->next_tree;
    my $rootn = $tree->get_root_node;
    # Collect all the strain ids
    for ($rootn->get_all_Descendents()) {
	$allids->{&get_id($_)}++ if $_->is_Leaf;
    }
    # Collect Clade data
    my $clades_h = {};
    my $clade_id = 0;
    # Is true if there is a conflict
    my $con;
    # Stores the ids of the conflicting clades ref to current
    my @is; # all the "i"s
    my @js; # all the "j"s
    &find_high_support($tree, $clades_h, \$clade_id);
    # Look for discordance

    for my $i (keys %$ref_clades) {
	for my $j (keys %$clades_h) {
	    if ( &conflict($ref_clades->{$i}, $clades_h->{$j}) ) {
		$con++;
		push @is, $i;
		push @js, $j;
	    }
	}
    }
    if ($con) {
	$conflict_count++;
	print "$file";
	print ": is in conflict ($con conflicts) with the reference tree" unless $short;
	print "\n";
	if ($verbose) {
	    print "\tclade conflicts (ref vs current):\n";
	    my $counter = 0;
	    for my $i (@is) {
		my $j = $js[$counter];
		print "\t\t(" . join(",", @{$ref_clades->{$i}}) . ") vs (" . join(",", @{$clades_h->{$j}}) . ")\n";
		$counter++;
	    }
	}
    }
}

exit if $short;
if ($conflict_count) {
    print "The number of trees which have clades (with support >= $min_sup) conflicting with the reference tree's clades is $conflict_count\n";
} else {
    print "No conflicts were found.\n"
}

#================Subroutines==================================================
sub collect_clades{
    # Collect clade data from the tree
    # All the inputs are references
    # Inputs: Tree, hash to store the clades and integer for clade numbering
    my ($tree, $hash, $i) = @_;
    # Get root
    my $root = $tree->get_root_node;
    # Go through all nodes
    for my $n ($root->get_all_Descendents()) {
	# Skip leaves
	next if $n->is_Leaf;
	# Get clade members, add to hash
	my $new++;
	my @clade = get_children($n);
	for (keys %$hash) {
	    my @now = @{ $hash->{$_} };
	    # Search already stored clades, whether it is already present
	    if (scalar(@now) == scalar(@clade)) {
		my %test;
		for (@now, @clade) {
		    $test{$_}++;
		}
		if (scalar(@now) == scalar(keys %test)) {
		    $new = 0;
		}
	    }
	}
	if ($new && scalar @clade > 1) {
	    # Increment the clade id
	    $$i++;
	    # Add clade to hash
	    $hash->{$$i} = \@clade;
	}
    }
}


sub create_tree{
    # Return a newick string
    # all the inputs are references to hashes
    # strain ids, clade id numbers, clade hash (id->array of members), clade support hash
    my ($names, $ids, $clade, $count) = @_;

    # If no more clades in ids then return the remaining strains
    if (scalar(keys %$ids) == 0) {
	return join(",", sort keys %$names);
    } else {
	# Group clades as "super" sets and subsets
	# (store the ids of subsets for the id of the superset)
	my %super;
	# Repeat untill all of them are chategorized
	while(keys %$ids) {
	    # Get the largest clade (sort based on number of members)
	    my ($now) = sort{scalar(@{ $clade->{$b} })<=>scalar(@{ $clade->{$a} })} keys %$ids;
	    delete $ids->{$now}; # Remove superclade
	    # Find subclusters
	    $super{$now} = {};
	    for (keys %$ids) {
		if (&is_subset($clade->{$_}, $clade->{$now})){
		    $super{$now}->{$_}++;
		    delete $ids->{$_};
		}
	    }
	}
	# Delete the strain names that are contained in the superclades
	# (they contain every id that are inside all the clades)
	for (keys %super) {
	    for (@{$clade->{$_}}) {
		delete $names->{$_};
	    }
	}
	# Recursive call for each superclade
	my $tree = join(",", map{# Create an anonymus hash to store the ids inside the super clade
	                         my $ns;
				 for (@{ $clade->{$_} }){
				     $ns->{$_}++
				 };
				 "(" . create_tree($ns,        # members of the clade
						   $super{$_}, # List of subclades
						   $clade,     # the two clade hash references
						   $count) . ")$count->{$_}"
			                                                    } sort{scalar(@{ $clade->{$b} })<=>scalar(@{ $clade->{$a} })} keys %super);
	# Append the remaining ids if there are
	return join(",", $tree, sort keys %$names);
    }
}

sub get_id{
    # Help get clean data from nodes
    my ($n) = @_;
    my $res = $n->id;
    $res =~ s/\s+$// if $res;
    return $res;
}

sub get_children{
    # Returns an array of all the leafs of the given node
    my ($node) = @_;
    my @children;
    # Examine all its descendents
    for ($node->get_all_Descendents()) {
	# Get the name of the leaves
	next unless $_->is_Leaf;
	my $id = &get_id($_);
	push @children, $id;
    }
    return sort @children;
}

sub find_high_support{
    my ($tree, $hash, $i) = @_;
    # Get root
    my $root = $tree->get_root_node;
    # Go through all nodes
    for my $n ($root->get_all_Descendents()) {
	# Skip leaves
	next if $n->is_Leaf;
	# Get support value if it exists
	my $value = &get_id($n);
	# If there is no min_sup then keep all
	#   else keep only well supported clades
	if (!$min_sup || ($value && $value >= $min_sup)) {
	    # Get clade members, add to hash
	    my $new++;
	    my @clade = get_children($n);
	    for (keys %$hash) {
		my @now = @{ $hash->{$_} };
		# Search already stored clades, whether it is already present
		if (scalar(@now) == scalar(@clade)) {
		    my %test;
		    for (@now, @clade) {
			$test{$_}++;
		    }
		    if (scalar(@now) == scalar(keys %test)) {
			$new = 0;
		    }
		}
	    }
	    if ($new && scalar @clade > 1) {
		# Increment the clade id
		$$i++;
		$hash->{$$i} = \@clade;
	    }
	}
    }
}


sub conflict{
    # Are the two clades in conflict
    # Inputs are references to arrays
    my ($c1, $c2) = @_;
    my ($u1, $u2, $s); # count of unique and shared elements
    for my $i (@{$c1}) {
	my $match;
	for my $j (@{$c2}) {
	    $match++ if $i eq $j;
	}
	# increment shared if it was found or unique if not
	$match ? $s++ : $u1++;
    }
    # If only one is defined then there is no conflict
    if ($u1 && $s) {
	for my $i (@{$c2}) {
	    my $match;
	    for my $j (@{$c1}) {
		$match++ if $i eq $j;
	    }
	    # increment shared if it was found or unique if not
	    $match ? $s++ : $u2++;
	}
	# If all three are defined then there is conflict
	return 1 if $u2;
    } 
    # No conflict
    return 0;
}

sub is_subset{
    # return true if a is subset of b
    # inputs are references to arrays
    my ($a, $b) = @_;
    my $found = 0;
    for my $i (@$a) {
	for my $j (@$b) {
	    $found++ if $i eq $j;
	}
    }
    # A is a subset of B if all the elements of A are found in B
    return $found == scalar @$a;
}
