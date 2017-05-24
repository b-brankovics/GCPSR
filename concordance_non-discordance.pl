#!/usr/bin/perl -w
use strict;
use Bio::Tree::Tree;
use Bio::TreeIO;

# Number the clades
my $clade_count = 0;
# Store the members of the clades and the support
my $clades_h = {}; # clades hash: id -> array (strain ids) ref
my $clades_c = {}; # clades count: id -> count/frequency
# Store the strain ids
my $allids = {};

# Deafult minimum support for clades to be kept
my $min_sup = 95;
# Deafult minimum count/frequency for clades to be kept
my $min_count = 2;

# Get parameters for min count and support
my ($sup_par) = grep{/^-min=(\d+(:?\.\d+)?)$/} @ARGV;
if ($sup_par) {
    $sup_par =~ /^-min=(\d+(:?\.\d+)?)$/;
    $min_sup = $1;
}
my ($count_par) = grep{/^-count=(\d+)$/} @ARGV;
if ($count_par) {
    $count_par =~ /^-count=(\d+)$/;
    $min_count = $1;
}

# Test for arguments that won't be processed
my @not_used = grep{$_ !~ /^-min=(\d+(:?\.\d+)?)$/} grep{$_ !~ /^-count=(\d+)$/} grep{$_ !~ /\.nwk$/} grep{$_ !~ /\.nex\.\S+\.t(re)$/} @ARGV;
if (@not_used) {
    for (@not_used) {
	if (/^(-)?-h(elp)?$/) {
	    last;
	} elsif ("-" eq $_) {
	    print{*STDERR} "$0 cannot read trees from STDIN (standard input).\n";
	}
	print{*STDERR} "'$_' is an incorrect argument\n";
    }
    die "Usage:\n\t$0 [-h | --help] [-min=<num>] [-count=<int>] tree1 tree2 ...\n" .
	"Description:\n\tA tool to implement the concordance and non-discordance analysis of GCPSR sensu Brankovics et al. 2017\n" .
	"Input:\t(tree1, tree2, ...)\n" .
	"\tMajority-rule concensus tree files in either newick or nexus format with bootstrap or Bayesian posterior probabilty support.\n" .
	"\tEach locus should have a separate tree and each files should contain only one tree.\n" .
	"Options:\n" .
	"\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
	"\t-min=<num>\n\t\tAll clades that have support values (BS or BPP) lower then <num> will be skipped. (Default: 95)\n" .
	"\t-count=<int>\n\t\tAll clades that are present in less then <int> majority-rule concensus trees will not be considered concordant. (Default: 2)" .
	"\n";
}
# Process the input files
for (grep{/\.nwk$/ || /\.nex\.\S+\.t(re)$/} @ARGV) {
    my $input;
    # Open trees or store minsupport value
    if (/\.nwk$/) {
	$input = new Bio::TreeIO(-file   => $_,
				 -format => "newick");
    } elsif (/\.nex\.\S+\.t(re)$/) {
	$input = new Bio::TreeIO(-file   => $_,
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
    &find_high_support($tree, $clades_h, $clades_c, \$clade_count);
}


# Remove Discordant clades
# Store final clades in a hash
my $nondis = {};
for my $i (keys %$clades_h) {
    # Is true if there is a conflict
    my $con;
    # Get support of the clade
    my $sup1 = $clades_c->{$i};
    next unless $sup1 >= $min_count; # Clade will not be kept
    for my $j (keys %$clades_h) {
	# Get support of the other clade
	my $sup2 = $clades_c->{$j};
	next unless $sup2 >= $min_count; # $sup1;
	# Test if there is a conflict
	if ( &conflict($clades_h->{$i}, $clades_h->{$j}) ) {
	    $con++;
	    last;
	}
    }
    # Add to the final hash
    $nondis->{$i} = $clades_h->{$i} unless $con;
}

$clades_h = $nondis;

my %clade_ids;
for (keys %$clades_h) {
    $clade_ids{$_}++;
}

# Print the final newick tree
print "(" . &create_tree($allids, \%clade_ids, $clades_h, $clades_c) . ");\n";

#================Subroutines==================================================
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
    my ($tree, $hash, $count, $i) = @_;
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
			$count->{$_}++;
		    }
		}
	    }
	    if ($new && scalar @clade > 1) {
		# Increment the clade id
		$$i++;
		$count->{$$i} = 1;
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
