#
# MOTree -- simple and fast data structures and parser for phylogenetic trees

package MOTree;
require Exporter;

use strict;

# A MOTree has the following entries:
# nNodes
# children: internal_id => list of children
# parents: internal_id => parent id or undef if root (usually 0)
# ids: internal_id => name of leaf or internal node (often a bootstrap value; optional)
# leaves: id => name for leaves
# branch_lengths: internal_id => branch length
# root: internal_id of the root
# Note -- MOTree::new does NOT check that branch lengths are valid
#
# In the below utilities, a node is just an internal integer id, not an actual object,
# so all the routines that operates on nodes are called like this:
# $tree->parent($node)
#
# Create a new tree from a string or a filehandle, e.g.:
# MOTree::new( newick => "((A:1,B:2)30,(C:3,D:4.1)40);" );
# MOTree::new( fh => \*STDIN );
# If reading from a file handle, returns undef on end-of-file, and
# each tree must end on its own line but a tree may be split across lines
sub new {
    my %params = @_;

    my $root = 0;
    my $nNodes = 1;
    my $children = [ [] ];
    my $parents = [ undef ];
    my $ids = [ undef ];
    my $branch_lengths = [ undef ];

    if (exists $params{fh}) {
	my $fh = $params{fh};
	my $string = "";
	while(my $line = <$fh>) {
	    $line =~ s/^\s+//;
	    $line =~ s/\s+$//;
	    $string .= $line;
	    last if $line =~ m/;$/;
	}
	return undef if ($string eq "");
	die "Not a newick tree: $string\n" unless $string =~ m/;$/;
	$params{newick} = $string;
    }

    if (exists $params{newick}) {
	my $string = $params{newick};
	$string =~ s/[\r\n]//g;
	return undef unless $string =~ m/;$/;
	$string =~ s/;$//;
	my @tokens = split /,/, $string;

	my @sofar = (); # the stack we are at inside the tree
	push @sofar, $root;
	my $firsttoken = 1;
	foreach my $token (@tokens) {
	    my $nDown = 0;
	
	    # go down
	    if ($token =~ m/^([\(]+)(.*)$/) {
		$nDown = length($1);
		$token = $2;
	    }
	    
	    # if not first, then had a comma, which implies up and then down
	    if ($firsttoken) {
		$firsttoken = 0;
	    } else {
		return undef if @sofar < 1;
		pop @sofar;
		$nDown++;
	    }
	    
	    while ($nDown-- > 0) {
		my $node = $nNodes++;
		$children->[$node] = [];
		my $parent = $sofar[-1];
		$parents->[$node] = $parent;
		push @{ $children->[$parent] }, $node;
		push @sofar, $node;
	    }
	    
	    my @pieces = split /[\)]/, $token, -1; # -1 means keep trailing empty fields
	    # We've removed leading left-parens (from first piece), commas (we're within a token),
	    # and these pieces are separated by right-parens
	    # Possibilities are:
	    # leading right-paren (implied up), which leads to recursion
	    # nodename:length
	    # bootstrap:length
	    # :length
	    # name       (no length)
	    # 
	    while (@pieces > 0) {
		my $piece = shift @pieces;
		if ($piece ne "") {
		    my @parts = split /:/, $piece;
		    my ($name, $length) = (undef,undef);
		    if (@parts == 1) {
			$name = $parts[0];
			$length = undef;
		    } elsif (@parts == 2) {
			($name,$length) = @parts;
			$length = undef if $length eq ""; # allow missing lengths
		    } else {
			return undef; # too many :
		    }
		    # note -- do not check if length is legit
		    my $node = $sofar[-1];
		    $ids->[$node] = $name if (defined $name);
		    $branch_lengths->[$node] = $length if defined $length;
		}
		if (@pieces > 0) {
		    return undef if @sofar==0;
		    pop @sofar;
		}
	    }
	}
	return undef unless @sofar == 1; # we should end at root!
    } elsif (scalar(keys %params) == 0) {
	# do nothing -- return an empty tree
    } else {
	die "Unknown paramters " . join(",", sort keys %params) . " to MOTree::new";
    }
    my $self =  { nNodes_ => $nNodes, children_ => $children, parents_ => $parents,
		  ids_ => $ids,
		  branch_lengths_ => $branch_lengths,
		  root_ => $root };
    bless $self;
    return $self;
}

sub get_root_node($) {
    my ($self) = @_;
    return $self->{root_};
}

sub set_root_node($$) {
    my ($self,$newroot) = @_;
    $self->{root_} = $newroot;
    $self->{parents_}[$newroot] = undef;
}

sub ancestor($$) {
    my ($self,$node) = @_;
    return $self->{parents_}[$node];
}

# get a node's id, usually a leaf name or a bootstrap value for an internal node
sub id($$) {
    my ($self,$node) = @_;
    return $self->{ids_}[$node];
}

sub set_id($$) {
    my ($self,$node,$newid) = @_;
    $self->{ids_}[$node] = $newid;
}

sub branch_length($$) {
    my ($self,$node) = @_;
    return $self->{branch_lengths_}[$node];
}

sub set_branch_length($$$) {
    my ($self,$node,$new_length) = @_;
    $self->{branch_lengths_}[$node] = $new_length;
}

# returns a list
sub children($$) {
    my ($self,$node) = @_;
    return @{ $self->{children_}[$node] };
}

sub nChildren($$) {
    my ($self,$node) = @_;
    return scalar(@{ $self->{children_}[$node] });
}

sub is_Leaf($$) {
    my ($self,$node) = @_;
    return scalar(@{ $self->{children_}[$node] }) == 0;
}

# Add a new child within specified node
sub newChild($$) {
    my ($self,$ancestor) = @_;
    my $newchild = $self->{nNodes_}++;
    push @{ $self->{children_}[$ancestor] }, $newchild;
    $self->{parents_}[$newchild] = $ancestor;
    $self->{children_}[$newchild] = [];
    return $newchild;
}

# delete node and have its parent link to the child(ren) instead,
# sum branch lengths of children, and leave internal ids (bootstraps) alone
# returns what used to be the node's ancestor
sub removeNode($$) {
    my ($self,$node) = @_;
    my $parents = $self->{parents_};
    my $children = $self->{children_};
    my $branch_lengths = $self->{branch_lengths_};

    my $ancestor = $parents->[$node];
    return if !defined $ancestor; # cannot remove root

    # Update children of ancestor and parents of children
    my @ancestorchildren = grep {$_ != $node} @{ $children->[$ancestor] };
    my @nodechildren = @{ $children->[$node] };
    foreach my $child (@nodechildren) {
	$parents->[$child] = $ancestor;
	push @ancestorchildren, $child;
    }
    $children->[$ancestor] = \@ancestorchildren;
    $children->[$node] = [];
    # so we can walk up, find the wrong root, and know this is a deleted node
    $parents->[$node] = undef;

    # adjust branch lengths
    my $oldlen = $branch_lengths->[$node];
    if (defined $oldlen) {
	foreach my $child (@nodechildren) {
	    $branch_lengths->[$child] += $oldlen if defined $branch_lengths->[$child];
	}
    }
    return $ancestor;
}

# id => leaf with that id
# note -- does a full scan not a hash lookup
# note -- checks that this is not a deleted node
sub findLeaf($$) {
    my ($self, $name) = @_;
    my $nNodes = $self->{nNodes_};
    my $ids = $self->{ids_};
    my $parents = $self->{parents_};
    my $children = $self->{children_};
    my $root = $self->{root_};
    for (my $node = 0; $node < $nNodes; $node++) {
	my $nodeid = $ids->[$node];
	if (defined $nodeid && $nodeid eq $name && $self->is_Leaf($node)) {
	    # check that we are connected to the root
	    my $ancestor = $node;
	    while(defined $parents->[$ancestor]) {
		$ancestor = $parents->[$ancestor];
	    }
	    return $node if $ancestor == $root;
	}
    }
    return undef;
}

# node -> reference to a list (includes the node as the first item and the root as the last item)
sub pathToRoot($$) {
    my ($self,$node) = @_;
    my $parents = $self->{parents_};
    my @path = ();
    do {
	push @path, $node;
	$node = $parents->[$node];
    } while (defined $node);
    return \@path;
}

# returns a reference to a list
# the optional second argument is a hash of internal_id => sortby value
sub depthfirst {
    my $tree = shift @_;
    die unless defined $tree;
    my $sorthash = undef;
    $sorthash = shift @_ if @_ > 0;
    my $list = $tree->all_descendents($tree->{root_}, $sorthash);
    unshift @$list, $tree->{root_};
    return $list;
}

# returns a reference to a list
# the optional second argument is a hash of internal_id => sortby value
sub all_descendents {
    my $self = shift @_;
    my $begin_node = shift @_;
    die unless defined $begin_node;
    my $sorthash = undef;
    $sorthash = shift @_ if @_ > 0;

    my @dfirst = ();
    my @work = ( $begin_node );
    my $children = $self->{children_};

    while (@work > 0) {
	my $node = shift @work;
	push @dfirst, $node unless $node == $begin_node;
	my @children = @{ $children->[$node] };
	if (defined $sorthash) {
	    @children = sort { $sorthash->{$a} <=> $sorthash->{$b} } @children;
	}
	unshift @work, @children;
    }
    return \@dfirst;
}

# the second argument is a hash of node => value for the nodes to stop at
sub all_descendents_limited_by($$$) {
    my ($tree,$begin_node,$limithash) = @_;
    my @dfirst = ();
    my @work = ( $begin_node );
    my $children = $tree->{children_};

    while (@work >0) {
	my $node = shift @work;
	push @dfirst, $node unless $node == $begin_node;
	if (!exists $limithash->{$node}) {
	    my @children = @{ $children->[$node]};
	    unshift @work, @children;
	}
    }
    return \@dfirst;
}

# the second argument is a hash of node => value for the nodes to stop at
sub all_above_limited_by($$$) {
    my ($tree,$begin_node,$limithash) = @_;
    my @out = ();
    my $ancestor = $tree->ancestor($begin_node);
    return [] if !defined $ancestor;

    my $children = $tree->{children_};
    my @work = ( $ancestor );
    my %visited = ( $begin_node => 1 );
    while (@work > 0) {
	my $node = shift @work;
	push @out, $node;
	$visited{$node} = 1;
	if (!exists $limithash->{$node}) {
	    my @children = @{ $children->[$node] };
	    unshift @work, grep { !exists $visited{$_} } @children;
	    push @work, $tree->ancestor($node) if defined $tree->ancestor($node);
	}
    }
    return \@out;
}

# format the tree in newick format, ending with a semicolon (but no newline)
# does not handle illegal characters in ids or non-numeric branch lengths
sub toNewick($) {
    my ($self) = @_;
    my $children = $self->{children_};
    my $ids = $self->{ids_};
    my $branch_lengths = $self->{branch_lengths_};

    my @pieces = ();
    my @work = ( [ $self->{root_}, 0] ); # a stack of [node, flag for entering (0) or exiting node (1)]
    my $depth = 0;
    while (@work > 0) {
	my ($node,$endflag) = @{ shift @work };
	my @children = @{ $children->[$node] };
	my $branch_length = $branch_lengths->[$node];
	my $id = $ids->[$node];
	my $parent = $self->ancestor($node);
	if (scalar(@children) == 0) {
	    push @pieces, "," if $children->[$parent][0] != $node;
	    push @pieces, $id if defined $id;
	    push @pieces, ":".$branch_length if defined $branch_length;
	} elsif ($endflag) {
	    push @pieces, ")";
	    push @pieces, $id if defined $id;
	    push @pieces, ":".$branch_length if defined $branch_length;
	} else {
	    push @pieces, "," if defined $parent && $children->[$parent][0] != $node;
	    push @pieces, "(";
	    unshift @work, [$node,1];
	    foreach my $child (reverse @children) {
		unshift @work, [$child,0];
	    }
	}
    }
    push @pieces, ";";
    return join("", @pieces);
}

# compute the maximum distance to root from every node; uses 1 if branch_length is not defined
# returns a reference to a hash of node => value
sub nodeDepth {
    my ($tree) = shift @_;
    my $useBranchLength = (@_ > 0) ? (shift @_) : 1;
    my $dfirst = $tree->depthfirst();
    my $nodeDepth = { $tree->get_root_node() => 0 };
    foreach my $node (@$dfirst) {
	my $selfd = $nodeDepth->{$node};
	die unless defined $selfd;
	foreach my $child ($tree->children($node)) {
	    my $down = $useBranchLength ? $tree->branch_length($child) : 1;
	    $down = 1 if !defined $down;
	    $nodeDepth->{$child} = $selfd + $down;
	}
    }
    return $nodeDepth;
}

# compute the maximum distance to child leaves from every node; uses 1 if branch_length is not defined
# returns a reference to a hash of node => value
sub nodeElevation {
    my ($self) = shift @_;
    my $useBranchLength = (@_>0) ? (shift @_) : 1;
    my $dfirst = $self->depthfirst();
    my %nodeElevation = ();
    my $children = $self->{children_};
    my $branch_lengths = $self->{branch_lengths_};

    foreach my $node (reverse @$dfirst) {
	my $childlist = $children->[$node];
	if(scalar(@$childlist) == 0) {
	    $nodeElevation{$node} = 0;
	} else {
	    my $max = -1e8;
	    foreach my $child (@$childlist) {
		die "No child elevation" unless defined $nodeElevation{$child};
		my $total = 0;
		if ($useBranchLength) {
		    my $childlen = $branch_lengths->[$child];
		    $total = $nodeElevation{$child} + (defined $childlen ? $childlen : 1);
		} else {
		    $total = $nodeElevation{$child} + 1;
		}
		$max = $total if $total > $max;
	    }
	    $nodeElevation{$node} = $max;
	}
    }
    return \%nodeElevation;
}


# Iterate from a node "up" the tree -- visit its parent, do a depthfirst traversal of its siblings,
# visit its uncle, do a depthfirst traversal of its cousins, visit uncle again with end=1, etc.
#
# Usage:
# my $iterator = $tree->up_iterator($start);
# while (my $where = $tree->next_up($iterator)) {
#     my ($node,$end,$elevation,$totbelow) = @$where;
# }
#
# Elevation/totbelow will be 0/1 for leaves (which are only entered, not exited)
# and will be set for internal nodes only when exiting them
#
# Totbelow counts leaves, not internal nodes
#
sub up_iterator($$) {
    my ($self,$node) = @_;
    return { stack => [ [$node,1] ], maxElevation => { $node => 0 }, totbelows => { $node => 0 } };
}

sub next_up($$) {
    my ($self,$iterator) = @_;
    my $stack = $iterator->{stack};
    my $maxElevation = $iterator->{maxElevation};
    my $totbelows = $iterator->{totbelows};
    return () if @$stack == 0;
    my ($node,$end) = @{ shift @$stack };

    my $childlist = $self->{children_}[$node];
    if (@$childlist == 0) {
	$maxElevation->{$node} = 0;
	$totbelows->{$node} = 1;
    } elsif ($end) {
	$maxElevation->{$node} = -1e8;
	$totbelows->{$node} = 0;
	foreach my $child (@$childlist) {
	    my $branchlen = $self->{branch_lengths_}[$child];
	    my $sum = $maxElevation->{$child} + (defined $branchlen ? $branchlen : 1);
	    $maxElevation->{$node} = $sum if $sum > $maxElevation->{$node};
	    $totbelows->{$node} += $totbelows->{$child};
	}
    } else {
	$maxElevation->{$node} = undef;
	$totbelows->{$node} = undef;
	unshift @$stack, [$node,1];
	foreach my $child (reverse @$childlist) {
	    unshift @$stack, [$child,0] unless exists $maxElevation->{$child};
	}
    }
    if (@$childlist == 0 || $end) {
	my $ancestor = $self->{parents_}[$node];
	if (defined $ancestor && !exists $maxElevation->{$ancestor}) {
	    push @$stack, [$ancestor,0];
	    $maxElevation->{$ancestor} = undef;
	}
    }
    return [ $node, $end, $maxElevation->{$node}, $totbelows->{$node} ];
}

sub distanceToRoot($$) {
    my ($tree,$node) = @_;
    die unless defined $node;
    my $root = $tree->get_root_node;
    my $branch_lengths = $tree->{branch_lengths_};
    my $parents = $tree->{parents_};
    my $depth;
    for ($depth = 0; defined($node) && $node != $root; $node = $parents->[$node]) {
	die "Undefined length for $node root $root name " . $tree->id($node) unless defined $branch_lengths->[$node];
	$depth += $branch_lengths->[$node];
    }
    return( $depth );
}

# Given a tree and a list of nodes (usually leaves), remove an singleton internal nodes
# that remain
# Only tested with nodes as leaves (but think it should work otherwise as well)
sub removeNodes($$) {
    my ($tree,$removeList) = @_;
    my $nodes = $tree->depthfirst();
    my %originalLeaf = map {$_ => 1} grep { $tree->is_Leaf($_) } @$nodes;
    my %deleted = ();
    foreach my $node (@$removeList) {
	$tree->removeNode($node);
	$deleted{$node} = 1;
    }
    my $changed = 0;
    my $root = $tree->get_root_node;
    do {
	$changed = 0;
	foreach my $node (reverse @{ $tree->depthfirst() }) {
	    if (!exists $deleted{$node} && !exists $originalLeaf{$node} && $tree->nChildren($node) < 2) {
		$tree->removeNode($node);
		$deleted{$node} = 1;
		$changed = 1;
	    }
	}
    } while($changed);
}

# returns a reference to a hash of node => number of leaves beneath it
sub countAllLeaves($) {
    my ($tree) = @_;
    my %nleaves = ();
    foreach my $node (reverse @{ $tree->depthfirst() }) {
	my @children = $tree->children($node);
	if (@children == 0) {
	    $nleaves{$node} = 1;
	} else {
	    my $total = 0;
	    foreach my $child (@children) {
		my $n = $nleaves{$child};
		die unless defined $n && $n > 0;
		$total += $n;
	    }
	    $nleaves{$node} = $total;
	}
    }
    return( \%nleaves );
}

# Compare two unrooted trees with the same list of unique leaf names
# Returns a hash of node1 => [node2,bDown2,nLeaves] for matching splits
#    bDown2 is 1 if node2's descendents match node1's,
#	 or 0 if the leaves NOT beneath node2 match node1's descendents
#    nLeaves is the total #leaves beneath node1
# Trivial (nLeaves==1 or nLeaves==all-1) splits are included, but not the root
# The Robinson-Foulds distance is then 2*(#nonmatching splits)
# and the topological %difference is nNonMatching/(nTotalLeaves-2)
sub Compare($$) {
    my ($t1,$t2) = @_;

    my $d1 = $t1->depthfirst();
    my @l1 = grep { $t1->is_Leaf($_) } @$d1;
    my %leaf1 = map { $t1->id($_) => $_ } @l1;
    die "Non-unique leaf names in 1st tree" unless scalar(keys %leaf1) == scalar(@l1);
    my $nDesc1 = $t1->countAllLeaves();

    my $d2 = $t2->depthfirst();
    my @l2 = grep { $t2->is_Leaf($_) } @$d2;
    my %leaf2 = map { $t2->id($_) => $_ } @l2;
    die "Non-unique leaf names in 2nd tree" unless scalar(keys %leaf2) == scalar(@l2);
    my $nDesc2 = $t2->countAllLeaves();

    # Note -- these includes leaves, but not root1 <=> root2, as this prevents traversals up
    # and would screw up all_above_limited_by
    my %match = (); # node1 => [node2,bDown2]
    my %matchReverse = (); # node2 => [node1,bDown2]
    my $root1 = $t1->get_root_node();
    my $root2 = $t2->get_root_node();

    # first, find all the greedy matches
    foreach my $node1 (reverse @$d1) {
	next if $node1 == $root1;
	my $n1 = $nDesc1->{$node1};
	if ($n1 == 1) {
	    my $id = $t1->id($node1);
	    my $node2 = $leaf2{$id};
	    die "id $id in 1st tree is not in second tree" unless defined $node2;
	    $match{$node1} = [$node2,1,$n1];
	    $matchReverse{$node2} = [$node1,1];
	    my $a1 = $t1->ancestor($node1);
	    my $a2 = $t2->ancestor($node2);
	} else {
	    # try to match via a matching child
	    my @children1 = $t1->children($node1);
	    die unless @children1 > 0;
	    my @match1 = grep {exists $match{$_}} @children1;
	    if (scalar(@children1) == scalar(@match1)) {
		# do they all match to corresponding nodes?
		my $child2 = $match{$children1[0]}[0];
		my $node2 = $t2->ancestor($child2);
		@match1 = grep { my $a2 = $t2->ancestor($match{$_}[0]);
				 defined($a2) && defined($node2) && $a2 == $node2; } @match1;
		if (scalar(@children1) == scalar(@match1)
		    && scalar(@match1) == scalar($t2->children($node2))) {
		    $match{$node1} = [$node2,1,$n1];
		    $matchReverse{$node2} = [$node1,1];
		    my $a1 = $t1->ancestor($node1);
		    my $a2 = $t2->ancestor($node2);
		    die "node1 $node1 node2 $node2 d1 $nDesc1->{$node1} d2 $nDesc2->{$node2}"
			unless $nDesc1->{$node1} == $nDesc2->{$node2};
		}
	    }
	}
    }

    # next, hash all the splits by the sum of a hashvalue for topmost matched descendents
    my %hashval1 = ();
    my %hashval2 = ();
    my $hashmax = 8 * scalar(@$d1) + 1;
    my $hashtot = 0;
    my %order = map { $_ => rand() } (0..($hashmax-1));
    my @values = sort {$order{$a} <=> $order{$b}} (0..($hashmax-1));
    use integer; # make sure no floating-point conversions happen

    while (my ($node1,$row) = each %match) {
	my $a1 = $t1->ancestor($node1);
	next if defined $a1 && exists $match{$a1}; # only consider maximal matching nodes
	my $node2 = $row->[0];
	$hashval1{$node1} = shift @values;
	$hashval2{$node2} = $hashval1{$node1};
	$hashtot += $hashval1{$node1};
	$hashtot -= $hashmax if $hashtot >= $hashmax;
    }
    foreach my $node1 (reverse @$d1) {
	next if $node1 == $root1 || exists $hashval1{$node1} || exists $match{$node1};
	my @children1 = $t1->children($node1);
	my $hashval = 0;
	foreach my $child1 (@children1) {
	    die unless exists $hashval1{$child1};
	    $hashval += $hashval1{$child1};
	}
	$hashval1{$node1} = $hashval % $hashmax;
    }
    foreach my $node2 (reverse @$d2) {
	next if $node2 == $root2 || exists $hashval2{$node2} || exists $matchReverse{$node2};
	my @children2 = $t2->children($node2);
	my $hashval = 0;
	foreach my $child2 (@children2) {
	    die unless exists $hashval2{$child2};
	    $hashval += $hashval2{$child2};
	}
	$hashval2{$node2} = $hashval % $hashmax;
    }
    my @hashed1 = map { [] } (0..($hashmax-1));
    while (my ($node,$val) = each %hashval1) {
	push @{ $hashed1[$val] }, $node;
	push @{ $hashed1[($hashmax+$hashtot-$val) % $hashmax] }, $node;
    }

    my @hashed2 = map { [] } (0..($hashmax-1));
    while (my ($node,$val) = each %hashval2) {
	push @{ $hashed2[$val] }, $node;
	push @{ $hashed2[($hashmax+$hashtot-$val) % $hashmax] }, $node;
    }

    my $nMatchThisRound;
    do {
	$nMatchThisRound = 0;
	foreach my $node1 (reverse @$d1) {
	    next if $node1 == $root1 || exists $match{$node1};
	    my $hashval = $hashval1{$node1};
	    my @desc1 = grep { exists $match{$_} }
	    @{ all_descendents_limited_by($t1,$node1,\%match) };
	    
	    my @matchOf1 = sort {$a <=> $b} map { $match{$_}[0] } @desc1;
	    my @candidate2 = @{ $hashed2[$hashval] };
	    foreach my $node2 (@candidate2) {
		next unless $nDesc1->{$node1} == $nDesc2->{$node2};
		my @desc2 = sort { $a <=> $b }
		grep { exists $matchReverse{$_} }
		@{ all_descendents_limited_by($t2,$node2,\%matchReverse) };
		
		if (scalar(@matchOf1) == scalar(@desc2)
		    && scalar( grep { $matchOf1[$_] != $desc2[$_] } (0..(scalar(@matchOf1)-1)) ) == 0) {
		    $match{$node1} = [$node2,1,$nDesc1->{$node1}];
		    $matchReverse{$node2} = [$node1,1];
		    $nMatchThisRound++;
		    last;
		}
	    }
	}

	# and now look for inverse matches
	foreach my $node1 (@$d1) {
	    next if $node1 == $root1 || exists $match{$node1};
	    my $hashval = $hashval1{$node1};
	    my @desc1 = grep { exists $match{$_} }
	    @{ all_descendents_limited_by($t1,$node1,\%match) };

	    my @matchOf1 = sort {$a <=> $b} map { $match{$_}[0] } @desc1;
	    my @candidate2 = @{ $hashed2[($hashmax+$hashtot-$hashval) % $hashmax] };
	    foreach my $node2 (@candidate2) {
		next unless $nDesc1->{$node1} == scalar(@l2) - $nDesc2->{$node2};
		my @up2 = sort {$a <=> $b} grep { exists $matchReverse{$_} }
		@{ all_above_limited_by($t2,$node2,\%matchReverse) };

		if (scalar(@matchOf1) == scalar(@up2)
		    && scalar( grep { $matchOf1[$_] != $up2[$_] } (0..(scalar(@matchOf1)-1)) ) == 0) {
		    $match{$node1} = [$node2,0,$nDesc1->{$node1}];
		    $matchReverse{$node2} = [$node1,0];
		    $nMatchThisRound++;
		    last;
		}
	    }
	}
    } while ($nMatchThisRound > 0);
    return \%match;
}

1;
