package Astral;
use strict;
use PifCosm_support_subs;
use File::Copy;

sub astral_tree {
    my $database = shift;
    my $path = shift;
    my $boot = shift; # to do bootstrapping, not included in code yet 
    my $anchor_gene = shift; # Taxa not in these gene trees will be excluded
    my @genes = @_; # Which gene trees to use
    my $dbh = PifCosm_support_subs::connect_to_database ( $database ); # connect to database
    if ($genes[0] eq 'all') { @genes = PifCosm_support_subs::get_gene_tables($dbh); } # if all genes should be used get them
    my $sth; #search statement handler
#    open INCLUDETAXA, '>', "XXXincludeTaxa.txt" or die "Could not open XXXincludeTaxa.txt: $!.\n";
    $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE gene LIKE ?"); # statment to search out genes
    my $TREEBENDER; # Variable for file or command line statement handler
    ### Should be put into function in PifCosm_support functions
    if ($anchor_gene ne 'n') { # If anchor genes given get the taxa that are in all their gene trees
	#$sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE gene LIKE ?");
	my @a_genes = split /,/, $anchor_gene; # separate individual genes
	for (my $i=0; $i < scalar @a_genes; ++$i) {
	    $sth->execute($a_genes[$i]);
	    my $tree = $sth->fetchrow_array(); # get tree, if several get first
	    if ($tree) { # if tree found
		if ($i==0) { # if first tree get tip labels
		    open $TREEBENDER, "| $External_program::treebender -t \",\n\" > XXXkeep_tips.txt";
		    print $TREEBENDER $tree;
		    close $TREEBENDER;
		}
		else { # if previous trees parsed, get tip labels that are also in those trees
		    open $TREEBENDER, "| $External_program::treebender -d file:XXXkeep_tips.txt -i | $External_program::treebender -t \",\n\" > XXXtemp.txt";
		    print $TREEBENDER $tree;
		    close $TREEBENDER;
		    move "XXXtemp.txt","XXXkeep_tips.txt";
		}
	    }
	}
    }
    ###########
    unlink "XXXtreesForAstral.trees"; # Clear tree file
    if ($anchor_gene ne 'n') { # If anchor genes given prepare to remove tips not in anchor gene trees
	open $TREEBENDER, "| $External_program::treebender -d file:XXXkeep_tips.txt -i >> XXXtreesForAstral.trees";
    }
    else { # Otherwise prepare to save trees as is
	open $TREEBENDER, '>>', "XXXtreesForAstral.trees" or die "Could not open XXXtreesForAstral.trees: $!.\n";
    }
    foreach my $gene (@genes) { # Action! Get the trees
	$sth->execute($gene);
	my $tree = $sth->fetchrow_array();
	print $TREEBENDER $tree;
    }
    close $TREEBENDER;
    $sth->finish();
    $astral_tree = `java -jar astral.4.10.12.jar -i XXXtreesForAstral.trees`; # Run astral and get tree
    unlink "XXXtreesForAstral.trees","XXXkeep_tips.txt";
    # Save Astral tree
    my $n_inserted_trees = $dbh->do("INSERT INTO alignment_groups (taxon,gene,tree,alignable) VALUES ('combined_astral','combined_astral',$astral_tree,2)") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
}
