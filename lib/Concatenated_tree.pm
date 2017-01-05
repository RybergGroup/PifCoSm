package Concatenated_tree;
use strict;
use PifCosm_support_subs;

sub concatenated_tree { # create a tree for the concatenated genes
    my $database = shift @_; # databse name
    my $path = shift @_; # path to tree building and tree manupilating programms
    my $final_tree_method = shift @_; # the method to use to tree
    my $n_threads = shift @_;
    my $store_boot_trees = shift @_;
    my $anchor_gene = shift @_; # gene that all taxa must have
    my @genes = @_; # genes to include in concatenated tree
    my $dbh = PifCosm_support_subs::connect_to_database ( $database ); # database handler
    if ($genes[0] eq 'all') { @genes = PifCosm_support_subs::get_gene_tables($dbh); } # if all genes should be used get them
    my $sth; #search statement handler
    { # make sure there are sequences for all genes
        my @temp;
        foreach (@genes) {
            $sth = $dbh->prepare("SELECT COUNT(taxon_name) FROM alignments WHERE $_\_sequence!='empty'") or die;
            $sth->execute();
            my $count= $sth->fetchrow_array();
            if ($count && $count>3) { push (@temp,$_); }
            else { print "$_ has fewer than four sequences and will be ignored.\n"; }
            $sth->finish();
        }
        @genes = @temp;
    }
    if (!@genes) { die "No genes to make concatenated tree from.\n"; }
    print "****************************************\n";
    print " Create tree for concatenated sequences\n";
    print "****************************************\n";
    print "Using genes:";
    foreach (@genes) { print " $_"; }
    print "\n";
    print "Getting sequences.\n";
    my $sequence_columns=''; # string for the gene columns
    my @sequence_lengths; # the length of each gene alignment
    my $length=0; # total length
    my $get_partitions = 'n';
    my $interleaved = 'n';
    my $bases_per_row = 5000;
    if ($final_tree_method =~ /_partition/) { $get_partitions = 'y'; }
    print "Creating ML tree.\n";
    my $method;
    my $method_append = '';
    if ($final_tree_method =~ /raxml/i) { $method = 'RAxML'; }
    elsif ($final_tree_method =~ /ExaML/i) {
        $method = 'ExaML';
        if ($final_tree_method =~ /RAx-br/i) { $method_append .= ',RAx-br'; }
    }
    elsif ($final_tree_method =~ /fasttree/i) { $method = 'fasttree'; $interleaved = 'y'; }
    else { die "Did not recognize any method for getting ML tree.\n"; }
    my $n_boot;
    if ($final_tree_method =~ /_bootstrap_([0-9]+)/) {
 	$n_boot = $1;
    }
    PifCosm_support_subs::print_alignment($database,'phylip',$get_partitions,'XXXtemp','taxon_name',$interleaved,$bases_per_row,$anchor_gene,@genes);
    if ($final_tree_method =~ /_partition/) { $get_partitions = 'XXXtemp.partitions.txt'; }
    my ($ml_tree,$ml_score,$boottrees_file) = PifCosm_support_subs::run_raxml ( 'XXXtemp.phy', $path, 'XXXtemp', $n_threads, $method . $method_append, $get_partitions,$n_boot);
    if ($ml_tree ne 'empty' && -e $ml_tree) {
	my $method = 'ML';
	if ($ml_tree =~ /parsimony/) { $method = 'MP'; }
	if (PifCosm_support_subs::read_tree_to_alignment_groups($ml_tree,$dbh,"combined","combined",$method) < 1) { print STDERR "WARNING!!! Was not able to store ML trees.\n"; }
	unlink $ml_tree;
    }
    else {
	print STDERR "WARNING!!! No ML tree found, so it cannot be saved.\n";
    }
    #unlink $get_partitions;
    #if ($ml_tree ne 'empty' and -e $ml_tree) {
    #    rename $ml_tree, 'XXX_ML.tree';
    #    $ml_tree = 'XXX_ML.tree';
        #if ($final_tree_method =~ /_bootstrap_([0-9]+)/) {
        #    my $n_boot = $1;
        #    print "Performing bootstrap analysis.\n";
        #    if ($method eq 'RAxML' or $method eq 'ExaML') {
        #        unlink glob "RAxML*.XXXtemp";
        #        if ($get_partitions eq 'n') {
        #            if ($n_threads > 1) {
        #                my @raxml = `${path}raxmlHPC-PTHREADS -x 12345 -p 54321 -# $n_boot -m GTRGAMMA -s XXXtemp.phy -n XXXtemp -T $n_threads`;
        #                @raxml = `${path}raxmlHPC-PTHREADS -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart`;
        #            }
        #            else {
        #                my @raxml = `${path}raxmlHPC -x 12345 -p 54321 -# $n_boot -m GTRGAMMA -s XXXtemp.phy -n XXXtemp`;
        #                @raxml = `${path}raxmlHPC -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart`;
        #            }
        #        }
        #        else {
        #            if ($n_threads > 1) {
        #                my @raxml = `${path}raxmlHPC-PTHREADS -x 12345 -p 54321 -# $n_boot -m GTRGAMMA -s XXXtemp.phy -q $get_partitions -n XXXtemp -T $n_threads`;
        #                @raxml = `${path}raxmlHPC-PTHREADS -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart -T $n_threads`;
        #            }
        #            else {
        #                my @raxml = `${path}raxmlHPC -x 12345 -p 54321 -# $n_boot -m GTRGAMMA -s XXXtemp.phy -q $get_partitions -n XXXtemp`;
        #                @raxml = `${path}raxmlHPC -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart`;
        #            }
        #        }
        #        rename "RAxML_bipartitions.XXXbipart",'XXX_ML.tree';
        #        unlink "XXXtemp.phy.reduced";
        #        $ml_tree = 'XXX_ML.tree';
	if ($store_boot_trees eq 'y') {
	    if ($boottrees_file ne 'empty' && -e $boottrees_file) {
		if (PifCosm_support_subs::read_tree_to_alignment_groups($boottrees_file,$dbh,"combined_boot","combined_boot","ML_boot") < 1) { print STDERR "WARNING!!! Was not able to store bootstraped trees.\n"; }
		#open BOOTTREES, $boottrees_file or die "Could not open $boottrees_file, that should contain bootstrap trees: $!.\n";
		#my $boot_trees = '';
		#while (my $temp = <BOOTTREES>) { $boot_trees .= $temp; }
		#close BOOTTREES or die;
		#my $sth = $dbh->prepare ("SELECT COUNT(taxon) FROM alignment_groups WHERE taxon='combined_boot' AND gene='combined_boot'") or die;
		#$sth->execute();
		#my $count = $sth->fetchrow_array();
		#$sth->finish();
		#my $string = $dbh->quote($boot_trees);
		#if ($count and $count >= 1) {
		#    my $updated = $dbh->do("UPDATE alignment_groups SET tree=$string WHERE taxon='combined_boot' and gene='combined_boot'") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
		#    if ($updated < 1) { print STDERR "WARNING!!! Was not able to store bootstrap trees.\n"; }
		#}
		#else {
		#    my $updated = $dbh->do("INSERT INTO alignment_groups (taxon,gene,tree,alignable) VALUES ('combined_boot','combined_boot',$string,2)") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
		#    if ($updated < 1) { print STDERR "WARNING!!! Was not able to store bootstraped trees.\n"; }
		#}
	    }
	    else { print STDERR "WARNING!!! No bootstrap trees found, so they cannot be saved.\n"; }
	}
	unlink $boottrees_file;
	#unlink glob "RAxML*.XXXbipart"
	#}
        #    elsif ($method eq 'fasttree') {
        #        system "${path}FastTree -quiet -nt -gtr -gamma -noml -boot $n_boot -intree $ml_tree > XXXtemp.tree";
        #        rename "XXXtemp.tree",'XXX_ML.tree';
        #        $ml_tree = 'XXX_ML.tree';
        #    }
        #    else { unlink "XXXtemp.phy"; print STDERR "Did not recognize method for bootstrap analysis. No bootstrap support generated.\n"; }
        #}
	#if ($ml_tree ne 'empty' && -e $ml_tree) {
	#    my $method = 'ML';
	#    if ($ml_tree =~ /parsimony/) { $method = 'MP'; }
	#    if (PifCosm_support_subs::read_tree_to_alignment_groups($ml_tree,$dbh,"combined","combined",$method) < 1) { print STDERR "WARNING!!! Was not able to store ML trees.\n"; }
	    #open MLTREE, "<$ml_tree" or die "Could not read ML tree ($ml_tree): $!.\n";
	    #$ml_tree = <MLTREE>;
	    #close MLTREE or die;
	    #unlink 'XXX_ML.tree';
	    #my $sth = $dbh->prepare ("SELECT COUNT(taxon) FROM alignment_groups WHERE taxon='combined' AND gene='combined'") or die;
	    #$sth->execute();
	    #my $count = $sth->fetchrow_array();
	    #$sth->finish();
	    #my $string = $dbh->quote($ml_tree);
	    #if ($count and $count >= 1) {
		#my $updated = $dbh->do("UPDATE alignment_groups SET tree=$string WHERE taxon='combined' and gene='combined'") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
		#if ($updated < 1) { print STDERR "WARNING!!! Was not able to store ML tree.\n"; }
	    #}
	    #else {
		#my $updated = $dbh->do("INSERT INTO alignment_groups (taxon,gene,tree,alignable) VALUES ('combined','combined',$string,2)") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
		#if ($updated < 1) { print STDERR "WARNING!!! Was not able to store ML tree.\n"; }
	    #}
	#    unlink $ml_tree;
	#}
	#else {
	#    print STDERR "WARNING!!! No ML tree found, so it cannot be saved.\n";
	#}
	unlink glob "$get_partitions*";
	unlink ("XXXtemp.phy","XXXtemp.phy.reduced");
	$dbh->disconnect();
}
1;
