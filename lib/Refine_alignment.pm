package Refine_alignment;
use strict;
use PifCosm_support_subs;

sub refine_alignment { # This module aligns the sequences similarily to SATe, it require aligned sequences as input 
    my $database = shift @_;
    my $path = shift @_; # path to treebuilding and treemanipulating programms
    my $dbh = PifCosm_support_subs::connect_to_database( $database );
    my $tree_method = shift @_;
    my $save_boot_trees = shift @_;
    my $stop_criterion = shift @_;
    my $max_size = shift @_; # max size of alignment groups
    my $use_guide_tree = shift @_; # if ML tree should be used as guide tree when aligning
    my $n_threads = shift @_;
    #my $n_genes = shift @_;
    my @align_genes = split /,/, shift @_;
    #for (my $i=0;$i<$n_genes;++$i) {
    #    $align_genes[$i] = shift @_;
    #}
    my @linked_genes = split /,/, shift @_;
    my $n_linked_genes = scalar @linked_genes;
    if ($align_genes[0] eq 'all') {
        @align_genes = PifCosm_support_subs::get_gene_tables($dbh);
        #$n_genes = scalar @align_genes;
    }
    my $split_switches = "longest_branch:n --split_stop t:2";
    $stop_criterion =~ /^([A-Za-z]+)_([0-9]+)/;
    my $stop_crit = $1;
    my $n_iterations = $2;
    print "Have been asked to refine alignments for the following genes:\n";
    foreach (@align_genes) { print "$_ "; }
    print "\n";
    { # Clean away genes that do not have any sequences
        my @temp;
        my $i=0;
        foreach (@align_genes) {
            my $count;
            my $sth=$dbh->prepare("SELECT COUNT(taxon_name) FROM alignments WHERE $_\_sequence!='empty'") or die;
            $sth->execute();
            $count=$sth->fetchrow_array();
            $sth->finish();
            if ($count && $count>0) { push(@temp,$_); }
            else {
                print "$_ have no sequences associated with it and will be ignored.\n";
            }
        }
        @align_genes = @temp;
    }
    for (my $i=0; $i<scalar @align_genes;++$i) {
        my $add = 'y';
        for (my $j=0; $j< scalar @linked_genes; ++$j) {
            if ($linked_genes[$j] =~ /$align_genes[$i]/) { $add='n'; last; }
        }
        if ($add eq 'y') { push (@linked_genes,$align_genes[$i]); }
    }
    for (my $i=0; $i < scalar @linked_genes; ++$i) {
        my @regions = split /;/,$linked_genes[$i];
        my $temp;
        for (my $j=0; $j<scalar @regions; ++$j) {
            my $add='n';
            foreach(@align_genes) { if ($_ eq $regions[$j]) { $add='y'; last;} }
            if ($add eq 'y') {
                if ($temp) { $temp .= ";$regions[$j]"; }
                else { $temp = $regions[$j]; }
            }
        }
        if ($temp) { $linked_genes[$i] = $temp; }
        else { $linked_genes[$i] = 'empty'; }
    }
    { # clean out linked gene groups that have no sequences that should be aligned
        my @temp;
        my $i=0;
        foreach (@linked_genes) { if ($_ ne 'empty') { $temp[$i++] = $_; } }
        @linked_genes = @temp;
    }
    $n_linked_genes = scalar @linked_genes;
    print "Will align the following groups. The same guide tree will be used for all genes within a group, but different trees will be used between groups: ";
    foreach (@linked_genes) { print "$_ "; }
    print "\n";
    if ($stop_crit eq 'max') { print "Will stop $n_iterations iterations after last improvment in likelihood score for each alignment group.\n"; }
    elsif ($stop_crit eq 'iterations') { print "Will run $n_iterations iterations for each alignment group.\n"; }
    for (my $i=0; $i<$n_linked_genes; ++$i) {
        print "Refining alignments for $linked_genes[$i].\n";
        my $tree;
        my @genes = split /;/,$linked_genes[$i];
        my $tree_gene;
        # if tree for combined, combined get it
        my $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE gene='combined' AND taxon='combined'") or die;
        $sth->execute();
        $tree=$sth->fetchrow_array();
        $sth->finish();
        my $j=0;
        # else if tree for first gene, combined get it
        while ((!$tree or !($tree=~ /\(.+\).*;/)) and $j < scalar @genes) {
            my $string = $dbh->quote($genes[$j]);
            $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE gene=$string AND taxon='combined'") or die;
            $sth->execute();
            $tree=$sth->fetchrow_array();
            $sth->finish();
            $tree_gene = $genes[$j];
            ++$j;
        }
        # switch names to taxon names
        if ($tree and $tree =~ /\(.+\).*;/ and !($tree =~ /Taxon_[1-9]+/)) {
            open TREE, ">XXXtree.tree" or die "Could not open XXXtree.tree: $!.\n";
            print TREE $tree;
            close TREE or die;
            my @tip_names = `${path}$External_program::treebender -t '\\n' < XXXtree.tree`;
            chomp(@tip_names);
            my $switch_names = '';
            foreach (@tip_names) {
                $_ =~ s/^'|'$//g;
                #print "$_\n";
                if ($_ =~ /^[A-Z_]{1,3}[0-9]+/) {
                    $sth = $dbh->prepare("SELECT taxon_name FROM alignments WHERE $tree_gene\_accno='$_'") or die;
                    $sth->execute();
                    my $taxon_name = $sth->fetchrow_array();
                    if ($taxon_name =~ /^Taxon_[0-9]+/) {
                        $switch_names .= "$_|$taxon_name,";
                    }
                    else { die "Could not find a taxon name for $_.\n"; }
                }
                elsif (!($_ =~ /^Taxon_[0-9]+/)) {
                    die "Cannot identify what type of taxon name '$_' is.\n";
                }
            }
            $switch_names =~ s/,$//;
            $tree=`${path}$External_program::treebender -c '$switch_names' < XXXtree.tree`;
            unlink 'XXXtree.tree';
        }
        my $interleaved = 'n';
        if ($tree_method eq 'fasttree') { $interleaved = 'y'; }
        PifCosm_support_subs::print_alignment($database,'phylip','y','XXXtemp','taxon_name',$interleaved,4000,'n',@genes); # print phylipfile
        open PHYLIP, "<XXXtemp.phy" or die "Could not read the alignment (XXXtemp.phy) for $linked_genes[$i]: $!\n";
        my $ntaxa;
        while (my $temp = <PHYLIP>) {
            if ($temp =~ /^([0-9]+)\s+[0-9]+/) {
                $ntaxa=$1;
                last;
            }
        }
        if ($ntaxa < 4) {
            print "To few taxa for $linked_genes[$i]. No refinment possible. Proceding to next group of linked genes.\n";
            unlink "XXXtemp.phy";
            next;
        }
        my $guideML;
        if ($tree and $tree=~ /\(.+\).*;/) {
            print "Found starting tree.\n";
            open TREE, ">XXXtree.tree" or die "Could not open XXXtree.tree: $!.\n";
            print TREE $tree;
            close TREE or die;
            my @tip_names = `${path}$External_program::treebender -t '\\n' < XXXtree.tree`;
            chomp(@tip_names);
            my %present_hash;
            while (@tip_names) { $present_hash{shift @tip_names} = 'n'; }
            my $add_taxa = 'n';
            open ALIGNMENT, "<XXXtemp.phy" or die "Could not open XXXtemp.phy: $!.\n";
            my $need_to_add = 'n';
            while (my $row = <ALIGNMENT>) {
                chomp ($row);
                if ($row =~ /(Taxon_[0-9]+)/) {
                    my $present='n';
                    my $taxon = $1;
                    foreach (keys %present_hash) {
                        my $temp = $_;
                        $temp =~ s/^'|'$//g;
                        if ($temp eq $taxon) { $present_hash{$_} = 'y'; $present = 'y'; }
                    }
                    if ($present eq 'n') { $need_to_add = 'y'; }
                }
            }
            close ALIGNMENT or die;
            my $drop_tips;
            foreach (keys %present_hash) {
                if ($present_hash{$_} eq 'n') {
                    if ($drop_tips) { $drop_tips .= ",$_"; }
                    else { $drop_tips = $_; }
                }
            }
            if ($drop_tips) {
                print "Removing superfluous taxa.\n";
                $drop_tips =~ s/'/\'/g;
                $drop_tips =~ s/"/\"/g;
                open TREE, "|${path}$External_program::treebender -d '$drop_tips' > XXXtree.tree" or die "Cannot pipe to treebender: $!\n";
                print TREE $tree;
                close TREE or die;
            }
            if ($need_to_add eq 'y') {
                print "Adding missing taxa.\n";
                unlink glob "RAxML_*.XXXtemp";
                my @raxml = `${path}$External_program::raxml -s XXXtemp.phy -m GTRGAMMA -n XXXtemp -p 54321-f p -t XXXtree.tree`;
                foreach (@raxml) {
                    if ($_ =~ /Final GAMMA-based Score of best tree (-{0,1}[0-9.]+)/) {
                        $guideML = $1;
                    }
                }
                unlink "XXXtree.tree";
                if (-e "RAxML_result.XXXtemp") {
                    rename "RAxML_result.XXXtemp", "XXXguide.tree";
                    unlink glob "RAxML_*.XXXtemp";
                }
                elsif (-e "RAxML_parsimonyTree.XXXtemp" )  {
                    rename "RAxML_parsimonyTree.XXXtemp", "XXXguide.tree";
                    unlink glob "RAxML_*.XXXtemp";
                }
                else { die "Could not find start tree (RAxML_result.XXXtemp) with added taxa.\n"; }
            }
            else { rename "XXXtree.tree", "XXXguide.tree"; }
            if (!$guideML or $tree_method eq 'fasttree') {
                print "Re-estimating branch lengths and getting ML score for guide tree.\n";
                unlink glob "RAxML_*.XXXtemp";
                if ($tree_method eq 'fasttree') {
                    system `${path}$External_program::fasttree -quiet -nosupport -log XXXlog.txt -nt -gtr -gamma -mllen -nome -intree XXXguide.tree XXXtemp.phy > XXXtemp.tree`;
                    rename "XXXtemp.tree","XXXguide.tree";
                    open LOGFILE, "<XXXlog.txt" or die "Could not open logfile to read ML score: $!.\n";
                    while (my $row = <LOGFILE>) {
                        if ($row =~ /Gamma20LogLk\s+(-{0,1}[0-9\.]+)/) { $guideML = $1; }
                    }
                    close LOGFILE or die;
                    unlink "XXXlog.txt";
                }
                else {
                    my @raxml=`${path}$External_program::raxml -f e -s XXXtemp.phy -m GTRGAMMA -n XXXtemp -t XXXguide.tree`; # optimize branchlengths on parsimony tree
                    foreach (@raxml) {
                        if ($_ =~ /Final GAMMA  likelihood: (-{0,1}[0-9\.]+)/) {
                            $guideML = $1;
                        }
                    }
                    if (-e "RAxML_result.XXXtemp") {
                        rename "RAxML_result.XXXtemp", "XXXguide.tree";
                        unlink glob "RAxML_*.XXXtemp";
                    }
                    else { die "Could not find start tree with new branch lengths (RAxML_result.XXXtemp).\n"; }
                }
                if ($guideML) { print "The likelihood of the guide tree is: $guideML\n"; }
                else { die "Could not find likelihood score for guide tree.\n" }
            }
        }
        else { # if no tree make one
            print "Did not find start tree... creating one.\n";
            my $partition_file;
            if (scalar @genes > 1) { $partition_file = 'XXXtemp.partitions.txt'; }
            else { $partition_file = 'n'; }
            my $guide_tree;
            ($guide_tree,$guideML) = PifCosm_support_subs::run_raxml("XXXtemp.phy", $path, "XXXtemp",$n_threads,$tree_method,$partition_file);
            if (-e $guide_tree and $guideML) {
                print "Created start tree. The likelihood is: $guideML\n";
                rename $guide_tree, "XXXguide.tree";
            }
            else { die "Could not create start tree. Quitting.\n"; }
        }
        unlink "XXXtemp.phy","XXXtemp.phy.reduced","XXXtemp.partitions.txt*";
        if (-e "XXXguide.tree") {
            my $i_since_best=0;
            my $iterations = 0;
            my $best_ML=$guideML;
            if (!$best_ML) { $best_ML = -9**9**9; }
            open TREE, "<XXXguide.tree" or die "Could not read XXXguide.tree: $!.\n";
            my $best_tree=<TREE>;
            close TREE or die;
            my $temp_tree="XXXguide.tree";
	    my $boot_file='empty';
            while (1) {
                print "Starting refinment itteration $iterations. It was $i_since_best itterations since the highest ML score.\n";
                ++$iterations;
                ++$i_since_best;
                my $temp_ML;
                for(my $i=0; $i<scalar @genes; ++$i) { # make alignment for each gene
                    $sth=$dbh->prepare("SELECT taxon_name FROM alignments WHERE $genes[$i]_sequence !='empty'") or die;
                    $sth->execute(); # get taxa for specific gene
                    my $n_taxa=0; # count the number of taxa
                    my $keep_tips=''; # string for the tips to keep
                    while (my $tip_name = $sth->fetchrow_array()) { # for each row returned from database query
                        $keep_tips .= "$tip_name,"; # store the taxon_name so we can keep it in the tree
                        ++$n_taxa; # count up the number of taxa
                    }
                    if ($n_taxa < 2) { last; } # if no or one sequence no need to align
                    $keep_tips =~ s/,$//; # get rid of last comma in keep tips, however I do not think this is nessesary
                    $split_switches = "longest_branch:n --split_stop t:2";
                    my $pruned_tree = `${path}$External_program::treebender -d $keep_tips -i < $temp_tree`;
                    &divide_and_align ($dbh,  $genes[$i], $max_size, $split_switches, 0, $use_guide_tree, $path, $pruned_tree);
                    rename "XXXalignment_0.fst", "XXXbuild_$genes[$i].fst";
                }
                my $partition_file;
                #my $alignment_file;
                open ALIGNMENT, "<XXXbuild_$genes[0].fst" or die "Could not open XXXbuild_$genes[0].fst: $!.\n";
                my $sequence_length=0; # count sequence length
                my $length_lable=0; # length of the taxon names
                my %new_seq; # hash to store new sequences
                while (my $row = <ALIGNMENT>) { # read first alignment row by row
                    chomp($row); # get rid of new line
                    if ($row =~ s/^>//) { # if sequence name remove leading >
                        $new_seq{$row}='empty'; # and initiate sequence in hash
                        if (length($row) > $length_lable) { $length_lable=length($row); } # store longest name
                        $sequence_length=0; # null the sequence length
                    }
                    else {
                        $row =~ /\s+/g; # get rid of white space
                        $sequence_length += length($row); # get length
                    }
                }
                close ALIGNMENT or die;
                open FASTA, "<XXXbuild_$genes[0].fst" or die "Could not open XXXbuild_$genes[0].fst: $!.\n"; # open alignment file for reading
                open PHYLIP, ">XXXconcatenated.phy" or die "Could not open XXXconcatenated.phy: $!.\n"; # open phylip file for writing
                print PHYLIP (scalar keys %new_seq) . " " . $sequence_length; # print the number of species and sequence length
                while (my $row = <FASTA>) { # for each row in the fasta file
                    chomp($row); # get rid of new line
                    if ($row =~ s/^>//) { # get rid of leading > in file name
                        print PHYLIP "\n$row" . (" " x (($length_lable+1)-length($row))); # print file name and make even seq rows
                    }
                    else {
                        $row =~ /\s+/g;
                        print PHYLIP uc($row); # print sequence
                    }
                }
                print PHYLIP "\n";
                close FASTA or die;
                close PHYLIP or die;
                if (scalar @genes > 1) { # if more than one gene
                    open PARTITIONFILE, ">XXXtemp.partitions.txt" or die "Could not open XXXtemp.partitions.txt: $!.\n";
                    print PARTITIONFILE "DNA, $genes[0] = 1-$sequence_length\n"; #write first partition to partition file
                    for (my $n=1; $n < scalar @genes; ++$n) { # for each gene
                        my $temp_sequence_length=0; # length of present gene
                        my $tip; # name of taxa the sequence belong to
                        open FASTA, "<XXXbuild_$genes[$n].fst" or die "Could not open XXXbuild_$genes[$n].fst: $!.\n";
                        while (my $row = <FASTA>) { # for each row in the alignment file
                            chomp($row); # get rid of new line
                            if ($row =~ s/^>//) { # get rid of leading >
                                $tip=$row; # store name so we know which taxa the sequence belong to
                                $new_seq{$tip}=''; # null the sequence
                                $temp_sequence_length=0; # null the sequence length for the gene
                                if (length($row) > $length_lable) { $length_lable=length($row); } # keep track of the longest taxon name
                            }
                            else {
                                $row =~ s/\s+//g; # get rid of white spaces
                                $new_seq{$tip} .= uc($row); # store sequence
                                $temp_sequence_length += length($row); # store length of sequence
                            }
                        }
                        close FASTA or die;
                        open INPHYLIP, "<XXXconcatenated.phy" or die "Could not open XXXconcatenated.phy: $!.\n";
                        open OUTPHYLIP, ">XXXtemp.phy" or die "Could not open XXXtemp.phy: $!.\n";
                        print OUTPHYLIP (scalar keys %new_seq) . " " . ($sequence_length+$temp_sequence_length) . "\n"; # print phylip header
                        while (my $row = <INPHYLIP>) { #for each row of the file with sequences from previous genes
                            if ($row =~ /[0-9]+\s+[0-9]+/) { next; } # if header row skip it
                            elsif ($row =~ /(\w+)\s+(.+)/) { # if containing name and sequence separated by white space 
                                print OUTPHYLIP $1 . (" " x (($length_lable+1)-length($1))); # print name and adjust for longest name
                                if ($new_seq{$1} eq 'empty') { # if no sequence is available for the gene
                                    print OUTPHYLIP $2 . "-" x $temp_sequence_length; # print gaps
                                }
                                else { # otherwise
                                    print OUTPHYLIP $2 . $new_seq{$1}; # print previous and present sequence
                                    $new_seq{$1} = 'empty'; # reset sequence
                                }
                                print OUTPHYLIP "\n";
                            }
                        }
                        close INPHYLIP or die;
                        foreach (keys %new_seq) { # for each taxon
                            if ($new_seq{$_} ne 'empty') { # if it has not already been added
                                print OUTPHYLIP $_ . (" " x (($length_lable+1)-length($_))) . ("-" x $sequence_length) . $new_seq{$_} . "\n"; # print name, gaps as long as previous sequences, and present sequence
                                $new_seq{$_} = 'empty';
                            }
                        }
                        close OUTPHYLIP or die;
                        rename "XXXtemp.phy", "XXXconcatenated.phy";
                        print PARTITIONFILE "DNA, $genes[$n] = " . ($sequence_length+1); # write partitioning
                        $sequence_length += $temp_sequence_length; # add sequence length
                        print PARTITIONFILE "-$sequence_length\n"; # finish printing partition
                    }
                    close PARTITIONFILE or die;
                    if ($tree_method =~ /fasttree/ and $sequence_length > 4000 and -e "XXXconcatenated.phy") {
                        unlink glob "XXXresidue_sequences_*.txt";
                        open SEQUENTIAL, "<XXXconcatenated.phy" or die "Could not open file XXXconcatenated.phy: $!.\n";
                        open INTERLEAVED, ">XXXresidue_sequences_0.txt" or die "Could not open XXXresidue_sequences_0.txt for writing: $!.\n";
                        while (my $row=<SEQUENTIAL>) {
                            if ($row =~ /(\w+\s+)(.+)/) {
                                print INTERLEAVED $1 . substr($2,0,4000) . "\n";
                                my $n_printed = 4000;
                                while (length($2)-$n_printed > 0) {
                                    open RESIDUES, ">>XXXresidue_sequences_$n_printed.txt" or die "Could not open temporary file XXXresidue_sequences_$n_printed.txt: $!.\n";
                                    print RESIDUES (' ' x length($1)) . substr($2,$n_printed,4000) . "\n";
                                    $n_printed += 4000;
                                    close RESIDUES;
                                }
                            }
                            elsif ($row =~ /[0-9]+\s+[0-9]+/) { print INTERLEAVED $row; }
                        }
                        close SEQUENTIAL or die;
                        close INTERLEAVED or die;
                        open OUTFILE, ">XXXconcatenated.phy" or die "Could not open file XXXconcatenated.phy: $!.\n";
                        foreach (glob "XXXresidue_sequences_*.txt") {
                            open SEQUENCES, "<$_" or die "Could not open temporary file $_: $!.\n";
                            while (my $row = <SEQUENCES>) {
                                print OUTFILE $row;
                            }
                            print OUTFILE "\n\n";
                            close SEQUENCES or die;
                        }
                        unlink glob "XXXresidue_sequences_*.txt";
                    }
                    $partition_file = 'XXXtemp.partitions.txt';
                }
                else {
                    $partition_file = 'n';
                    #$alignment_file = "XXXbuild_$genes[0].fst";
                }
                unlink $temp_tree;
                undef %new_seq;
                if (-e "XXXconcatenated.phy") {
		    my $n_boot;
		    if ($tree_method=~ s/_bootstrap_([0-9]+)//) { $n_boot = $1; } # if doing bootstraps
                    ($temp_tree,$temp_ML,$boot_file) = PifCosm_support_subs::run_raxml("XXXconcatenated.phy", $path, "XXXtemp",$n_threads,$tree_method,$partition_file,$n_boot);
                    unlink glob "XXXconcatenated.phy*";
                    unlink glob "XXXtemp.partitions.txt*";
                }
                else {
		    unlink glob "XXXconcatenated.phy*";
		    unlink glob "XXXtemp.partitions.txt*";
		    die "Could not create concatenated alignment.\n";
		}
                if ($temp_tree eq 'empty' or !(-e $temp_tree)) {
                    die "Could not construct tree.\n";
                }
                # compare tree to previous trees
                if ((!$best_ML and $temp_ML) or ($temp_ML and ($temp_ML > $best_ML))) {
                    print "    Found better tree ($temp_ML compared to $best_ML)!\n";
                    $best_ML = $temp_ML;
                    open TREE, "<$temp_tree" or die "Could not read ML tree: $!.\n";
                    $best_tree = <TREE>;
                    close TREE or die;
                    foreach (@genes) { # save alignments
                        my ($ntaxa,$sequence_length) = PifCosm_support_subs::read_fasta_to_alignments_table("XXXbuild_$_.fst",$dbh,$_,'taxon_name');
                        print "    Saved $ntaxa aligned sequences, $sequence_length bp long, to $_ alignment.\n";
                    }
                    $i_since_best=0;
                }
                else { print "    New likelihood not better than previous ($temp_ML compared to $best_ML).\n"; }
                unlink glob "XXXbuild_*.fst"; #remove alignment files
                # chack if stop criteria reached
                if ($stop_crit eq 'max') {
                    if ($i_since_best >= $n_iterations) {
                        last;
                    }
                }
                elsif ($stop_crit eq 'iterations') {
                    if ($iterations >= $n_iterations) {
                        last;
                    }
                }
                else {
                    last;
                }
            }
            if (-e $temp_tree) { unlink $temp_tree; }
            if ($best_tree =~ /^\(.+\).*;/) { # save best tree
                print "Saving tree for $linked_genes[$i] in table alignment_groups.\n";
                #print "$best_tree\n";
		my $method='ML';
		if ($best_tree =~ /parsimony/) { $method = 'MP'; }
		PifCosm_support_subs::read_tree_to_alignment_groups($best_tree,$dbh,'combined',$linked_genes[$i],$method);
                my $string = $dbh->quote($best_tree);
                my $updates = $dbh->do("UPDATE alignment_groups SET tree=$string WHERE taxon='combined' AND gene='$linked_genes[$i]'") or die;
                if ( $updates != 1 ) {
                    $updates = $dbh->do("INSERT INTO alignment_groups (taxon,gene,tree) VALUES ('combined','$linked_genes[$i]',$string)") or die;
                }
                print "Updated $updates rows.\n";
            }
	    if ($save_boot_trees eq 'y') {
		if (-e $boot_file) {
		    PifCosm_support_subs::read_tree_to_alignment_groups($boot_file,$dbh,'combined_boot',"$linked_genes[$i]_boot","ML_boot"); # save boot trees to database
		    unlink $boot_file;
		}
	    }
	    else { unlink $boot_file; }
        }
        else { die "No start tree.\n"; }
    }
    $dbh->disconnect();
}

sub divide_and_align {
    no warnings 'recursion';
    my $dbh = shift @_;
    my $gene = shift @_;
    my $max_size = shift @_;
    my $split_switches = shift @_;
    my $number = shift @_; # number to give different alignment files different names
    $number++; #make sure it is higher than previous run
    my $use_guide_tree = shift @_;
    my $path = shift @_;
    my $tree = \$_[0];
    my @groups = `echo "$$tree" | ${path}$External_program::treebender --split $split_switches`;
    my $n_aligned=0;
    for (my $i=0; $i < scalar @groups; $i++) {
        chomp ($groups[$i]);
        my $temp = `echo "$groups[$i]" | ${path}$External_program::treebender -t`;
        chomp ($temp); # get rid of new line
        my @tip_labels = split /,/, $temp; # separate the tip lables
        if (scalar @tip_labels < $max_size) {
            open ALIGNMENT, ">XXXtemp.fst" or die "Could not open XXXtemp.fst: $!.\n"; # open alignment file
            for (my $n=1; $n <= scalar @tip_labels; ++$n) { # For each number that the taxa should be replased with in the mafft guide tree
                my $string = $dbh->quote($tip_labels[$n-1]);
                my $sth = $dbh->prepare("SELECT $gene\_sequence FROM alignments WHERE taxon_name=$string") or die;
                $sth->execute();
                my $sequence = $sth->fetchrow_array(); #get sequence
                $sth->finish();
                print ALIGNMENT ">$tip_labels[$n-1]\n$sequence\n"; # print sequence
                $groups[$i]=~ s/$tip_labels[$n-1]/$n/; # replace taxon name with number in guide tree
            }
            close ALIGNMENT or die;
            if (scalar @tip_labels > 1) {
                my $mafft_switches = '--quiet';
                if ($use_guide_tree eq 'y' and scalar @tip_labels > 3) {
                    my $mafft_tree=PifCosm_support_subs::newick_to_mafft($groups[$i],1); # get guide tree in mafft format
                    open GUIDETREE, ">XXXmafft.tree" or die "Could not open XXXmafft.tree: $!.\n";
                    print GUIDETREE $mafft_tree; # print guide tree to file
                    close GUIDETREE or die;
                    $mafft_switches .= " --treein XXXmafft.tree";
                }
                system "${path}$External_program::mafft $mafft_switches XXXtemp.fst > XXXalignment_$number.fst"; # align the group
                if ($use_guide_tree eq 'y' and scalar @tip_labels > 3) { unlink "XXXmafft.tree"; }
            }
            else { rename "XXXtemp.fst", "XXXalignment_$number.fst"; } # if only one sequence no need to align
        }
        else {
            &divide_and_align($dbh,$gene,$max_size,$split_switches,$number,$use_guide_tree,$path,$groups[$i]);
        }
        if ($i>0) {
            if (-e "XXXbuilt_$number.fst") {
                if (scalar @tip_labels < 2) {
                    system "${path}$External_program::mafft --quiet --add XXXalignment_$number.fst XXXbuilt_$number.fst > XXXtemp.fst";
                    rename "XXXtemp.fst", "XXXbuilt_$number.fst";
                    unlink "XXXalignment_$number.fst";
                }
                elsif ($n_aligned < 2) {
                    system "${path}$External_program::mafft --quiet --add XXXbuilt_$number.fst XXXalignment_$number.fst > XXXtemp.fst";
                    rename "XXXtemp.fst", "XXXbuilt_$number.fst";
                    unlink "XXXalignment_$number.fst";
                }
                else {
                    if ($n_aligned >= scalar @tip_labels ) {
                        #system "${path}mafft --quiet --addprofile XXXalignment_$number.fst XXXbuilt_$number.fst > XXXtemp.fst"; # add sequences to previous alignment
                        system "${path}$External_program::muscle -quiet -profile -in1 XXXbuilt_$number.fst -in2 XXXalignment_$number.fst -out XXXtemp.fst"; # add sequences to previous alignment
                    }
                    else { system "${path}$External_program::muscle -quiet -profil -in1 XXXalignment_$number.fst -in2 XXXbuilt_$number.fst -out XXXtemp.fst"; } # add sequences to previous alignment
                    #else { system "${path}mafft --quiet --addprofil XXXbuilt_$number.fst XXXalignment_$number.fst > XXXtemp.fst"; } # add sequences to previous alignment
                }
                rename "XXXtemp.fst", "XXXbuilt_$number.fst";
                unlink "XXXalignment_$number.fst";
                $n_aligned += scalar @tip_labels;
            }
            else { die "File error!!! Could not find XXXbuilt.fst.\n" }
        }
        else {
            rename "XXXalignment_$number.fst", "XXXbuilt_$number.fst";
            $n_aligned += scalar @tip_labels;
        }
    }
    rename "XXXbuilt_$number.fst", "XXXalignment_" . ($number-1) . ".fst";
}

1;
