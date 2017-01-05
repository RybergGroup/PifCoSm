package Align_sequences;
use strict;
use PifCosm_support_subs;

#############################################################################
### Subroutine to align the separate genes and produce phylogenetic trees ###
#############################################################################

sub align_sequences {
    # General variables
    my $database = shift @_; #variable to store database name
    my $path = shift @_; # variable with the path to alignment programs
    my $n_threads = shift @_; # store the number of threads to use, 0 means no multi threading
    my $create_tree = shift @_; # y if tree should be created, else n
    my $tree_method = shift @_; #'Light,RAx-br';
    my $store_boot = shift @_; # if 'y' store bootstraped trees
    my $rm_outliers = shift @_; # if 'y' outliers will be removed
    my @linked_genes = @_; # array to store linked genes 
    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$ARGV[0]","","") or die "Couldn't connect to database: " . DBI->errstr;
    my $sth;
    my @tables = PifCosm_support_subs::get_gene_tables($dbh);
    my %groups_no; # hash to store the number of alignment groups there are for each gene
    my %groups_taxa; # hash to store the alignment groups
    # get the number of and alignment groups for each gene
    foreach (@tables) {
        $sth = $dbh->prepare("SELECT taxon,alignable FROM alignment_groups WHERE gene='$_'") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        while ( my @row = $sth->fetchrow_array() ) {
            if ($row[0] ne "combined" and $row[1]) {
                if ($groups_taxa{$_}) { $groups_taxa{$_} .= "\|$row[0],$row[1]"; }
                else { $groups_taxa{$_} = "$row[0],$row[1]"; }
                ++$groups_no{$_};
            }
        }
    }
    my @rest_groups; # array to store alignment groups we were not able to align
    my @ordered_keys = PifCosm_support_subs::keys_in_order(%groups_no);
    foreach (@ordered_keys) {
        print "*" x (12+length($_)) . "\n";
        print " Aligning $_.\n";
        print "*" x (12+length($_)) . "\n";
        my @groups = split /\|/, $groups_taxa{$_}; # get individual groups
        for (my $i=0; $i < scalar @groups; ++$i) {
            my @entries = split /,/, $groups[$i]; # get taxon and if it is alignable
            if ($entries[1] > 0) { # if alignable
                &align($dbh,$path,$_,$entries[0],$n_threads); # align
            }
            else { # if not alignable
                print "Sequences $_ for $entries[0] cannot be easily aligned, getting guide tree from linked region.\n";
                my ($guide_gene) = &get_guide_tree_gene($dbh,$_,$entries[0],@linked_genes); # get guide tree
                if ($guide_gene) { # if we got a guide tree
                    print "Will use $guide_gene for guide tree.\n";
                    my $ntaxa = &align_using_guide_tree ($dbh,$path,$_,$entries[0],$tree_method,$guide_gene,$n_threads); # align using guide tree
                    if (!$ntaxa or $ntaxa < 1) { # if no alignment
                        print "Unable to create alignment using guide tree. Will return to $_ for $entries[0] later.\n";
                        push(@rest_groups, "$_ | $entries[0]"); # save for later
                    }
                }
                else { # if no guide tree
                    print "Non of the genes had an apropriate tree. Will return to $_ for $entries[0] later.\n";
                    push(@rest_groups, "$_ | $entries[0]"); # save it for later
                }
            }
        }
    }
    ### Try to align unaligned groups again
    foreach (@rest_groups) {
        my ($gene,$taxon) = split / \| /, $_;
        print "*" x (24+length($gene)+length($taxon)) . "\n";
        print " Revisiting $gene, aligning $taxon.\n";
        print "*" x (24+length($gene)+length($taxon)) . "\n";
        print "Trying to align using guide tree.\n";
        my ($guide_gene) = &get_guide_tree_gene($dbh,$gene,$taxon,@linked_genes); # see if we find a guide tree
        my $ntaxa=0;
        if ($guide_gene) { # if we find a guide tree
            print "Found guide tree ($guide_gene).\n";
            $ntaxa = &align_using_guide_tree ($dbh,$path,$gene,$taxon,$tree_method,$guide_gene,$n_threads); # align using guide tree
        }
        if ($ntaxa < 1) { # if no guide tree or alignment with guide tree failed
            print "No alignment using guide tree, will align ordinarily.\n";
            &align($dbh,$path,$gene,$taxon,$n_threads); # align without guide tree
        }
    }

    ### Merge alignments ###
    foreach (@ordered_keys) {
        print "*" x (26+length($_)) . "\n";
        print " Merging alignments for $_.\n";
        print "*" x (26+length($_)) . "\n";
        my @groups = split /\|/, $groups_taxa{$_}; # get the separate groups
        my $n_taxa=0;
        if (scalar @groups > 1) {
            my @aligned_groups; # array to store the taxa for which there are sequences
            # Get the number of taxa for which there are sequences
            my $taxon_tree;
            for (my $i=0; $i < scalar @groups; ++$i) { # for each alignment group
                my @entries = split /,/, $groups[$i]; # get the taxon name
                $entries[0] =~ s/_/ /g;
                $sth = $dbh->prepare("SELECT DISTINCT taxon_string FROM gb_data WHERE (taxon_string LIKE '%; $entries[0];%' OR taxon_string LIKE '$entries[0];%' OR taxon_string LIKE '%; $entries[0]' OR taxon_string='$entries[0]')")
                    or die "Could not prepare statement: " . $dbh->errstr;
                $sth->execute();
                my $taxon_string = $sth->fetchrow_array();
                if (!($taxon_string =~ /\w/)) { die "Could not find the taxonomic hierarchy for $entries[0]. Quitting.\n"; }
                my @taxon_array = split /; /, $taxon_string;
                while ($taxon_array[-1] and $taxon_array[-1] ne $entries[0]) { pop @taxon_array; }
                my $flag;
                if (@taxon_array and (scalar @taxon_array) > 0) { $taxon_tree = &insert_taxon_in_taxon_tree($taxon_tree,@taxon_array); }
            }
            #print "$taxon_tree;\n";
            for (my $i=0; $i < scalar @groups; ++$i) { # for each alignment group
                my @entries = split /,/, $groups[$i]; # get the taxon name
                if ($entries[0] eq "combined") { next; }
                my $temp=$entries[0];
                $temp =~ s/_/ /g;
                $sth = $dbh->prepare("SELECT COUNT(alignments.$_\_accno) FROM alignments INNER JOIN gb_data ON alignments.$_\_accno=gb_data.accno WHERE alignments.$_\_sequence!='empty' AND (gb_data.taxon_string LIKE '%; $temp;%' OR gb_data.taxon_string LIKE '$temp;%' OR gb_data.taxon_string LIKE '%; $temp' OR gb_data.taxon_string='$temp')")
                    or die "Could not prepare statement: " . $dbh->errstr;
                $sth->execute();
                my $ntaxa = $sth->fetchrow_array(); # get number of seq
                $sth->finish();
                if ($ntaxa and $ntaxa > 0) { push (@aligned_groups,$entries[0]); } # store the taxon name
                else { print "Missing alignment for $temp. $_ sequences of this taxon will therefore be concidered singletons.\n"; }
            }
            # See if it all adds up, or if there are singletons
            {
                $sth = $dbh->prepare("SELECT COUNT($_\_accno) FROM alignments WHERE $_\_accno!='empty' AND $_\_sequence='empty'") or die "Could not prepare statement: " . $dbh->errstr;
                $sth->execute();
                my $ntaxa = $sth->fetchrow_array(); # get the number of sequences that should had been included but are not
                $sth->finish();
                if ($ntaxa and $ntaxa > 0) { # if there are sequences missing
                    push (@aligned_groups,'singletons'); # save it
                    print "Found $ntaxa singletons that will be aligned to the rest of the taxa.\n";
                }
                else { print "No sigletons.\n"; }
            }
            if ($taxon_tree) { $n_taxa=&merge_alignments_according_to_tree($dbh,$path,$taxon_tree,$_,$n_threads); }
            else { die "Could not merge alignments for $_ since no taxonomic hierarchy was available.\n"; }
            print "Mearged $n_taxa sequences into one alignment for $_.\n";
        }
        elsif (scalar @groups == 1) {
            print "Only one alignment group, no need to merge.\n";
            $sth = $dbh->prepare("SELECT COUNT($_\_accno) FROM alignments WHERE $_\_accno != 'empty' AND $_\_sequence = 'empty'") or die;
            $sth->execute();
            my $n_singletons = $sth->fetchrow_array();
            $sth->finish();
            if ($n_singletons > 0) {
                print "Found $n_singletons singletons. They will be added to the alignment.\n";
                $sth = $dbh->prepare("SELECT alignments.$_\_accno,$_.sequence FROM alignments INNER JOIN $_ ON alignments.$_\_accno=$_.accno WHERE alignments.$_\_accno != 'empty' AND alignments.$_\_sequence = 'empty'") or die;
                my $n_temp = PifCosm_support_subs::print_fasta("XXXsingletons.fst",$sth);
                if ($n_singletons != $n_temp) { unlink "XXXsingletons.txt"; die "Could not write a proper fasta file for singletons.\n"; }
                $sth = $dbh->prepare("SELECT $_\_accno,$_\_sequence FROM alignments WHERE $_\_accno != 'empty' AND $_\_sequence != 'empty'") or die;
                my $n_aligned = PifCosm_support_subs::print_fasta("XXXaligned_$_.fst",$sth);
                if (-e "XXXsingletons.fst" and -e "XXXaligned_$_.fst") {
                    if ($n_threads > 0) { # if using threads
                        my $threads = $n_threads;
                        if ($n_threads > 1) { --$threads; } # only give extra threads
                            system "${path}mafft --quiet --thread $threads --add XXXsingletons.fst XXXaligned_$_.fst > XXXtemp_temp.fst"; # align
                    }
                    else { system "${path}mafft --quiet --add XXXsingletons.fst XXXaligned_$_.fst > XXXtemp_temp.fst"; } # in no threads align ordinarily
                    unlink ("XXXsingletons.fst","XXXaligned_$_.fst");
                    if (-e "XXXtemp_temp.fst") {
                        my ($ntaxa,$seqlength) = PifCosm_support_subs::read_fasta_to_alignments_table("XXXtemp_temp.fst",$dbh,$_,'accno'); # read alignment to database
                        unlink "XXXtemp_temp.fst";
                    }
                    else { die "Failed to add singletons to alignment.\n"; }
                }
                else { unlink ("XXXsingletons.fst","XXXaligned_$_.fst"); die "Failed to write apropriate fasta files.\n"; }
            }
            $sth = $dbh->prepare("SELECT COUNT($_\_accno) FROM alignments WHERE $_\_sequence != 'empty'") or die;
            $sth->execute();
            $n_taxa = $sth->fetchrow_array();
            $sth->finish();
            print "Aligned $n_taxa sequences for $_.\n";
        }
        else { die "Could not find what groups have been aligned for $_. Quitting.\n"; }
        if ($n_taxa > 3 and $create_tree eq "y") {
            $sth = $dbh->prepare("SELECT COUNT($_\_accno),LENGTH($_\_sequence) FROM alignments WHERE $_\_sequence!='empty' GROUP BY LENGTH($_\_sequence)") or die "Could not prepare statement: " . $dbh->errstr;
            $sth->execute();
            my $seqlength;
            my $n = 0;
            while (my @row=$sth->fetchrow_array()) {
                if ($row[0]>$n) {
                    $n = $row[0];
                    $seqlength = $row[1];
                }
            }
            $sth->finish();
            if ( $n != $n_taxa ) { die "Something wrong with the alignment for $_. Quitting.\n"; }
            # get sequences from database
            $sth = $dbh->prepare("SELECT $_\_accno,$_\_sequence FROM alignments WHERE $_\_sequence!='empty'") or die "Could not prepare statement: " . $dbh->errstr;
            print "Write phylip alignment ($n_taxa taxa, $seqlength bp).\n";
            my $miss = PifCosm_support_subs::print_phylip ( "XXXtemp.phy",$sth,$n_taxa,$seqlength); # write phylip alignment
            if ($miss == 0) { # if right number of sequences written
                print "Create tree for combined data.\n";
                my $sth2 = $dbh->prepare("SELECT COUNT(taxon) FROM alignment_groups WHERE taxon='combined' AND gene='$_'");
                $sth2->execute();
                my $number = $sth2->fetchrow_array();
                if ($number!=1) {
                    my $string = $dbh->quote($_);
                    $dbh->do("INSERT INTO alignment_groups (taxon,gene) VALUES ('combined',$string)") or die "Could not insert value in alignment_groups: " . $dbh->errstr; # make room for tree in database
                }
                &create_and_save_tree("XXXtemp.phy",$path,$dbh,$_,"combined",$n_threads,$tree_method,$store_boot); # create and save a tree for the gene
                if ($rm_outliers eq 'y') { &remove_outliers ($dbh, $path, $_, 'combined'); }
            }
            else { print "Faild to write a correct phylip file. Missed with $miss taxa.\n"; }
            unlink "XXXtemp.phy","XXXtemp.phy.reduced"; # remove temporary file
        }
        else {
            print "Less than four taxa, skip creating tree.\n";
        }
    }
    $dbh->disconnect();
}
sub merge_alignments_according_to_tree {
    my $dbh = shift @_;
    my $path = shift @_;
    my $taxon_tree = shift @_;
    my $gene = shift @_;
    my $n_threads = shift @_;
    my $sth;
    my $ntaxa=0;
    print "Merging separate $gene sequence alignments.\n";
    while ($taxon_tree) {
        if ($taxon_tree =~ /(.*?)\(([^\(\)]*)\)\[([^\[\]]*)\](.*)/) {
            my $left = $1;
            my @taxa = split /,/,$2;
            my $parent = $3;
            my $right = $4;
            my %number_taxa;
            my $query_parent = $parent;
            $query_parent =~ s/'/''/g;
            print "Merging ";
            foreach (@taxa) { print $_ , ', '; }
            print "into $parent.\n";
            for ( my $i=0; $i< scalar @taxa; ++$i) {
                my $taxa = $taxa[$i];
                $taxa =~ s/'/''/g;
                $sth = $dbh->prepare("SELECT alignments.$gene\_accno,alignments.$gene\_sequence FROM alignments INNER JOIN gb_data ON alignments.$gene\_accno=gb_data.accno WHERE alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $taxa;%' OR gb_data.taxon_string LIKE '$taxa;%' OR gb_data.taxon_string LIKE '%; $taxa' OR gb_data.taxon_string='$taxa')")
                    or die "Could not prepare statement: " . $dbh->errstr;
                $taxa[$i] =~ s/\W/_/g;
                $number_taxa{$taxa[$i]} = PifCosm_support_subs::print_fasta("XXX$gene$taxa[$i].fst",$sth); # write sequences for taxa
            }
            my @keys;
            if (%number_taxa) {
                @keys = PifCosm_support_subs::keys_in_order(%number_taxa);
            }
            else { unlink glob "XXX$gene*.fst"; die "Could not find taxa within $parent.\n"; }
            unlink "XXXbuild.fst";
            my $all_singletons = 'y';
            foreach ( @keys ) { if ($number_taxa{$_} > 1) { $all_singletons = 'n'; last;} }
            if ( $all_singletons ne 'y') {
                for (my $i=0; $i < scalar @keys; ++$i) { # for all taxa to be merged
                    if ($number_taxa{$keys[$i]} < 1) { # if no sequences in taxon delete file
                        unlink "XXX$gene$keys[$i].fst";
                        next;
                    }
                    if ($i==0 or !(-e "XXXbuild.fst")) { # if it is the first taxon to merge move it directly to the merging file
                        if (-e "XXX$gene$keys[$i].fst") { rename "XXX$gene$keys[$i].fst", "XXXbuild.fst"; }
                    }
                    elsif (-e "XXX$gene$keys[$i].fst") { # otherwise add the taxon to the merging file
                        #system "${path}mafft --addprofile XXX$gene$keys[$i].fst XXXbuild.fst > XXXtemp_temp.fst"; # add sequences to previous alignment
                        system "${path}muscle -quiet -profile -in1 XXXbuild.fst -in2 XXX$gene$keys[$i].fst -out XXXtemp_temp.fst"; # add sequences to previous alignment
                        rename "XXXtemp_temp.fst", "XXXbuild.fst";
                        unlink "XXX$gene$keys[$i].fst";
                    }
                    else { print "Lacking one or both files to merge. Proceeding reluctantly.\n"; }
                }
                $sth = $dbh->prepare("SELECT alignments.$gene\_accno,$gene.sequence FROM alignments INNER JOIN gb_data ON alignments.$gene\_accno=gb_data.accno INNER JOIN $gene ON alignments.$gene\_accno=$gene.accno WHERE alignments.$gene\_accno!='empty' AND alignments.$gene\_sequence='empty' AND (gb_data.taxon_string LIKE '%; $query_parent;%' OR gb_data.taxon_string LIKE '$query_parent;%' OR gb_data.taxon_string LIKE '%; $query_parent' OR gb_data.taxon_string='$query_parent')") # get singleton sequences
                    or die "Could not prepare statement: " . $dbh->errstr;
                my $n_taxa = PifCosm_support_subs::print_fasta("XXX${gene}_singletons.fst",$sth); # write sequences for taxa
                if ($n_taxa > 0) {
                    if (-e "XXXbuild.fst") {
                        #print "Found singletons for $parent. Adding them to alignment.\n";
                        if ($n_threads > 0) { # if using threads
                            my $threads = $n_threads;
                            if ($n_threads > 1) { --$threads; } # only give extra threads
                                system "${path}mafft --quiet --thread $threads --add XXX${gene}_singletons.fst XXXbuild.fst > XXXtemp_temp.fst"; # align
                        }
                        else { system "${path}mafft --quiet --add XXX${gene}_singletons.fst XXXbuild.fst > XXXtemp_temp.fst"; } # in no threads align ordinarily
                        rename "XXXtemp_temp.fst","XXXbuild.fst"; # move new alignment to the file buing built on
                    }
                    else {
                        if ($n_taxa > 1) {
                            print "WARNING!!! No previous alignment, for $query_parent, to add singletons too. Aligning singletons.\n";
                            if ($n_threads > 0) {
                                my $threads = $n_threads;
                                if ($n_threads > 1) { --$threads; }
                                print "Using mafft, pthread version.\n";
                                system "${path}mafft --auto --quiet --thread $threads XXX${gene}_singletons.fst > XXXbuild.fst";
                            }
                            else {
                                print "Using mafft, serial version.\n";
                                system "${path}mafft --auto --quiet XXX${gene}_singletons.fst > XXXbuild.fst";
                            }
                        }
                    }
                }
                unlink "XXX${gene}_singletons.fst";
            }
            else {
               unlink glob "XXX$gene*.fst"; # one seq per taxa at most, no need to keep them in separate files
               $sth = $dbh->prepare("SELECT alignments.$gene\_accno,$gene.sequence FROM alignments INNER JOIN gb_data ON alignments.$gene\_accno=gb_data.accno INNER JOIN $gene ON alignments.$gene\_accno=$gene.accno WHERE alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $query_parent;%' OR gb_data.taxon_string LIKE '$query_parent;%' OR gb_data.taxon_string LIKE '%; $query_parent' OR gb_data.taxon_string='$query_parent')") # get all the sequences for the heigher taxa
                    or die "Could not prepare statement: " . $dbh->errstr;
                my $parent_no_space = $parent;
                $parent_no_space =~ s/\W/_/g;
                my $n_taxa = PifCosm_support_subs::print_fasta("XXX${gene}$parent_no_space.fst",$sth); # write sequences for taxa
                if ($n_taxa > 0) {
                    if ($n_taxa > 1) {
                        print "All taxa within $parent are singletons. Aligning them.\n";
                        if ($n_threads > 0) {
                            my $threads = $n_threads;
                            if ($n_threads > 1) { --$threads; }
                            print "Using mafft, pthread version.\n";
                            system "${path}mafft --auto --quiet --thread $threads XXX${gene}$parent_no_space.fst > XXXbuild.fst";
                        }
                        else {
                            print "Using mafft, serial version.\n";
                            system "${path}mafft --auto --quiet XXX${gene}$parent_no_space.fst > XXXbuild.fst";
                        }
                        unlink "XXX${gene}$parent_no_space.fst";
                    }
                    else { rename "XXX${gene}$parent_no_space.fst","XXXbuild.fst"; }
                }
                else {
                    print "WARNING!!! No taxa in $parent.\n";
                }
            }
            if (-e "XXXbuild.fst") {
                my $seqlength;
                ($ntaxa,$seqlength) = PifCosm_support_subs::read_fasta_to_alignments_table("XXXbuild.fst",$dbh,$_,'accno'); # read alignment to database
                unlink "XXXbuild.fst";
            }
            else { print "Failed to merge alignments for $parent ($gene).\n"; }
            $taxon_tree = $left . $parent . $right;
        }
        else { last; }
    }
    return $ntaxa;
}
sub insert_taxon_in_taxon_tree {
    my $taxon_tree = shift @_;
    my $taxon;
    my $left = '';
    my $right = '';
    if (!$taxon_tree) {
        $taxon = shift @_;
        if (scalar @_ > 0) {
            $left = "((";
            $right = ")[$taxon])";
        }
        else { $left = "($taxon)"; }
    }
    else { $left = $taxon_tree; }
    while (@_) {
        $taxon = shift @_;
        if ($left =~ /(.*)(\)\[$taxon\].*)/) {
            $left = $1;
            $right = $2 . $right;
        }
        elsif ($left =~ /(.*)(\b$taxon\b)(.*)/) {
            if (scalar @_ > 0) { $left = "$1()[$taxon]"; }
        }
        else {
            if (scalar @_ > 0) {
                if ($left =~ /\($/) {
                    $left .= "(";
                    $right = ")[$taxon]$right";
                }
                else {
                    $left .= ",(";
                    $right = ")[$taxon]$right";
                }
            }
            else {
                if ($left =~ /\($/) {
                    $left .=  $taxon;
                }
                else {
                    $left .= ",$taxon";
                }
            }
        }
    }
    return $left . $right;
}

sub align {
    my $dbh = shift @_; # database handler
    my $path = shift @_; # variable to store the path to mafft
    my $gene = shift @_; # gene to align (column in table)
    my $taxon = shift @_; # what taxon to select (from taxon string in gb_data)
    my $n_threads = shift @_; # how many threads to use if using multi threads (otherwise 0)
    my $n_taxa=0; # to store number of taxa
    my $seq_length; # to store the sequence length
    print "Creating fasta file for $taxon, $gene.\n";
    # get the sequences
    my $sth = $dbh->prepare("SELECT $gene.accno,$gene.sequence FROM $gene INNER JOIN alignments ON $gene.accno=alignments.$gene\_accno INNER JOIN gb_data ON $gene.accno=gb_data.accno WHERE alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')")
        or die "Could not prepare statement: " . $dbh->errstr;
    $n_taxa = PifCosm_support_subs::print_fasta("XXXtemp_fastafile.fst",$sth); # make fasta file and get number of sequences
    if ($n_taxa < 2) { # if less than two sequences, no point in aligning it
        print "Only one sequence of $gene for $taxon. Copying it to alignments.\n";
        $sth = $dbh->prepare("SELECT $gene.accno,$gene.sequence FROM $gene INNER JOIN alignments ON $gene.accno=alignments.$gene\_accno INNER JOIN gb_data ON $gene.accno=gb_data.accno WHERE alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')")
        or die "Could not prepare statement: " . $dbh->errstr; # get accno and sequence
        $sth->execute();
        my @row = $sth->fetchrow_array(); # store accno and sequence
        $sth->finish();
        if ($row[1] and $row[0]) {
            my $update = $dbh->do("UPDATE alignments SET $gene\_sequence = '$row[1]' WHERE $gene\_accno='$row[0]'") or die "Could not update sequence in alignments: " . $dbh->errstr; # include sequence in alignments table
        }
        unlink "XXXtemp_fastafile.fst";
        return; # return nothing
    }
    if ($n_threads > 0) { # if pthreads
        my $threads = $n_threads;
        if ($n_threads > 1) { --$threads; } # mafft need one thread to coordinate the others
        print "Aligning $n_taxa sequences using mafft, pthread version.\n";
        system "${path}mafft --auto --quiet --thread $threads XXXtemp_fastafile.fst > XXXtemp_aligned.fst"; # align
    }
    else { # if not multi-threads
        print "Aligning $n_taxa sequences using mafft.\n";
        system "${path}mafft --auto --quiet XXXtemp_fastafile.fst > XXXtemp_aligned.fst"; # align
    }
    unlink "XXXtemp_fastafile.fst"; # remove file with unaligned sequences
    print "Reading aligned sequences to database.\n";
    ($n_taxa,$seq_length) = PifCosm_support_subs::read_fasta_to_alignments_table("XXXtemp_aligned.fst",$dbh,$gene,'accno'); # read the aligned sequences to the database
    unlink "XXXtemp_aligned.fst";
    return ($n_taxa); # return the number of sequences
}

# create an alignment and tree using a guide tree given a database handler, gene, taxon, linking tree, and linking gene
sub align_using_guide_tree {
    my $dbh = shift @_; # database handler
    my $path = shift @_; # path to programms
    my $gene = shift @_; # gene to align (column)
    my $taxon = shift @_; # taxon to align (gb_data taxon string)
    my $tree_method = shift @_; # the tree to use as guide tree
    my $guide_gene = shift @_; # the gene to use for guide tree
    my $n_threads = shift @_; # number of threads to use (0 if just serial)
    my $no_guide_tree = 'n'; # flag if no guide tree was possible
    my $guide_tree = 'empty';
    $tree_method =~ s/_bootstrap_[0-9]*//;
    ### Get sequences to make guide tree
    my $sth = $dbh->prepare("SELECT COUNT(alignments.$guide_gene\_accno),LENGTH(alignments.$guide_gene\_sequence) FROM alignments INNER JOIN gb_data ON alignments.$guide_gene\_accno=gb_data.accno WHERE alignments.$guide_gene\_sequence!='empty' AND alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')")
        or die "Could not prepare statement: " . $dbh->errstr;
    my ($n_taxa,$seq_length) = $sth->fetchrow_array();
    if ($n_taxa > 4) {
        print "Making guide tree.\n";
        $sth = $dbh->prepare("SELECT alignments.$guide_gene\_accno,alignments.$guide_gene\_sequence FROM alignments INNER JOIN gb_data ON alignments.$guide_gene\_accno=gb_data.accno WHERE alignments.$guide_gene\_sequence!='empty' AND alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')")
        or die "Could not prepare statement: " . $dbh->errstr;
        my $miss = PifCosm_support_subs::print_phylip("XXXtemp.phy",$sth,$n_taxa,$seq_length); # write phylip file
        {    my $MLscore;
            ($guide_tree,$MLscore) = PifCosm_support_subs::run_raxml("XXXtemp.phy",$path,"$guide_gene$taxon",$n_threads,$tree_method,'n',); # create tree
        }
        unlink ("XXXtemp.phy","XXXtemp.phy.reduced");
        if ($guide_tree eq 'empty' or $guide_tree =~ /RAxML_parsimonyTree/) { # if parsimony tree or no tree
            $no_guide_tree = 'y'; # flag as no suitable guide tree
            if ($guide_tree ne 'empty') { unlink $guide_tree; } # if parsimony tree remove parsimony tree file
        } # if parsimony tree or no tree, flag as no suitable guide tree
        
    }
    else { $no_guide_tree = 'y'; } # flag that it should not be aligned using guide tree
    print "Getting sequences from $gene that corrspond to sequences in $guide_gene.\n";
    # get the sequences that overlap between the sequences to be aligned and the sequences that will be used for the guide tree
    $sth = $dbh->prepare("SELECT alignments.$guide_gene\_accno,$gene.accno,$gene.sequence FROM $gene INNER JOIN alignments ON $gene.accno=alignments.$gene\_accno INNER JOIN gb_data ON $gene.accno=gb_data.accno WHERE alignments.$guide_gene\_sequences!='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')")
        or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    my @keep_tips; # the sequences that will be keept in the guide tree
    print "Writing fasta file.\n";
    open FASTAFILE, ">XXXtemp_guided.fst" or die "Could not open XXXtemp_guided.fst: $!.\n";
    my $n_guided=0; # number of sequences
    while (my @row = $sth->fetchrow_array()) {
        push (@keep_tips,$row[0]); # get the guide accno
        print FASTAFILE ">$row[1]\n$row[2]\n"; # print the fasta file with sequences to align
        ++$n_guided; # count the number of sequences
    }
    close FASTAFILE or die;
    $sth->finish();
    if ($n_guided != $n_taxa) { $no_guide_tree = 'y'; } # flag that it should not be aligned using guide tree
    if ($n_guided < 4) { $no_guide_tree = 'y'; } # flag that it should not be aligned using guide tree
    if ($no_guide_tree eq 'y') { # if no suitable guide tree
        print "No suitable guide tree.\n";
        return 0; # no alignment with guide tree, return 0
    }
    print "Preparing guide tree with " . scalar @keep_tips . " tips.\n";
    # replace the taxon labels in the tree with numbers for mafft
    my $bender_string;
    for (my $j=0; $j < scalar @keep_tips; ++$j) { # for each tip
        my $order = $j+1; # numbers from 1+
        if ($bender_string) { $bender_string .= ",$keep_tips[$j]|$order"; } # make string for treebender for swithing tip names
        else { $bender_string = "'$keep_tips[$j]|$order"; } # if first accno
    }
    $bender_string .= "'"; # end benderstring with '
    open TREEBENDER, "${path}treebender -c $bender_string < $guide_tree |"; # pipe tree into perl
    $guide_tree = <TREEBENDER>; # get tree
    close TREEBENDER or die "Could not close treebender properly: $!";
    my ($scale,$comparisons) = PifCosm_support_subs::comp_distance($dbh,$path,$guide_gene, $gene, $taxon, 10); # find value to rescale branch lengths
    if (!$scale) { # if no scale was found
        print "Could not compare sequence distances between $guide_gene and $gene. Will not rescale branches in guide tree.\n";
        $scale=1; # do not rescale
    }
    else {
        print "Compared sequence pairs of $guide_gene and $gene for $comparisons taxa. Will multiply the branches in the guide tree with $scale.\n";
    }
    $guide_tree = PifCosm_support_subs::newick_to_mafft($guide_tree,$scale); # reformat guide tree to mafft format
    open TREEFILE, ">XXXtemp_mafft.tree" or die "Could not open XXXtemp_mafft.tree: $!.\n"; # open file handler for mafft guide tree
    print TREEFILE $guide_tree; # print mafft guide tree
    close TREEFILE or die;
    print "Aligning sequences present in guide tree.\n";
    system "${path}mafft --auto --quiet --treein XXXtemp_mafft.tree XXXtemp_guided.fst > XXXtemp.aligned.fst"; # align sequences using guide tree
    print "Getting sequences not in guide tree.\n";
    # get sequences that did not overlap with the guide gene
    $sth = $dbh->prepare("SELECT $gene.accno,$gene.sequence FROM $gene INNER JOIN alignments ON $gene.accno=alignments.$gene\_accno INNER JOIN gb_data ON $gene.accno=gb_data.accno WHERE alignments.$guide_gene\_sequence='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')")
        or die "Could not prepare statement: " . $dbh->errstr;
    my $n_added = PifCosm_support_subs::print_fasta("XXXtemp.fst",$sth); # print them to fasta file
    print "Aligning sequences to guided alignment.\n";
    system "${path}mafft --auto --quiet --add XXXtemp.fst XXXtemp.aligned.fst > XXXtemp_all.fst"; # add them to alignment using mafft
    print "Reading aligned sequences to database.\n";
    my ($ntaxa,$seqlength) = PifCosm_support_subs::read_fasta_to_alignments_table("XXXtemp_all.fst",$dbh,$gene,'accno'); # read the aligned sequences to the database
    if ($ntaxa != $n_guided+$n_added) { # if it doesn't add up
        print "WARNING!!! " . (($n_guided+$n_added)-$ntaxa) . " fewer sequences were read into the database than written to be aligned.\n";
    }
    return ($ntaxa); # return the number of taxa aligned
}

# returns linked tree and gene given database handler, gene, taxon, and linked gene array
sub get_guide_tree_gene {
    my $dbh = shift @_; # database handler
    my $gene = shift @_; # gene to get guide tree for
    my $taxon = shift @_; # taxon that should be aligned
    #rest of @_ are linked genes
    my @possibilities; # to store possible guide genes
    #my $guide_tree; # to store guide tree
    my $guide_gene; # to store the guide gene
    # get the full taxon string
    my $query_taxon = $taxon;
    $query_taxon =~ s/'/''/g;
    my $sth = $dbh->prepare("SELECT DISTINCT taxon_string FROM gb_data WHERE (taxon_string LIKE '%; $query_taxon;%' OR taxon_string LIKE '$query_taxon;%' OR taxon_string LIKE '%; $query_taxon' OR taxon_string='$query_taxon')")  or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    my $temp = $sth->fetchrow_array(); # save full taxon string
    $sth->finish();
    my @hirarchy;
    if ($temp){ @hirarchy = split /; /, $temp; } # get the hierarchy
    unshift (@hirarchy,'combined'); # add combined to the hierarchy
    my $taxa_no; # the position in the hierarchy the taxon is
    for (my $j=0; $j < scalar @hirarchy; ++$j) { # for each taxon in the hierarchy
        if ($hirarchy[$j] eq $taxon) { $taxa_no = $j; } # save the position if it match the taxon we are interested in
    }
    for (my $j=0; $j < scalar @_; ++$j) { # for each gene that is linked
        if ($_[$j] =~ $gene) { # if our gene of interest is among the linked genes
            @possibilities = split ';', $_[$j]; # save the linked genes as possible guide genes
        }
    }
    if (@possibilities) { # if there are linked genes
        print "Found linked regions: ";
        for (my $j=0; $j < scalar @possibilities; ++$j) { # print each linked gene
            if ($possibilities[$j] eq $gene ) { next; }
            print $possibilities[$j];
            if ($j < scalar @possibilities -1) { print ', '; }
        }
        print ".\n";
        for (my $j=0; $j < scalar @possibilities; ++$j) { # for each possible guide gene
            my $level=-1;
            if ($possibilities[$j] eq $gene ) { next; } # can not guide itself, try next
            else {
                # get the possible trees for the gene
                my $sth = $dbh->prepare("SELECT COUNT(alignments.$possibilities[$j]\_accno),LENGTH(alignments.$possibilities[$j]\_sequence) FROM alignments INNER JOIN gb_data ON alignments.$possibilities[$j]\_accno=gb_data.accno WHERE alignments.$possibilities[$j]\_sequence!='empty' AND alignments.$gene\_accno!='empty' AND (gb_data.taxon_string LIKE '%; $query_taxon;%' OR gb_data.taxon_string LIKE '$query_taxon;%' OR gb_data.taxon_string LIKE '%; $query_taxon' OR gb_data.taxon_string='$query_taxon')")
                    or die "Could not prepare statement: " . $dbh->errstr;
                $sth->execute();
                my @row = $sth->fetchrow_array();
                $sth->finish();
                if ($row[0] > 4) { $guide_gene = $possibilities[$j]; }
            }
        }
    }
    return ($guide_gene); # return guide tree and guide gene
}

sub create_and_save_tree { # sub to create and save tree
    my $phylip_file = shift @_; # phylip file
    my $path = shift @_; # path for treebuilding program
    my $dbh = shift @_; # databse handler
    my $gene = shift @_; # gene that is treed
    my $taxon = shift @_; # taxon that is treed
    my $n_threads = shift @_; # number of threads if pthreads, else = 0
    my $tree_method = shift @_; # what method to use to tree
    my $store_boot = shift @_; # if 'y' store bootstrap trees
    my $tree_file; # to store tree file name
    my $boottree_file;
    my $n_boot;
    if ($tree_method =~ s/_bootstrap_([0-9]+)//) { $n_boot = $1; }
    if ($tree_method =~ /(Light)|(RAx)|(fasttree)/) { # if OK method
        my $MLscore;
        ($tree_file,$MLscore,$boottree_file) = PifCosm_support_subs::run_raxml($phylip_file,$path,"$gene$taxon",$n_threads,$tree_method,'n',$n_boot); # create tree
    }
    else {
        print "No recognized method for phylogenetic analysis.\n";
        $tree_file ='empty';
    }
    if ($tree_file eq 'empty') { print "Failed to create tree for $taxon, $gene.\n"; } # if no tree
    elsif ($tree_file =~ /RAxML_parsimonyTree/) { # if parsimony tree
        print "Only parsimony tree for $taxon, $gene ($tree_file), no ML tree.\n";
        PifCosm_support_subs::read_tree_to_alignment_groups($tree_file,$dbh,$gene,$taxon,"Parsimony"); # save tree in database
	unlink $tree_file;
    }
    else { # if not parsimony it is ML
        print "Tree created for $taxon, $gene ($tree_file).\n";
        PifCosm_support_subs::read_tree_to_alignment_groups($tree_file,$dbh,$gene,$taxon,"ML"); # save tree in database
	unlink $tree_file;
    }
    if ($store_boot eq 'y') {
	if ($boottree_file ne 'empty' && -e $boottree_file) {
	    print "Bootstraped trees created for $taxon, $gene ($boottree_file).\n";
	    PifCosm_support_subs::read_tree_to_alignment_groups($boottree_file,$dbh,"$gene\_boot","$taxon\_boot","ML_boot"); # save trees in database
	    unlink $boottree_file;
	}
	else { print "Failed to create bootstrap trees for $taxon, $gene.\n"; }
    }
    if ( -e $boottree_file) { unlink $boottree_file; } # get rid of boot tree file
    if ( -e $tree_file ) { unlink $tree_file; } # get rid of tree file
}

#sub read_tree_to_alignment_groups { # saves tree to database
#    my $file = shift @_; # tree file
#    my $dbh = shift @_; # database handler
#    my $gene = shift @_; # gene to save it under
#    my $taxon = shift @_; # taxon to save it under
#    my $method = shift @_; # method to annotate it with
#    open TREEFILE, "<$file" or die "Could not open tree file: $!.\n"; # open and
#    my $tree;
#    while (my $row = <TREEFILE>) { $tree .= $row; } # read tree
#    close TREEFILE or die;
    # uppdate database
#    my $string = $dbh->quote($tree);
#    my $update = $dbh->do("UPDATE alignment_groups SET tree=$string,tree_method='$method' WHERE gene='$gene' AND taxon='$taxon'") or die "Could not update database: " . $dbh->errstr;
#    return $update; # return how many rows were uppdated
#}

sub remove_outliers { # sub that remove sequences that are decendents of "extremly" long branches (outliers)
    my $dbh = shift @_; # get database handler 
    my $path = shift @_; # path to clustertree and treebender
    my $gene = shift @_; # get the gene the tree is for
    my $taxa = shift @_; # get taxon the tree is for
    # get tree from database
    my $string_taxa = $dbh->quote($taxa);
    my $string_gene = $dbh->quote($gene);
    my $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE taxon=$string_taxa AND gene=$string_gene") or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    my $tree = $sth->fetchrow_array(); # save tree
    $sth->finish();
    print "Identifying outlier sequences for $gene ($taxa).\n";
    if ($tree and $tree ne 'empty' and $tree =~ /;/) { # if we got a tree
        open TREEFILE, ">XXXtemp_tree_file.tree" or die "Could not open XXXtemp_tree_file.tree: $!.\n";
        print TREEFILE $tree; #print tree to file
        close TREEFILE or die;
        my $cut_off = &get_outlier_length("XXXtemp_tree_file.tree", $path); # get the branchlength cut off from outliers according to gamma distribution
        my @remove_seq = `${path}treebender --cluster long_branch --cut_off $cut_off < XXXtemp_tree_file.tree`; # get groups separated by branches longer than cut-off
        my %largest_cluster; # hash to store number of taxa per cluster
	if ($remove_seq[0] =~ /^### tree/) { shift @remove_seq; }
        for (my $i=0; $i< scalar @remove_seq; ++$i) { # for each cluster
            chomp($remove_seq[$i]); # get rid of line breaks
            $remove_seq[$i] =~ s/\s+$//; # get rid of trailing white spaces
            $remove_seq[$i] =~ s/^\s+//; # get rid of leading white spaces
            my @temp = split /\s+/,$remove_seq[$i]; # separate different tips
            $largest_cluster{$i} = scalar @temp; # count sequences in cluster
            $remove_seq[$i] = join ",", @temp; # join with , between
        }
        my @min_first = PifCosm_support_subs::keys_in_order(%largest_cluster); # get the clusters ordered by how many sequences they have
        #foreach (@min_first) { print "$_ - $remove_seq[$_]\n"; }
        my $drop_tips = '';
        for (my $i=1; $i < scalar @min_first; ++$i) { # add all but largest cluster to sequences to remove
            if ($i==1) { $drop_tips = $remove_seq[$min_first[$i]]; } # make string with tips to remove
            else { $drop_tips .= ",$remove_seq[$min_first[$i]]"; } # as above
        }
        @remove_seq = split /,/, $drop_tips; # get all sequences separately
        print "Found " . (scalar @remove_seq) . " sequences nested in clades delimited by a branch longer than $cut_off. These sequences will be removed from the alignment.\n";
        if (scalar @remove_seq > 0 ) { # if sequences that should be removed
            print "Removing:";
            my $changes_accnos = 0; # save the number of accnos that has been removed
            foreach (@remove_seq) {
                $changes_accnos += PifCosm_support_subs::remove_accno_from_alignments($dbh,$gene,$_);
                print " $_...";
            } # remove each accno from the alignments table
            $tree = `${path}treebender -d '$drop_tips' < XXXtemp_tree_file.tree`; # remove sequences from tree
            unlink "XXXtemp_tree_file.tree"; # delete tree file
            my $string = $dbh->quote($tree);
            my $changes = $dbh->do("UPDATE alignment_groups SET tree=$string WHERE taxon='$taxa' AND gene='$gene'"); # save pruned tree in database
            print "\nDone cleaning.\n";
            return $changes_accnos; # return the number of accnos that were deleted
        }
        else { unlink "XXXtemp_tree_file.tree"; return 0; } # if no sequences to remove, no need to do anything but delete the tree file
    }
    else { print "No tree recovered for $gene of $taxa. No sequences pruned.\n"; return 0; } # if no tree was recovered print warning and return 0
}

sub get_outlier_length { # get branchlength that is unlikely
    my $treefile = shift @_; # get tree file
    my $path  = shift @_; # path to treebender
    chomp(my @data = `${path}treebender -a \\\\n < $treefile`); # get branchlengths from tree
    @data=sort by_number (@data); # sort the branchlengths
    sub by_number { $a <=> $b }
    while ($data[0] =~ /[^0-9\.eE-]/ or $data[0] eq '') {shift @data;} # if there are non numeric values, remove them
    my $avarage=0; # to get the mean
    foreach (@data) { $avarage+=$_;} # first sum up branch lengths
    $avarage /= scalar @data; # then divide by number of values
    my ($shape,$scale) = &optimize_gamma(0,3,@data); # get shape and scale ML estimates of a gamma distribution given the observed branch lengths
    my $cut_off; # branch length that is smaller than the highest non likely branch length but higher than the highest likely branch length
    if ($shape > 0 && $scale > 0) {
        my $i = scalar @data -1; # the number of the longest branch
        while (&pgamma($shape,$scale,$data[$i])>1-(0.05/($i+1))) { --$i;} # while the P value for the branch is higher that the Holm–Bonferroni adjusted alpha of 0.05 (P cut of)
        if ($i == scalar @data -1) { # if the longest branch was likely
            $cut_off = ($data[$i]+1); # the cut off is longer than this branch
        }
        else { # if longest branch not unlikely
            $cut_off = (($data[$i]+$data[$i+1])/2); # cut off is avarage of the highest non likely branch length and the highest likely branch length
        }
    }
    else {
        print STDERR "ERROR IN GAMMA ESTIMATION!!! NO SEQUENCES EXCLUDED!!!\n";
        $cut_off = ($data[scalar @data-1]+1);
    }
    return $cut_off;
}

sub optimize_gamma { # optimization of parameters of gamma distribution using a golden section search (see eg. Yang Z, 2006, Computational Molecular Evolution, OXFORD UNIVERSITY PRESS
                     # since the ML value of the scale parameter can be calculated given the shape parameter, optimization can be done using a univariate method
    my @shape;
    $shape[0] = shift @_; # set min boundary of for shape
    $shape[3] = shift @_; # set max boundary of for shape
    if ($shape[0] < 0.00001) { $shape[0] = 0.00001 }; # can not handle zero values
    if (!$shape[3]) { $shape[3] = 10; } # if we have no max value set it to 10
    $shape[1] = $shape[3]-($shape[3]-$shape[0])*0.6180; # set intermediet values based on the golden cut
    $shape[2] = $shape[0]+($shape[3]-$shape[0])*0.6180; # set intermediet values based on the golden cut

    my $avarage=0; # get avarage by
    foreach (@_) { $avarage+=$_;} # summing the values and
    $avarage /= scalar @_; # dividing them by the number of values

    my @log_lik; # the log likelihood of the boundary values
    $log_lik[0] = &dgamma($shape[0],$avarage/$shape[0],@_); # get start values
    $log_lik[1] = &dgamma($shape[1],$avarage/$shape[1],@_); # get start values
    $log_lik[2] = &dgamma($shape[2],$avarage/$shape[2],@_); # get start values
    $log_lik[3] = &dgamma($shape[3],$avarage/$shape[3],@_); # get start values
    while (1) { # until we reach stop criterium
        if ($log_lik[1]<$log_lik[2]) { # if lower intermediate value is worse than the higher
            $log_lik[0] = $log_lik[1]; # we can move the lower boundary up
            $shape[0] = $shape[1]; # both likelihood and value
        }
        else { # if higher value worse
            $log_lik[3]=$log_lik[2]; # we can move higher bound down
            $shape[3] = $shape[2]; # both likelihood and value
        }
        $shape[1] = $shape[3]-($shape[3]-$shape[0])*0.6180; # get new lower intermediate value
        $shape[2] = $shape[0]+($shape[3]-$shape[0])*0.6180; # get new higher intermediate value
        $log_lik[1] = &dgamma($shape[1],$avarage/$shape[1],@_); # get log like for lower
        $log_lik[2] = &dgamma($shape[2],$avarage/$shape[2],@_); # get log like for higher
        if (($shape[0]/$shape[3])>0.9999) { return ($shape[1],$avarage/$shape[1]); } # if the boundaries are really close return the lower intermediate shape, and calculate scale parameter and return
    }
}

sub dgamma { # calculate the log of the probability dencity
    my $shape = shift @_; # get the shape
    my $scale = shift @_; # and scale
    my $log_lik = 0; # null log like
    foreach(@_) { # for each value
        if (!$_ or $_ =~ /[^0-9\.eE-]/ or $_ eq '' or $_ < 0) { next;} # if it is not numerical or negative skip it
        $log_lik += (log($_)*($shape-1))-($shape*log($scale))-log(&gamma_function($shape))-($_/$scale); # add log likelihood
    }
    return $log_lik; # return total log like
}

sub pgamma { # calculate the culumative likelihood
    my $shape = shift @_; # get shape
    my $scale = shift @_; # and scale
    my $x = shift @_; # get value to check
    return &low_incomp_gamma($shape,$x/$scale)/&gamma_function($shape); # calculate P and return
}

sub gamma_function { # Lanczos approximation of gamma function
    my $z = shift @_;
    $z -= 1;
    my @p = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7);
    my $g = 7;
    my $x = $p[0];
    foreach my $i (1 .. $g+1) {
        $x += $p[$i]/($z+$i);
    }
    my $t = $z+$g+0.5;
    return sqrt(2*3.14159265359) * $t ** ($z+0.5) * exp(-$t) * $x;
}

sub low_incomp_gamma { # approximate the integral of the gamma function from 0 to value
    my $s = shift @_; # parameter
    my $x = shift @_; # value to integrate to
    my $i_gamma = 0; # incomplete gamma value
    my $i = 0; # counter
    while (1) { # do until satisfactory precition
        my $temp = $i_gamma; # store temporarily
        $i_gamma += ($x**$i)/&gamma_function(($s+$i+1)); # sum up incomplete gamma value
        if ($temp/$i_gamma > 0.9999) { last; } # if sufficient precition, stop adding and move on
        ++$i; # add up counter
    }
    return ($x**$s)*&gamma_function($s)*exp(-$x)*$i_gamma; # finnish calculation and return
}

1;
