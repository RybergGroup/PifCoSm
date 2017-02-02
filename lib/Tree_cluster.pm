package Tree_cluster;
use strict;
use PifCosm_support_subs;

##################################################
### Subroutine to cluster monophyletic species ###
##################################################

sub tree_cluster {
    my $database = shift @_; # get the database
    my $path = shift @_; # path to RAxML and phylomand programs
    my $species_name = shift @_; # if sequences of monophyletic species should be clustered
    my $seq_sim = shift @_; # should monophyletic groups united by short branch lengths (single link) be clustered
    my %cut_off = split(/,/, shift @_); # cut off value for branch length clustering
    my $min_length = shift @_; # get the minimum length of sequences to cluster
    my $n_threads = shift @_; # get the number of threads (0 if pthread programs are not gone be used)
    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr;
    my @tables = keys %cut_off;
    if ($tables[0] =~ /^all$/i) {
        my $cut_off = $cut_off{$tables[0]};
        undef %cut_off;
        @tables = PifCosm_support_subs::get_gene_tables($dbh); # get gene tables
        foreach (@tables) { $cut_off{$_} = $cut_off; }
    }
    foreach (@tables) { $cut_off{$_} = 1-$cut_off{$_}; }
    my $sth = $dbh->prepare("SELECT DISTINCT taxon_string FROM gb_data ORDER BY taxon_string") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    my @taxon_strings; # to store distinct taxon anotations
    my $i = 0;
    while (my $ref = $sth->fetchrow_hashref()) {
        $taxon_strings[$i++] = $ref->{taxon_string};
    }
    $sth->finish();
    print "There are " . ($i-1) . " taxa (mostly at genus level).\n";
    # cluster for each gene separately
    for (my $j = 0; $j < scalar @tables; ++$j) {
        if (PifCosm_support_subs::table_present($dbh,$tables[$j]) eq 'n') { next; }
        my $n_uppdate = 0;
        my $n_taxa_without_sequence=0;
        print '-' x (length($tables[$j])+2) . "\n $tables[$j] \n" . '-' x (length($tables[$j])+2) . "\n";
        # get number of accnos
        $sth = $dbh->prepare("SELECT COUNT(accno) FROM $tables[$j] WHERE cluster != 'empty'");
        $sth->execute();
        my $count = $sth->fetchrow_array();
        # cluster within each taxon separately
        for ( my $k = 0; $k < scalar @taxon_strings; ++$k ) {
            my $taxon = $dbh->quote($taxon_strings[$k]);
            if ($count > 0) { # if there are previous clusters
                $sth = $dbh->prepare("SELECT $tables[$j].accno,$tables[$j].sequence FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno = gb_data.accno WHERE gb_data.taxon_string = $taxon and $tables[$j].cluster = 'lead' AND LENGTH($tables[$j].sequence) >= $min_length")
                    or die "Couldn't prepare statement: " . $dbh->errstr; # only get the 'best' sequence of each cluster of the gene for the taxon
            }
            else { # if no previous clusters
                $sth = $dbh->prepare("SELECT $tables[$j].accno,$tables[$j].sequence FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno = gb_data.accno WHERE gb_data.taxon_string = $taxon AND LENGTH($tables[$j].sequence) >= $min_length")
                    or die "Couldn't prepare statement: " . $dbh->errstr; # get all accnos and sequences of the gene for the taxon
            }
            $sth->execute();
            my $counter = 0; # to count the number of sequences
            open FASTAFILE, ">XXXunaligned.fst" or die "Couldn't open XXXunaligned.fst: $!.\n";
            while (my $seq_ref = $sth->fetchrow_hashref()) {
                print FASTAFILE ">" . $seq_ref->{accno} . "\n" . $seq_ref->{sequence} . "\n";
                ++$counter;
            }
            close FASTAFILE or die;
            $sth->finish();
            if ($counter > 0) { print "--- Clustering $counter sequences in $taxon_strings[$k] ---\n"; }
            if ($counter > 3) { # if more than three sequences align and create phylogeny
                # align
                if ($n_threads > 0) { `${path}$External_program::mafft --auto --quiet --thread $n_threads XXXunaligned.fst > XXXaligned.fst`; } # if using threads use PTHREAD version of mafft
                else { `${path}$External_program::mafft --auto --quiet XXXunaligned.fst > XXXaligned.fst`; } # otherwise not
                open ALIGNED, "<XXXaligned.fst" or die "Couldn't open XXXaligned.fst: $!.\n"; # file handle to read aligned sequences
                open PHYLIP, ">XXXaligned.phy" or die "Couldn't open XXXaligned.phy: $!.\n"; # file handle to print sequences in phylip format
                my $accno; # to store accno
                my $sequence; # to store sequence
                $i = 0; # counter
                # convert alignment to phylip format
                while (my $infile = <ALIGNED>) { # read aligned sequences row by row
                    chomp ($infile); # get rid of line breaks
                    if ($infile =~ s/^>//) { # if sequence name, remove leading >
                        if ($sequence) { # if sequence have been read
                            if ($i == 0) { # if first sequence
                                print PHYLIP $counter . " " . length($sequence) . "\n"; # number of taxa and number of bases
                                ++$i; # flag that it is not the first sequence any more
                            }
                            print PHYLIP $accno . ' ' . $sequence . "\n"; # output accno and sequence
                        }
                        $accno = $infile; # get the new accno
                        undef $sequence; # empty sequence
                    }
                    else { $sequence .= $infile; } # if it is not name it has to be sequence
                }
                print PHYLIP $accno . ' ' . $sequence . "\n"; # print last sequence
                close ALIGNED or die;
                close PHYLIP or die;
                # make phylogeny
                if (-e "RAxML_info.XXXtemp") { die "Can not run raxmlHPC since there already are files with the ending XXXtemp.\n"; } # check if files will block running RAxML
                my $raxml = `${path}$External_program::raxml -y -s XXXaligned.phy -m GTRGAMMA -n XXXtemp -p 12345`; # run parsimonator
                if ( -e "RAxML_parsimonyTree.XXXtemp" ) { # if tree file from parsimonator
                    if ( $seq_sim eq 'y') {
                        `${path}$External_program::raxml -f e -s XXXaligned.phy -m GTRGAMMA -n XXXoptim -t RAxML_parsimonyTree.XXXtemp`; # optimize branchlengths on parsimony tree
                        if (-e "RAxML_result.XXXoptim") {
                            print "Clustering based on similarity.\n";
                            system "${path}$External_program::treebender --cluster branch_length --cut_off $cut_off{$tables[$j]} < RAxML_result.XXXoptim > XXXclusters.txt"; # single link clustering based on branch lengths
                            my $n_clusters=&process_clusters($dbh,'XXXclusters.txt',$tables[$j]); # select 'best' sequence in cluster and make changes in database
                            print "$n_clusters clusters.\n";
                            # get sequences that are still lead in clusters
                            my $taxon = $dbh->quote($taxon_strings[$k]);
                            $sth = $dbh->prepare("SELECT $tables[$j].accno FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno = gb_data.accno WHERE gb_data.taxon_string = $taxon and $tables[$j].cluster = 'lead' AND LENGTH($tables[$j].sequence) >= $min_length")
                                or die "Couldn't prepare statement: " . $dbh->errstr; # only get the 'best' sequence of each cluster of the gene for the taxon
                            $sth->execute();
                            my $keep = ''; # string of sequences to keep, for tree bender
                            #my $n_to_keep=0;
                            while (my $temp = $sth->fetchrow_array() ) { $keep .= "$temp,"; } #++$n_to_keep;} # add each lead sequence
                            $keep =~ s/,$//; # remove trailing ,
                            system "${path}$External_program::treebender -d $keep -i < RAxML_result.XXXoptim > RAxML_parsimonyTree.XXXtemp"; # remove non-leade sequences from tree
                            unlink "XXXclusters.txt","RAxML_info.XXXoptim","RAxML_result.XXXoptim","RAxML_log.XXXoptim","RAxML_binaryModelParameters.XXXoptim","XXXaligned.phy.reduced"; # clean up
                        }
                        else { print STDERR "WARNING!!! No tree with branch lengths optimized from RAxML. Can not cluster based on similarity in tree. Continuing reluctantly...\n"; }
                    }
                    if ($species_name eq 'y') {# cluster sequences based on taxon annotation
                        print "Clustering based on species annotation.\n";
                        system "${path}$External_program::treebender --cluster database:$database,gb_data,accno,species < RAxML_parsimonyTree.XXXtemp > XXXclusters.txt"; # use it to cluster monophyletic species
                        my $n_clusters=&process_clusters($dbh,'XXXclusters.txt',$tables[$j]); # select 'best' sequence in cluster and make changes in database
                        print "$n_clusters clusters.\n";
                        unlink "XXXclusters.txt"; # clean up
                    }
                }
                else {
                    print STDERR "WARNING!!! RAxML failed, no tree, no clustering. Continuing...\n";
                }
                unlink "XXXaligned.fst", "XXXaligned.phy", "XXXaligned.phy.reduced", "RAxML_info.XXXtemp", "RAxML_parsimonyTree.XXXtemp"; # clean up
            }
            elsif ($counter > 0) { # if less than three but more than 0 sequences
                if ($count > 0) { # if previous clusters
                    my $taxon = $dbh->quote($taxon_strings[$k]);
                    $sth = $dbh->prepare("SELECT $tables[$j].accno,$tables[$j].sequence,gb_data.species FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno = gb_data.accno WHERE gb_data.taxon_string = $taxon AND $tables[$j].cluster = 'lead'  AND LENGTH($tables[$j].sequence) >= $min_length")
                        or die "Couldn't prepare statement: " . $dbh->errstr;
                }
                else { # if no previous cluster
                    my $taxon = $dbh->quote($taxon_strings[$k]);
                    $sth = $dbh->prepare("SELECT $tables[$j].accno,$tables[$j].sequence,gb_data.species FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno = gb_data.accno WHERE gb_data.taxon_string = $taxon AND LENGTH($tables[$j].sequence) >= $min_length")
                        or die "Couldn't prepare statement: " . $dbh->errstr;
                }
                $sth->execute();
                my @accnos; # array for accnos
                my @sequences; # array for sequences
                my @species; # array for corresponding species names
                my $i=0; # counter
                # get the values
                while (my @row = $sth->fetchrow_array()) {
                    $accnos[$i] = $row[0];
                    $sequences[$i] = $row[1];
                    $species[$i] = $row[2];
                    ++$i;
                }
                my @same_species; # array to store if it is the same species as another sequence
                if ($i != $counter) { die "Missmatch in number of sequences read for $taxon_strings[$k], $tables[$j]. Quitting.\n"; }
                if ($species_name eq 'y') { print "Cluster seq with same species annotation.\n"; }
                if ($seq_sim eq 'y') { print "Cluster based on pairwise sequence similarity.\n"; }
                my $n_cells;
                if ($i == 2) { $n_cells = 1; }
                elsif ($i == 3) { $n_cells = 3; }
                else { $n_cells = 0; }
                for ($i=0; $i < $n_cells; ++$i) { $same_species[$i] = 0; } # assume they are not the same species
                for ($i=0; $i < scalar @accnos; ++$i) {           # check each pair if they are the same
                    for (my $j=$i+1; $j < scalar @accnos; ++$j) { #
                        my $cluster_flag = 0;                     # flag if they should be clustered
                        if ($species_name eq 'y') {
                            if ($species[$i] eq $species[$j]) {       # if so
                                ++$cluster_flag;
                            }
                        }
                        if ( $seq_sim eq 'y') {
                            open SEQPAIR, ">XXXseq_pair.fst" or die "Could not open XXXseq_pair.fst: $!.\n";
                            print SEQPAIR ">$accnos[$i]\n$sequences[$i]\n>$accnos[$j]\n$sequences[$j]\n";
                            chomp (my $JCdistance = `${path}$External_program::pairalign -j < XXXseq_pair.fst`);
                            unlink "XXXseq_pair.fst";
                            $JCdistance =~ s/[^0-9\.-]//g;
                            if ($JCdistance < $cut_off{$tables[$j]}) { ++$cluster_flag; }
                        }
                        if ($cluster_flag > 0) {
                            if ($i==0 && $j==1) { $same_species[0] = 1; } # flag em
                            elsif ($i==0 && $j==2) { $same_species[1] = 1; }
                            elsif ($i==1 && $j==2) { $same_species[2] = 1; }
                        }
                    }
                }
                # make uppdates for each accno
#                for ($i=0; $i < scalar @accnos; ++$i) {
                open CLUSTERFILE, ">XXXclusters.txt" or die "Could not open temp file XXXclusters.txt: $!\n";
                if ($n_cells == 1) {
                    if ($same_species[0] == 1) { print CLUSTERFILE "$accnos[0] $accnos[1]\n"; }
                    else { print CLUSTERFILE "$accnos[0]\n$accnos[1]\n"; }
                }
                elsif ($n_cells == 3) {
                    if ($same_species[0] == 0 && $same_species[1] == 0 && $same_species[2] == 0) { print CLUSTERFILE "$accnos[0]\n$accnos[1]\n$accnos[2]\n"; }
                    elsif ($same_species[0] == 1 && $same_species[1] == 0 && $same_species[2] == 0) { print CLUSTERFILE "$accnos[0] $accnos[1]\n$accnos[2]\n"; }
                    elsif ($same_species[0] == 0 && $same_species[1] == 1 && $same_species[2] == 0) { print CLUSTERFILE "$accnos[0] $accnos[2]\n$accnos[1]\n"; }
                    elsif ($same_species[0] == 0 && $same_species[1] == 0 && $same_species[2] == 1) { print CLUSTERFILE "$accnos[0]\n$accnos[1] $accnos[2]\n"; }
                    else { print CLUSTERFILE "$accnos[0] $accnos[1] $accnos[2]\n"; }
                }
                else { foreach (@accnos) { print CLUSTERFILE "$_\n"; } }
                close CLUSTERFILE or die;
                my $n_clusters=&process_clusters($dbh,'XXXclusters.txt',$tables[$j]); # select 'best' sequence in cluster and make changes in database
                print "$n_clusters clusters.\n";
                unlink 'XXXclusters.txt';
 #               }
                $sth->finish();
            }
            else { ++$n_taxa_without_sequence; } #
            unlink "XXXunaligned.fst";
        }
        print "--------------------------\n$n_taxa_without_sequence taxa had no sequence associated with them for $tables[$j].\n";
    }
    $dbh->disconnect();
}
sub process_clusters {
    my $dbh = shift @_;
    my $cluster_file = shift @_; # file with each cluster on a separate row and accnos in cluster separated by space
    my $table = shift @_; # gene that is being clustered
    my $n_cluster = 0; # count number of clusters
    my $sth; # search string handler
    my $i = 0; #counter
    open CLUSTERS, "<$cluster_file" or die "Couldn't open $cluster_file: $!.\n"; # file handeler for clustree output
    while (my $infile = <CLUSTERS>) { # read clustertree file line by line
        chomp ($infile);
        $infile =~ s/\[.*\]//g; # remove anything between square brachets
        if ($infile =~ /[A-Za-z0-9]+/) { #if matching letters or numbers
            $infile =~ s/(^\s+)|(\s+$)//g; # substitute leading and trailing white spaces
            $infile =~ s/\s+/,/g; # substitute whitespace for comma
            PifCosm_support_subs::process_cluster_array($dbh,$infile,$table);
            ++$n_cluster;
        }
    }
    close CLUSTERS or die;
    return $n_cluster;
}
1;
