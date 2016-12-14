package CV_gene_linker;
use strict;
use PifCosm_support_subs;
use Matching;

sub CV_gene_linker {
    my $database = shift @_; # get database
    my $comp_columns= shift @_; # get the columns in gb_data to compare
    my @comp_type = split (/,/, shift @_); # get the type of scoring to use for each column
    my @comp_score = split (/,/, shift @_); # get the score for each column
    my $dbh = PifCosm_support_subs::connect_to_database($database); # connect to database
    my @tables = PifCosm_support_subs::get_gene_tables ($dbh); # get the names of the gene tables
    print "Will link the following genes: ";
    foreach(@tables) { print "$_ "; } # print what genes will be used
    print "\n";
    my $query_genes = ''; # create string for what columns to create
    foreach (@tables) {
        $query_genes .= ", $_\_accno TEXT DEFAULT 'empty', $_\_sequence TEXT DEFAULT 'empty'"; # add column for accno and sequence for each gene
    }
    $dbh->do("CREATE TABLE alignments (taxon_name TEXT$query_genes , PRIMARY KEY (taxon_name))") or die "Could not create table: " . $dbh->errstr; # create alignments table
    # start by defining a taxa for each cluster in first gene
    my $sth=$dbh->prepare("SELECT accno FROM $tables[0] WHERE cluster='lead'"); # get the lead sequence of each cluster
    $sth->execute();
    my $number=1; # set the number for unique names of taxa
    while (my $accno=$sth->fetchrow_array()) {
        $dbh->do("INSERT INTO alignments (taxon_name,$tables[0]_accno) VALUES ('Taxon_$number','$accno')"); # insert taxa for each cluster in database
        ++$number; # increase the number for taxon name
    }
    $sth->finish();
    # go through the rest of the genes one by one and link or add taxa
    for (my $i=1; $i< scalar @tables; ++$i) { 
        print "Linking sequences of $tables[$i] to taxa.\n";
        my $columns=''; # string for the columns with accnos of previously added genes
        my %cluster_names; # hash to store the accno of the lead of each cluster 
        my @graph; # array to store the edges of the graph for linking
        for (my $j=0; $j<$i; ++$j) {
            $columns.= ",$tables[$j]_accno"; # add accno column for each gene preceding the one to be linked in this itteration
        }
        $sth=$dbh->prepare("SELECT taxon_name$columns FROM alignments"); # get the taxon name and accnos for each gene
        $sth->execute();
        while (my $row = $sth->fetchrow_hashref()) { # for each taxon
            my $taxon=$row->{'taxon_name'}; # save the taxon name
            my %taxon_accno; # hash to store the info to compare for linkeage of clusters
            for (my $j=0; $j<$i; ++$j) { # for each previous gene
                 my $accno=$row->{"$tables[$j]_accno"}; # get the accno of the lead sequence of that cluster
                 if ($accno eq 'empty') { next; } # if no accno for the gene go to next gene
                 my $sth2 = $dbh->prepare("SELECT $comp_columns FROM gb_data INNER JOIN $tables[$j] on gb_data.accno=$tables[$j].accno WHERE $tables[$j].accno='$accno' OR $tables[$j].cluster='$accno'"); # get the info for each sequence of the taxon
                 $sth2->execute();
                 while (my @row=$sth2->fetchrow_array()) {
                     if (!$taxon_accno{$accno}) { # if accno not checked before
                         $taxon_accno{$accno}=\@row; # save the info for the accno
                     }
                 }
                 $sth2->finish;
            }
            my $sth2=$dbh->prepare("SELECT gb_data.accno,$tables[$i].cluster,$comp_columns FROM gb_data INNER JOIN $tables[$i] on gb_data.accno=$tables[$i].accno WHERE $tables[$i].cluster!='empty'"); # get the info to compare for each sequence in the gene to link in this iteration
            $sth2->execute();
            my %edges; # hash to store edges between the taxa and the clusters to add
            while (my @row = $sth2->fetchrow_array()) {
                my $lead; # store lead sequence for cluster
                if ($row[1] eq 'lead') { $lead=$row[0]; } # if it is the lead sequence store the accno
                else { $lead=$row[1]; } # if not the lead sequence store the accno for the lead sequence
                ++$cluster_names{$lead}; # flag that we have read a sequence for the cluster
                foreach my $accno (keys %taxon_accno) { # For each accno belonging to the taxa
                    if (!$edges{$lead}) { $edges{$lead}=[]; } # initiate the edge
                    for (my $j=0; $j< scalar @{$taxon_accno{$accno}}; ++$j) { # for each column to compare
                        if ($taxon_accno{$accno}->[$j] eq $row[$j+2]) { # if it is the same for the taxon and cluster
                            if ($comp_type[$j] eq 'o' ) { $edges{$lead}->[$j]= $comp_score[$j]; } # if it should be scored based on occurence save score
                            elsif ($comp_type[$j] eq 'a' ) { $edges{$lead}->[$j]+= $comp_score[$j]; } # if it should be scored additatively based on occurences add score
                            elsif ($comp_type[$j] eq 'e' && !$edges{$lead}->[$j]) {$edges{$lead}->[$j] = $comp_score[$j]; } # if it should be based on only occurence of the same annotation add score
                            elsif ($comp_type[$j] eq 'h' && !$edges{$lead}->[$j]) {$edges{$lead}->[$j] = $comp_score[$j]; } # if exclusive occurence should get full score
                        }
                        else {
                            if ($comp_type[$j] eq 'h' && $edges{$lead}->[$j] && $edges{$lead}->[$j]==$comp_score[$j]) { $edges{$lead}->[$j] /= 2; } # if unexclusive occurence should get half score
                            elsif ($comp_type[$j] eq 'e' && $edges{$lead}->[$j] && $edges{$lead}->[$j]==$comp_score[$j]) { $edges{$lead}->[$j] = 0; } # if unexclusive occurence should get zero score
                        }
                    }
                    
                }
            }
            $sth2->finish;
            foreach my $accno (keys %edges) { # add each edges to the graph
                my $sum=0;
                foreach (@{$edges{$accno}}) { if ($_) { $sum += $_; } } # add the scores for each column
                if ($sum>0) { push (@graph,[$taxon,$accno,$sum]); } # add the edge to graph
            }
        }
        my %links =  Graph::Matching::max_weight_matching(\@graph); # get the matching taxa and clusters
        foreach (keys %links) { # go through links and only save links where a Taxon name is in the key
            if (!($_ =~ /^Taxon_[0-9]+$/)) { # key is not a taxon name
                if (!$links{$links{$_}}) { # if there is no corresponding link with the taxon name as key
                    $links{$links{$_}}=$_; # save the link with the taxon name as key
                }
                elsif ($_ ne $links{$links{$_}}) { # if taxon that is linked to is also linked to oter cluster print a warning
                    print STDERR "WARNING!!! Uncompatible linkege made, will be resolved arbitrarily.\n";
                }
                delete $links{$_}; # delete the link where the key is not a taxon name
            }
        }
        foreach (keys %links) { # for link
            if ($_ =~ /^Taxon_[0-9]+$/) { # should always be the case but check
                $dbh->do("UPDATE alignments SET $tables[$i]_accno='$links{$_}' WHERE taxon_name='$_'"); # add the linked cluster to the taxa
                $cluster_names{$links{$_}}=0; # flag that the cluster have been added
            }
        }
        foreach my $lead (keys %cluster_names) { # check if all cluster have been added
            if ($cluster_names{$lead}>0) { # if not added
               $dbh->do("INSERT INTO alignments (taxon_name,$tables[$i]_accno) VALUES ('Taxon_$number','$lead')"); # add cluster to new taxon
               ++$number; # increase the number for taxon name
            }
        }
        print "There are now " . ($number-1) . " taxa.\n";
    }   
}
1; 
