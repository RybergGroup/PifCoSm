package Gene_linker;
use strict;
use PifCosm_support_subs;

################################################################
### Subroutine to create table for alignments and link genes ###
################################################################
sub gene_linker {
    my $database = shift @_; # get database
    my $min_length = shift @_; # get minimi sequence length to include
    my $max_cluster = shift @_; # y if it is enough that things are clustered in one gene to be considered the same taxa
                                # n if it enogh that things are in separate clusters for one gene to be considered separate taxa
    my $dbh = PifCosm_support_subs::connect_to_database($database);
    my @tables = PifCosm_support_subs::get_gene_tables ($dbh);
    print "Will link the following genes: ";
    foreach(@tables) { print "$_ "; }
    print "\n";
    my $query_genes = ''; # create string for what columns to create
    foreach (@tables) {
        $query_genes .= ", $_\_accno TEXT DEFAULT 'empty', $_\_sequence TEXT DEFAULT 'empty'";
    }
    # create table for alignments
    $dbh->do("CREATE TABLE alignments (taxon_name TEXT$query_genes , PRIMARY KEY (taxon_name))") or die "Could not create table: " . $dbh->errstr;

    my %numbers_gene; # to store the number of sequences for each gene
    # get the number of sequences to include for each gene
    for (my $i=0; $i < scalar @tables; ++$i) {
        my $sth = $dbh->prepare("SELECT COUNT(accno) FROM $tables[$i] WHERE cluster='lead' AND LENGTH(sequence)>$min_length") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        $numbers_gene{$tables[$i]} = $sth->fetchrow_array();
        $sth->finish();
    }
    # sort the tables based on the number of associated sequences
    @tables = PifCosm_support_subs::keys_in_order(%numbers_gene);#sort in_order keys %numbers_gene;

    my $taxon_counter=1; # counter to give unique name to each taxon
    # get sequences for each gene and sort them into taxa
    for (my $i=0; $i < scalar @tables; ++$i) {
        print "Looking for new taxa among the $tables[$i] sequences.\n";
        # get the accnos already in the alignment table
        my $sth = $dbh->prepare("SELECT $tables[$i]_accno FROM alignments WHERE $tables[$i]_accno!='empty'") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        my @present_accnos; # accnos already present
        while (my $accno = $sth->fetchrow_array()) {
            push (@present_accnos, $accno);
        }
        $sth->finish();
        # Get the sequences for the gene to include in the alignment
        my $sth_i = $dbh->prepare("SELECT accno FROM $tables[$i] WHERE cluster='lead' AND LENGTH(sequence)>$min_length") or die "Could not prepare statement: " . $dbh->errstr;
        $sth_i->execute();
        my @add_accnos; # accnos to add
        while (my $accno = $sth_i->fetchrow_array()) {
            my $add = 'y'; # assume the accno should be added
            my $updated; # to store how many rows in the database that were effected
            foreach (@present_accnos) {
                if ($accno eq $_) {
                    $add = 'n'; # if accno already in table flag it
                    last;       # and stop looking
                }
            }
            if ($add eq 'y') { # if the sequence should be added
                my @matching_gene_accnos; # to store the sequences that should be linked between genes
                $matching_gene_accnos[$i] = $accno; # store present accno for this gene
                my @conflicting_gene_accnos; # array to store accnos that are in conflict with the linking (already linked to other accnos)
                # Step I, match accnos between tables based on individuals
                for (my $j=0; $j<scalar @tables; ++$j) { # for each gene
                    if ($j == $i) { next; } # if present gene skip to next
                    my %corresponding = &switch_accno_by_ind($dbh, $tables[$j],$accno); # get accno or cluster lead from other gene that are from same individual
                    if ($corresponding{$accno}) { # if matching accno found
                        $matching_gene_accnos[$j] = $corresponding{$accno}; # store the matching accno
                        # get accno annotation in present table for corresponding gene in other table
                        $sth = $dbh->prepare("SELECT $tables[$j]_accno,$tables[$i]_accno,taxon_name FROM alignments WHERE $tables[$j]_accno=\'$corresponding{$accno}\'") or die "Could not prepare statement: " . $dbh->errstr;
                        $sth->execute();
                        my @duplicate_accnos; # to store posible conflicting accnos
                        while (my @temp = $sth->fetchrow_array()) { # get row by row
                            push (@duplicate_accnos, join ',', @temp); # store potential conflicting accnos
                        }
                        $sth->finish();
                        $conflicting_gene_accnos[$j] = join " \| ", @duplicate_accnos; # store all potential conflicting accnos
                    }
                    else { $matching_gene_accnos[$j] = 'empty'; } # if no corresponding accno set it to 'empty'
                }
                # Step II, if any of the potential conflicting taxa do not have a sequence for the gene join that taxa
                my $join_taxon;
                if ($max_cluster && $max_cluster eq 'y') {
                    my @merge; # taxa to merge
                    for (my $j=0; $j<scalar @conflicting_gene_accnos; ++$j) { # for each gene
                        my @duplicates;
                        if ($conflicting_gene_accnos[$j]) { @duplicates= split / \| /, $conflicting_gene_accnos[$j]; }# get each duplicate
                        foreach my $dup (@duplicates) {             # foreach duplicate
                            my @entrys = split /,/,$dup;      # get conflicting accnos
                            my $f='n'; # assume taxon not previously present
                            foreach my $taxa (@merge) {
                                if ($taxa eq $entrys[2]) { $f='y'; last; } # flag if present
                            }
                            if ($f eq 'n') { push (@merge, $entrys[2]); };   # add taxon to merge if not previously present
                        }
                    }
                    if (@merge) {
                        for (my $j=0; $j<scalar @tables; ++$j) { # for each gene
                            my %values;
                            if ($matching_gene_accnos[$j] && $matching_gene_accnos[$j] ne 'empty') {
                                $sth= $dbh->prepare("SELECT LENGTH($tables[$j].sequence),gb_data.proportion_N FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno=gb_data.accno WHERE $tables[$j].accno='$matching_gene_accnos[$j]'") or die;
                                $sth->execute();
                                my @entries = $sth->fetchrow_array();
                                $values{$matching_gene_accnos[$j]} = $entries[0]*(1-$entries[1]);
                                $sth->finish();
                            }
                            foreach my $taxa (@merge) {
                                $sth=$dbh->prepare("SELECT $tables[$j]_accno FROM alignments WHERE taxon_name='$taxa'") or die;
                                $sth->execute();
                                my $alt_accno=$sth->fetchrow_array();
                                $sth->finish();
                                if ( $alt_accno ne 'empty' && !$values{$alt_accno} ) {
                                    $sth= $dbh->prepare("SELECT LENGTH($tables[$j].sequence),gb_data.proportion_N FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno=gb_data.accno WHERE $tables[$j].accno='$alt_accno'") or die;
                                    $sth->execute();
                                    my @entries = $sth->fetchrow_array();
                                    $sth->finish();
                                    $values{$alt_accno} = $entries[0]*(1-$entries[1]);
                                }
                            }
                            if (%values) {
                                my $maximum = 0;
                                my $keeper = 'empty';
                                foreach my $key (keys %values) {
                                    if ($values{$key} > $maximum) {
                                        $maximum = $values{$key};
                                        $keeper = $key;
                                    }
                                }
                                $matching_gene_accnos[$j] = $keeper;
                            }
                            else { $matching_gene_accnos[$j] = 'empty'; }
                        }
                        my $low_number_taxa;
                        foreach my $taxa (@merge) {
                            my $number=$taxa;
                            $number=~s/Taxon_//;
                            if (!$low_number_taxa || $number > $low_number_taxa) { $low_number_taxa = $number; }
                        }
                        $join_taxon="Taxon_$low_number_taxa";
                        foreach my $taxa (@merge) {
                            if ($taxa && $taxa ne $join_taxon) {
                                my $deleted=$dbh->do("DELETE FROM alignments WHERE taxon_name='$taxa'") or die;
                                print "Deleted $taxa ($deleted) duplicate taxa.\n";
                            }
                        }
                    }
                }
                else {
                    for (my $j=0; $j<scalar @conflicting_gene_accnos; ++$j) { # for each gene
                        if ($conflicting_gene_accnos[$j]) { # if there is conflict on wich accnos should be joined
                            my @duplicates = split / \| /, $conflicting_gene_accnos[$j]; # get each duplicate
                            foreach (@duplicates) {             # foreach duplicate
                                my @entrys = split /,/,$_;      # get conflicting accnos
                                if ( $entrys[1] eq 'empty' ) {   # if there is no conflict (as Darth Vader said)
                                    $join_taxon = $entrys[2];   # join the taxon of matching gene (or the dark side as Darth would say)
                                }
                            }
                        }
                    }
                    # Step III, update database
                    for (my $j=0; $j<scalar @tables; ++$j) { # for each table
                        if ($join_taxon) { # if we have a taxa to join
                            # get the accno for the table and taxa
                            $sth = $dbh->prepare("SELECT $tables[$j]_accno FROM alignments WHERE taxon_name=\'$join_taxon\'") or die "Could not prepare statement: " . $dbh->errstr;
                            $sth->execute();
                            my $alt_link = $sth->fetchrow_array(); # this is an alternative match to what we found before based on individuals
                            $sth->finish();
                            if ( $alt_link eq $matching_gene_accnos[$j] ) {next;} # if it is the same accno as we found before, there is no need to do anything
                            elsif ($alt_link ne 'empty' and (!$matching_gene_accnos[$j] or $matching_gene_accnos[$j] eq 'empty')) { # if we did not have an earlier matching sequence
                                $matching_gene_accnos[$j] = $alt_link; # we keep this
                                next;                                  # and move on
                            }
                            elsif ( $matching_gene_accnos[$j] ne 'empty' and $alt_link ne 'empty') { # damn, we have a conflict and need to resolve it
                                # get data for the alternative accno
                                $sth= $dbh->prepare("SELECT LENGTH($tables[$j].sequence),gb_data.proportion_N FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno=gb_data.accno WHERE $tables[$j].accno=\'$alt_link\'")
                                    or die "Could not prepare statement: " . $dbh->errstr;
                                $sth->execute();
                                my @temp = $sth->fetchrow_array();
                                my $alt_value = $temp[0]*(1-$temp[1]); # value to compare for sequence already in alignments
                                $sth->finish();
                                # get data for the match based on individual
                                $sth= $dbh->prepare("SELECT LENGTH($tables[$j].sequence),gb_data.proportion_N FROM $tables[$j] INNER JOIN gb_data ON $tables[$j].accno=gb_data.accno WHERE $tables[$j].accno=\'$matching_gene_accnos[$j]\'")
                                    or die "Could not prepare statement: " . $dbh->errstr;
                                $sth->execute();
                                @temp = $sth->fetchrow_array();
                                my $pres_value = $temp[0]*(1-$temp[1]); # value to compare for sequence linkt to the sequence to add in this iteration
                                $sth->finish();
                                if ($alt_value >= $pres_value) { $matching_gene_accnos[$j] = $alt_link; next; } # if the alternative is better get it
                            }
                        }

                        if ($conflicting_gene_accnos[$j]) { # if we have potential conflicts
                            my @duplicates = split / \| /, $conflicting_gene_accnos[$j]; # get separate entrys
                            my @duplicate_accnos; # to store accnos
                            my @duplicate_taxon;  # and taxon name
                            foreach (@duplicates) { # for each potential conflict
                                my @entrys = split /,/,$_; # get data
                                if ($join_taxon) { # if we have a taxon for the accno we are checking
                                    if ($entrys[2] ne $join_taxon ) { # if the potential conflicting taxa is not of the same taxa we have a conflict
                                        push (@duplicate_accnos,$entrys[1]); # save accno
                                        push (@duplicate_taxon,$entrys[2]);  # and taxon
                                    }
                                }
                                else { # if we do not have a taxon for the present sequence we have a conflict (since the accno we want already represents other taxa)
                                    push (@duplicate_accnos,$entrys[1]); # save accno
                                    push (@duplicate_taxon,$entrys[2]);  # and taxon
                                }
                            }
                            # get the individual for the taxon we have
                            $sth = $dbh->prepare("SELECT individual FROM gb_data WHERE accno=\'$accno\'") or die "Could not prepare statement: " . $dbh->errstr;
                            $sth->execute();
                            my $individual = $sth->fetchrow_array(); # individual of the present taxon
                            $sth->finish();
                            my @compeating_individual; # to store individuals of taxa competing for the match
                            # get the individual of each duplicate
                            foreach (@duplicate_accnos) {
                                if ($_ eq 'empty') { push(@compeating_individual,'empty'); next; }
                                $sth = $dbh->prepare("SELECT individual FROM gb_data WHERE accno=\'$_\'") or die "Could not prepare statement: " . $dbh->errstr;
                                $sth->execute();
                                while (my $temp = $sth->fetchrow_array()) {
                                    push(@compeating_individual,$temp);
                                }
                                $sth->finish();
                            }
                            # get the individual of the matching sequence
                            $sth = $dbh->prepare("SELECT individual FROM gb_data WHERE accno=\'$matching_gene_accnos[$j]\'") or die "Could not prepare statement: " . $dbh->errstr;
                            $sth->execute();
                            my $corr_individual = $sth->fetchrow_array();
                            $sth->finish();
                            if ($corr_individual eq $individual) { # if the matching sequence is from the same individual as the accno being incorporated they should be in the same taxon
                                for (my $k=0; $k < scalar @duplicate_taxon; ++$k) { # try to find an alternative for the competing taxa
                                    my $alternative; # to store the alternative
                                    if ($compeating_individual[$k] ne 'empty') { # if we have an individual annotation for the competing accno
                                        $alternative = &find_alternative_accno( $dbh, $tables[$j], $compeating_individual[$k], $min_length, $matching_gene_accnos[$j] ); # find alternative
                                    }
                                    if ($alternative) { # if an alternative was found update table
                                        $updated = $dbh->do("UPDATE alignments SET $tables[$j]_accno=\'$alternative\' WHERE taxon_name=\'$duplicate_taxon[$k]\'") or die "Could not uppdate alignments: " . $dbh->errstr;
                                    }
                                    else { # otherwise remove the conflicting accno from the competing taxa
                                        my $columns_string=''; # make a query string to get the accnos for each gene for the taxon
                                        foreach (@tables) {
                                            if ($_ ne $tables[$j]) { $columns_string .= "$_\_accno,"; }
                                        }
                                        $columns_string=~ s/,$//; # remove trailing ','
                                        # get the accnos
                                        $sth = $dbh->prepare("SELECT $columns_string FROM alignments WHERE taxon_name=\'$duplicate_taxon[$k]\'") or die "Could not prepare statement: " . $dbh->errstr;
                                        $sth->execute();
                                        my @columns = $sth->fetchrow_array();
                                        $sth->finish();
                                        my $empty_flag = 'y'; # assume that every column is empty
                                        foreach ( @columns ) {
                                            if ( $_ ne 'empty' ) {
                                                $empty_flag = 'n'; # if not empty flag
                                                last;
                                            }
                                        }
                                        if ( $empty_flag eq 'n' ) { # if at least one column is not empty set the accno to empty
                                            $updated = $dbh->do("UPDATE alignments SET $tables[$j]_accno=\'empty\' WHERE taxon_name=\'$duplicate_taxon[$k]\'") or die "Could not uppdate alignments: " . $dbh->errstr;
                                        }
                                        else { # if the duplicating sequence is the only for the taxon, the taxon is redundant so it is removed
                                            $updated = $dbh->do("DELETE FROM alignments WHERE taxon_name=\'$duplicate_taxon[$k]\'") or die "Could not delete row from alignments: " . $dbh->errstr;
                                        }
                                    }
                                }
                            }
                            else { # if the individual do not match the accno being incorporated find an alternative for the match for this accno
                                my $alternative = &find_alternative_accno( $dbh, $tables[$j], $individual, $min_length, $matching_gene_accnos[$j] ); # get alternative sequence
                                if ($alternative) { # if alternative found
                                    $matching_gene_accnos[$j] = $alternative; # save it
                                }
                                else {
                                    $matching_gene_accnos[$j] = 'empty'; # null the match
                                }
                            }
                        }
                    }
                }
                if ($join_taxon) { # if we have a taxon to join
                    my $update_string = ''; # create a string for updating the table
                    for (my $j=0; $j < scalar @tables; ++$j) {
                        if ($matching_gene_accnos[$j] and $matching_gene_accnos[$j] ne 'empty') {
                            $update_string .= "$tables[$j]_accno=\'$matching_gene_accnos[$j]\',"; # add each column and its value
                        }
                    }
                    $update_string =~ s/,$//; # remove trailing ','
                    # do the uppdate
                    $updated = $dbh->do("UPDATE alignments SET $update_string WHERE taxon_name=\'$join_taxon\'") or die "Could not update $update_string in alignments for $join_taxon: " . $dbh->errstr;
                }
                else { # if we do not have a taxa to join create a new one
                    my $taxa = "Taxon_$taxon_counter";  # create unique name
                    my $column_string = "(taxon_name,"; # create a string with the columns to name
                    my $values_string = "(\'$taxa\',";  # create a string with values for each column
                    ++$taxon_counter; # prepare for next taxa
                    for (my $j=0; $j < scalar @tables; ++$j) { # for each gene
                        if ($matching_gene_accnos[$j] and $matching_gene_accnos[$j] ne 'empty') { # if there is a matching sequence
                            $column_string .= "$tables[$j]_accno,"; # add the gene to the collumns to update
                            $values_string .= "\'$matching_gene_accnos[$j]\',"; # and the accno to put in
                        }
                    }
                    $column_string =~ s/,$/\)/; # remove trailing ','
                    $values_string =~ s/,$/\)/; # -"-
                    $updated = $dbh->do("INSERT INTO alignments $column_string VALUES $values_string") or die "Could not insert new values: " . $dbh->errstr; # insert the taxa
                    if (($taxon_counter-1)%100 == 0) { print "Inserted " . ($taxon_counter-1) . " taxa.\n"; }
                }
            }
        }
        $sth_i->finish();
    }
}
### Subroutine to check for matching sequences in specified table for an array of sequences
sub switch_accno_by_ind {
    my $dbh = shift @_; # get database handler
    my $table = shift @_; # get gene
    my @start_accnos = @_; # get accnos
    my %return_hash; # stor results to return
    for (my $i=0; $i<scalar @start_accnos; ++$i) { # for each accno
        # get the individual annotation
        my $sth = $dbh->prepare("SELECT individual FROM gb_data WHERE accno=\'$start_accnos[$i]\'");
        $sth->execute();
        my $indv = $sth->fetchrow_array();
        $sth->finish();
        # Get data fror sequences of the same individual in other table
        $sth = $dbh->prepare("SELECT $table.accno,$table.cluster,LENGTH($table.sequence),gb_data.proportion_N FROM $table INNER JOIN gb_data ON $table.accno=gb_data.accno WHERE gb_data.individual=\'$indv\'");
        $sth->execute();
        my $accno; # to store accno
        my $lead_flag = 'n'; # if it is the lead in a cluster
        my $value = 0; # to store sequence length times proportion non N
        while ( my @row = $sth->fetchrow_array()) { # for each accno
            my $temp_value = $row[2]*(1-$row[3]); # get sequence length times proportion non N
            if ($row[1] eq 'lead') { # if sequence is lead
                $accno = $row[0]; # get the accno
                last;             # and quit
            }
            elsif ($temp_value > $value and $row[1] ne 'empty') { # if better than competing sequences and in a cluster
                $value = $temp_value; # get value
                $accno = $row[1];     # and store which sequence it points to
            }
        }
        $sth->finish();
        if ($accno) { # if a matching sequence was found
            my $flag = 'n'; # assume accno has not been added to return hash
            foreach (keys %return_hash) { # for each entry in return hash
                if ($accno eq $return_hash{$_} ) { # if the accno has been added to one of the sequences being checked
                    $flag = 'y'; # flag it
                    last; # stop looking
                }
            }
            if ($flag eq 'n') { $return_hash{$start_accnos[$i]} = $accno; } # if not added add accno as match
        }
    }
    return %return_hash; # return
}

### subroutine to find alternative matching accno
sub find_alternative_accno {
    my $dbh = shift @_; # get the database handler
    my $table = shift @_; # get the gene column for which an alternative is needed
    my $individual = shift @_; # get the individual of sequence to get match for
    my $min_length = shift @_; # get the minimum length of viable sequence
    my $not_this = shift @_; # get the accno that it is an alternative for
    my $sth = $dbh->prepare("SELECT $table.accno,LENGTH($table.sequence),gb_data.proportion_N FROM $table INNER JOIN gb_data ON $table.accno=gb_data.accno WHERE gb_data.individual=\'$individual\' AND LENGTH($table.sequence) >= $min_length")
        or die "Could not prepare statement: " . $dbh->errstr;  # get potential alternatives
    $sth->execute();
    my @posibilities;
    my @comp_values;
    while (my @temp = $sth->fetchrow_array()) {
        push(@posibilities, $temp[0]); # save the alternatives
        push(@comp_values, $temp[1]*(1-$temp[2])); # save how good they are (sequence length times proportion non Ns)
    }
    $sth->finish();
    foreach (@posibilities) { # for each possible alternative
        if ($_ eq $not_this) { # if it is the sequence it should be an alternative for
            $_ = 'empty';      # it is not an alternative, so set it to 'empty'
            next;              # move on
        }
        # See if the altenative accno is already in the alignment table
        $sth = $dbh->prepare("SELECT $table\_accno FROM alignments WHERE $table\_accno=\'$_\'") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        my $temp = $sth->fetchrow_array();

        if ($temp) {         # if so
            #print "$temp\n";
            $_ = 'empty';    # it is not an alternative
        }
    }
    my $alternative; # for the value to return
    my $value=0;     # value to compare for the best sequence
    for (my $k=0; $k<scalar @posibilities; ++$k) { # for each potential alternative
        if ($posibilities[$k] eq 'empty') { next; } # if not an alternative, move on
        else {
            if ($comp_values[$k] > $value) {       # if it is the best alternative
                $value = $comp_values[$k];         # update best value
                $alternative = $posibilities[$k];  # store the accno as the best alternative
            }
        }
    }
    return $alternative; # return the best alternative
}
1;
