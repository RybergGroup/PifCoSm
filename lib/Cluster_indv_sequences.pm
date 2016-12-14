package Cluster_indv_sequences;
use strict;
use PifCosm_support_subs;

##########################################################
### Subroutine to cluster individuals within each gene ###
##########################################################

sub cluster_indv_sequences {
    my $database = shift @_;    # variable for database name
    my $min_length = shift @_;  # minimum length of sequences to include in clustering, shorter sequences will have the cluster 'empty'
    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr;
    my @tables = PifCosm_support_subs::get_gene_tables($dbh); # the tables with the separate gene sequences
    foreach my $table (@tables) {
        print "Clustering individuals for $table.\n";
        my $sth = $dbh->prepare("SELECT gb_data.individual,$table.accno,$table.cluster FROM $table INNER JOIN gb_data ON $table.accno=gb_data.accno WHERE gb_data.individual!='empty' AND LENGTH($table.sequence)>=$min_length ORDER BY gb_data.individual"); # get the 'individual' annotation of each accno in the table
        $sth->execute();
        my $previous = 'empty'; # variable for the individual annotation in previous iteration (below), start at 'empty'
        my $n_ind=0; # counter of how many accnos of the same individual
#        my %accnos_data; # store cluster and non N length for each accno belonging to the same individual
        my $cluster;
        my $update=0; # store how many entrys have been uppdated
        # Cluster successively for each accno
        while (my @row = $sth->fetchrow_array()) {
            if ($row[0] eq $previous or $previous eq 'empty') {
                if(!$cluster) { $cluster=$row[1]; }
                else { $cluster .= ",$row[1]"; }
            }
            elsif ($cluster) {
                $update+=PifCosm_support_subs::process_cluster_array($dbh,$cluster,$table);
                $cluster = $row[1];
            }
            else { print STDERR "Warning!!! Something went wrong! No cluster! Individual = $row[0]; accno = $row[1]; gene = $table.\n"; }
#            if ($row[0] eq $previous) { # if the same individual as previous accno
#                ++$n_ind;
#            }
#            elsif ($n_ind > 0) { # if new individual and several individuals in previous individual, cluster them
#                my $lead; # store the "best" accno in the cluster
#                my $value=0; # store sequence length times proportion not n bases to define best accno
#                foreach(keys %accnos_data) { # for each previous accno for the individual
#                    my @temp = split /\|/, $accnos_data{$_}; # get lead and modified seq. length
#                    if ($temp[1] > $value) { # if modified seq. length longer than previous seq. lengths
#                        $value=$temp[1]; # store modified seq. length
#                        $lead=$_;        # store lead accno
#                    }
#                }
#                my @temp = split /\|/, $accnos_data{$lead}; # get value for "best" accno
#                if ($temp[0] eq 'empty') { # if the lead is 'empty'
#                    $update += $dbh->do("UPDATE $table SET cluster='lead' WHERE accno='$lead'"); # set accno to lead in database
#                }
#                elsif ($temp[0] ne 'lead') { # if "best" accnos cluster is not 'empty' and not lead
#                    $lead = $temp[0]; # lead should be the lead of the cluster the accno is in
#                }
#                foreach(keys %accnos_data) { # time to uppdate database
#                    if ($_ eq $lead) { next; } # if the sequence is lead there is no need to update
#                    else {
#                        $update += $dbh->do("UPDATE $table SET cluster='$lead' WHERE accno='$_'"); # point the sequences to the cluster sequence
#                        my @temp = split /\|/, $accnos_data{$_};
#                        if ( $temp[0] ne 'lead' and  $temp[0] ne 'empty' ) {
#                            $update += $dbh->do("UPDATE $table SET cluster='$lead' WHERE cluster='$temp[0]'");
#                        }
#                        else { $update += $dbh->do("UPDATE $table SET cluster='$lead' WHERE cluster='$_'"); }  # if other sequences were clustered to the accno, point them to new lead
#                    }
#                }
#                undef %accnos_data; # null the hasch to prepare for next individual
#                $n_ind = 0; # null number of accnos in individual
#            }
#            elsif (%accnos_data) { # if only one accno for individual
#                my @keys = keys %accnos_data;
#                my @temp = split /\|/, $accnos_data{$keys[0]}; # get cluster for accno
#                if ($temp[0] eq 'empty') { # if not previously clustered
#                    $update += $dbh->do("UPDATE $table SET cluster='lead' WHERE accno='$keys[0]'"); # its the lead in it's own cluster
#                }
#                undef %accnos_data; # null the hasch to prepare for next individual
#                $n_ind = 0; # null number of accnos in individual
#            }
#            my $temp_accno;
#            if ($row[2] eq 'lead' or $row[2] eq 'empty') { $temp_accno=$row[1]; }
#            else { $temp_accno = $row[2]; }
#            my $value_sth = $dbh->prepare("SELECT LENGTH($table.sequence),gb_data.proportion_N FROM $table INNER JOIN gb_data ON $table.accno=gb_data.accno WHERE gb_data.accno='$temp_accno'") or die;
#            $value_sth->execute();
#            my @values_array = $value_sth->fetchrow_array();
#            $accnos_data{$row[1]} = $row[2] . '|' . $values_array[0]*(1-$values_array[1]); # save info for accno
            $previous = $row[0]; # save individual 'name'
#            $value_sth->finish();
        }
        $sth->finish();
        if ($cluster) { $update+=PifCosm_support_subs::process_cluster_array($dbh,$cluster,$table); }
        # repeat above procedure for last individual
#        if ($n_ind > 0) {
#            my $lead;
#            my $value=0;
#            foreach(keys %accnos_data) {
#                my @temp = split /\|/, $accnos_data{$_};
#                if ($temp[1] > $value) {
#                    $value=$temp[1];
#                    $lead=$_;
#                }
#            }
#            my @temp = split /\|/, $accnos_data{$lead};
#            if ($temp[0] eq 'empty') {
#                $update += $dbh->do("UPDATE $table SET cluster=\'lead\' WHERE accno=\'$lead\'");
#            }
#            elsif ($temp[0] ne 'lead') {
#                $lead = $temp[0];
#            }
#            foreach(keys %accnos_data) {
#                if ($_ eq $lead) { next; }
#                else {
#                    $update += $dbh->do("UPDATE $table SET cluster=\'$lead\' WHERE accno=\'$_\'");
#                    my @temp = split /\|/, $accnos_data{$_};
#                    if ( $temp[0] ne 'lead' and  $temp[0] ne 'empty' ) {
#                        $update += $dbh->do("UPDATE $table SET cluster=\'$lead\' WHERE cluster=\'$temp[0]\'");
#                    }
#                    else { $update += $dbh->do("UPDATE $table SET cluster=\'$lead\' WHERE cluster=\'$_\'"); }  # if other sequences were clustered to the accno, point them to new lead
#                }
#            }
#            undef %accnos_data;
#            $n_ind = 0;
#        }
#        elsif (%accnos_data) {
#            my @keys = keys %accnos_data;
#            my @temp = split /\|/, $accnos_data{$keys[0]};
#            if ($temp[0] eq 'empty') {
#                $update += $dbh->do("UPDATE $table SET cluster=\'lead\' WHERE accno=\'$keys[0]\'");
#            }
#            undef %accnos_data;
#            $n_ind = 0;
#        }
        print "Changed $update rows for $table.\n";
    }
    $dbh->disconnect();
}
1;
