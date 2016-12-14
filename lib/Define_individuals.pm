package Define_individuals;
use strict;
use PifCosm_support_subs;

########################################
### Subroutine to define individuals ###
########################################

sub define_individuals {
    my $database = shift @_;
    my $n_columns = shift @_; # the number of columns used to define individuals
    my @columns; # the columns used to define individuals
    my $i; # counter
    for ($i=0; $i<$n_columns; ++$i) {
        $columns[$i] = shift @_; # get column names
    }
    my $conditional = shift @_; # a column that has to be the same for individuals
    my $ind_number=1; # number to designate different individuals

    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr;
    my $indv_column_flag = PifCosm_support_subs::column_present($dbh,'gb_data','individual'); #'n'; # flag if individual column present
    # if no column for individuals, make one
    if ($indv_column_flag eq 'n') {
        my $updates = $dbh->do("ALTER TABLE gb_data ADD individual TEXT DEFAULT 'empty'") or die "Couldn't uppdate database: " . $dbh->errstr;
        print "Added the column individual to gb_data.\n";
    }
    else { # otherwise
        print "Column for individuals already present in table gb_data. Uppdating it.\n";
        my $sth = $dbh->prepare("SELECT individual FROM gb_data WHERE individual!='empty'") or die "Couldn't prepare statement: " . $dbh->errstr; # get the names of the individuals
        $sth->execute();
        while (my $indv = $sth->fetchrow_array()) {
            if ( $indv =~ s/Ind_// ) {
                if ($indv > $ind_number) { $ind_number = $indv; } # get the highest number designating an individual
            }
        }
        $sth->finish();
        print "The highest numbered individual is Ind_$ind_number.\n";
        ++$ind_number; # prepare for new individual
    }
    # for each column that defines an individual
    for ($i=0; $i < scalar @columns; ++$i) {
        print "Defining individuals based on $columns[$i].\n";
        my $sth;
        if ($conditional) {
            $sth = $dbh->prepare("SELECT $columns[$i],individual,$conditional FROM gb_data WHERE $columns[$i]!='empty' and $conditional!='empty' ORDER BY $conditional,$columns[$i]")
                or die "Couldn't prepare statement: " . $dbh->errstr; # get data for column
        }
        else {
            $sth = $dbh->prepare("SELECT $columns[$i],individual FROM gb_data WHERE $columns[$i]!='empty' ORDER BY $columns[$i]")
                or die "Couldn't prepare statement: " . $dbh->errstr; # get data for column
        }
        $sth->execute();
        my $n_equal=0; # keep track of number of accnos belonging to the individual
        my @previous= ('empty','empty'); # keep track of values in previous iteration
        my $individual; # store individual name
        while (my @rows = $sth->fetchrow_array()) { # for each column entry
            if (($conditional and ($rows[0] eq $previous[0] and $rows[2] eq $previous[1])) or (!$conditional and ($rows[0] eq $previous[0]))) { # if column and conditional annotation are the same as in previous iteration
                ++$n_equal; # count up the numbers of entrys with identical annotation
                if ( $rows[1] ne "empty" ) { # if an individual is defined
                    if (!$individual) {
                        # if matching previous defined individual use this name
                        $individual = $rows[1];
                    }
                    elsif ( $rows[1] ne $individual ) {
                        # if not same as individual of previous accno something is suspect and is flaged by putting individual as empty
                        $individual = 'empty';
                    }
                }
            }
            elsif ($n_equal > 0 ) { # if not same as previous iteration and more than one accno
                if (!$individual) { # if we do not have an individual annotation yet
                    # create new individual
                    $individual = "Ind_$ind_number";
                    ++$ind_number; # prepare for next individual
                }
                $previous[0]=~ s/'/''/g;
                $previous[1]=~ s/'/''/g;
                if ($conditional) {
                    my $updates = $dbh->do("UPDATE gb_data SET individual=\'$individual\' WHERE $columns[$i]=\'$previous[0]\' AND $conditional=\'$previous[1]\' AND individual=\'empty\'") or print "$individual $columns[$i] $previous[0] $conditional $previous[1]\n"; # put individual annotation ino database
                }
                else {
                    my $updates = $dbh->do("UPDATE gb_data SET individual=\'$individual\' WHERE $columns[$i]=\'$previous[0]\' AND individual=\'empty\'") or print "$individual $columns[$i] $previous[0]\n"; # put individual annotation ino database
                }
                undef $individual; # prepare for next individual
                $n_equal = 0; # prepare for next individual
            }
            if ($conditional) { @previous = ($rows[0],$rows[2]); } # store column and conditional annotation to next iteration
            else { @previous = ($rows[0],'empty'); }
        }
        $sth->finish();
        print "Highest number for an individual is ". ($ind_number-1) . ".\n";
    }
    # set individual for remaining accnos, that did not group with others
    print "Assigning singleton sequences to individuals.\n";
    my $sth = $dbh->prepare("SELECT accno FROM gb_data WHERE individual='empty'") or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    while (my $accno = $sth->fetchrow_array()) {
        my $updates = $dbh->do("UPDATE gb_data SET individual=\'Ind_$ind_number\' WHERE accno=\'$accno\'");
        ++$ind_number;
    }
    $sth->finish();
    $dbh->disconnect();
    print "Highest number for an individual is ". ($ind_number-1) . ".\n";
}
1;
