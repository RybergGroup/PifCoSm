###############################################################
### Subroutine to parse genbank data into a SQLite database ###
###############################################################
use strict;
package Gb_parser;
use PifCosm_support_subs;

sub gb_parser {
    my $database = shift @_; # The first argument will be used as database name
    my $gbfile = shift @_ ; # The second argument will be used as name on file with gen bank data
    my $max_seq_length = shift @_; # The thired argument will be used to limit the length of
                                   # sequences includd in the database
    # Variables to store separate fields for each GenBank entry
    my $sequence = 'empty';
    my $accession = 'empty';
    my $species = 'empty';
    my $voucher = 'empty';
    my $isolate = 'empty';
    my $strain = 'empty';
    my $country = 'empty';
    my $collected_by = 'empty';
    my $note = 'empty';
    my $identified_by = 'empty';
    my $environmental_sample = 'No';
    my $host = 'empty';
    my $lat = 'NULL';
    my $lon = 'NULL';
    my $sub_species = 'empty';
    my $variety = 'empty';
    my $taxon_string = 'empty';
    my $reference = 'empty';
    # Flags to tell wich field is being pased for fields longer than a row
    my $field_flag = "start";
    my $note_flag = 'n';
    my $reference_flag = 'n';
    # Variable to count occurences of sequences longer than the max
    my $over_max_counter = 0;
    # Variable to count number of included rows
    my $update = 0;
    # Create table for gb data if there isn't one
    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr;
    my $gb_data_flag = 'n'; # Assume ther is no table for the gb data
    # check if there is a table for gb data
    my $sth = $dbh->prepare("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name" ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    while (my $table = $sth->fetchrow_array()) {
        if ($table eq "gb_data") { $gb_data_flag = 'y'; } # if there is a table named gb_data flag it
    }
    $sth->finish();
    # if no table for gb data was found create one
    if ( $gb_data_flag eq 'n' ) {
            $dbh->do("CREATE TABLE gb_data (accno TEXT PRIMARY KEY, species TEXT DEFAULT 'empty', sub_species TEXT DEFAULT 'empty', variety TEXT DEFAULT 'empty', voucher TEXT DEFAULT 'empty', isolate TEXT DEFAULT 'empty', strain TEXT DEFAULT 'empty', taxon_string TEXT DEFAULT 'empty', country TEXT DEFAULT 'empty', collected_by TEXT DEFAULT 'empty', identified_by TEXT DEFAULT 'empty', environmental_sample TEXT DEFAULT 'No', host TEXT DEFAULT 'empty', latitude REAL, longitude REAL, reference TEXT DEFAULT 'empty', sequence TEXT DEFAULT 'empty', proportion_N REAL)") or die "Could not create table gb_data: " . $dbh->errstr;
    }

    # start to parse gb data file one row at the time
    open INFILE, $gbfile or die "Could not open $gbfile: $!.\n"; # open gb data file
    print "Inserted... ";
    $| = 1; # Flush buffert after each print
    while ( my $infile = <INFILE> ) { # read next row
        # If at end of entery print entery
        if ( $infile =~ /^\/\// ) { # if at end of a entery (marked by //) insert the data into the database
            if (length($sequence) <= $max_seq_length) { # only include sequences shorter than the given max length
                    $reference =~ s/\`//g;
                    $reference = $dbh->quote($reference);
                    $identified_by =~ s/\`//g;
                    $identified_by = $dbh->quote($identified_by);
                    $collected_by =~ s/\`//g;
                    $collected_by = $dbh->quote($collected_by);
                    $country =~ s/\`//g;
                    $country = $dbh->quote($country);
                    $species =~ s/\`//g;
                    $species = $dbh->quote($species);
                    $sub_species =~ s/\`//g;
                    $sub_species = $dbh->quote($sub_species);
                    $variety =~ s/\`//g;
                    $variety = $dbh->quote($variety);
                    $voucher =~ s/\`//g;
                    $voucher = $dbh->quote($voucher);
                    $isolate =~ s/\`//g;
                    $isolate = $dbh->quote($isolate);
                    $strain =~ s/\`//g;
                    $strain = $dbh->quote($strain);
                    $host =~ s/\`//g;
                    $host = $dbh->quote($host);
                    $taxon_string =~ s/\`//g;
                    $taxon_string = $dbh->quote($taxon_string);
                    $sequence = uc($sequence);
                    my $n = $sequence;
                    $n =~ s/[^N]//g;
                    my $proportion_n = length($n)/length($sequence);

                    if (!$lat) {$lat = 'NULL';}
                    if (!$lon) {$lon = 'NULL';}
                    $update += $dbh->do("INSERT INTO gb_data (accno,species,sub_species,variety,voucher,isolate,strain,taxon_string,country,collected_by,identified_by,environmental_sample,host,latitude,longitude,reference,sequence,proportion_N) values ('$accession',$species,$sub_species,$variety,$voucher,$isolate,$strain,$taxon_string,$country,$collected_by,$identified_by,'$environmental_sample',$host,$lat,$lon,$reference,'$sequence',$proportion_n)") or die "Could not insert into gb_data: " . $dbh->errstr;
                    if ($update % 10000 == 0) {
                        print "$update seq... ";
                        #print "'$accession','$species','$sub_species','$variety','$voucher','$isolate','$strain','$country','$collected_by','$identified_by','$environmental_sample','$host',$lat,$lon,'$reference','$sequence'\n";
                    }
            }
            else { ++$over_max_counter; } # count the number of sequences longer than max length
            # reset variables for each field to default
            $sequence = 'empty';
            $accession = 'empty';
            $species = 'empty';
            $taxon_string = 'empty';
            $voucher = 'empty';
            $isolate = 'empty';
            $strain = 'empty';
            $country = 'empty';
            $collected_by = 'empty';
            $note = 'empty';
            $identified_by = 'empty';
            $environmental_sample = 'No';
            $host = 'empty';
            $field_flag = "start";
            $note_flag = 'n';
            $lat = 'NULL';
            $lon = 'NULL';
            $sub_species = 'empty';
            $variety = 'empty';
            $reference = 'empty';
            $reference_flag = 'n';
            next; # done on this row, move on
        }
        # Otherwise see if we enter new field
        # if so, note this and move on to next row
        elsif ( $infile =~ /^ORIGIN/ ) {
            $field_flag = 'origin';
            undef $sequence;
            next;
        }
        elsif ( $infile =~ /^ACCESSION\s+(\w+)/ ) {
            $accession = $1;
            next;
        }
        elsif ( $infile =~ /^FEATURES/ ) {
            $field_flag = 'features';
            next;
        }
        elsif ( $infile =~ /^\s*ORGANISM\s+(.+)/ ) {
            $species = $1;
	    $species =~ s/ var\.|f\.|ssp. \w+$//;
            $field_flag = 'organism';
            next;
        }
        elsif ( $infile =~ /^REFERENCE/ ) {
            $field_flag = 'reference';
            if ( $infile =~ /^REFERENCE\s+([0-9]+)/ ) {
                if ( $1 == 1 ) {
                    $field_flag = 'reference1';
                }
                else { $reference_flag = 'n'; }
            }
            next;
        }
        # Parse depending on what field we are in
        if ($field_flag eq 'origin') { # if in origin field
            $infile =~ s/[^a-zA-Z]//g; # delete everything that is not letters
            $sequence .= $infile;      # and add it to the sequence
        }
        elsif ( $field_flag eq 'organism' ) { # if in organism field
            chomp($infile);                   # remove end of line
            $infile =~ s/^\s+//;              # remove leading white spaces
            $infile =~ s/\s+$//;              # remove trailing white spaces
            if ($infile =~ /^[\w ]+;|(\.$)/) {
                if ($infile =~ s/\.$//) {         # if removing dot at the end
                    $field_flag = 'empty';        # reset the field flag
                }
                if ($taxon_string eq 'empty') { $taxon_string = $infile; } # if nothing has been read to the taxon string it is set to the infile
                else { $taxon_string .= " $infile"; } # if there is something in the taxon_string, add infile to it
            }
            else { print STDERR "'$infile' for $accession does not pass test for inclusion in taxon string. May cause problems in downstream modules.\n"; }
        }
        elsif ( $field_flag eq 'reference1') {    # if we are in the reference field
            if ($infile =~ /\s+AUTHORS\s+(.+)/) { # if we are reading authors
                $reference = $1;                  # get whats there (after AUTHOR)
                $reference_flag = 'y';            # flag that we got something
            }
            elsif ($infile =~ /\s+TITLE\s+(.+)/) { # Similatily get the title
                $reference .= " $1";
                $reference_flag = 'y';
            }
            elsif ($infile =~ /\s+JOURNAL\s+(.+)/) { # and journal
                $reference .= ". $1";
                $reference_flag = 'y';
            }
            elsif ($reference_flag eq 'y' and $infile =~ /\s+(.+)/) { #if we got something add it to the reference annotation
                $reference .= " $1";
            }
            chomp($reference); # rempve any end of line
        }
        elsif ( $field_flag eq 'features' ) {                 # if in the features field, for the following annotations
            if ( $infile =~ /\/specimen_voucher="([^"]+)/ ) { # get whats after = and between the " signs
                $voucher = $1;                                # store it in apropriate variable
                chomp($voucher);                              # get rid of new line character
            }
            elsif ( $infile =~ /\/strain="([^"]+)/ ) {        # -"-
                $strain = $1;
                chomp ($strain);
            }
            elsif ( $infile =~ /\/isolate="([^"]+)/ ) {       # -"-
                $isolate = $1;
                chomp ($isolate);
            }
            elsif ( $infile =~ /\/country="([^"]+)/ ) {       # -"-
                $country = $1;
                chomp ($country);
            }
            elsif ( $infile =~ /\/collected_by="([^"]+)/ ) {  # -"-
                $collected_by = $1;
                chomp ($collected_by);
            }
            elsif ( $infile =~ /\/identified_by="([^"]+)/ ) { # -"-
                $identified_by = $1;
                chomp ($identified_by);
            }
            elsif ( $infile =~ /\/environmental_sample/ ) {   # -"-
                $environmental_sample = 'Yes';
            }
            elsif ( $infile =~ /\/host="([^"]+)/ ) {          # -"-
                $host = $1;
                chomp ($host);
            }
            elsif ( $infile =~ /\/sub_species="([^"]+)/ ) {   # -"-
                $sub_species = $1;
                chomp ($sub_species);
            }
            elsif ( $infile =~ /\/variety="([^"]+)/ ) {       # -"-
                $variety = $1;
                chomp ($variety);
            }
            elsif ( $infile =~ /\/lat_lon="([^"]+)/ ) {       # if long lat annotation, take care of different ways to annotate
                my @temp = split /\s+/,$1;                    # separate the different parts
                if ( !($temp[0] =~ /[^0-9\.]/) and !($temp[2] =~ /[^0-9\.]/) and # if first and third value are numbers
                     ($temp[1] eq 'N' or $temp[1] eq 'n' or $temp[1] eq 'S' or $temp[1] eq 's') and # and second north-south
                     ($temp[3] eq 'E' or $temp[3] eq 'e' or $temp[3] eq 'W' or $temp[3] eq 'w') ) { # and forth east-west
                    $lat = $temp[0]; # use first as latitude
                    if ($temp[1] eq 'S' or $temp[1] eq 's') { $lat *= -1; } # if south make it negative
                    $lon = $temp[2]; # use third as longitude
                    if ($temp[3] eq 'W' or $temp[3] eq 'w') { $lat *= -1; } # if west make it negative
                }
            }
            elsif ( $infile =~ /\/note="([^"]+)/ and $infile eq 'empty') { # see above
                $note = $1;
                chomp ($note);
            }
            elsif ( $note_flag eq 'y' ) {  # if in note field
                if ( $infile =~ /(.*)"/) { # if end "
                    $note .= " $1";        # add to note
                    $note_flag = 'n';      # and note that the note has ended
                }
                else { $note .= $infile; } # else just add to note
                chomp($note);              # remove newline char
            }
        }
    }
    $| = 0; # reset default buffering
    close INFILE or die;
    $dbh->disconnect(); # disconnect from database
    print "$update seq. \n";
    print "$over_max_counter sequences were longer than $max_seq_length and not included in $database.\n"; # output how many sequences were longer than max length
}
1;
