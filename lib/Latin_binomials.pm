package Latin_binomials;
use strict;
use PifCosm_support_subs;

################################################################################
### Subroutines to get 'latin binomials' for the taxa in the alignment table ###
################################################################################
sub create_species_names_column {
    my $database = shift @_; # name of database
    my $dbh = PifCosm_support_subs::connect_to_database($database); # get database handler
    my @gene_names = PifCosm_support_subs::get_gene_tables($dbh); # get the gene tables
    # create a column for the names
    my $add_column = PifCosm_support_subs::column_present($dbh,'alignments','species_name');
    if ($add_column eq 'n') { my $update = $dbh->do("ALTER TABLE alignments ADD species_name TEXT DEFAULT 'empty'") or die "Could not add column species_name: " . $dbh->errstr; }
    &insert_new_species_names ($dbh,@gene_names); # get 'latin binomials' for each taxon
}

sub insert_new_species_names {
    my $dbh = shift @_; # get database handler
    my @gene_names = @_; # get the gene names
    my $gene_columns=''; # for the query string
    my $update=0;
    #print "Finding appropriate species names.\n";
    foreach (@gene_names) { $gene_columns .= ",$_\_accno"; } # add each column to the query string
    # get the taxon names and accnos for each taxon
    my $sth = $dbh->prepare("SELECT taxon_name$gene_columns FROM alignments WHERE species_name='empty'") or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    while (my @row = $sth->fetchrow_array()) { # for each taxon
        my %names; # hash to count the number of times each name occure
        my @accnos; # array for the accnos
        for (my $i=1; $i < scalar @row; ++$i) { # for each gene
            if ($row[$i] ne 'empty') { # if there is an accno
                my $already_checked = 'n'; # assume the accno has not been checked before
                foreach (@accnos) {
                     if ($_ eq $row[$i]) { $already_checked = 'y'; } # flag if accno checked before
                }
                if ($already_checked eq 'y') { next; } # if accno already checked skip to next
                push (@accnos, $row[$i]); # else, add accno
                # get species name for accno
                my $name_sth = $dbh->prepare("SELECT species FROM gb_data WHERE accno='$row[$i]'") or die "Could not prepare statement: " . $dbh->errstr;
                $name_sth->execute();
                my $string = $name_sth->fetchrow_array();
                ++$names{$string}; # count up the number of times the name has occured
                $name_sth->finish();
            }
        }
        my $species; # for the species name
        my $count=0; # to store highest occurens of a name
        foreach (keys %names) { # for each posible name
            if ($names{$_} > $count) { # if it is the most frequent name for the taxa
                $count = $names{$_};   # save the number of times it occured
                $species = $_;         # save the name
            }
        }
        # remove problematic characters
        $species =~ s/\s/_/g;
        $species =~ s/\(.*\)//g;
        $species =~ s/[-,:;'"\.\\\/]/_/g;
        $species =~ s/__/_/g;
        # get previous taxon annotated as the same species
        my $species_sth = $dbh->prepare("SELECT species_name FROM alignments WHERE species_name LIKE '$species%'")  or die "Could not prepare statement: " . $dbh->errstr;
        $species_sth->execute();
        $count=0; # to store the highest number added to taxon name
        while (my $species_no = $species_sth->fetchrow_array()) { # for each taxon
            $species_no =~ /_([0-9]+)$/; # get the number added to make it distinct
            if ($1 > $count) { $count = $1; } # if the number is the highest so far store it
        }
        $species_sth->finish();
        ++$count; # increase one from the highest number
        # update taxon with the species name
        my $string = $dbh->quote("$species\_$count");
        $update += $dbh->do("UPDATE alignments SET species_name=$string WHERE taxon_name='$row[0]'") or die "Could not update species name in alignments for $row[0]: " . $dbh->errstr;
        #print "Updated $update rows in alignments for $row[0].\n";
    }
    $sth->finish();
    print "Gave $update taxa 'latin binomials' in alignments.\n";
}
1;
