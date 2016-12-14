###########################################################################################################
### Subroutine to create table of which taxa to use as minimal alignment groups for the different genes ###
###########################################################################################################

package Create_alignment_groups;
use strict;
use PifCosm_support_subs;

sub create_alignment_groups {
    my $database = shift @_; # database name
    my $min = shift @_; # if to define on lowest or highest level
    my $dbh = PifCosm_support_subs::connect_to_database($database); # get database handler
    my @genes = PifCosm_support_subs::get_gene_tables($dbh); # get all gene names
    my $n_groups = 0;
    my $taxon_level; # if to lock in first or last entry in array
    if ($min eq 'y') { $taxon_level = -1; } # if lowest level should be used look in last entry
    else { $taxon_level = 0; } # if highest look in first entry
    #my @taxa;
    #my @hirarchy;
    my @tables = PifCosm_support_subs::get_tables($dbh); # get all table names
    foreach(@tables) {
        if ($_ eq 'alignment_groups') { # if alignment_groups present
            $dbh->do("DROP TABLE alignment_groups"); # remove it
            last;
        }
    }
    $dbh->do("CREATE TABLE alignment_groups (gene TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', tree TEXT DEFAULT 'empty', tree_method TEXT DEFAULT 'empty',alignable INTEGER, PRIMARY KEY (gene, taxon))") or die "Could not create table gb_data: " . $dbh->errstr; # create alignment_group table
    foreach (@genes) { # for each gene
        my @taxa; # the alignment group taxa for this gene
        my @hirarchy; # the taxonomic hirarchy for this gene
        my $sth = $dbh->prepare("SELECT DISTINCT taxon_string FROM gb_data INNER JOIN $_ ON gb_data.accno=$_.accno WHERE $_.cluster='lead'"); # get distinct taxonomic string
        $sth->execute();
        while (my $taxon_string = $sth->fetchrow_array() ) { # for each distinct taxonomic string
            my @temp = split /;/,$taxon_string; # separate taxonomic levels
            $temp[$taxon_level] =~ s/^\s+//; # remove leading white space in the level that should be used
            $temp[$taxon_level] =~ s/\s+$//; # remove trailing white space in the level that should be used
            $temp[$taxon_level] =~ s/\s/_/g; # remove white space in the level that should be used            
            if ($min eq 'y') {
                while ($temp[$taxon_level] =~ /^[a-z]/) { # if using lowest level, remove lowest level if it begins with lowercase letter
                    pop @temp;
                    $temp[$taxon_level] =~ s/^\s+//; # clean white spaces
                    $temp[$taxon_level] =~ s/\s+$//;
                    $temp[$taxon_level] =~ s/\s/_/g;
                }
            }
            my $new = 'y'; # assume taxon is new
            foreach my $rank (@hirarchy) { if ($rank =~ /\b$temp[$taxon_level]\b/ ) { $new = 'n'; last;} } # look if taxon is already in @hirarchy
            if ($new eq 'y') { # if not just higher level in same hirarchy as previous taxa
                foreach my $old (@taxa) { if ($old =~ /\b$temp[$taxon_level]\b/ ) { $new = 'n'; last;} } # may be in taxa avoiding lowest level began with lowercase letter
                if ($new eq 'y') { # if not among taxa and not higher taxa in hirarchy
                    my $i = 0; # counter for old hirarche
                    my $j = 0; # counter for new cleaned hirarchy
                    my @temp_taxa = @taxa; # save old array of taxa to align
                    undef @taxa; # make room for cleaned array of taxa to align
                    my @temp_hirarchy = @hirarchy; # save old taxon string array
                    undef @hirarchy; # make room for new taxon string array
                    for ($i = 0; $i < scalar @temp_taxa; ++$i) { # for each old taxon
                        if (!($taxon_string =~ /\b$temp_taxa[$i]\b/)) { # if the old taxon is not just higher taxa in the hirarchy of this taxa
                            $taxa[$j]=$temp_taxa[$i];               # keep it
                            $hirarchy[$j]=$temp_hirarchy[$i];       # and its hirarchy
                            ++$j;                                   # and move on
                        }                                           # else do not keep it
                    }
                    push (@taxa,$temp[$taxon_level]);               # add new taxa to the cleaned old once
                    push (@hirarchy,$taxon_string);                 # same for hirarchy
                }
            }
        }
        #my $k=0;
        for (my $i=0; $i < scalar @taxa; ++$i) {
            my $string_taxa = $dbh->quote($taxa[$i]);
            $n_groups += $dbh->do("INSERT INTO alignment_groups (gene, taxon, alignable) VALUES ('$_',$string_taxa,1)");
        }
    }
    return $n_groups;
}
1;
