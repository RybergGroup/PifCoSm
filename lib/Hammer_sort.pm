package Hammer_sort;
use strict;
use PifCosm_support_subs;

#################################################################
### Sub routines to sort genes to their own tables using HMMER ###
#################################################################
sub hammer_sort {
    my $database = shift @_; #use the first argument as the name of the fasta file 
    my $path = shift @_; # path to hmmscan
    my $hmmdatabase = shift @_; #takes the name of the database from the second argument
    my $min_sequence_length_out = shift @_; # min length of sequence to include in gene table
    my $e_value_cut_off = shift @_; # cut off for e-value
    my $n_cut_off = shift @_; # cut off for how large fraction of n's we accept
    my $print_non_match = shift @_;
    my %intergenes = @_; # defines genes situated between other genes
    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr; # database connection
    my $sequence="empty"; # variable for sequence
    my %hmmresults; # hash to store the result from the HMMER search
    my $n_included_seg=0; # to store the number of sequences that are identified as one or several specific genes
    my $n_gene_entries=0; # number of separate enteries made to gene tables
    my $n_revcomp = 0; # number of sequences that has been reversed and complemented
    # get sequences to sort 
    print "The following genes are defined to be between other genes\n";
    foreach (keys %intergenes) { print "$_ is between $intergenes{$_}.\n"; }
    my $sth = $dbh->prepare("SELECT accno FROM gb_data WHERE proportion_N<$n_cut_off AND LENGTH(sequence)>$min_sequence_length_out") or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    my @accnos; # array to store acnos
    while ($accnos[scalar @accnos] = $sth->fetchrow_array()) { } # get each accno
    $sth->finish();
    # Check with gene each accno belong to
    my $i=0;
    if ($print_non_match eq 'y') {
        open NONMATCH, ">seq_not_matching_genes.txt" or die "Could not open seq_not_matching_genes.txt: $!.\n";
    }
    for ($i=0; $i<scalar @accnos; ++$i) {
        if (!$accnos[$i]) { next; }
        # get sequence for accno
        $sth = $dbh->prepare("SELECT sequence FROM gb_data WHERE accno='$accnos[$i]'") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        $sequence = $sth->fetchrow_array();
        $sth->finish();
        # if sequence is found
        if ($sequence and $sequence ne "empty") {
            %hmmresults = &runhammer($path,$accnos[$i],$sequence,$e_value_cut_off,$hmmdatabase); # get results from HMMER
            #If no match try reverse and complement sequence
            if ($hmmresults{"empty"} and $hmmresults{"empty"} eq "empty") {
                my $temp = &reverse_and_compliment($sequence);
                %hmmresults = &runhammer($path,$accnos[$i],$temp,$e_value_cut_off,$hmmdatabase);
                #if the reverse and compliment sequence gives a result save the sequence in rev. comp.
                if (!$hmmresults{"empty"} or $hmmresults{"empty"} ne "empty") {
                    $sequence = $temp;
                    $n_revcomp += $dbh->do("UPDATE gb_data SET sequence='$sequence' WHERE accno='$accnos[$i]'") or die "Could not performe update: " . $dbh->errstr;;
                }
            }
            if (%hmmresults) { # if any results from HMMER
                my $result_flag = 'f'; # flag to tell if any meaningful result
                foreach (keys %hmmresults) { # for each entry
                    if ($_ ne 'empty' and $hmmresults{$_} =~ /;/) { $result_flag = 't'; } # if key not 'empty' and there is a ; separating the different enteries flag as ok
                }
                if ( $result_flag eq 't' ) { # if ok process results
                    # interprete HMMER results and split and put sequence fragments in appropriate tables (if any)
                    my ($temp1,$temp2) = &sort_and_output_fasta($accnos[$i], $sequence, $dbh, %intergenes, '|||', %hmmresults);
                    $n_included_seg += $temp1; # keep count
                    $n_gene_entries += $temp2; # keep count
                }
                elsif ($print_non_match eq 'y') { print NONMATCH "$accnos[$i]\n"; }
            }
        }
        # if no sequence!
        else { print STDERR "WARNING!!! No sequence for $accnos[$i].\n"; }
        if (($i+1) % 10000 == 0) { print $i+1 . " accnos processed, $n_revcomp were reversed and complemented, $n_included_seg were identified as belonging to one or more specific genes, and $n_gene_entries specific gene enteries have been made.\n"; }
    }
    $dbh->disconnect(); #finish
    if ($print_non_match eq 'y') {
        close NONMATCH or die;
    }
    print $i+1 . " accnos processed, $n_included_seg were identified as belonging to one or more specific genes, and $n_gene_entries specific gene enteries have been made.\n";
}

### Subrutines associated with sorting using HAMMER ###

#Reverses and complement a sequence inputted as a string
sub reverse_and_compliment {
    my @sequence = reverse(split //, $_[0]); #reverse sequence and split to array
    for (my $i=0; $i < scalar @sequence; ++$i) { # complement each base, only capital letters in output
        if ($sequence[$i] eq 'a' or $sequence[$i] eq 'A') { $sequence[$i] = 'T'; }
        elsif ($sequence[$i] eq 't' or $sequence[$i] eq 'T') { $sequence[$i] = 'A'; }
        elsif ($sequence[$i] eq 'g' or $sequence[$i] eq 'G') { $sequence[$i] = 'C'; }
        elsif ($sequence[$i] eq 'c' or $sequence[$i] eq 'C') { $sequence[$i] = 'G'; }
        elsif ($sequence[$i] eq 'y' or $sequence[$i] eq 'Y') { $sequence[$i] = 'R'; }
        elsif ($sequence[$i] eq 'r' or $sequence[$i] eq 'R') { $sequence[$i] = 'Y'; }
        elsif ($sequence[$i] eq 'w' or $sequence[$i] eq 'W') { $sequence[$i] = 'W'; }
        elsif ($sequence[$i] eq 's' or $sequence[$i] eq 'S') { $sequence[$i] = 'S'; }
        elsif ($sequence[$i] eq 'k' or $sequence[$i] eq 'K') { $sequence[$i] = 'M'; }
        elsif ($sequence[$i] eq 'm' or $sequence[$i] eq 'M') { $sequence[$i] = 'K'; }
        elsif ($sequence[$i] eq 'b' or $sequence[$i] eq 'B') { $sequence[$i] = 'V'; }
        elsif ($sequence[$i] eq 'd' or $sequence[$i] eq 'D') { $sequence[$i] = 'H'; }
        elsif ($sequence[$i] eq 'h' or $sequence[$i] eq 'H') { $sequence[$i] = 'D'; }
        elsif ($sequence[$i] eq 'v' or $sequence[$i] eq 'V') { $sequence[$i] = 'B'; }
        elsif ($sequence[$i] eq 'n' or $sequence[$i] eq 'N') { $sequence[$i] = 'N'; }
    }
    return join "", @sequence; # join sequence to a string and return it

}

# sub to create table for gene
sub test_and_create_table {
    my $dbh = shift @_;
    my $lc_gene = shift @_;
    # get all table names
    my @tables = PifCosm_support_subs::get_gene_tables($dbh);
    #my $sth = $dbh->prepare("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name" ) or die "Couldn't prepare statement: " . $dbh->errstr;
    #$sth->execute();
    my $tableflag = 'n'; # assume the table is not present
    #while (my $temp = $sth->fetchrow_array() ) { 
    foreach(@tables) {
        if ($lc_gene eq $_) { # if table present
            $tableflag = 'y';    # note it
            last;                # and stop looking
        }
    }
    # if table not already present create it
    if ( $tableflag eq 'n' ) {
        my $update = $dbh->do("CREATE TABLE $lc_gene (accno TEXT PRIMARY KEY, sequence TEXT DEFAULT 'empty', cluster TEXT DEFAULT 'empty')") or die "Could not create table: " . $dbh->errstr;
        print "Creating table $lc_gene.\n"; # let the world know
    }
}

#Divides the sequence and output the parts to apropriate fasta files
sub sort_and_output_fasta {
    my $sequencename = shift @_;     #The first item should be the sequence name
    my $sequence = shift @_;         #The second item should be the sequence
    my $dbh = shift @_;              #The thired item should be the database handeler
    my %intergenes;                  #hash defining genes that are between other genes
    while ($_[0] ne '|||') {         #the hash is separated from the next by |||
        my $key = shift @_;
        $intergenes{$key} = shift @_; # get the values
    }
    if ($_[0] eq '|||') { shift @_; } #remove separator
    my %hmmresults=@_;               #the rest is the hammer results in a hash
    my @keys = keys %hmmresults;
    my $printflag = "not printed";   # flag to see if the sequence had a significant match to any gene
    my %gene;                        # hash to store start, middle, and end position for each gene
    my %e_value_gene;                # hash to store the e-value for each gene
    #Go through each matched model, if two hits to the same position in the same gene set the one with highest e-value to "empty"
    for (my $i=0; $i < (scalar @keys)-1; ++$i) {
        if (!$keys[$i] or !$hmmresults{$keys[$i]} or $hmmresults{$keys[$i]} eq "empty" or $keys[$i] eq 'empty') { next; }
        my $gene_i;
        my $position_i;
        if ($keys[$i] =~ /^([A-Za-z0-9]+)[_\.]([A-Za-z0-9]+)[_\.]([A-Za-z0-9]+)/) { # get the info about the HMM that was matched from the hash keys
            $gene_i = $1; # the first value is the gene
            $position_i = $2; #get the gene and the position
            $gene{$1} = "0;0;0";           #make sure we have values for start middle and end of gene and null all positions
        }
        else { next; }
        #Go throug subsequent models
        for (my $j=$i+1; $j < scalar @keys; ++$j) {
            if (!$keys[$i] or !$hmmresults{$keys[$i]} or $hmmresults{$keys[$j]} eq "empty" or $keys[$i] eq 'empty') { next; }
            if ($keys[$j] =~ /^([A-Za-z0-9]+)[_\.]([A-Za-z0-9]+)[_\.]([A-Za-z0-9]+)/) { # do the same for next gene/key
                my $gene_j = $1;
                my $position_j = $2;
                $gene{$1} = "0;0;0";
                if ( $gene_i eq $gene_j and ( $position_i eq $position_j or ($position_i =~ /(whole)|(full)/ or $position_j =~ /(whole)|(full)/) ) ) { #if the same gene and position in gene
                    my @temp_i = split /;/, $hmmresults{$keys[$i]}; #get the e-value for the first model
                    my @temp_j = split /;/, $hmmresults{$keys[$j]}; #get the e-value for the second model
                    # ignore worse matches for the same position
                    if (($temp_j[2] and !$temp_i[2]) or ($temp_j[2] and $temp_i[2] and $temp_j[2] < $temp_i[2])) {
                        $hmmresults{$keys[$i]} = "empty";           #if the first model have higher e-value than the second it is set as "empty"
                    }
                    else { $hmmresults{$keys[$j]} = "empty"; }      #othervise we set the second to "empty"
                }
            }
            else { next; }
        }
    }
    #get the start, midle, and end position for each gene
    foreach (@keys) { #keys %hmmresults) {
        my $key = $_; # get the key
        if ( !$key or !$hmmresults{$key} or $hmmresults{$key} eq "empty" or $key eq 'empty' ) { next; } # if entry empty proceed
        if ($key =~ /^([A-Za-z0-9]+)[_\.]([A-Za-z0-9]+)[_\.]([A-Za-z0-9]+)/) { # match gene and position and separate identifier
            my $gene=$1;     # get gene
            my $position=$2; # get position
            if (!$gene{$gene}) { $gene{$gene} = "0;0;0"; } # this should not hapen but apers to do happen
            my @temp = split /;/, $gene{$gene}; # get the positions from previous matches to the gene, we should just have one match to each region (se above)
            if (!$gene{$gene}) { print "$gene\n"; }
            my @basepairposition = split /;/, $hmmresults{$key}; # get basepair position
            if ($position =~ /(start)|(beginning)|(whole)|(full)/ and $basepairposition[0] and $basepairposition[0] =~ /[0-9]/) { $temp[0] = $basepairposition[0]; } # if position include beginig of gene get start position
            if ($position =~ /(middle)|(between)/ and $basepairposition[0] and $basepairposition[0] =~ /[0-9]/) { $temp[1] = $basepairposition[0]; } #if position is in the middle get middle position
            if ($position =~ /(stop)|(end)|(whole)|(full)/ and $basepairposition[1] and $basepairposition[0] =~ /[0-9]/) { $temp[2] = $basepairposition[1]; } # if position include end of gene get end position
            if (!$e_value_gene{$gene}) { $e_value_gene{$gene} = $basepairposition[2]; } #get e-value for the gene
            elsif ($e_value_gene{$gene} > $basepairposition[2]) { $e_value_gene{$gene} = $basepairposition[2]; } # if new e-value lower than previous get it
            $gene{$gene} = join ";", @temp; # save the new matches
        }
    }

    #Check regions that are defined to be between genes
    foreach (keys %intergenes) { # for each gene defined to be between two other genes
        my $key = $_; # get intergene
        my @temp = split /;/, $intergenes{$key}; # get the genes it is between
        my $start = 0; # variable for start position
        my $end = 0; # variable for end position
        my ($trash1,$trash2); # variables for values that will not be used

        if ( $gene{$temp[0]} ) { # if we have a match for the leading gene
            #if ($gene{$temp[0]} ne undef) { # and values for that gene
                ($trash1,$trash2,$start) = split /;/,$gene{$temp[0]}; # get the end position of that gene
                if ($start > 0) { # if a real end position
                    ++$start; # set start of intergene to one more than the end
                    if (!$e_value_gene{$key}) { $e_value_gene{$key} = $e_value_gene{$temp[0]}; } # set e-value for intergene
                    elsif ($e_value_gene{$key} > $e_value_gene{$temp[0]}) { $e_value_gene{$key} = $e_value_gene{$temp[0]}; } # if better evalue is available get it
                }
            #}
        }
        if ( $gene{$temp[1]} ) { # if we have the trailing gene
            #if ($gene{$temp[1]} ne undef) { # and it is defined
                ($end,$trash1,$trash2) = split /;/,$gene{$temp[1]}; # get its start position
                if ($end > 0) { # if start position is a positive value
                   --$end;      # put end of intergene to one before start of trailing gene
                   if (!$e_value_gene{$key}) { $e_value_gene{$key} = $e_value_gene{$temp[0]}; } # get e-value
                   elsif ($e_value_gene{$key} > $e_value_gene{$temp[0]}) { $e_value_gene{$key} = $e_value_gene{$temp[0]}; } # if better e-value available get it
                }
            #}
        }
        if ($start > 0 or $end > 0) {
            $gene{$key} = "$start;0;$end"; # if we have start or end of intergene save it
        }
    }
    #go through each gene and save in database
    my $update = 0;
    foreach(keys %gene) {
        my $key = $_; # get key
        if ( !$key or !$gene{$key} or !$gene{$key} or $key eq "empty" or $gene{$key} eq "empty" ) { next; } # if no value skip it
        my ($start, $middle, $end) = split /;/, $gene{$key}; # get the start, middle, and end positions
        #if (!$start) { print "$gene{$key}\n"; }
        if ($start > 0 and $end > 0 and $end > $start) { # if end after start include gene
            my $lc_gene = lc($_); # use lower case for uniformity
            &test_and_create_table ( $dbh, $lc_gene ); # se if there is a table for the gene, else make one
            my $sub_sequence=substr($sequence,$start,$end-$start); # pars out gene sequence
            $update += $dbh->do("INSERT INTO $lc_gene (accno,sequence) VALUES ('$sequencename','$sub_sequence')") or die "Could not insert: " . $dbh->errstr; # save it in database
            $printflag = "printed"; # flag the sequence as printed
        }
        elsif ($start > 0 and $end == 0 and $start < length($sequence)) { # if start but no end, and start before end of sequence
            #print "Start found at: $start.\n";
            my $lc_gene = lc($_);
            &test_and_create_table ( $dbh, $lc_gene );
            my $sub_sequence=substr($sequence,$start);
            $update += $dbh->do("INSERT INTO $lc_gene (accno,sequence) VALUES ('$sequencename','$sub_sequence')") or die "Could not insert: " . $dbh->errstr;
            $printflag = "printed";
        }
        elsif ($end > 0 and $start == 0 ) { # if end but no start
            #print "End found at: $end.\n";
            my $lc_gene = lc($_);
            &test_and_create_table ( $dbh, $lc_gene );
            my $sub_sequence=substr($sequence,0,$end);
            $update += $dbh->do("INSERT INTO $lc_gene (accno,sequence) VALUES ('$sequencename','$sub_sequence')") or die "Could not insert: " . $dbh->errstr;
            $printflag = "printed";
        }
        elsif ($middle != 0 and $end == 0 and $start == 0 ) { # if only middle, no need to parse sequence
            #print "Recognized a part of the gene starting at: $middle.\n";
            my $lc_gene = lc($_);
            &test_and_create_table ( $dbh, $lc_gene );
            $update += $dbh->do("INSERT INTO $lc_gene (accno,sequence) VALUES ('$sequencename','$sequence')") or die "Could not insert: " . $dbh->errstr;
            $printflag = "printed";
        }
    }

    if ( $printflag eq "not printed" ) { #if no model was matched print it to garbage file
            return (0,$update);
    }
    else { return (1,$update); }

}

sub runhammer { # takes "sequence name (not starting with >", sequence, and hmm model database as arguments
    my $hmmdatabase = pop @_;        # hmm model database, should be last in argument string
    my $path = shift @_;             # path to hmmscan
    chomp(my $query = shift @_);     # sequence name, no > and no white space
    chomp(my $sequence = shift @_);  # sequence
    my $e_value_cut_off = shift @_;  # maximum e-value to report
    my %start_end; # hash to store start and end position (on the sequence) for the match to the different models

    # write file with query sequence
    open QUERYFILE, ">query.fst" or die;
    print QUERYFILE ">$query\n$sequence";
    close QUERYFILE or die;

    # Start HMMER search and read output
    my @hmmoutput = `${path}$External_program::hmmscan -E $e_value_cut_off $hmmdatabase query.fst`;
    unlink "query.fst" or die; # remove query sequence file
    if (!($hmmoutput[0] =~ /hmmscan :: search sequence/)) { die "Error in execution of hmmscan\n"; }

    # parse HMMER output
    for (my $i=0; $i < scalar @hmmoutput; ++$i) {
        if ($hmmoutput[$i] =~ /\[No hits detected that satisfy reporting thresholds\]/) {
            # print "returning (\"empty\",\"empty\")\n";
            return ("empty","empty");
        }
        #Get the starting and end positions
        elsif ($hmmoutput[$i] =~ /^>>\s+(\S+)/) { # Found hmm name
            my $temp = $1; # get hmm name
            my $evalue = 1; # set e-value to improbable value
            for ($i += 3; $i < scalar @hmmoutput; ++$i) { # continue three rows down
                if ($hmmoutput[$i] =~ /^\s+Alignments/) { last; } # stop parsing of entery when hitting alignments
                if (!$hmmoutput[$i] or $hmmoutput[$i] =~ /^\s*$/) { next; } # if blank row continue to next
                my @string = split /\s+/, $hmmoutput[$i]; # othervise it should be data about the hit
                if ($string[6] < $evalue) { # if best hit
                    $evalue = $string[6]; # save e-value
                    my $start = $string[10]-$string[7]+1; # start should be at start of hmm not start of hit
                    if ($start < 1) { $start = 1; } # deletitions may mess this up
                    $start_end{$temp} = $start . ";" . $string[11] . ";" . $evalue; # get start, end, and e-value for the hit
                }
            }
        }
    }
    return %start_end; # return data for each hit
}
1;
