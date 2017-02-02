package PifCosm_support_subs;
use strict;

### sort subroutine
#sub by_number { $a <=> $b }
### Subroutines that many other subroutinse are dependent on ###
sub get_tables {
    my $dbh = shift @_;
    my $sth = $dbh->prepare("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name" ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    my @tables;
    while (my $table = $sth->fetchrow_array()) {
        push (@tables, $table);
    }
    $sth->finish();
    return @tables
}
sub get_gene_tables {
    my $dbh = shift @_;
    my $sth = $dbh->prepare("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name" ) or die "Couldn't prepare statement: " . $dbh->errstr;
    $sth->execute();
    my @tables;
    while (my $table = $sth->fetchrow_array()) {
        if ($table ne "alignments" and $table ne "gb_data" and $table ne "alignment_groups") { push (@tables, $table); }
    }
    $sth->finish();
    return @tables
}
sub table_present {
    my $dbh = shift @_;
    my $table = shift @_;
    my @candidates = &get_tables( $dbh );
    foreach (@candidates) {
        if ($_ eq $table) { return 'y'; }
    }
    return 'n';
}

sub exist_in_array {
    my $string = shift @_;
    foreach (@_) {
	if ($string eq $_) { return 'y' }
    }
    return 'n';
}

sub connect_to_database {
    my $database = shift @_;
    return (DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr);
}
sub column_present {
    my $dbh = shift @_;
    my $table = shift @_;
    my $column = shift @_;
    my $sth = $dbh->prepare("PRAGMA TABLE_INFO($table)"); # SQLite specific?
    $sth->execute();
    my $indv_column_flag = 'n'; # flag if individual column present
    while (my $ref = $sth->fetchrow_hashref()) {
         if ($ref->{name} eq $column) { $indv_column_flag = 'y'; }
    }
    $sth->finish();
    return $indv_column_flag;
}
sub get_columns {
    my $dbh = shift @_;
    my $table = shift @_;
    my @return;
    my $sth = $dbh->prepare("PRAGMA TABLE_INFO($table)"); # SQLite specific?
    $sth->execute();
    while (my $ref = $sth->fetchrow_hashref()) {
        push (@return,$ref->{name});
    }
    $sth->finish();
    return @return;
}
sub date {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my @week_day = ('Sun','Mon','Tue','Wed','Thu','Fri','Sat');
    my @month = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
    $year += 1900; # get year after Christ
    if ($min < 10) { $min = "0$min"; }
    if ($sec < 10) { $sec = "0$sec"; }
    return "$week_day[$wday] $month[$mon] $mday $hour:$min:$sec $year";
}
sub get_cut_off_string {
    my %cut_offs = split /,/, shift @_;
    my @genes = @_;
    my $return_string = '';
    if ($cut_offs{'all'}) {
        foreach my $gene (@genes) {
            $return_string .= "$gene,$cut_offs{'all'},";
        }
    }
    else {
        foreach my $gene (@genes) {
            if ($cut_offs{$gene}) {
                $return_string .= "$gene,$cut_offs{$gene},";
            }
            elsif ($cut_offs{'default'}) {
                $return_string .= "$gene,$cut_offs{'default'},";
            }
            else { die "Could not find cut-off value for $gene!!!\n"; }
        }
    }
    $return_string =~ s/,$//;
    return $return_string;
}
my %hash;
sub by_value { $hash{$b} <=> $hash{$a}; }
sub keys_in_order {
    %hash = @_;
    return sort by_value keys %hash;
}

sub check_taxon_string_consistency {
    my $database = shift @_;
    my $non_taxa = shift @_;
    my $output = shift @_;
    my $dbh = &connect_to_database( $database );
    my $sth = $dbh->prepare("SELECT DISTINCT taxon_string FROM gb_data") or die;
    $sth->execute();
    my %taxon_hash;
    while (my $taxon_string = $sth->fetchrow_array()) {
        my @taxon_array = split /; /, $taxon_string;
	for (my $i=0; $i<scalar @taxon_array; ++$i) {
            my $parent;
            if ($i == 0) { $parent = 'o'; }
            else { $parent = $taxon_array[$i-1]; }
	    if (!$taxon_hash{$taxon_array[$i]}) { $taxon_hash{$taxon_array[$i]} = [ {$parent,1}, [] ]; }
	    else { $taxon_hash{$taxon_array[$i]}->[0]->{$parent} += 1; }
        }
    }
    $sth->finish();
    if ($output eq 'p') { print "Conflict in classification of taxa\n" }
    my %uppdate_hash;
    foreach my $taxon (keys %taxon_hash) {
	@{$taxon_hash{$taxon}->[1]} = keys_in_order(%{$taxon_hash{$taxon}->[0]});
	if ( $taxon_hash{$taxon}->[1]->[0] eq $taxon ) {
	    shift @{$taxon_hash{$taxon}->[1]};
	    if (scalar @{$taxon_hash{$taxon}->[1]} < 1) { $taxon_hash{$taxon}->[1]->[0] = 'o'; }
	}
        if (scalar @{$taxon_hash{$taxon}->[1]} > 1) {
	    if ($output eq 'p') {
                print "Conflict in classification of $taxon could be either ";
		for (my $j=0; $j< scalar @{$taxon_hash{$taxon}->[1]}; ++$j) {
                    print $taxon_hash{$taxon}->[1]->[$j];
		    if ($j == scalar @{$taxon_hash{$taxon}->[1]} -2) { print ", or "; }
		    elsif ($j < scalar @{$taxon_hash{$taxon}->[1]} -1) { print ", "; }
                    
                }
                print "; ordered by number of occurences\n";
            }
        }
    }
    if ($output eq 'r') {
	print "Making new taxon_strings\n";
	$sth = $dbh->prepare("SELECT DISTINCT taxon_string FROM gb_data") or die;
	$sth->execute();
    	while (my $taxon_string = $sth->fetchrow_array()) {
	    my @taxon_array = split /; /, $taxon_string;
	    my @ignor_taxa = split /;/, $non_taxa;
	    my $new_string;
	    my $taxon = $taxon_array[-1];
            my $i=-1;
	    while ($taxon and $taxon ne 'o' and &exist_in_array($taxon,@ignor_taxa) eq 'y') {
		print "Ignoring $taxon in $taxon_string\n";
		if (!$new_string) { $new_string = $taxon; }
                else { $new_string = $taxon . "; " . $new_string; }
		$taxon = $taxon_array[--$i];
	    }
	    while ($taxon and $taxon ne 'o') {
		#print $taxon , ' ' , $taxon_array[-1] , "\n";
		if (!$new_string) { $new_string = $taxon; }
		else { $new_string = $taxon . "; " . $new_string; }
		$taxon = $taxon_hash{$taxon}->[1]->[0];
	    }
	    if ($new_string and $new_string ne $taxon_string) {
		$new_string =~ s/'/''/g;
		$taxon_string =~ s/'/''/g;
		my $update = $dbh->do("UPDATE gb_data SET taxon_string = '$new_string' WHERE taxon_string = '$taxon_string'");
		if ($update and $update > 0) { print "$new_string  -- replaces -- $taxon_string\n"; }
		else { print STDERR "WARNING!!! Could not replace '$taxon_string' with '$new_string'.\n"; }
	    }
	}
    }
    $dbh->disconnect();
}

sub prepare_search_for_taxon_in_taxon_string {
    my $dbh = shift @_;
    my $taxon = shift @_;
    my $columns = shift @_;
    my $table = shift @_;
    $taxon =~ s/'/''/g;
    return $dbh->prepare("SELECT $columns FROM $table WHERE (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon')");# or die;
}

################################################################
sub run_raxml { # sub to run RAxML
    my $phylipfile = shift @_; # sequence file
    my $path  = shift @_; # path to treebuilding programms
    my $run_name = shift @_; # name to give the run
    my $n_threads = shift @_; # number of threads if pthreads, =0 if serial
    my @methods = split /,/, shift @_; # the methods to use
    my $partition_file = shift @_;
    my $bootstraps = shift @_; # undef or zero if no bootstrap, else the number of bootstrap replicates
    my $ml_score;
    #my @return_array;
    my $ml_tree = 'empty';
    if (!$methods[0]) {
        print STDERR "No method for phylogenetic analysis given. No tree produced.\n";
        return 'empty',$ml_score,'empty';
    }
    my $additional = '';
    if ( $methods[0] eq 'RAxML' and $partition_file ne 'n') { $additional .= " -q $partition_file"; }
    my @raxml;
    unlink glob "RAxML*.$run_name*"; # remove previous files with the same name
    if ($methods[0] eq 'ExaML') { # if runing ExaML
	if ($n_threads > 0) {
	    @raxml=`${path}$External_program::raxml -y -s $phylipfile -n $run_name -m GTRGAMMA -p 12345`; # make starting tree if runing EXaml
	    #foreach (glob "*") { print "$_\n"; }
	    if ( -e "RAxML_parsimonyTree.$run_name" ) { # if parsimony tree exists
		rename "RAxML_parsimonyTree.$run_name","RAxML_parsimonyTree.$run_name.0";
		unlink glob "RAxML*.$run_name"; # remove all files with the same run name
		my @pars_out;
		if ($partition_file ne 'n') { @pars_out = `$External_program::parsexaml -s $phylipfile -m DNA -n $phylipfile -q $partition_file`; }
		else { @pars_out = `$External_program::parsexaml -s $phylipfile -m DNA -n $phylipfile`; }
		@raxml=`$External_program::openmpi $n_threads ${path}$External_program::examl -s $phylipfile.binary -m GAMMA -n $run_name -t RAxML_parsimonyTree.$run_name.0`; # run RAxML
		foreach (@raxml) {
		    if ($_ =~ /Likelihood\s+:\s+(-{0,1}[0-9.]+)/) {
			$ml_score = $1;
		    }
		}
		if ( -e "ExaML_result.$run_name" ) { # if a tree file is found
		    #unlink "RAxML_parsimonyTree.$run_name.0"; # delete parsimony tree
		    #unlink glob "ExaML_binaryCheckpoint.$run_name\_*"; # remove checkpoint files
		    #unlink "ExaML_info.$run_name"; # remove info file
		    #unlink "ExaML_modelFile.$run_name"; # remove log file
		    #unlink "$phylipfile.reduced";
		    rename "ExaML_result.$run_name", "XXX_ML.$run_name.tree"; # Change name so it is easy to delete all generated files
		    $ml_tree = "XXX_ML.$run_name.tree";
		    #unlink (glob "RAxML_*.$run_name*");
		    #unlink (glob "ExaML_*.$run_name*");
		    #return "RAxML_result.$run_name",$ml_score; # return tree file
		}
		elsif ($methods[1] and $methods[1] eq 'RAx-br') { # if ExaML fail to produce a tree
		    print STDERR "No ExaML tree. Doing ML on parsimony topology.\n";
		    unlink glob "RAxML*.$run_name"; # remove files with same run name
		    undef @raxml;
		    @raxml=`${path}$External_program::raxml -f e -s $phylipfile -m GTRGAMMA -n $run_name -t RAxML_parsimonyTree.$run_name.0$additional`; # optimize branchlengths on parsimony tree
		    foreach (@raxml) {
			if ($_ =~ /Final GAMMA  likelihood: (-{0,1}[0-9\.]+)/) {
			    $ml_score = $1;
			}
		    }
		    rename "RAxML_result.$run_name","XXX_ML.$run_name.tree";
		    $ml_tree = "XXX_ML.$run_name.tree";
		#    unlink glob "RAxML_*.$run_name*"; # remove file for model parameters
		#    unlink "$phylipfile.reduced";
		}
		#if ( $ml_tree eq 'empty' && -e "RAxML_result.$run_name" ) { # if tree file produced
		#    rename "RAxML_result.$run_name", "XXX_ML.$run_name.tree";
                #    $ml_tree = "XXX_ML.$run_name.tree";
		    #unlink "RAxML_parsimonyTree.$run_name.0"; # remove parsimony tree file
		    #unlink "RAxML_info.$run_name"; # remove info file
		    #unlink "RAxML_log.$run_name"; # remove log file
		#    unlink "RAxML_*.$run_name*";
		#    unlink "$phylipfile.reduced";
		    #return "RAxML_result.$run_name",$ml_score; # return tree file
		#}
		#else
	       if ($ml_tree eq 'empty') { # if RAxML failed totaly return parsimony tree
		    print STDERR "No ML tree, returning parsimony tree.\n";
		    foreach (@raxml) { print $_; }
		    #return "RAxML_parsimonyTree.$run_name.0",$ml_score;
		    rename "RAxML_parsimonyTree.$run_name.0","XXXparsimony.$run_name.tree";
		    $ml_tree = "XXXparsimony.$run_name.tree";
		    #unlink "RAxML_*.$run_name*";
		}
		unlink "$phylipfile.reduced";
		unlink (glob "RAxML_*.$run_name*");
	    	unlink (glob "ExaML_*.$run_name*");
	    }
	    else {
		print STDERR "No starting tree generated for ExaML.\n";
		foreach (@raxml) { print $_; }
		#return 'empty',$ml_score; # if no parsimony tree return empty
	    }
	}
	else {
	    print STDERR "Examl requires multiple cores.\n";
	    #return 'empty',$ml_score; # if numbers of threads not given return empty
	}
    }
    elsif ($methods[0] eq 'RAxML') {
        my @raxml=`${path}$External_program::raxml -f d -s $phylipfile -m GTRGAMMA -p 54321 -n $run_name$additional`;
        foreach (@raxml) {
            if ($_ =~ /Final GAMMA-based Score of best tree (-{0,1}[0-9\.]+)/) {
            #if ($_ =~ /Likelihood\s+:\s+(-{0,1}[0-9.]+)/) {
                $ml_score = $1;
            }
        }
        if ( -e "RAxML_result.$run_name" ) { # if tree file produced
	    rename "RAxML_result.$run_name", "XXX_ML.$run_name.tree";
	    $ml_tree = "XXX_ML.$run_name.tree";
	    unlink "RAxML_*.$run_name*";
            #unlink "RAxML_info.$run_name"; # remove info file
            #unlink "RAxML_log.$run_name"; # remove log file
            unlink "$phylipfile.reduced";
            #return "RAxML_result.$run_name",$ml_score; # return tree file
        }
        else {
            print STDERR "No tree generated.\n";
            foreach (@raxml) { print $_; }
            #return 'empty',$ml_score;
        }
    }
    elsif ($methods[0] eq 'fasttree') {
        system "${path}$External_program::fasttree -quiet -nt -gtr -gamma -log XXXlogfile.txt < $phylipfile > XXX_ML.$run_name.tree";
        if (-e "XXXlogfile.txt") {
            open LOGFILE, "<XXXlogfile.txt" or die "Could not open logfile.\n";
            while (my $row=<LOGFILE>) {
                if ($row =~ /Gamma20LogLk\s+(-{0,1}[0-9\.]+)/) {
                    $ml_score = $1;
                }
            }
            close LOGFILE or die;
            unlink "XXXlogfile.txt";
        }
        if (-e "XXX_ML.$run_name.tree") { $ml_tree = "XXX_ML.$run_name.tree"; }
        else {
            print STDERR "No tree generated.\n";
            #return 'empty',$ml_score;
        }
    }
    else {
        print STDERR "No recognized method for the phylogenetic analysis given. No tree produced.\n";
        return 'empty',$ml_score,'empty';
    }
############################
# save ML tree in $ml_tree
    my $boottree_file = 'empty';
    if ($bootstraps) {
	print "Performing bootstrap analysis.\n";
	if ($methods[0] eq 'RAxML' or $methods[0] eq 'ExaML') {
	    unlink glob "RAxML*.XXXtemp";
	    if ($partition_file eq 'n') {
		if ($n_threads > 1) {
		    my @raxml = `${path}$External_program::raxmlPTHREADS -x 12345 -p 54321 -# $bootstraps -m GTRGAMMA -s $phylipfile -n XXXtemp -T $n_threads`;
		    if ($ml_tree ne 'empty' and -e $ml_tree) { @raxml = `${path}$External_program::raxmlPTHREADS -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart`; }
		    else { print STDERR "No ML tree for $run_name to draw support values on for $run_name.\n"; }
		}
		else {
		    my @raxml = `${path}$External_program::raxml -x 12345 -p 54321 -# $bootstraps -m GTRGAMMA -s $phylipfile -n XXXtemp`;
		    if ($ml_tree ne 'empty' and -e $ml_tree) { @raxml = `${path}$External_program::raxml -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart`; }
		    else { print STDERR "No ML tree for $run_name to draw support values on for $run_name.\n"; }
		}
	    }
	    else {
		if ($n_threads > 1) {
		    my @raxml = `${path}$External_program::raxmlPTHREADS -x 12345 -p 54321 -# $bootstraps -m GTRGAMMA -s $phylipfile -q $partition_file -n XXXtemp -T $n_threads`;
		    if ($ml_tree ne 'empty' and -e $ml_tree) { @raxml = `${path}$External_program::raxmlPTHREADS -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart -T $n_threads`; }
		    else { print STDERR "No ML tree for $run_name to draw support values on for $run_name.\n"; }
		}
		else {
		    my @raxml = `${path}$External_program::raxml -x 12345 -p 54321 -# $bootstraps -m GTRGAMMA -s $phylipfile -q $partition_file -n XXXtemp`;
		    if ($ml_tree ne 'empty' and -e $ml_tree) { @raxml = `${path}$External_program::raxml -f b -t $ml_tree -z RAxML_bootstrap.XXXtemp -m GTRGAMMA -n XXXbipart`; }
		    else { print STDERR "No ML tree for $run_name to draw support values on for $run_name.\n"; }
		}
	    }
	    if (-e "RAxML_bipartitions.XXXbipart") {
		if (!$ml_tree) { $ml_tree = 'XXX_ML.tree'; }
		rename "RAxML_bipartitions.XXXbipart",$ml_tree;
	    }
	    unlink "$phylipfile.reduced";
	    if (-e 'RAxML_bootstrap.XXXtemp') {
		rename 'RAxML_bootstrap.XXXtemp', "XXX_boot.$run_name.trees";
		$boottree_file = "XXX_boot.$run_name.trees";
	    }
	    else {
		print STDERR "Failed to do bootstrap in RAxML. Continuing reluctantly...\n";
	    }
	    unlink glob "RAxML*.XXXtemp";
	    unlink glob "RAxML*.XXXbipart"
	}
	elsif ($methods[0] eq 'fasttree') {
	    system "${path}$External_program::fasttree -quiet -nt -gtr -gamma -noml -boot $bootstraps -intree $ml_tree < $phylipfile > XXXtemp.tree";
	    if (-e 'XXXtemp.tree') {
		rename "XXXtemp.tree","XXX_ML.$run_name.tree";
		$ml_tree = "XXX_ML.$run_name.tree";
	    }
	    else {
		print STDERR "Failed to do bootstrap in FastTree. Continuing reluctantly...\n";
	    }
	}
	else { print STDERR "Did not recognize method for bootstrap analysis. No bootstrap support generated.\n"; }
    }
#############################
    return $ml_tree,$ml_score,$boottree_file;
}

# Print a relaxed phylip file (argument 1) from database. Need a prepared statement (argument 2) with accno in the first column and
# sequence in the second column returned, the number of taxa that will be printed and the length of the alignment. Returns the
# difference between the printed taxa and the number of taxa that should be printed.
sub print_phylip {
    my $filename=shift @_; # file to write to
    my $sth = shift @_; # search handler with accno and sequence
    my $n_taxa = shift @_; # number of taxa that will be written
    my $seq_length = shift @_; # sequence length that will be writen
    my $printed_taxa=0; # number of taxa that has been printed to file
    my $row_length = 5000; # the number of bases to have on each row
    $sth->execute();
    open PHYLIPFILE, ">$filename" or die "Could not open $filename: $!.\n"; # open file handler
    print PHYLIPFILE "$n_taxa $seq_length\n"; # print number of taxa and alignment length
    unlink glob "XXXtemp_sequences_*.txt";
    while ( my @row = $sth->fetchrow_array()) { # for each sequence
        ++$printed_taxa; # count number of taxa
        if ($row[1] eq 'empty') { $row[1] = '-' x $seq_length; }
        print PHYLIPFILE "$row[0] " . substr($row[1],0,$row_length) . "\n"; # print accno and sequence
        my $n_printed = $row_length;
        while ($seq_length-$n_printed > 0) {
            open RESIDUAL, ">>XXXtemp_sequences_$n_printed.txt" or die "Could not open temporary file XXXtemp_sequences_$n_printed.txt: $!.\n";
            print RESIDUAL substr($row[1],$n_printed,$row_length) . "\n";
            $n_printed += $row_length;
            close RESIDUAL or die;
        }
    }
    $sth->finish();
    close PHYLIPFILE or die;
    my $mearged = $row_length;

    foreach (glob "XXXtemp_sequences_*.txt") {
        open PHYLIPFILE, ">>$filename" or die "Could not open $filename: $!.\n"; # open file handler to append text
        open RESIDUAL, "<$_" or die "Could not open temporary file $_: $!.\n";
        while (my $row = <RESIDUAL>) {
            print PHYLIPFILE $row;
        }
        close RESIDUAL or die;
        close PHYLIPFILE or die;
    }
    unlink glob "XXXtemp_sequences_*.txt";
    return $printed_taxa-$n_taxa; # return the difference between how many sequences were printed and how many should be printed
}

# Print a fastafile (argument 1) from database. Need a prepared statement (argument 2) with accno in the first column and
# sequence in the second column returned. Returns number of printed taxa
sub print_fasta {
    my $filename=shift @_; # file to write to
    my $sth = shift @_; # search handler with accno and sequences
    my $n_taxa = 0; # number of taxa printed
    $sth->execute();
    open FASTAFILE, ">$filename" or die "Could not open $filename: $!.\n"; # oprn file handler
    while ( my @row = $sth->fetchrow_array()) { # for each sequence
        print FASTAFILE ">$row[0]\n$row[1]\n"; # print accno and sequence
        ++$n_taxa; # count number of taxa printed
    }
    $sth->finish();
    close FASTAFILE or die;
    return $n_taxa; # return how many taxa were printed
}

#sub add_to_mafft_mearge_files {
#    my $fastafile = shift;
#    my $mearge_file = shift;
#    my $previous_seqs = shift;
#    my $sth = shift;
#    my $n_taxa = 0;
#    $sth->execute();
#    open FASTAFILE, ">>$fastafile" or die "Could not open $fastafile: $!.\n"; # oprn file handler
#    while ( my @row = $sth->fetchrow_array()) { # for each sequence
#        print FASTAFILE ">$row[0]\n$row[1]\n"; # print accno and sequence
#        ++$n_taxa; # count number of taxa printed
#    }
#    $sth->finish();
#    close FASTAFILE or die;
#    if ($n_taxa > 1) {
#	open SEQUENCES, ">>", $mearge_file or die "Could not open $mearge_file: $!.\n";
#	for (my $i=1; $i <= $n_taxa; ++$i) {
#	    if ($i == 1) { print SEQUENCES $i+$previous_seqs; }
#	    else { print SEQUENCES ' ', $i+$previous_seqs; }
#	}
#	print SEQUENCES "\n";
#	close SEQUENCES;
#    }
#    return $n_taxa+$previous_seqs;
#}

# Put sequences in fasta file (argument 1) into the alignments table in a database (connected database string is argument 2)
# for given gene (argument 3). Return number of taxa and length of last sequence
sub read_fasta_to_alignments_table {
    my $alignment_file = shift @_; # file to read
    my $dbh = shift @_; # database handler
    my $gene = shift @_; # gene (column) to include the sequences under
    my $name_column = shift @_;
    if ($name_column eq 'accno') { $name_column = "$gene\_accno"; }
    my $n_taxa = 0; # number of taxa included
    open ALIGNMENT, "<$alignment_file" or die "Could not open $alignment_file: $!.\n"; # open file handler
    my $sequence; # to store sequence
    my $accno; # to store accno
    while ( my $infile = <ALIGNMENT> ) { # read infile row by row
        chomp($infile); # remove line breaks
        if ( $infile =~ s/^>// ) { # if sequence name remove > and
            ++$n_taxa; # count up the number of taxa read
            if ($accno) { # if a sequence have been read put sequence and accno into database
                &insert_sequence($dbh,$gene,$name_column,$accno,$sequence);
            }
            $accno = $infile; # set the accno to the new sequence name
            undef $sequence; # reset sequence
        }
        else { $sequence .= $infile; } # if not a sequence name it is sequence
    }
    if ($accno) { &insert_sequence($dbh,$gene,$name_column,$accno,$sequence); }
    return ($n_taxa,length($sequence)); # return the number of sequences and the length of the last sequence
}

sub insert_sequence {
    my $dbh = shift @_; # database handler
    my $gene = shift @_; # gene (column) to include the sequences under
    my $name_column = shift @_;
    my $accno = shift @_;
    my $sequence = shift @_;
    if ($sequence and $sequence =~ /[^-]/) {
        $sequence =~ s/\s//g;
        $dbh->do("UPDATE alignments SET $gene\_sequence=\'$sequence\' WHERE $name_column=\'$accno\'") or die "Could not update alignments: " .$dbh->errstr;
    }
    elsif ($accno) {
        my $sth = $dbh->prepare("SELECT $gene\_accno FROM alignments WHERE $name_column=\'$accno\'") or die;
        $sth->execute();
        my $remove_accno = $sth->fetchrow_array();
        if ($remove_accno) {
            &remove_accno_from_alignments($dbh,$gene,$remove_accno);
        }
    }

}

sub comp_distance { # subroutine to compare the distances between pairs of sequences from different genes
    my $dbh = shift @_; # database handler
    my $path = shift @_; # path to pairalign
    my $gene1 = shift @_; # first gene to compare
    my $gene2 = shift @_; # second gene to compare
    my $taxon = shift @_; # taxon to compare for
    my $n_comp = shift @_; # number of sequences to compare
    # get sequences
    my $sth = $dbh->prepare("SELECT $gene1.accno,$gene1.sequence,$gene2.accno,$gene2.sequence FROM alignments INNER JOIN $gene1 ON alignments.$gene1\_accno=$gene1.accno INNER JOIN $gene2 ON $gene2.accno=alignments.$gene2\_accno INNER JOIN gb_data ON alignments.$gene2\_accno=gb_data.accno WHERE $gene1\_sequence!='empty' AND (gb_data.taxon_string LIKE '%; $taxon;%' OR gb_data.taxon_string LIKE '$taxon;%' OR gb_data.taxon_string LIKE '%; $taxon' OR gb_data.taxon_string='$taxon') ORDER BY RANDOM() LIMIT $n_comp") or die "Could not prepare sequence query: " . $dbh->errstr;
    $sth->execute();
    my $fasta1=''; # sequences of first gene in fasta format
    my $fasta2=''; # sequences of second gene in fasta format
    my $n=0; # number of comparisons made
    while (my @row= $sth->fetchrow_array()) { # get sequences
        $fasta1 .= ">$row[0]\n$row[1]\n";
        $fasta2 .= ">$row[2]\n$row[3]\n";
    }
    # write sequences for first gene
    open FASTAFILE, ">XXXtemp_pair.fst" or die "Could not open XXXtemp_pair.fst: $!.\n";
    print FASTAFILE $fasta1;
    close FASTAFILE or die;
    # get distances between sequences for first gene
    my @distances1 = `${path}$External_program::pairalign -j < XXXtemp_pair.fst`;
    # write sequences for second gene
    open FASTAFILE, ">XXXtemp_pair.fst" or die "Could not open XXXtemp_pair.fst: $!.\n";
    print FASTAFILE $fasta2;
    close FASTAFILE or die;
    # get distances between sequences for second gene
    my @distances2 = `${path}$External_program::pairalign -j < XXXtemp_pair.fst`;
    unlink "XXXtemp_pair.fst"; # delete sequence file
    if (scalar @distances1 != scalar @distances2) { return (undef,undef); } # if diferent number of comparisons within each gene return undef
    else { # if OK
        my $sum=0; # to count sum
        for (my $i=0; $i < scalar @distances1; ++$i) { # for each pairwise difference
            if ($distances1[$i]<=0) { ++$n; next; } # count number of identical sequences, can not divide by 0
            $sum += $distances2[$i]/$distances1[$i]; # otherwise sum up the fraction
        }
        if (((scalar @distances1)-$n)<=0) { return (undef,undef); } # if no comparisons left return undel
        return ($sum/((scalar @distances1)-$n),(scalar @distances1)-$n); # else return the mean fraction and number of comparisons
    }
}

sub remove_accno_from_alignments { # sub to remove accnos from alignments table
    my $dbh = shift @_; # get database
    my $gene = shift @_; # get gene to now which column to look in
    my $accno = shift @_; # get accno to remove
    my $changes = 0; # to save number of changes made
    my $string_accno = $dbh->quote($accno);
    my $sth = $dbh->prepare("SELECT taxon_name FROM alignments WHERE ${gene}_accno=$string_accno"); # get taxon name
    $sth->execute();
    my $taxon = $sth->fetchrow_array(); # get taxon name for sequence
    $sth->finish();
    my @columns = &get_gene_tables($dbh); # get all the gene columns
    my $column_string = ''; # initiate search string for columns
    foreach(@columns) { if ($_ ne $gene) { $column_string .= "${_}_accno,"; } } # add each column to search string
    $column_string =~ s/,$//; # remove trailing ,
    # get accno data for each column of the taxon
    my $string_taxon = $dbh->quote($taxon);
    $sth = $dbh->prepare("SELECT $column_string FROM alignments WHERE taxon_name=$string_taxon");
    $sth->execute();
    my @temp = $sth->fetchrow_array(); # save the output
    $sth->finish();
    my $delete_taxa = 'y'; # assume there is no other accnos for the taxon
    foreach (@temp) { if ($_ ne 'empty') { $delete_taxa = 'n'; last; } } # if there is accno for another gene, flag it
    if ($delete_taxa eq 'y') { # if no other gene for the taxon
        $changes = $dbh->do("DELETE FROM alignments WHERE taxon_name=$string_taxon"); # delete taxon
    }
    else { # if other gene is available for the accno
        $changes = $dbh->do("UPDATE alignments SET ${gene}_accno='empty',${gene}_sequence='empty' WHERE ${gene}_accno='$accno'"); # set accno and sequence to 'empty' for the accno
    }
    return $changes; # return how many changes have been made
}

#############################################
### Subroutine to print alignment to file ###
#############################################

sub print_alignment {
    my $database = shift @_; # database name
    my $alignment_format = shift @_; # format to print the alignment in phylip, fasta, or nexus.
    my $get_partitions = shift @_; # indicates if a file (or block for nexus) should be printed with the start and end positions for each gene
    my $file_stem = shift @_; # the stem of the output file name
    my $column_name = shift @_; # column to use as name for sequences
    my $interleaved = shift @_; # if the interleaved format should be used
    my $row_length = shift @_; # the number of bases on one row in phylip format
    my $anchor_gene = shift @_; # a gene that need to be present in all taxa
    my @gene_alignments = @_; # the genes to include in the alignment file
    my $alignment_file; # alignment file name
    my $partition_file; # name for partition file
    if ($alignment_format eq 'phylip') { # if phylip format
        $alignment_file = "$file_stem.phy"; # phylip file name
        $partition_file = "$file_stem.partitions.txt"; # partition file name
    }
    elsif ($alignment_format eq 'fasta') { # if fasta format
        $alignment_file = "$file_stem.fst"; # fasta file name
        $partition_file = "$file_stem.partitions.txt"; # partition file name
        if ($interleaved eq 'y') {
            print STDERR "WARNING!!! Cannot print fastaformat interleaved. Will print ordinary fastaformat.\n";
            $interleaved = 'n';
        }
    }
    elsif ($alignment_format eq 'nexus') { # if nexus format
        $alignment_file = "$file_stem.nex"; # nexus file name, partitions will be writen to Mr Bayes block
    }
    print "    Alignment will be printed to $alignment_file.\n";
    my $sth; # string for database queries
    my $sequence_columns=''; # query string for which columns to get
    my @sequence_lengths; # array for the length of each region
    my $length=0; # total sequence length
    my $dbh = DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr; # conect to database
    if ($gene_alignments[0] eq 'all') { # Get all genes if asked
        @gene_alignments = &get_gene_tables($dbh);
    }
    foreach (@gene_alignments) { # Get the length of each gene, assume all sequences are aligned to same length
        $sth = $dbh->prepare("SELECT LENGTH($_\_sequence) FROM alignments WHERE $_\_sequence!='empty'") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        my $temp = $sth->fetchrow_array();
        if (!$temp) { $temp=0; }
        if ($temp < 1) { print "    No sequences for $_. It will therefore not be printed to alignment or partitions.\n"; }
        push (@sequence_lengths, $temp ); # save length of each gene
        $length += $temp; # sum up length of all genes
        $sth->finish();
    }
    foreach (@gene_alignments) {
        $sequence_columns .= ",$_\_sequence"; # add column to query
    }
    my $condition = &sql_where_genes($dbh,$anchor_gene,join(',',@gene_alignments),'AND'); # condition to get sequence in query
    # get number of taxa and length of taxon name
    $sth = $dbh->prepare("SELECT COUNT($column_name),MAX(LENGTH($column_name)) FROM alignments WHERE $condition") or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    my ($n_taxa,$max_length) = $sth->fetchrow_array(); # store the result
    $sth->finish();
    # get taxon name and sequences
    $sth = $dbh->prepare("SELECT $column_name$sequence_columns FROM alignments WHERE $condition") or die "Could not prepare statement: " . $dbh->errstr;
    $sth->execute();
    open FILE, ">$alignment_file" or die "Could not open $alignment_file: $!.\n"; # open alignment file
    # print head of file if not fasta file
    if ($alignment_format eq 'phylip') { print FILE $n_taxa . ' ' . $length . "\n"; } # if phylip print number of sequences and alignmnet length
    elsif ($alignment_format eq 'nexus') {
        print FILE "#NEXUS\nBEGIN DATA;\nDIMENSIONS  NTAX=$n_taxa NCHAR=$length;\nFORMAT DATATYPE=DNA MISSING=N GAP=-";
        print FILE " INTERLEAVE=YES";
        print FILE ";\n\nMATRIX\n"; } # if nexus print header
    if ($interleaved eq 'y') { unlink glob "XXXtemp_alignment_*.txt"; }

    # print sequences
    while (my @row = $sth->fetchrow_array()) { # for each taxon
        if ($interleaved eq 'y') {
            my $n_printed = 0; # number of bases or genes printed
            my $total; # total number of bases or genes
            if ($row_length eq 'genes') { $total = scalar @gene_alignments; }
            else { $total = $length; }
            my $residue = 'empty'; # bases left to print from previous gene
            my $i=0; # counter for which column to print
            while ($residue ne 'empty' or $i < scalar @row) {
                open ALIGNMENT, ">>", "XXXtemp_alignment_$n_printed.txt" or die "Could not open temporary file XXXtemp_alignment_$n_printed.txt: $!.\n";
                if ($i == 0 or $alignment_format eq 'nexus') { # print sequence name if first row for sequence or if nexus
                    print ALIGNMENT "$row[0]" . ' ' x ($max_length - length($row[0]) + 1);
                    if ($i == 0) { ++$i; }
                }
                else { # else add space
                    print ALIGNMENT ' ' x ($max_length+1);
                }
                if ($row_length eq 'genes') { # if printing one gene per row
                    while ($sequence_lengths[$i-1] < 1 && $i < scalar @row) { # while no bases for a gene proceed to next gene
                        ++$n_printed;
                        ++$i;
                    }
                    if ($row[$i] eq 'empty') { $row[$i] = '-' x $sequence_lengths[$i-1]; } # if no bases set sequence to all gaps
                    print ALIGNMENT uc($row[$i]); # print sequence
                    ++$n_printed;
                    ++$i;
                }
                else {
                    my $bp_left = $row_length;
                    if ($residue ne 'empty') {
                        print ALIGNMENT substr(uc($residue),0,$row_length);
                        $n_printed += length(substr($residue,0,$row_length));
                        $bp_left -= length(substr($residue,0,$row_length));
                        if (length($residue) > $row_length) {
                            $residue = substr($residue,$row_length);
                        }
                        else { $residue = 'empty'; }
                    }
                    while ( $bp_left > 0 and $i < scalar @row) {
                        if ($row[$i] eq 'empty') { $row[$i] = '-' x $sequence_lengths[$i-1]; }
                        print ALIGNMENT substr(uc($row[$i]),0,$bp_left);
                        $n_printed += length(substr($row[$i],0,$bp_left));
                        if (length($row[$i]) > $bp_left) {
                            $residue = substr($row[$i],$bp_left);
                        }
                        $bp_left -= length(substr($row[$i],0,$bp_left));
                        ++$i;
                    }
                }
                print ALIGNMENT "\n";
                close ALIGNMENT or die;
            }
        }
        else {
            for (my $i=0; $i<scalar @row; ++$i) { # for each entry
                if ($i == 0) { # the first entry is the taxon name
                    if ($alignment_format eq 'phylip' or $alignment_format eq 'nexus') { # phylip and nexus have name and sequence on same row
                        print FILE "$row[$i]" . ' ' x ($max_length - length($row[$i]) + 1); # se to that the sequences start at the same position and
                    }                                                                       # that there is at leased one spase between the name and sequence
                    elsif ($alignment_format eq 'fasta') { # fasta look different
                        print FILE ">$row[$i]\n";
                    }
                }
                else { # the rest of the entries are sequences
                    if ( $row[$i] ne 'empty') { # if sequence present
                        print FILE uc($row[$i]);   # print it
                    }
                    else {
                        print FILE '-' x $sequence_lengths[$i-1]; # if no sequence for present gene print gaps instead
                    }
                }
            }
            print FILE "\n"; # new line after sequence in all three formats
        }
    }
    $sth->finish();
    if ($interleaved eq 'y') {
        my @files = glob "XXXtemp_alignment_*.txt";
        foreach (@files) {
            open ALIGNMENT, "<$_" or die "Could not open temporary file $_: $!.\n";
            while (my $row=<ALIGNMENT>) {
                print FILE $row;
            }
            close ALIGNMENT;
            print FILE "\n\n";
        }
        unlink glob "XXXtemp_alignment_*.txt";
    }
    if ($alignment_format eq 'nexus') { print FILE ";\n\nEND;\n"; } # if nexus end DATA block
    close FILE or die;
    if ($get_partitions eq 'y') { # if we should write partitions
        if ($alignment_format eq 'nexus') { # if nexus print charset in a Mr Bayes block in the alignment file
            open PARTITION,">>$alignment_file" or die "Could not open $alignment_file: $!.\n"; # open alignment file
            print PARTITION "\nBEGIN MRBAYES;\n\n"; # print header
            print "    Partitions will be printed in MrBayes block in alignment file.\n";
        }
        else { # if not nexus create new file
            open PARTITION,">$partition_file" or die "Could not open $partition_file: $!.\n";
            print "    Partitions will be printed to $partition_file.\n";
        }
        my $i=0; # counter
        my $start = 1; # first sequence start at first base
        foreach (@gene_alignments) { # for each gene
            if ($sequence_lengths[$i] > 0) {
                if ($alignment_format eq 'nexus') { print PARTITION "    charset ";} # for nexus
                else { print PARTITION "DNA, "; } # otherwise
                print PARTITION "$_ = $start-" . ($start+$sequence_lengths[$i]-1); # print gene name, and start and end base position
                if ($alignment_format eq 'nexus') { print PARTITION ";"; } # if nexus end row with ;
                print PARTITION "\n"; # new row in every format
                $start += $sequence_lengths[$i]; # get start of next partition
            }
            ++$i; # count up which partition to print next
        }
        if ($alignment_format eq 'nexus') { print PARTITION "\nEND;\n"; } # if nexus end Mr Bayes block
        close PARTITION or die;
    }
    $dbh->disconnect();
}

############################
sub sql_where_genes {
    my $dbh = shift @_;
    my $anchor_gene = shift @_; # a gene that need to be present in all taxa
    my @gene_alignments = split /,/, shift @_; # the genes to include in the alignment file
    my $and_or = shift @_;
    if (!$and_or or !($and_or =~ /^and$/i)) { $and_or = 'OR'; }
    if ($gene_alignments[0] eq 'all') { # Get all genes if asked
        @gene_alignments = &get_gene_tables($dbh);
    }    
    my $condition; # condition to get sequence in query
    foreach (@gene_alignments) {
        if (!$condition) {
            $condition = "$_\_sequence != 'empty'"; # add column not equal to empty to query
        }
        else {
            $condition .= " OR $_\_sequence != 'empty'"; # OR since it's enough that the taxon has one sequence
        }
    }
    if ($anchor_gene ne 'n') {
        $condition = '(' . $condition;
        $condition .= ") AND ";
        my @temp_genes = split /,/, $anchor_gene;
        if (scalar @temp_genes > 1) {
            $condition .= "($temp_genes[0]\_sequence != 'empty'";
            for (my $g=1; $g < scalar @temp_genes; ++$g) {
                $condition .= " $and_or $temp_genes[$g]\_sequence != 'empty'";
            }
            $condition .= ")";
        }
        else {
            $condition .= "$temp_genes[0]\_sequence != 'empty'"
        }
    }
    return $condition;
}

##########################################
### Subroutine to print trees to files ###
##########################################

sub print_trees {
    my $database = shift @_; # database name
    my $tree_file_stem = shift @_; # file name stem
    my $taxon = shift @_; # taxon to print for
    my @gene_trees = @_; # genes for which trees should be printed
    my $dbh = DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr; # connect to database
    if ($gene_trees[0] eq 'all') { @gene_trees = &get_gene_tables($dbh); } # if printing trees for all genes, get all gene names
    foreach (@gene_trees) { # for each gene that tree should be printed for
        # get tree
        my $string_taxon = $dbh->quote($taxon);
        my $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE gene = '$_' AND taxon = $string_taxon") or die "Could not prepare statement: " . $dbh->errstr;
        $sth->execute();
        open TREEFILE, ">$tree_file_stem\_$_.tree" or die "Could not open $tree_file_stem\_$_.tree: $!.\n"; # open tree file
        my $printed = 'n'; # flag
        while (my $row = $sth->fetchrow_array()) { # for each recieved tree, but it should only be one
            if (!($row =~ /\(/)) { next; } # if there is no parenthesis there is no tree, move on
            print TREEFILE "$row\n"; # else print tree to file
            $printed = 'y'; # flag that it was printed
        }
        $sth->finish();
        close TREEFILE or die;
        if ($printed eq 'n') { # if no tree was printed
            print "    No tree printed for $_\n"; # print warning
            unlink "$tree_file_stem\_$_.tree"; # and remove file
        }
        else {
            print "    " . ucfirst($_) . " tree for $taxon printed to $tree_file_stem\_$_.tree.\n";
        }
    }
    $dbh->disconnect();
}

####################################################################################
### Subroutine to get the largest block of genes with more than X taxa in common ###
####################################################################################

sub find_lardgest_block {
    my $database = shift @_; # database name
    my $condition = shift @_; # number of taxa that has to be shared between genes for them to be clustered to the same block
    my $dbh = &connect_to_database($database); # connect to database
    my @tables = &get_gene_tables($dbh); # get the genes that are present
    my @matrix; # matrix to store if the gene pair share enough taxa
    for (my $i=0; $i < scalar @tables -1; ++$i) { # for each gene except the last
        $matrix[$i] = ['n']; # initiate matrix
        for (my $j=$i+1; $j < scalar @tables; ++$j) { # for each succesive gene get number of shared taxa
            my $sth = $dbh->prepare("SELECT COUNT(taxon_name) FROM alignments WHERE $tables[$i]\_sequence!='empty' AND $tables[$j]\_sequence!='empty'") or die;
            $sth->execute();
            my $n = $sth->fetchrow_array(); # get the number of taxa
            if ($n >= $condition) { $matrix[$i][$j] = 'y'; } # if number of taxa is greater than cut off, mark it
            else { $matrix[$i][$j] = 'n'; } # if not mark that
        }
    }
    $dbh->disconnect(); # done with database
    my @clusters; # to store separate gene blocks
    my $counter=0; # count the lowest number of a new block
    my @to_check = (0..$#tables); # the genes that are left to check (number in tables array)
    while (@to_check) { # while there are still genes to check
        my @check; # genes to check for specific block
        $check[0] = shift @to_check; # get the first from remaining to check
        while (@check) { # while there are still genes to check for linking genes in that block
            my $present = shift @check; # the gene to check at present
            if ($clusters[$counter]) { # if cluster already initiated
                if (!($clusters[$counter] =~ / $present;/)) { # and the present gene is not present
                    $clusters[$counter] .= " $present;"; # add the gene
                }
            }
            else { $clusters[$counter] = " $present;";} # otherwise initiate cluster
            for (my $i=0; $i < scalar @tables; ++$i) { # for each table
                if ($i<$present) { # just half a matrix, so if iteration is lower than the number of present gene
                    if ($matrix[$i][$present] eq 'y') { # find value in the row of the taxa comparing with, if it is yes
                        my $flag='y'; # flag if the gene to add has aready been checked for linking genes
                        my @temp; # array to temporary store genes that are not linked to the present
                        foreach (@to_check) { # for each gene left to check
                            if ($_ != $i) { push (@temp, $_); } # only keep genes that are not linked to the present gene
                            else { $flag = 'n'; } # if gene linked to the present has not already been checked
                        }
                        @to_check = @temp; # get the genes that are not linked to the present
                        if ($flag eq 'n') { push (@check,$i); } # if the gene that is linked to the present has not already been checked, check it in conjection with this block
                    }
                }
                elsif ($i>$present) { # if the iteration number is higher than the number for the present gene
                    if ($matrix[$present][$i] eq 'y') { # get the value from the row of present, otherwise the same as above
                        my $flag='y';
                        my @temp;
                        foreach (@to_check) {
                            if ($_ != $i) { push (@temp, $_); }
                            else { $flag = 'n'; }
                        }
                        @to_check = @temp;
                        if ($flag eq 'n') { push (@check,$i); }
                    }
                }
            }
        }
        ++$counter; # all the genes linked to the present block has been added, prepare for next block
    }
    my $max=0; # max number of genes for a block
    my $best=0; # the block with most genes
    for (my $i=0; $i<scalar @clusters; ++$i) { # for each block
        my @temp = split /;/, $clusters[$i]; # separate the genes
        if (scalar @temp > $max) { # if it is the block with most genes
            $max = scalar @temp; # save the number of genes
            $best=$i; # and the block
        }
    }
    $clusters[$best] =~ s/ //g; # remove spaces
    my @return_array = split /;/, $clusters[$best]; # separate the gens
    foreach (@return_array) { $_ = $tables[$_]; } # get the names of the genes
    return @return_array; # return the array
}

###################################################
### Subroutine to change tip names in tree file ###
###################################################

sub change_names {
    my $database = shift @_;
    my $path  = shift @_; # path to treebender
    my $tree = shift @_;
    my $file_stem = shift @_;
    my $file = shift @_;
    my $dbh = &connect_to_database ($database);
    my $sth;
    if ($file eq 'n') {
        my ($gene,$taxon) = split /,/,$tree;
        print "Getting tree for taxon $taxon, gene $gene from $database.\n";
        my $string_taxon = $dbh->quote($taxon);
        my $string_gene = $dbh->quote($gene);
        $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE gene=$string_gene AND taxon=$string_taxon") or die;
        $sth->execute();
        my $temp = $sth->fetchrow_array();
        $sth->finish;
        if ($temp) {
            open TREE,">XXXtemp.tree" or die "Could not open XXXtemp.tree: $!.\n";
            print TREE $temp;
            close TREE or die;
        }
        else { die "Could not find tree for gene \"$string_gene\" and taxon \"$string_taxon\". Quitting...\n"; }
        $tree = 'XXXtemp.tree';
    }
    else { print "Getting tree from $tree.\n"; }
    if (-e $tree) {
        my @tip_lables = `${path}$External_program::treebender -t '\\n' < $tree`;
        chomp (@tip_lables);
        if (scalar @tip_lables > 1) {
            my $n_accno=0;
            my $n_taxon=0;
            my $switch_statement='';
            foreach (@tip_lables) {
                if ($_ =~ /Taxon_/) { ++$n_taxon; }
                elsif ($_ =~ /^[A-Z_]{2,3}[0-9]+/) { ++$n_accno; }
            }
            if ($n_accno == $n_taxon) {
                print "Could not decide what type of tip lables the tree has so could not switch tip names.\n";
            }
            elsif ($n_accno > $n_taxon) {
                foreach (@tip_lables) {
                    my $string = $dbh->quote($_);
                    $sth=$dbh->prepare("SELECT accno,species FROM gb_data WHERE accno=$string") or die;
                    $sth->execute();
                    my ($accno,$species) = $sth->fetchrow_array();
                    $species =~ s/[ \'\"\(\)\[\]\{\}\:\;]/_/g;
                    $switch_statement .= "$_|$accno\_$species,";
                }
            }
            elsif ($n_accno < $n_taxon) {
                foreach (@tip_lables) {
                    my $sql_ = $_;
                    $sql_ =~ s/(^')|('$)//g;
                    $sql_ = $dbh->quote($sql_);
                    $sth=$dbh->prepare("SELECT species_name FROM alignments WHERE taxon_name=$sql_") or die;
                    $sth->execute();
                    my $species = $sth->fetchrow_array();
                    if (!$species) { $species = $_; }
                    $species =~ s/[ \'\"\(\)\[\]\{\}\:\;]/_/g;
                    $switch_statement .= "$_|$species,";
                }
            }
            if ($switch_statement =~ s/,$//) {
                system "${path}$External_program::treebender -c \"$switch_statement\" < $tree > $file_stem.switched_names.tree";
                if (-e 'XXXtemp.tree') { unlink 'XXXtemp.tree'; }
                if (-e "$file_stem.switched_names.tree") { return "$file_stem.switched_names.tree"; }
                else { print "No tree file with new names created.\n"; return 'No'; }
            }
            else {
                print "Could not construct a statement for switching tip names. Sorry.\n";
                if (-e 'XXXtemp.tree') { unlink 'XXXtemp.tree'; }
                return 'No';
           }
        }
    }
    else { print "WARNING!!! No tree file named $tree. Cannot switch tip names.\n"; return 'No'; }
}

sub get_distinct_entries {
    my $database = shift @_;
    my $out_stem = shift @_;
    my $dbh = &connect_to_database ( $database );
    foreach (@_) { # treat remaining arguments as column names
        print "Getting distinct entries for column $_.\n";
        my $sth = $dbh->prepare("SELECT DISTINCT($_) FROM gb_data WHERE $_!='empty'") or die;
        $sth->execute();
        open OUTPUT, ">$out_stem\_$_.txt" or die "Could not write to $out_stem\_$_.txt: $!.\n";
        while (my @row=$sth->fetchrow_array()) {
            for (my $i=0; $i<scalar @row; ++$i) {
                if ($i==0) { print OUTPUT $row[$i]; }
                else { print OUTPUT "\t$row[$i]"; }
            }
            print OUTPUT "\n";
        }
        $sth->finish();
        close OUTPUT or die;
        print "Output printed to $out_stem\_$_.txt.\n";
    }
    $dbh->disconnect();
}

sub change_entries {
    my $database = shift @_;
    my $rm_entries = shift @_;
    my $dbh = &connect_to_database ($database);
    my $updates=0;
    foreach (@_) { # treat rest of arguments as file names
        if ($rm_entries ne 'y') { print "Making changes given in $_.\n"; }
        else { print "Removing sequences matching entries given in $_.\n" }
        my $local_updates=0;
        my $i=0;
        open INFILE, "<$_" or die "Could not open $_: $!.\n";
        my @columns;
        while (my $row=<INFILE>) {
            chomp($row);
            if (!$row or !($row =~ /.+/)) { next; }
            if ($i == 0) {
                my @temp = split /\t/, $row;
                foreach (@temp) {
                    if ($_ and $_ =~ /[a-zA-Z0-9_]+/) { push (@columns,$_); }
                }
            }
            else {
                my $col=0;
                my @temp = split /\t/, $row;
                if ($rm_entries ne 'y') {
                    for (my $j=0; $j<((scalar @temp) -1); $j+=2) {
                        my $string1 = $dbh->quote($temp[$j+1]);
                        my $string2 = $dbh->quote($temp[$j]);
                        my $n=$dbh->do("UPDATE gb_data SET $columns[$col]=$string1 WHERE $columns[$col]=$string2") or die;
                        ++$col;
                        if ($n <1) { print STDERR "WARNING no updates were made from '$temp[$j]' to '$temp[$j+1]'.\n"; }
                        else { $local_updates += $n; }
                    }
                }
                else {
                    my $string='';
                    for (my $j=0; $j<scalar @columns; ++$j) {
                        if ($j>=scalar @temp) { last; }
                        if ($j == 0) { $string="$columns[$j]=" . $dbh->quote($temp[$j]); }
                        else { $string=" AND $columns[$j]=" . $dbh->quote($temp[$j]); }
                    }
                    my $n=$dbh->do("DELETE FROM gb_data WHERE $string") or die;
                    if ($n <1) { print STDERR "WARNING sequences were deleted where $string.\n"; }
                    else { $local_updates += $n; }
                }
            }
            ++$i;
        }
        close INFILE or die;
        if ($rm_entries ne 'y') { print "Made $local_updates updates.\n"; }
        else { print "Deleted $local_updates rows.\n"; }
        $updates += $local_updates;
    }
    $dbh->disconnect();
    if ($rm_entries ne 'y') { print "In total $updates updates were done.\n"; }
    else { print "In total $updates rows were deleted.\n"; }
}

### the two folowing subroutines are translations from the ruby script newic2fasta.rb available at http://mafft.cbrc.jp/alignment/software/treein.html
sub resolve { # Not sure this routine works properly, should not be used since treebender give fully resolved trees
    my $tree = shift @_;
    while (1) {
        $tree =~ s/,([0-9]+):(\-{0,1}[0-9\.]+),([0-9]+):(\-{0,1}[0-9\.]+)/XXX/;
        my $hit1 = $1;
        my $hit2 = $2;
        my $hit3 = $3;
        my $hit4 = $4;
        if (!($tree =~ /XXX/)) { last; }
        my $poshit = @-;
        my $i = $poshit;
        my $height = 0;
        while ($i >= 0) {
            if ($height == 0 && index($tree,$i) eq '(') {last;}
            if (index($tree,$i) eq ')') {
                ++$height;
            }
            elsif (index($tree,$i) eq '(') {
                --$height;
            }
            --$i;
        }
        my $poskakko = $i;
        my $zenhan = substr($tree,0,$poskakko);
        if ($poskakko == -1) { $zenhan = ""; }
        my $treelen = length($tree);
        $tree = $zenhan . "(" . substr($tree,$poskakko+1,$treelen);
        $tree =~ s/XXX/,$hit1:$hit2\):0,$hit3:$hit4/;
    }
    return $tree;
}

sub newick_to_mafft {
    my $tree = shift @_;
    my $scale = shift @_;
    $tree =~ s/_.*?:/:/g;
    $tree =~ s/[0-9]\.[0-9]*e-[0-9][0-9]/0/g;
    $tree =~ s/\[.*?\]//g;
    $tree =~ s/ //g;
    my $mafft_tree = '';
    my @memi = (-1,-1);
    my @leni = (-1,-1);
    while ($tree =~ /\(/ ) {
        $tree = &resolve( $tree );
        $tree =~ s/\(([0-9]+):(\-{0,1}[0-9\.]+),([0-9]+):(\-{0,1}[0-9\.]+)\)/XXX/;
        $memi[0] = $1;
        $leni[0] = $2 * $scale;
        $memi[1] = $3;
        $leni[1] = $4 * $scale;

        if ($leni[0] > 10 || $leni[1] > 10) {
            print STDERR "\n";
            print STDERR "Please check the scale of branch length!\n";
            print STDERR "The unit of branch lengths must be 'substitution/site'\n";
            print STDERR "If the unit is 'substition' in your tree, please\n";
            print STDERR "use the scale argument,\n";
            print STDERR "% newick2mafft scale in > out\n";
            print STDERR "where scale = 1/(alignment length)\n";
            print STDERR "\n";
            return 'empty';
        }
        if ($memi[1] < $memi[0]) {
            @memi = reverse (@memi);
            @leni = reverse (@leni);
        }
        $tree =~ s/XXX/$memi[0]/;
        $mafft_tree .= sprintf( "%5d %5d %10.5f %10.5f", $memi[0], $memi[1], $leni[0], $leni[1] );
        $mafft_tree .= "\n";
   }
   return $mafft_tree;
}

sub read_tree_to_alignment_groups { # saves tree to database
    my $file = shift @_; # tree file
    my $dbh = shift @_; # database handler
    my $gene = shift @_; # gene to save it under
    my $taxon = shift @_; # taxon to save it under
    my $method = shift @_; # method to annotate it with
    open TREEFILE, "<$file" or die "Could not open tree file $file: $!.\n"; # open and
    my $tree='';
    while (my $row = <TREEFILE>) { $tree .= $row; } # read tree
    close TREEFILE or die;
    # check if entry for taxon and gene is present
    my $sth = $dbh->prepare ("SELECT COUNT(taxon) FROM alignment_groups WHERE taxon='$taxon' AND gene='$gene'") or die;
    $sth->execute();
    my $count = $sth->fetchrow_array();
    $sth->finish();
    # uppdate database
    my $string = $dbh->quote($tree);
    #print STDERR "Tree to save from $file (method: $method, gene: $gene, and taxon: $taxon):\n$string\n";
    if ($count and $count >= 1) {
	my $update = $dbh->do("UPDATE alignment_groups SET tree=$string,tree_method='$method' WHERE gene='$gene' AND taxon='$taxon'") or die "Could not update database: " . $dbh->errstr;
	return $update; # return how many rows were uppdated
    }
    else {
	my $update = $dbh->do("INSERT INTO alignment_groups (taxon,gene,tree,alignable) VALUES ('$taxon','$gene',$string,2)") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
	return $update; # return how many rows were inserted
    }
}


sub process_cluster_array { # this should be moved to Support functions and be implemented in the tree clustering as well
    my $dbh = shift @_;
    my @accnos = split(/,/,shift @_); # file with each cluster on a separate row and accnos in cluster separated by space
    my $table = shift @_; # gene that is being clustered
    my $n_uppdate = 0; # count number of rows that has been uppdated
    my $sth; # search string handler
    my $i = 0; #counter
    #print "processing.\n";
    if (scalar @accnos > 1) { # if more than one accno we have a cluster
        my $longest = $accnos[0]; # set first accno as the one with longest sequence
        my $max = 0; # set the max observed sequence length to zero
        my %values;
        my %clusters;
        foreach (@accnos) {
            $sth = $dbh->prepare("SELECT cluster FROM $table WHERE accno='$_'") or die;
            $sth->execute();
            my $cluster = $sth->fetchrow_array();
            $sth->finish();
            my $query_accno;
            if ($cluster eq 'lead' or $cluster eq 'empty') { $query_accno = $_; }
            else { $query_accno = $cluster; }
            $sth = $dbh->prepare("SELECT LENGTH($table.sequence),gb_data.proportion_N FROM $table INNER JOIN gb_data on $table.accno=gb_data.accno WHERE $table.accno='$query_accno'")
                or die "Couldn't prepare statement: " . $dbh->errstr; #get sequence length and proportion of N's
            $sth->execute();
            my ($length,$proportion_N);
            ($length,$proportion_N) = $sth->fetchrow_array();
            $clusters{$_} = $cluster;
            $sth->finish();
            $values{$_} = $length*(1-$proportion_N);
        }
        foreach (keys %values) { # get individual row
            if ( $values{$_} > $max ) { # if length of the sequence times poroportion of N's in sequence is > than any previous value
                $max = $values{$_};     # save value
                $longest = $_;  # and accno
            }
        }
        undef %values;
        # get the cluster annotation for the 'best' sequence
        my $lead = $clusters{$longest}; # save the cluster annotation
        $sth->finish();
        if ($lead eq 'empty') { # if no previous cluster
            $lead = $longest; # the leade is the 'best' sequence
            $n_uppdate += $dbh->do("UPDATE $table SET cluster=\'lead\' WHERE accno=\'$longest\'")
                or die "Couldn't uppdate database: " . $dbh->errstr; # set 'best' sequence to lead in it's cluster
        }
        elsif ($lead eq 'lead') { $lead = $longest; } # if the 'best' sequence is already the lead of a cluster set lead to the 'best' sequence
        my $condition; # reuse condition
        for ($i=0; $i < scalar @accnos; ++$i) { # for each accno in cluster
            if ($accnos[$i] eq $lead) { next; } # if it is the lead, skip to next
            elsif ($condition) { $condition .= " or accno=\'$accnos[$i]\'"; } # if it is not the first begin with 'or'
            else { $condition = "accno=\'$accnos[$i]\'"; } # if it is the first to add just add
            if ($i == (scalar @accnos)-1 or !($i % 100)) {
                $n_uppdate += $dbh->do("UPDATE $table SET cluster='$lead' WHERE $condition") or die "Couldn't uppdate database: " . $dbh->errstr;
                undef $condition;
            }
        }
        for ($i=0; $i < scalar @accnos; ++$i) { # for each accno in cluster
            my $query_accno;
            if ($clusters{$accnos[$i]} eq 'lead' or $clusters{$accnos[$i]} eq 'empty') { $query_accno = $accnos[$i]; }
            else { $query_accno = $clusters{$accnos[$i]}; }
            if ($condition) { $condition .= " or cluster='$query_accno'"; } # if it is not the first begin with 'or'
            else { $condition = "cluster='$query_accno'"; } # if it is the first to add just add
            if ($i == (scalar @accnos)-1 or !($i % 100)) {
                $n_uppdate += $dbh->do("UPDATE $table SET cluster='$lead' WHERE $condition") or die "Couldn't uppdate database: " . $dbh->errstr;
                undef $condition;
            }
        }
        undef $condition; # reuse condition
    }
    else { # if only one sequence in cluster
        # get the cluster
        $sth = $dbh->prepare("SELECT cluster FROM $table WHERE accno=\'$accnos[0]\'") or die "Couldn't prepare statement: " . $dbh->errstr;
        $sth->execute();
        my $cluster = $sth->fetchrow_array();
        $sth->finish();
        # if not previously included in a cluster make it lead it its own cluster
        if ($cluster eq 'empty') {
            $n_uppdate = $dbh->do("UPDATE $table SET cluster=\'lead\' WHERE accno=\'$accnos[0]\'") or die "Couldn't uppdate database: " . $dbh->errstr;
        }
    }
    return $n_uppdate;
}

sub stats {
    my $database = shift @_;
    my $stats_type = shift @_;
    my $detailed=shift @_;
    my $path = shift @_;
    my $anchor_gene = shift @_; # a gene that need to be present in all taxa
    if ($stats_type eq 'alignment') { &alignment_stats($database,$detailed,$path,$anchor_gene,@_); }
    elsif ($stats_type eq 'cluster') { &cluster_stats($database,$detailed,$path,@_); }
    else { print STDERR "Do not recognice '$stats_type' as type of statistics. No statistics produced.\n"; }
}

sub cluster_stats{
    my $database = shift @_;
    my $detailed=shift @_;
    my $path = shift @_;
    my @genes = @_;
    my $dbh = PifCosm_support_subs::connect_to_database($database);
    if ($genes[0] eq 'all') { @genes=&get_gene_tables($dbh); }
    foreach my  $gene (@genes) {
        print "### Clusters for $gene ###\n";
        my $query = "SELECT $gene.accno";
        my @lead_accno;
        if ($detailed eq 'y') { $query .= ",gb_data.species"; }
        $query .= " FROM $gene ";
        if ($detailed eq 'y') { $query .= "INNER JOIN gb_data ON $gene.accno=gb_data.accno "; }
        $query .= "WHERE $gene.cluster='lead'";
        my $sth = $dbh->prepare($query) or die;
        $sth->execute();
        while ( my @row = $sth->fetchrow_array) { push(@lead_accno,@row); }
        $sth->finish();
        while (@lead_accno) {
            my $accno=shift @lead_accno;
            my $species;
            if ($detailed eq 'y') { $species=shift @lead_accno; }
            print "$accno ";
            if ($species) { print "$species "; }
            print "- ";
            $query = "SELECT $gene.accno";
            if ($detailed eq 'y') { $query .= ",gb_data.species"; }
            $query .= " FROM $gene ";
            if ($detailed eq 'y') { $query .= "INNER JOIN gb_data ON $gene.accno=gb_data.accno "; }
            $query .= "WHERE $gene.cluster='$accno'";
            my $sth = $dbh->prepare($query) or die;
            $sth->execute();
            while ( my @row = $sth->fetchrow_array) {
                print "$row[0]";
                if ($detailed eq 'y') { print " $row[1]"; }
                print ", ";
            }             
            $sth->finish();
            print "\n";
        }
    }
}

sub alignment_stats { ### Sub to calculate some statistics for the alignments
    my $database = shift @_;
    my $detailed=shift @_;
    my $path = shift @_;
    my $anchor_gene = shift @_; # a gene that need to be present in all taxa
    my @genes = @_;
    my $dbh = PifCosm_support_subs::connect_to_database($database);
    if ($genes[0] eq 'all') { @genes=&get_gene_tables($dbh); }
    my $condition;
    if (&table_present( $dbh, "alignments") eq 'y') {
        my $columns='taxon_name';
        my $latin_name = &column_present($dbh, "alignments","species_name");
        my $n=0;
        if ($latin_name eq 'y') {
            $columns.= ',species_name';
        }
        foreach my $gene (@genes) {
            $columns .= ",$gene\_accno,$gene\_sequence";
	}
	$condition = &sql_where_genes ( $dbh, $anchor_gene, join(',', @genes), 'AND' );
        my $sth = $dbh->prepare("SELECT $columns FROM alignments WHERE $condition") or die;
        $sth->execute();
        my %seq_lengths;
        my %seq_present;
        my %seq_lead;
        my %seq_aligned_length;
        my %seq_not_lead;
        my $n_taxa=0;
        my $n_taxa_all=0;
        my $gene_taxa_string='';
        while (my @row = $sth->fetchrow_array) {
            ++$n_taxa;
            my $taxon_name=shift @row;
            my $species_name;
            if ($latin_name eq 'y') { $species_name=shift @row; }
            if ($detailed eq 'y') {
                print "$taxon_name";
                if ($species_name) { print " | $species_name"; }
                print "   ";
            }
            my $seq_length=0;
            my $j=0;
            my $all_genes='y';
            my $non_lead='n';
            my %taxon_string;
            for (my $i=0; $i<scalar @row; ++$i) {
                if ($i % 2 == 0) {
                    if ($detailed eq 'y') { print " $genes[$j]:"; }
                    if ($row[$i] ne 'empty') {
                        $gene_taxa_string.="$j,";
                        my $row_sth=$dbh->prepare("SELECT cluster,taxon_string FROM $genes[$j] INNER JOIN gb_data ON $genes[$j].accno=gb_data.accno WHERE gb_data.accno='$row[$i]'") or die;
                        $row_sth->execute();
                        my @data=$row_sth->fetchrow_array();
                        if ($data[0] eq 'lead') {
                            ++$seq_lead{$genes[$j]};
                        }
                        else { ++$seq_not_lead{$genes[$j]}; }
                        if ($detailed eq 'y') { print " $row[$i] ($data[0])"; }
                        $taxon_string{$data[1]}++;
                        ++$seq_present{$genes[$j]};
                    }
                    else {
                        if ($detailed eq 'y') { print ' No'; }
                        $all_genes='n';
                    }
                    ++$j;
                }
                else {
                    if ($row[$i] ne 'empty') {
                        $seq_aligned_length{$genes[$j-1]}=length($row[$i]);
                        $row[$i]=~ s/-//g;
                        my $unaligned_length = length($row[$i]);
                        $seq_lengths{$genes[$j-1]}+= $unaligned_length;
                        if ($detailed eq 'y') { print " " . $unaligned_length . "bp"; }
                    }
                    if ($detailed eq 'y') { print ";"; }
                }
            }
            $gene_taxa_string =~ s/,$/|/;
            my @keys = sort { length $a <=> length $b } keys %taxon_string;
            if (scalar @keys > 1) {
                my $n_taxa=1;
                for (my $i=0; $i<scalar @keys-1; ++$i) {
                    my $present='n';
                    for (my $j=$i+1; $j<scalar @keys; ++$j) {
                        if ($keys[$j] =~ /$keys[$i]/) { $present='y'; last; }
                    }
                    if ($present eq 'n') { ++$n_taxa; }
                }
                if ($n_taxa > 1) { print " [WARNING, more than one taxon present]"; }
            }
            if ($detailed eq 'y') { print "\n"; }
            if ($all_genes eq 'y') { ++$n_taxa_all; }
        }
        if ($detailed eq 'y') { print "*******************************************\n"; }
        my $total=0;
        my $total_seq_length=0;
        foreach (keys %seq_aligned_length) { $total_seq_length+=$seq_aligned_length{$_}; }
        print "Number of taxa = $n_taxa\nAlignment length = $total_seq_length\n";
        foreach (keys %seq_present) {
             print "    $_ is present for $seq_present{$_} taxa.\n";
             $total+=$seq_present{$_};
        }
        print "    $n_taxa_all taxa have all of the genes.\n";
        print "    Gene covarage is " . ($total/($n_taxa * scalar @genes)) . "\n";
        print "    The nucleotide covarage for each gene is:\n";
        $total=0;
        foreach (keys %seq_lengths) {
            print "    $_ " . ($seq_lengths{$_}/($n_taxa * $seq_aligned_length{$_})) . "\n";
            $total+=$seq_lengths{$_};
        }
        print "    The overall nucleotide covarage is " . ($total/($total_seq_length*$n_taxa)) . "\n";
        $total=0;
        print "    The proportion non-lead sequences for each gene is:\n";
        foreach (keys %seq_not_lead) {
            print "    $_ " . ($seq_not_lead{$_}/$seq_present{$_}) . "\n";
            $total+=$seq_not_lead{$_};
        }
        print "The total number of non-lead sequences is $total\n";
        my @decisiveness=`${path}$External_program::contree -D "$gene_taxa_string" -i 1000`;
        print @decisiveness;
    }
    else { print STDERR "WARNING!!! No alignments table present, no statistics created.\n"; }

}

1;
