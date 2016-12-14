package CD_hit_clust;
use strict;
use PifCosm_support_subs;

##################################################
### Subroutine to cluster monophyletic species ###
##################################################

sub CD_hit_clust {
    my $database = shift @_; # get the database
    my $path = shift @_; # path to RAxML and phylomand programs
    my %cut_off = split(/,/, shift @_); # cut off value for branch length clustering
    my $min_length = shift @_; # get the minimum length of sequences to cluster
    my $n_treads = shift @_; # get number of threads
    my $dbh = PifCosm_support_subs::connect_to_database($database); #DBI->connect("dbi:SQLite:dbname=$database","","") or die "Couldn't connect to database: " . DBI->errstr;
    my @tables = keys %cut_off;
    if ($tables[0] =~ /^all$/i) { 
        my $cut_off = $cut_off{$tables[0]};
        undef %cut_off;
        @tables = PifCosm_support_subs::get_gene_tables($dbh); # get gene tables
        foreach (@tables) { $cut_off{$_} = $cut_off; }
    }
    # cluster for each gene separately
    for (my $j = 0; $j < scalar @tables; ++$j) {
        if (PifCosm_support_subs::table_present($dbh,$tables[$j]) eq 'n') { next; }
        print "Clustering $tables[$j] using the cut off $cut_off{$tables[$j]}.\n";
        my $n_uppdate = 0;
        print $tables[$j] . "\n";
        # get number of accnos
        my $sth= $dbh->prepare("SELECT COUNT(accno) FROM $tables[$j] WHERE cluster = 'lead'") or die;
        $sth->execute();
        my $count= $sth->fetchrow_array();
        $sth->finish();
        if ($count && $count >0) {
            $sth = $dbh->prepare("SELECT accno,sequence FROM $tables[$j] WHERE cluster = 'lead' AND LENGTH(sequence) >= $min_length") or die;
        }
        else {
            $sth = $dbh->prepare("SELECT accno,sequence FROM $tables[$j] WHERE LENGTH(sequence) >= $min_length") or die;
        }
        my $n_seq = PifCosm_support_subs::print_fasta("XXX$tables[$j].fst",$sth);
        if ($n_treads > 1) { system "${path}cdhit-est -i XXX$tables[$j].fst -o XXX$tables[$j] -c $cut_off{$tables[$j]} -T $n_treads > XXX$tables[$j].out"; }
        else { system "${path}cdhit-est -i XXX$tables[$j].fst -o XXX$tables[$j] -c $cut_off{$tables[$j]} > XXX$tables[$j].out"; }
        open CLUSTERS, "<XXX$tables[$j].clstr" or die "Could not open XXX$tables[$j].clstr: $!.\n";
        my $cluster;
        while (my $row = <CLUSTERS>) {
            if ($row =~ /^>/ and $cluster) {
                PifCosm_support_subs::process_cluster_array($dbh,$cluster,$tables[$j]);
                undef $cluster;
            }
            elsif ($row =~ /^[0-9]+\s+[0-9A-Za-z,\s]+>(\w+)/) {
                if (!$cluster) { $cluster = $1; }
                else { $cluster .= ",$1"; }
            }
        }
        if ($cluster) { PifCosm_support_subs::process_cluster_array($dbh,$cluster,$tables[$j]); }
        unlink glob ("XXX$tables[$j]*");
    }
    $dbh->disconnect();
}
1;
