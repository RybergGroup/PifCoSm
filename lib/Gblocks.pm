package Gblocks;
use strict;
use PifCosm_support_subs;

sub Gblocks {
    my $database = shift @_;
    my $path = shift @_; # path to Gblocks
    my $b0 = shift @_;
    my $b1 = shift @_;
    my $b2 = shift @_;
    my $b3 = shift @_;
    my $b4 = shift @_;
    my $b5 = shift @_;
    my $min_length = shift @_;
    my @genes = @_;
    if ($genes[0] eq 'all') { # Get all genes if asked
        my $dbh = PifCosm_support_subs::connect_to_database($database);
        @genes = PifCosm_support_subs::get_gene_tables($dbh);
        $dbh->disconnect;
    }
    foreach (@genes) {
       print "***************\nRunning Gblocks for $_.\n";
       my $dbh = PifCosm_support_subs::connect_to_database($database);
       my $sth = $dbh->prepare("SELECT $_\_accno,$_\_sequence FROM alignments WHERE $_\_sequence != 'empty'") or die;
       my $n_seq = PifCosm_support_subs::print_fasta("XXXalignment.fst",$sth);
       $sth->finish();
       if (!$n_seq || $n_seq < 2) { print "Less than two sequences for $_. Ignoring it.\n"; unlink "XXXalignment.fst"; next; }
       my $switches ='';
       if ($b0 > 2) { $switches .= " -b0=$b0"; }
       if ($b1 > 0) { $switches .= " -b1=$b1"; }
       if ($b2 > 0) { $switches .= " -b2=$b2"; }
       if ($b3 > 0) { $switches .= " -b3=$b3"; }
       if ($b4 > 1) { $switches .= " -b4=$b4"; }
       if ($b5 eq 'n' or $b5 eq 'h' or $b5 eq 'a') { $switches .= " -b5=$b5"; }
       my @gblocks = `${path}$External_program::gblocks XXXalignment.fst -t=d$switches`;
       PifCosm_support_subs::read_fasta_to_alignments_table ("XXXalignment.fst-gb",$dbh,$_,'accno');
       unlink glob "XXXalignment.fst*";
       $sth = $dbh->prepare("SELECT $_\_accno,$_\_sequence FROM alignments WHERE $_\_sequence != 'empty'") or die;
       $sth->execute();
       my $updates=0;
       while (my @seq = $sth->fetchrow_array()) {
           if (length $seq[1] < $min_length or !($seq[1] =~ /[^\s-]/)) { $updates += PifCosm_support_subs::remove_accno_from_alignments($dbh,$_,$seq[0]); }
       }
       $sth->finish();
       $sth = $dbh->prepare("SELECT COUNT($_\_accno) FROM alignments WHERE $_\_sequence != 'empty'") or die;
       $sth->execute();
       my $count = $sth->fetchrow_array();
       foreach (@gblocks) { if ($_ =~ /^Original alignment:/ or $_ =~ /^Gblocks alignment:/) { print $_; } }
       print "$updates sequences were deleted from the alignment for $_ due to that they were to short or only contained gaped sites.\n";
       if ($count < 1) { print STDERR "WARNING!!! No sequences left for $_!!!\n"; }
       $sth->finish();
       $dbh->disconnect();
    }
    print "***************\n";
}
1;
