package Exclude_rogues;
use strict;
use PifCosm_support_subs;

sub exclude_rogues {
    my $database = shift @_;
    my $path = shift @_; # path to RogueNaRok, RAxML, and treebender
    my $print_taxa = shift @_;
    my $output_tree = shift @_;
    my $remove_from_trees = shift @_;
    my $remove_from_alignment = shift @_;
    my $name_stem = shift @_;
    my $dbh = PifCosm_support_subs::connect_to_database ($database);
    my $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE taxon='combined_boot' AND gene='combined_boot'") or die;
    $sth->execute();
    open BOOTTREE, ">XXXboot.trees" or die "Could not open temporary file XXXboot.trees: $!.\n";
    print BOOTTREE $sth->fetchrow_array();
    close BOOTTREE or die;
    unlink glob "*.XXXtempRogue";
    my $RogueNaRok = `${path}RogueNaRok -i XXXboot.trees -n XXXtempRogue`;
    undef $RogueNaRok;
    unlink "XXXboot.trees";
    open ROGUES, "<RogueNaRok_droppedRogues.XXXtempRogue" or die "Could not open temporary file RogueNaRok_droppedRogues.XXXtempRogue: $!.\n";
    my @rogues;
    my $i = 0;
    while (my $infile = <ROGUES>) {
        my @temp = split /\t/, $infile;
        if ($temp[2] ne 'taxon' and $temp[2] ne 'NA') {
            $rogues[$i++] = $temp[2];
        }
    }
    close ROGUES or die;
    unlink glob "*.XXXtempRogue";
    if (scalar @rogues > 0) {
        if ($print_taxa eq 'y') {
            print "The following " . scalar @rogues . " rogue taxa were identified by RogueNaRok:\n";
            for ($i = 0; $i < scalar @rogues; ++$i) {
                if ($i>0) { print ", "; }
                print $rogues[$i];
            }
        }
        print ".\n";
        if ($remove_from_alignment eq 'y') {
            my $deleted = 0;
            foreach (@rogues) {
                $deleted += $dbh->do("DELETE FROM alignments WHERE taxon_name = '$_'") or die;
            }
            print "$deleted rogue taxa were removed from the alignment.\n";
        }
        if ($output_tree eq 'y' or $remove_from_trees eq 'y') {
            my $drop_tips = join ",", @rogues;
            $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE taxon='combined' AND gene='combined'") or die;
            $sth->execute();
            open MLTREE, ">XXXml.tree" or die "Could not open temporary file XXXml.tree: $!.\n";
            print MLTREE $sth->fetchrow_array();
            $sth->finish();
            close MLTREE or die;
            system "${path}treebender -d $drop_tips < XXXml.tree > XXXml.norogue.tree";
            unlink "XXXml.tree";
            $sth = $dbh->prepare("SELECT tree FROM alignment_groups WHERE taxon='combined_boot' AND gene='combined_boot'") or die;
            $sth->execute();
            unlink "XXXboot.trees";
            my @boot_trees = split /\n/,$sth->fetchrow_array();
            foreach (@boot_trees) {
                open BOOTTREE, ">XXXtemp.trees" or die "Could not open temporary file XXXboot.trees: $!.\n";
                print BOOTTREE $_;
                system "${path}treebender -d $drop_tips < XXXtemp.trees >> XXXboot.trees";
                close BOOTTREE or die;
            }
            unlink "XXXtemp.trees";
            unlink glob "*.XXXsupport";
            my $raxml = `${path}raxmlHPC -f b -t XXXml.norogue.tree -z XXXboot.trees -m GTRGAMMA -n XXXsupport`;
            undef $raxml;
            unlink glob "XXXschematic.alignment*";
            unlink "XXXml.norogue.tree";
            if ($remove_from_trees eq 'y') {
                open MLTREE, "<RAxML_bipartitions.XXXsupport" or die "Could not open temporary file RAxML_bipartitions.XXXsupport: $!.\n";
                my $mltree = <MLTREE>;
                close MLTREE or die;
                $mltree = $dbh->quote($mltree);
                my $updated = $dbh->do("UPDATE alignment_groups SET tree=$mltree WHERE taxon='combined' and gene='combined'") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
                if ($updated < 1) { print STDERR "WARNING!!! Was not able to store ML tree.\n"; }
                else { print "Updated tree for taxon combined and gene combined, with tree without rogue taxa.\n"; }
                open BOOTTREE, "<XXXboot.trees" or die "Could not open temporary file XXXboot.trees: $!.\n";
                my $boot_trees = <BOOTTREE>;
                close BOOTTREE or die;
                $updated=0;
                $boot_trees = $dbh->quote($boot_trees);
                $updated = $dbh->do("UPDATE alignment_groups SET tree=$boot_trees WHERE taxon='combined_boot' and gene='combined_boot'") or die "Could not insert value in alignment_groups: " . $dbh->errstr;
                if ($updated < 1) { print STDERR "WARNING!!! Was not able to store bootstrap trees.\n"; }
                else { print "Updated trees for taxon combined_boot and gene combined_boot, with trees without rogue taxa.\n"; }
            }
            if ($output_tree eq 'y') {
                rename "RAxML_bipartitions.XXXsupport", "$name_stem.no_rogue.tree" or die "Could not find tree without rogue taxa for output.\n";
                print "Tree without rogue taxa saved in $name_stem.no_rogue.tree.\n";
            }
            unlink glob "*.XXXsupport";
            unlink "XXXboot.trees";
        }
    }
    else { print "No rogues were identified.\n"; }
}

1;
