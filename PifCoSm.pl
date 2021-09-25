#!/usr/bin/env perl

use warnings;
use strict;
use DBI;
use FindBin ();
use lib "$FindBin::Bin/lib";
use PifCosm_support_subs;
use Create_alignment_groups;
use Gb_parser;
use Hammer_sort;
use Cluster_indv_sequences;
use Define_individuals;
use Tree_cluster;
use CD_hit_clust;
use Gene_linker;
use CV_gene_linker;
use Latin_binomials;
use Align_sequences;
use Refine_alignment;
use Gblocks;
use Concatenated_tree;
use Exclude_rogues;
use External_program;

### Usage string ###
my $usage = "\
USAGE: [perl] PifCoSm.pl batch_file.txt\
       [perl] PifCoSm.pl [options] (-o|--modules) module[,module..]\
       [perl] PifCoSm.pl (-h|--help)\
";

my $help = "\
MODULES and OPTIONS:\
\
gb_parser - will parse a text file in full GenBank format and store the\
            information in a table called gb_data in the database.\
        -d/--database\
        -g/--gb_file\
        -m/--max_seq_length\
\
change_entries - will change the annotation gb_data table.\
        -d/--database\
        -H/--changes_files\
        -ch/--change_option\
\
define_individuals - will define individuals based on one or several columns.\
        -d/--database\
        -v/--individual_columns\
        -V/--individual_condition\
\
gene_parser - will try to match each sequence to hmms of different genes and\
              sort out the different genes to separate tables.\
        -d/--database\
        -M/--hmm_database\
        -Z/--print_non_match\
        -l/--min_seq_length_gene\
        -e/--e_value_cut_off\
        -n/--n_cut_off\
        -i/--inter_genes\
        -P/--path\
\
cluster_individuals - will cluster the sequences of each gene that belong to\
                      the same individual as defined by define_individuals.\
        -d/--database\
        -l/--min_seq_length_gene\
\
tree_clustering - will cluster sequences of each gene based on trees built for\
                  each distinct taxon string.\
        -d/--database\
        -C/--tree_clust_type\
        -s/--sim_cut_off\
        -l/--min_seq_length_gene\
        -T/--threads\
        -P/--path\
\
CD_hit_cluster - will cluster sequences based on similarity using CDhit-est.\
        -d/--database\
        -s/--sim_cut_off\
        -l/--min_seq_length_gene\
        -T/--threads\
        -P/--path\
\
similarity_cluster - will cluster sequences of each gene based on pairwise\
                     sequence hamming distances.\
        -d/--database\
        -s/--sim_cut_off\
        -l/--min_seq_length_gene\
        -T/--threads\
        -P/--path\
\
alignment_groups - will define groups to align for each gene based on the\
                   taxon string annotation.\
        -d/--database\
        -G/--group_on_genera\
\
link_genes - will create the table alignments for the aligned sequences.\
        -d/--database\
        -L/--maximize_linking\
        -l/--min_seq_length_gene\
\
multimotu - will link genes using graph theory based on Chesters and Vogler\
            (2013).\
        -d/--database\
        -cc/--comp_columns\
        -ck/--comp_kind\
        -cs/--comp_score\
\
latin_names - will try to find appropriate species names for the ‘taxa’\
              defined by linked_genes.\
        -d/--database\
\
align_sequences - will align the sequences that were linked by link_genes.\
        -d/--database\
        -T/--threads\
        -E/--create_gene_trees\
        -p/--phylo_method\
        -B/--store_boot_trees\
        -Y/--remove_outliers\
        -k/--linked_genes\
        -P/--path\
\
refine_alignments - will refine the alignments created by align_sequences\
                    using an algorithm based on Liu et al (2009 and 2012).\
        -d/--database\
        -p/--phylo_method\
        -B/--store_boot_trees\
        -r/--stop_criterion\
        -ms/--max_alignment_group_size\
        -ug/--use_guide_tree\
        -T/--threads\
        -u/--use_genes\
        -k/--linked_genes\
        -P/--path\
\
gblocks - will remove poorly aligned regions using Gblocks.\
        -d/--database\
        -b0\
        -b1\
        -b2\
        -b3\
        -b4\
        -b5\
        -u/--use_genes\
        -l/--min_seq_length_gene\
        -P/--path\
\
catenated_tree - will produce a tree from a catenated alignment of the given\
                 genes, possibly partitioning the separate genes.\
        -d/--database\
        -c/--cat_phylo_method\
        -T/--threads\
        -B/--store_boot_trees\
        -u/--use_genes\
        -Q/--anchor_gene\
        -P/--path\
\
exclude_rogues - will identify rogue taxa using RogueNaRok.\
        -d/--database\
        -N/--out_file_name\
        -w/--print_rogues\
        -x/--output_no_rogues_tree\
        -z/--rm_rogues_from_ml_tree\
        -q/--rm_rogues_from_alignment\
        -P/--path\
\
get_distinct_entries - will output each unique entry in one or several columns\
                       in the table gb_data.\
        -d/--database\
        -U/--column_name\
        -N/--out_file_name\
\
get_shared_values - will output values in a given column that are also present\
                    in another given column in the table gb_data.\
        -d/--database\
        -U/--column_name\
        -N/--out_file_name\
\
get_alignments - will output a concatenated sequence alignment for the given\
                ‘genes’ in either, fasta, phylip, or nexus format.\
        -d/--database\
        -a/--alignment_format\
        -D/--interleaved\
        -W/ --bases_per_row\
        -A/--partition_file\
        -N/--out_file_name\
        -S/--sequence_name_column\
        -u/--use_genes\
        -Q/--anchor_gene\
\
get_trees - will output one or more trees produced by PifCoSm as stored in\
            PifCoSm.\
        -d/--database\
        -N/--out_file_name\
        -t/--get_taxon\
        -u/--use_genes\
        -O/--only_show\
\
switch_tip_names - will switch the tip names of a tree (and only one tree)\
                   from accession numbers to species names as given in gb_data\
                   or from ‘taxon’ names to species names as given by the\
                   latin_names module.\
        -d/--database\
        -I/--tree_infile\
        -R/--tree_file/--gene_tree\
        -N/--out_file_name\
        -P/--path\
\
get_otu_taxonomy - will give a tab separated list of the taxonomic hierarchy\
                   for each OTU for as long as it agrees for all sequences for\
                   the given genes.\
        -d/--database\
        -u/--use_genes\
        -Q/--anchor_gene\
        -S/--sequence_name_column\
\
stats - will present some statistics about either the clustering ('cluster')\
        or the alignments ('alignment') for given genes.\
        -d/--database\
        --stats_type\
        -u/--use_genes\
        -Q/--anchor_gene (not used for cluster stats)\
        -P/--path\
";

### Global variables ###
my $taxon_tree=0;
my $max_seq_length = 10000; # variable for maximum sequence length
my $batch_file; # variable for batch file name
my $database; # database name
my $gb_data_file = "sequence.gb"; # name of file with gb_data (full genbank format)
my @modules; # array with which modules to ru and in which order
my %backup;
my $hmmdatabase; # name of the HAMMER model database
my $print_non_match = 'n'; # print the accnos that are not sorted to gene in file
my $min_sequence_length_out = 50; # min length of sequences to include in the separate gene databases
my $e_value_cut_off = 0.0000001; # e value cut off for including a sequence in a gene table
my $n_cut_off = 0.1; # cut off for how large fraction of n's we accept for inclusion in gene tables
my $sim_cut_off = 'all,0.99'; # similarity cut of for sequences to cluster togather
my $sim_tree_clust = 'y'; # should tree similarity clust be performed (if tree clust is performed)
my $species_tree_clust = 'y'; # should monophyletic species be clustered if tree clust is performed
my %intergenes; # hash to define genes situated between two genes (hmms) e.g. biven like ("ITS1", "nucSSU;nuc58","ITS2","nuc58;nucLSU","IGS","nucLSU;nucSSU");
my @individual_criteria; #= ("voucher", "isolate", "strain");
my $individual_condition; #= "species";
my $min_group = 'y';
my $create_gene_tree = 'y';
my $tree_method = 'RAxML';
my $rm_outliers = 'y';
my $final_tree_method = 'RAxML_partitioned_bootstrap_100';
my @linked_genes; # array to store linked genese eg. 'nuclsu;nucssu;nuc58;its1;its2' or 'mtlsu;mtssu'
my $alignment_format = 'phylip';
my $get_partitions = 'y';
my @gene_alignments = ('all');
my $anchor_gene = 'n';
my @changes_files;
my $taxon = 'combined';
my $align_seq_name_col = 'taxon_name';
my $out_file_stem = 'out';
my $input_tree_file = 'n';
my $max_size = 200;
my $use_guide_tree = 'n';
my $stop_criterion = "iterations_3";
my $tree_file = 'combined,combined';
my $n_threads = 0;
my $path = '';
my $store_boot_trees = 'n';
my @column_name;
my $interleaved = 'n';
my $bases_per_row = 1000;
my $print_batch_file = 'n';
my $only_show = 'n';
my $maximum_linking = 'n'; # sets if taxa should be defined by linkage at a single gene (the alternative is that they are formed to include all (not really) lead sequences)
# variables for multimotu
my $comp_columns='species';
my $comp_type='h';
my $comp_score='2';

# Variables for gblocks
my $b0 = 0;
my $b1 = 0;
my $b2 = 0;
my $b3 = 0;
my $b4 = 0;
my $b5 = 'n';
# Variables for the exclude_rogues module
my $print_taxa = 'y';
my $output_tree = 'n';
my $remove_from_trees = 'n';
my $remove_from_alignment = 'n';

# variables for change_entries
my $change_option = 'c';

my $stats_type = 'alignment';
########################

######################
### Main of script ###
######################
### Parsing arguments ###
{
    if (scalar @ARGV == 0) {
        print $usage;
        exit;
    }
    my $i;
    for ($i=0; $i < scalar @ARGV; ++$i) {
        if ($ARGV[$i] eq "-h" or $ARGV[$i] eq "--help") {
           print $usage;
           print $help;
           exit;
        }
        if ($ARGV[$i] eq "-f" or $ARGV[$i] eq "--batch_file" or scalar @ARGV == 1) {
           if (@ARGV == 1) { $batch_file = $ARGV[0]; }
           elsif (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $batch_file = $ARGV[++$i]; }
           &pars_batch_file($batch_file);
        }
    }
    for ($i=0; $i < scalar @ARGV; ++$i) {
        if ($ARGV[$i] eq "-f" or $ARGV[$i] eq "--batch_file") {
           if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $batch_file = $ARGV[++$i]; }
           next;
        }
        elsif ($ARGV[$i] eq "-d" or $ARGV[$i] eq "--database") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $database = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-g" or $ARGV[$i] eq "--gb_file") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $gb_data_file = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-o" or $ARGV[$i] eq "--modules") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { @modules = split /,/, $ARGV[++$i]; }
            else { die "-o / --modules need a ',' separated string with wich modules to run, e.g. -o gb_parser,gene_parser.\n"; }
        }
        elsif ($ARGV[$i] eq "-M" or $ARGV[$i] eq "--hmm_database") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $hmmdatabase = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-Z" or $ARGV[$i] eq "--print_non_match") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $print_non_match = $ARGV[++$i]; }
            if ($print_non_match =~ /^(y|Y)/) { $print_non_match = 'y'; }
            else { $print_non_match = 'n'; }
        }
        elsif ($ARGV[$i] eq "-v" or $ARGV[$i] eq "--individual_columns") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { @individual_criteria = split /,/, $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-V" or $ARGV[$i] eq "--individual_condition") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $individual_condition = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-i" or $ARGV[$i] eq "--inter_genes") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { %intergenes = split /,/, $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-b" or $ARGV[$i] eq "--backup") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                my @temp = split /,/,$ARGV[++$i];
                foreach (@temp) { $backup{$_} = 'y'; }
            }
            else { $backup{"all"} = 'y'; }
        }
        elsif ($ARGV[$i] eq "-e" or $ARGV[$i] eq "--e_value_cut_off") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $e_value_cut_off = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-n" or $ARGV[$i] eq "--n_cut_off") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $n_cut_off = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-s" or $ARGV[$i] eq "--sim_cut_off") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $sim_cut_off = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-m" or $ARGV[$i] eq "--max_seq_length") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $max_seq_length = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-l" or $ARGV[$i] eq "--min_seq_length_gene") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $min_sequence_length_out = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-C" or $ARGV[$i] eq "--tree_clust_type") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^b/) {
                    $species_tree_clust = 'y';
                    $sim_tree_clust = 'y';
                }
                elsif ($ARGV[$i] =~ /^sp/) {
                    $species_tree_clust = 'y';
                    $sim_tree_clust = 'n';
                }
                elsif ($ARGV[$i] =~ /^si/) {
                    $species_tree_clust = 'n';
                    $sim_tree_clust = 'y';
                }
            }
        }
        elsif ($ARGV[$i] eq "-G" or $ARGV[$i] eq "--group_on_genera") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $min_group = $ARGV[++$i]; }
            if ($min_group =~ /^(y|Y)/) { $min_group = 'y'; }
            else { $min_group = 'n'; }
        }
        elsif ($ARGV[$i] eq "-cc" or $ARGV[$i] eq "--comp_columns") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $comp_columns = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-ck" or $ARGV[$i] eq "--comp_kind") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $comp_type = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-cs" or $ARGV[$i] eq "--comp_score") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $comp_score = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-p" or $ARGV[$i] eq "--phylo_method") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $tree_method = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-c" or $ARGV[$i] eq "--concat_phylo_method") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $final_tree_method = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-r" or $ARGV[$i] eq "--stop_criterion") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $stop_criterion = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-ms" or $ARGV[$i] eq "--max_alignment_group_size") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $max_size = $ARGV[++$i]; }
            if ($max_size =~ /[^0-9]/) { die "-ms/--max_alignment_group_size must be a integer value.\n"; }
        }
        elsif ($ARGV[$i] eq "-ug" or $ARGV[$i] eq "--use_guide_tree") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $use_guide_tree = $ARGV[++$i]; }
            if ($use_guide_tree =~ /^n/i) { $use_guide_tree = 'n'; }
            elsif ($use_guide_tree =~ /^n/i) { $use_guide_tree = 'y'; }
            else { die "-ug/--use_guide_tree must be followed by y/n (yes or no).\n"; }
        }
        elsif ($ARGV[$i] eq "-a" or $ARGV[$i] eq "--alignment_format") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $alignment_format = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-D" or $ARGV[$i] eq "--interleaved") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ( $ARGV[$i] =~ /^(y|Y)/ ) { $interleaved = 'y'; }
                elsif ( $ARGV[$i] =~ /^(n|N)/ ) { $interleaved = 'n'; }
                else { die "Do not recognize argument given after $ARGV[$i-1], it should be yes or no.\n"; }
            }
            else { $interleaved = 'y'; }
        }
        elsif ($ARGV[$i] eq "-W" or $ARGV[$i] eq "--bases_per_row") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { 
                $bases_per_row = $ARGV[++$i];
                if ($bases_per_row =~ /^gene/i) {
                    $bases_per_row = 'genes';
                }
                elsif ($bases_per_row =~ /[^0-9]/) { die "$ARGV[$i-1] may only be followed by the key word genes or a integer number. Try again.\n"; }
            }
            else { die "$ARGV[$i] should be folowed by the number of bases (e.g. $ARGV[$i] 100) or the key word genes (for each gene in separate blocks).\n"; }
        }
        elsif ($ARGV[$i] eq "-A" or $ARGV[$i] eq "--partition_file") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $get_partitions = $ARGV[++$i]; }
            if ($get_partitions =~ /^(y|Y)/) { $get_partitions = 'y'; }
            else { $get_partitions = 'n'; }
        }
        elsif ($ARGV[$i] eq "-u" or $ARGV[$i] eq "--use_genes") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { @gene_alignments = split /,/,$ARGV[++$i]; }
            else { die "--use_genes/-u require a comma separated string of genes, or the key word all, e.g. rpb1,rpb2,nuclsu,nuc58.\n"; }
        }
        elsif ($ARGV[$i] eq "-Q" or $ARGV[$i] eq "--anchor_gene") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $anchor_gene = $ARGV[++$i]; $anchor_gene =~ s/\s//g; }
            else { $anchor_gene = 'n'; }
        }
        elsif ($ARGV[$i] eq "-L" or $ARGV[$i] eq "--maximize_linking") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ( $ARGV[$i] =~ /^(y|Y)/ ) { $maximum_linking = 'y'; }
                elsif ( $ARGV[$i] =~ /^(n|N)/ ) { $maximum_linking = 'n'; }
                else { die "Do not recognize argument given after $ARGV[$i-1], it should be yes or no.\n"; }
            }
            else { $maximum_linking = 'y'; }
        }
        elsif ($ARGV[$i] eq "-t" or $ARGV[$i] eq "--get_taxon") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $taxon = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-S" or $ARGV[$i] eq "--sequence_name_column") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $align_seq_name_col = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-N" or $ARGV[$i] eq "--out_file_name") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $out_file_stem = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-H" or $ARGV[$i] eq "--changes_files") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { @changes_files = split /,/,$ARGV[++$i]; }
            else { die "--changes_files/-H require a file name or a comma separated string of file names e.g. changes_species.txt,changes_voucher.txt.\n"; }
        }
        elsif ($ARGV[$i] eq "-R" or $ARGV[$i] eq "--tree_file" or $ARGV[$i] eq "--gene_tree") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $tree_file = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-I" or $ARGV[$i] eq "--tree_infile") {
            $input_tree_file = 'y';
        }
        elsif ($ARGV[$i] eq "-k" or $ARGV[$i] eq "--linked_genes") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { @linked_genes = split /,/,$ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-B" or $ARGV[$i] eq "--store_boot_trees") {
            $store_boot_trees='y';
        }
        elsif ($ARGV[$i] eq "-P" or $ARGV[$i] eq "--path") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $path = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-U" or $ARGV[$i] eq "--column_name") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { @column_name = split /,/,$ARGV[++$i]; }
            else { die "--column_name/-U require a column name or comma separated string of column names, e.g. species,voucher.\n"; }
        }
        elsif ($ARGV[$i] eq "-O" or $ARGV[$i] eq "--only_show") { 
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $only_show = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $only_show = 'n'; }
            }
            else { $only_show = 'y'; }
        }
        elsif ($ARGV[$i] eq "-E" or $ARGV[$i] eq "--create_gene_trees") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $create_gene_tree = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $create_gene_tree = 'n'; }
            }
            else { $create_gene_tree = 'y'; }
        }
        elsif ($ARGV[$i] eq "-Y" or $ARGV[$i] eq "--remove_outliers") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $rm_outliers = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $rm_outliers = 'n'; }
            }
            else { $rm_outliers = 'y'; }
        }
        elsif ($ARGV[$i] eq "-ch" or $ARGV[$i] eq "--change_option") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^remove$/i || $ARGV[$i] =~ /^rm$/i) { $change_option = 'r'; }
                elsif ($ARGV[$i] =~ /^change$/i) { $change_option = 'c'; }
		elsif ($ARGV[$i] =~ /^update$/i) { $change_option = 'u'; }
		else { die "Unknown option $ARGV[$i] to -ch/--change_option.\n"; }
            }
            else { $rm_outliers = 'y'; }
        }
        elsif ($ARGV[$i] eq "-w" or $ARGV[$i] eq "--print_rogues") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $print_taxa = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $print_taxa = 'n'; }
                else { die "Invalid argument given for -x/--print_rogues, the alternatives are y/yes/n/no.\n"; }
            }
            else { $print_taxa = 'y'; }
        }
        elsif ($ARGV[$i] eq "-x" or $ARGV[$i] eq "--output_no_rogues_tree") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $output_tree = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $output_tree = 'n'; }
                else { die "Invalid argument given for -q/--output_no_rogues_tree, the alternatives are y/yes/n/no.\n"; }
            }
            else { $output_tree = 'y'; }
        }   
        elsif ($ARGV[$i] eq "-z" or $ARGV[$i] eq "--rm_rogues_from_ml_tree") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $remove_from_trees = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $remove_from_trees = 'n'; }
                else { die "Invalid argument given for -z/--rm_rogues_from_ml_tree, the alternatives are y/yes/n/no.\n"; }
            }
            else { $remove_from_trees = 'y'; }
        }
        elsif ($ARGV[$i] eq "-q" or $ARGV[$i] eq "--rm_rogues_from_alignment") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $remove_from_alignment = 'y'; }
                elsif ($ARGV[$i] =~ /^n/i) { $remove_from_alignment = 'n'; }
                else { die "Invalid argument given for -q/--rm_rogues_from_alignment, the alternatives are y/yes/n/no.\n"; }
            }
            else { $remove_from_alignment = 'y'; }
        }
        elsif ($ARGV[$i] eq "-b0") {
            if (($i+1) < scalar @ARGV and $ARGV[$i+1] =~ /^[0-9]+/) {
                ++$i;
                $b0 = $ARGV[$i];
            }
            else { die "-b0 should be followed by a integer value.\n"; }
        }
        elsif ($ARGV[$i] eq "-b1") {
            if (($i+1) < scalar @ARGV and $ARGV[$i+1] =~ /^[0-9]+/) {
                ++$i;
                $b1 = $ARGV[$i];
            }
            else { die "-b1 should be followed by a integer value.\n"; }
        }
        elsif ($ARGV[$i] eq "-b2") {            
            if (($i+1) < scalar @ARGV and $ARGV[$i+1] =~ /^[0-9]+/) {
                ++$i;
                $b2 = $ARGV[$i];
            }
            else { die "-b2 should be followed by a integer value.\n"; }
        }
        elsif ($ARGV[$i] eq "-b3") {            
            if (($i+1) < scalar @ARGV and $ARGV[$i+1] =~ /^[0-9]+/) {
                ++$i;
                $b3 = $ARGV[$i];
            }
            else { die "-b3 should be followed by a integer value.\n"; }
        }
        elsif ($ARGV[$i] eq "-b4") {            
            if (($i+1) < scalar @ARGV and $ARGV[$i+1] =~ /^[0-9]+/) {
                ++$i;
                $b4 = $ARGV[$i];
            }
            else { die "-b4 should be followed by a integer value.\n"; }
        }
        elsif ($ARGV[$i] eq "-b5") {            
            if (($i+1) < scalar @ARGV and $ARGV[$i+1] =~ /^(n|h|a)$/) {
                ++$i;
                $b5 = $ARGV[$i];
            }
            else { die "-b5 should be followed by n,h, or a.\n"; }
        }
        elsif ($ARGV[$i] eq "-T" or $ARGV[$i] eq "--threads") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) { $n_threads = $ARGV[++$i]; }
        }
        elsif ($ARGV[$i] eq "-X" or $ARGV[$i] eq "--print_batch_file") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                if ($ARGV[$i] =~ /^y/i) { $print_batch_file = 'y'; }
            }
            else { $print_batch_file = 'y'; }
        }
        elsif ($ARGV[$i] eq "--stats_type") {
            if (($i+1) < scalar @ARGV and !($ARGV[$i+1] =~ /^-/)) {
                ++$i;
                $stats_type=$ARGV[$i];
                if ($stats_type ne 'alignment' and $stats_type ne 'cluster') { die "Illeagal value for --stats_type.\n"; }
            }
            else {
                print STDERR "--stats_type should be followed by either 'alignment' or 'cluster' to give what type of statistics to present. No value given. Will produce 'alignment' stats.\n";
            }
        }
        else { # if non of above switches
            # if two arguments interpret them as batch file followed by GenBank data file
            if (scalar @ARGV == 2 and !($ARGV[0] =~ /^-/) and !($ARGV[1] =~ /^-/)) {
                $batch_file = $ARGV[0];
                $gb_data_file = $ARGV[1];
            }
            # if one argument interprete it as a batch file
            elsif (scalar @ARGV == 1 and !($ARGV[0] =~ /^-/)) {
                $batch_file = $ARGV[0];
            }
            # if more arguments something is wrong
            else { die "Do not recognize $ARGV[$i].\n"; }
        }
    }
}
### Tests ###
{
    if (!$database) { die "A database name must be provided by -d/--database or if using a batch file database.\n"; }
    if (@modules) {
        my %run_module;
        my $gb_data_check = 'n';
        for (my $i=0; $i < scalar @modules; ++$i) {
            if ($modules[$i] eq "gb_parser") {
                $run_module{"gb_parser"} = 'y';
                if ($max_seq_length =~ /[^0-9]/g) { die "Invalid max sequence length (--max_seq_length). Will not be able to parse GenBank file. Quitting.\n"; }
                elsif (!$database) { die "No database name available. Will not be able to parse GenBank file. Quitting.\n"; }
                elsif (!(-e $gb_data_file)) { die "Could not find file for GenBank data ($gb_data_file). Quitting.\n"; }
            }
            else {
                if (!$run_module{"gb_parser"} and $gb_data_check eq 'n') {
                    my $dbh = PifCosm_support_subs::connect_to_database ($database);
                    if (PifCosm_support_subs::table_present($dbh,'gb_data') eq 'n') { die "There must be a table gb_data in $database or you must run the module gb_parser before any other module.\n"; }
                    $dbh->disconnect();
                    $gb_data_check = 'y';
                }
		if ($modules[$i] eq "check_taxonomy") { $run_module{"check_taxonomy"} = 'y'; }
                elsif ($modules[$i] eq "change_entries") {
                    $run_module{"change_entries"} = 'y';
                    if (@changes_files) {
                        foreach (@changes_files) { if (!(-e $_)) { die "Could not find the file $_ for the change_entries module. Quitting.\n"; } }
                    }
                    else { die "change_entries module require file/files with information on what fields to change and in what column, see manual.\n"; }
                }
                elsif ($modules[$i] eq "define_individuals") {
                    $run_module{"define_individuals"} = 'y';
                    if (@individual_criteria) {
                        my @test_columns = @individual_criteria;
                        if ($individual_condition) { push (@test_columns, $individual_condition); }
                        my $checked = 'n';
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        if ($gb_data_check eq 'y' or ($run_module{"gb_parser"} and PifCosm_support_subs::table_present($dbh,'gb_data') eq 'y')) {
                            foreach(@test_columns) { if (PifCosm_support_subs::column_present($dbh,'gb_data',$_) eq 'n') { die "The column $_ is not present in the table gb_data in the database $database. Use other column to define individual (-v and -V).\n"; } }
                            $checked = 'y';
                        }
                        $dbh->disconnect();
                        if ($run_module{"gb_parser"} and $checked eq 'n') {
                            foreach(@test_columns) {
                                if (!($_ eq 'accno' or $_ eq 'species'  or $_ eq 'sub_species' or $_ eq 'variety' or $_ eq 'voucher' or $_ eq 'isolate' or $_ eq 'strain' or $_ eq 'taxon_string' or $_ eq 'country' or $_ eq 'collected_by' or $_ eq 'identified_by' or $_ eq 'environmental_sample' or $_ eq 'host' or $_ eq 'latitude' or $_ eq 'longitude' or $_ eq 'reference' or $_ eq 'sequence' or $_ eq 'proportion_N')) {
                                    die "The column $_ will not be created in the table gb_data by gb_parser. Use other column to define individual (-v).\n";
                                }
                            }
                        }
                    }
                    else { die "Could not find the names of the columns to use for the define_individuals module (given by -v). Quitting.\n" }
                }
                elsif ($modules[$i] eq "gene_parser") {
                    $run_module{"gene_parser"} = 'y';
                    if (!(-e "$hmmdatabase")) { die "Could not find the HAMMER database $hmmdatabase needed for the gene_parser module.\n"; }
                    if ($e_value_cut_off =~ /[^e0-9\.-]/) { die "The e value cut off value (-e/--e_value_cut_off) must be a real number and may not contain other characters.\n"; }
                    if ($n_cut_off =~ /[^e0-9\.-]/) { die "The cut of value for proportion N's in a sequence (-n/--n_cut_off) must be a real number and may not contain other characters.\n"; }
                    if (!($print_non_match eq 'y' or $print_non_match eq 'n')) { die "The flag for printing accnos of sequence that lack match to gene to a file (-Z/--print_non_match) must be either y or n.\n"; }
                    if (%intergenes) {
                        foreach (keys %intergenes) {
                            if (!$intergenes{$_} or !($intergenes{$_} =~ /;/)) {
                                print STDERR "Intergenes (-i/--inter_genes) must be given as 'gene,leding_gene;trailing_gene'";
                                if ($intergenes{$_}) { die ", not '$_,$intergenes{$_}'.\n"; }
                                else { die ".\n"; }
                            }
                        }
                    }
                    if ($min_sequence_length_out =~ /[^0-9]/g) { die "$min_sequence_length_out is an invalid value for minimum sequence length for sorting into genes (--min_seq_length_gene). Will not be able to sort genes. Quitting.\n"; }
                }
                elsif ($modules[$i] eq "cluster_individuals") {
                    $run_module{"cluster_individuals"} = 'y';
                    my $dbh = PifCosm_support_subs::connect_to_database ($database);
                    if (!$run_module{"define_individuals"}) {
                        if ($gb_data_check eq 'y' or ($run_module{"gb_parser"} and PifCosm_support_subs::table_present($dbh,'gb_data') eq 'y')) {
                            if (PifCosm_support_subs::column_present($dbh,'gb_data','individual') eq 'n') { die "Could not find the column individual in table gb_data, you need to run the module define_individuals, before you can cluster_individuals.\n"; }
                        }
                        else {
                            die "You need to run the module define_individuals before you can cluster individuals.\n";
                        }
                    }
                    if ($min_sequence_length_out =~ /[^0-9]/g) { die "$min_sequence_length_out is an invalid value for minimum sequence length for clustering (--min_seq_length_gene).\n"; } 
                }
                elsif ($modules[$i] eq "tree_clustering" or $modules[$i] eq "similarity_cluster" or $modules[$i] eq "CD_hit_cluster") {
                    $run_module{"tree_clustering"} = 'y';
                    $run_module{"similarity_cluster"} = 'y';
                    $run_module{"CD_hit_cluster"} = 'y';
                    my $dbh = PifCosm_support_subs::connect_to_database ($database);
                    my %test_cutoff=split /,/,$sim_cut_off;
                    foreach (keys %test_cutoff) {
                        if (!(($_ =~ /^all$/i) or ($_ =~ /^default$/i)) and !$run_module{"gene_parser"} and PifCosm_support_subs::table_present($dbh,$_) eq 'n') { die "Do not recognize gene $_ for similarity cut off value (-s/--sim_cut_off).\n"; }
                        if ($test_cutoff{$_} =~ /[^0-9\.]/) { die "$test_cutoff{$_} is an invalid similarity cut off value (-s/--sim_cut_off) for gene $_.\n"; }
                    }
                    if ($min_sequence_length_out =~ /[^0-9]/g) { die "$min_sequence_length_out is an invalid value for minimum sequence length for clustering (--min_seq_length_gene).\n"; }
                    if ($n_threads =~ /[^0-9]/g) { die "$n_threads is a invalid value for the numberof threads (-T/--threads).\n"; }
                }
                elsif ($modules[$i] eq "alignment_groups") {
                    $run_module{"alignment_groups"} = 'y';
                    if ( !($min_group eq 'y' or $min_group eq 'n')) { die "Error in what groups to align (-G/--group_on_genera).\n"; }
                }
                elsif ($modules[$i] eq "link_genes") {
                    $run_module{"link_genes"} = 'y';
                    if ($min_sequence_length_out =~ /[^0-9]/g) { die "$min_sequence_length_out is an invalid value for minimum sequence length for link_genes module (--min_seq_length_gene).\n"; }
                }
                elsif ($modules[$i] eq "multimotu") {
                    $run_module{"multimotu"} = 'y';
                    my $dbh = PifCosm_support_subs::connect_to_database ($database);
                    my @comp_column = split /,/,$comp_columns;
                    foreach (@comp_column) {
                        if (PifCosm_support_subs::column_present($dbh,'gb_data',$_) eq 'n') { die "Could not find the column $_ in table gb_data. Give other column with -cc/--comp_columns.\n"; }
                        if ($_ eq 'accno') { $_ = 'gb_data.accno'; }
                    }
                    $comp_columns = join ',', @comp_column;
                    my @comp_type = split /,/,$comp_type;
                    if ($comp_type=~ /[^oaeh,]/) { die "Do not recognize the kind of comparison (-ck/--comp_kind). a, o, e, and h are the only valid options. Each kind should be separated by a comma and there should be no white space.\n"; }
                    my @comp_score = split /,/,$comp_score;
                    if ($comp_score=~/[^0-9,\.]/) { die "The comparison score (-cs/--comp_score) may only contain numbers separated by commas and no white space.\n"; }
                    if (scalar @comp_column > scalar @comp_type || scalar @comp_column > scalar @comp_score) {
                        die "Each column that is compared (-cc/--comp_columns) by multimotu need to have a type (-ck/--comp_kind) of comparison and a score (-cs/--comp_score).\n";
                    }
                }
                elsif ($modules[$i] eq "latin_names" or $modules[$i] eq "align_sequences" or $modules[$i] eq "refine_alignments" or $modules[$i] eq "concatenated_tree" or $modules[$i] eq "get_alignments" or $modules[$i] eq "gblocks") {
                    if (!($run_module{"link_genes"} or $run_module{"multimotu"})) {
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        if (PifCosm_support_subs::table_present($dbh,'alignments') eq 'n') { die "There must be a table alignments in $database or you must run the module link_genes before any of the module latin_names, align_sequences, refine_alignments, concatenated_tree, get_otu_taxonomy, and get_alignments.\n"; }
                        $dbh->disconnect();
                    }
                    if ($modules[$i] eq "align_sequences") {
                        $run_module{"align_sequences"} = 'y';
                        if ($n_threads =~ /[^0-9]/g) { die "$n_threads is a invalid value for the numberof threads (-T/--threads).\n"; }
                        if (!($store_boot_trees eq 'y' or $store_boot_trees eq 'n')) { die "$store_boot_trees is not a valid option for if bootstrap trees should be stored or not.\n"; }
			my $test_method = $tree_method;
			$test_method =~ s/_bootstrap_[0-9]*//;
                        if (!($test_method eq "ExaML" or $test_method eq "ExaML,RAx-br" or $test_method eq "RAxML" or $test_method eq "fasttree")) { die "$tree_method is not a recognized phylogenetic method (-p/--phylo_method).\n"; }
                        if ($rm_outliers eq 'y' and $create_gene_tree eq 'n') { die "Must create gene trees (-E/--create_gene_trees set to yes) to define which outliers to remove (-Y/--remove_outliers).\n"; }
                        foreach (@linked_genes) { if ($_ =~ /[^a-zA-Z0-9_;]/) { die "The linked gene annotation $_ contain illegal characters.\n"; }}
                    }
                    elsif ($modules[$i] eq "refine_alignments") {
                        $run_module{"refine_alignments"} = 'y';
                        if ($n_threads =~ /[^0-9]/g) { die "$n_threads is a invalid value for the numberof threads (-T/--threads).\n"; }
                        if (!($store_boot_trees eq 'y' or $store_boot_trees eq 'n')) { die "$store_boot_trees is not a valid option for if bootstrap trees should be stored or not.\n"; }
                        if (!($tree_method eq "ExaML" or $tree_method eq "ExaML,RAx-br" or $tree_method eq "RAxML" or $tree_method eq "fasttree")) { die "$tree_method is not a recognized phylogenetic method (-p/--phylo_method).\n"; }
                        foreach (@linked_genes) { if ($_ =~ /[^a-zA-Z0-9_;]/) { die "The linked gene annotation $_ contain illegal characters.\n"; }}
                        if (!($stop_criterion =~ /^iterations_[0-9]+/ or $stop_criterion =~ /^max_[0-9]+/)) { die "$stop_criterion is not a valid stop criterion when refining the alignments (-r/--stop_criterion).\n"; }
                        if (!($gene_alignments[0] eq 'all' or $gene_alignments[0] =~ /^max_block_[0-9]+/)) {
                            if (!$run_module{"gene_parser"}) {
                                my $dbh = PifCosm_support_subs::connect_to_database ($database);
                                foreach (@gene_alignments) { if (PifCosm_support_subs::table_present($dbh,$_) ne 'y') { die "The gene $_ is not present in $database and can not be used, remove or give other gene to use (-u/--use_genes).\n"; }}
                                $dbh->disconnect();
                            }
                        }
                    }
                    elsif ($modules[$i] eq "gblocks") {
                        $run_module{"gblocks"} = 'y';
                        if (!$run_module{"align_sequences"}) {
                            my $dbh = PifCosm_support_subs::connect_to_database ($database);
                            if (PifCosm_support_subs::table_present($dbh,"alignments") ne 'y') { die "\n"}
                            $dbh->disconnect();
                        }
                        if ($b0 =~ /[^0-9]/) { die "Minimum Length Of An Initial Block must be an integer. $b0 is not a valid number, give other number by -b0.\n"; }
                        if ($b1 =~ /[^0-9]/) { die "Minimum Number Of Sequences For A Conserved Position must be an integer. $b1 is not a valid number, give other number by -b1.\n"; }
                        if ($b2 =~ /[^0-9]/) { die "Minimum Number Of Sequences For A Flank Position must be an integer. $b2 is not a valid number, give other number by -b2.\n"; }
                        if ($b3 =~ /[^0-9]/) { die "Maximum Number Of Contiguous Nonconserved Positions must be an integer. $b3 is not a valid number, give other number by -b3.\n"; }
                        if ($b4 =~ /[^0-9]/) { die "Minimum Length Of A Block must be an integer. $b0 is not a valid number, give other number by -b4.\n"; }
                        if ($b5 =~ /[^nha]/) { die "Allowed Gap Positions must be n (none), h (with half), or a (all). $b5 is not a valid option, give other option by -b5.\n"; }
                        if ($min_sequence_length_out =~ /[^0-9]/) { die "$min_sequence_length_out is not a valid minimum sequence length. The minimum sequence length, given by -l/--min_seq_length_gene, must be an integer number.\n"; }
                        if (!($gene_alignments[0] eq 'all' or $gene_alignments[0] =~ /^max_block_[0-9]+/)) {
                            if (!$run_module{"gene_parser"}) {
                                my $dbh = PifCosm_support_subs::connect_to_database ($database);
                                foreach (@gene_alignments) { if (PifCosm_support_subs::table_present($dbh,$_) ne 'y') { die "The gene $_ is not present in $database and can not be used, remove or give other gene to use (-u/--use_genes).\n"; }}
                                $dbh->disconnect();
                            }
                        }
                    }
                    elsif ($modules[$i] eq "concatenated_tree") {
                        $run_module{"concatenated_tree"} = 'y';
                        if (!($final_tree_method =~ /(ExaML)|(ExaML,RAx-br)|(RAxML)|(fasttree)/i)) {die "Could not find viable phylogenetic method in '$final_tree_method', give other method (-c/--concat_phylo_method).\n"; }
                        if ($n_threads =~ /[^0-9]/g) { die "$n_threads is a invalid value for the numberof threads (-T/--threads).\n"; }
                        if (!($store_boot_trees eq 'y' or $store_boot_trees eq 'n')) { die "$store_boot_trees is not a valid option for if bootstrap trees should be stored or not.\n"; }
                        if (!($gene_alignments[0] eq 'all' or $gene_alignments[0] =~ /^max_block_[0-9]+/)) {
                            if (!$run_module{"gene_parser"}) {
                                my $dbh = PifCosm_support_subs::connect_to_database ($database);
                                foreach (@gene_alignments) { if (PifCosm_support_subs::table_present($dbh,$_) ne 'y') { die "The gene $_ is not present in $database and can not be used, remove or give other gene to use (-u/--use_genes).\n"; }}
                                $dbh->disconnect();
                            }
                        }
                        if (!$run_module{"link_genes"}) {
                            my $dbh = PifCosm_support_subs::connect_to_database ($database);
                            if ($anchor_gene ne 'n') { foreach (split /,/, $anchor_gene) { if (PifCosm_support_subs::column_present ($dbh,"alignments","$_\_sequence") ne 'y') { die "Could not find $_\_sequence as column name in table alignments in $database.\n" } } }
                            $dbh->disconnect();
                        }
                    }
                    elsif ($modules[$i] eq "get_alignments") {
                        $run_module{"get_alignments"} = 'y';
                        if (!($alignment_format eq 'phylip' or $alignment_format eq 'nexus' or $alignment_format eq 'fasta')) { die "Do not recognize $alignment_format as a valid alignment format (-a/--alignment_format).\n" }
                        if (!($get_partitions eq 'y' or $get_partitions eq 'n')) { die "Can not determin if $get_partitions means that partition file/block should be created or not (-A/--partition_file).\n"; }
                        if ($out_file_stem =~ /\W/) { die "$out_file_stem is not valid as part of a file name (-N/--out_file_name).\n"; }
                        if (!$run_module{"link_genes"}) {
                            my $dbh = PifCosm_support_subs::connect_to_database ($database);
                            if (PifCosm_support_subs::column_present ($dbh,"alignments",$align_seq_name_col) ne 'y') { die "Could not find $align_seq_name_col as column name in table alignments in $database (-S/--sequence_name_column).\n" }
                            if ($anchor_gene ne 'n') { foreach (split /,/, $anchor_gene) { if (PifCosm_support_subs::column_present ($dbh,"alignments","$_\_sequence") ne 'y') { die "Could not find $_\_sequence as column name in table alignments in $database.\n" } } }
                            $dbh->disconnect();
                        }
                        if (!($gene_alignments[0] eq 'all' or $gene_alignments[0] =~ /^max_block_[0-9]+/)) {
                            if (!$run_module{"gene_parser"}) {
                                my $dbh = PifCosm_support_subs::connect_to_database ($database);
                                foreach (@gene_alignments) { if (PifCosm_support_subs::table_present($dbh,$_) ne 'y') { die "The gene $_ is not present in $database and can not be used, remove or give other gene to use (-u/--use_genes).\n"; }}
                                $dbh->disconnect();
                            }
                        }
                    }
		    if ($modules[$i] eq "get_otu_taxonomy") {
			$run_module{"get_otu_taxonomy"} = 'y';
			if (!$run_module{"link_genes"}) {
                            my $dbh = PifCosm_support_subs::connect_to_database ($database);
                            if (PifCosm_support_subs::column_present ($dbh,"alignments",$align_seq_name_col) ne 'y') { die "Could not find $align_seq_name_col as column name in table alignments in $database (-S/--sequence_name_column).\n" }
                            if ($anchor_gene ne 'n') { foreach (split /,/, $anchor_gene) { if (PifCosm_support_subs::column_present ($dbh,"alignments","$_\_sequence") ne 'y') { die "Could not find $_\_sequence as column name in table alignments in $database.\n" } } }
                            $dbh->disconnect();
                        }
			if (!($gene_alignments[0] eq 'all' or $gene_alignments[0] =~ /^max_block_[0-9]+/)) {
                            if (!$run_module{"gene_parser"}) {
                                my $dbh = PifCosm_support_subs::connect_to_database ($database);
                                foreach (@gene_alignments) { if (PifCosm_support_subs::table_present($dbh,$_) ne 'y') { die "The gene $_ is not present in $database and can not be used, remove or give other gene to use (-u/--use_genes).\n"; }}
                                $dbh->disconnect();
                            }
                        }
		    }
                }
                elsif ($modules[$i] eq "exclude_rogues") {
                    $run_module{"exclude_rogues"} = 'y';
                    if (!($run_module{"similarity_cluster"} or $run_module{"alignment_groups"})) {
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        if (PifCosm_support_subs::table_present ($dbh,"alignment_groups") ne 'y') { die "The module exclude_rogues need bootstrap trees, but the alignment_groups table where these should be stored could not be found. Run the similarity_cluster or alignment_groups modules first, followed by concatenated_tree.\n"; }
                        $dbh->disconnect();
                    }
                    if (!$run_module{"concatenated_tree"}) {
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        my $sth = $dbh->prepare("SELECT gene FROM alignment_groups WHERE taxon='combined_boot'") or die;
                        $sth->execute();
                        my $gene = $sth->fetchrow_array();
                        if (!$gene or $gene ne 'combined_boot') { die "Could not find bootstrap trees. The exclude_rogues module need bootstrap trees, run the catenate_tree module first and save bootstrap trees.\n"; }
                        $sth->finish();
                        $dbh->disconnect();
                    }
                    elsif ($store_boot_trees eq 'n' or !($final_tree_method =~ /bootstrap/)) {
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        my $sth = $dbh->prepare("SELECT gene FROM alignment_groups WHERE taxon='combined_boot'") or die;
                        $sth->execute();
                        my $gene = $sth->fetchrow_array();
                        if (!$gene or $gene ne 'combined_boot') { die "Could not find bootstrap trees. The exclude_rogues module need bootstrap trees. You need to save the bootstrap trees when running the catenate_tree module.\n"; }
                    }
                }
                elsif ($modules[$i] eq "get_distinct_entries" || $modules[$i] eq "get_shared_values") {
                    #$run_module{"get_distinct_entries"} = 'y';
                    if (!$run_module{"gb_parser"}) {
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        foreach (@column_name) { if (PifCosm_support_subs::column_present ($dbh,"gb_data",$_) ne 'y') { die "Column $_ not present in table gb_data in $database. Try other column (-U/--column_name).\n" } }
                        $dbh->disconnect();
                    }
                    if ($out_file_stem =~ /\W/) { die "$out_file_stem is not valid as part of a file name (-N/--out_file_name).\n"; }
                }
                elsif ($modules[$i] eq "get_trees") {
                    $run_module{"get_trees"} = 'y';
                    if (!$run_module{"similarity_cluster"} and !$run_module{"alignment_groups"}) {
                        my $dbh = PifCosm_support_subs::connect_to_database ($database);
                        if (PifCosm_support_subs::table_present($dbh,"alignment_groups") ne 'y') { die "Could not find the table alignment_groups so can not get trees (module get_trees).\n Run the similarity_cluster or alignment_groups cluster first, and at least one module creating trees.\n"; }
                    }
                    if ($out_file_stem =~ /\W/) { die "$out_file_stem is not valid as part of a file name (-N/--out_file_name).\n"; }
                    if (!$taxon and $taxon =~ /\W/) { die "The taxon $taxon contains illegal characters. Try another taxon name (-t/--get_taxon).\n"; }
                    if (!($gene_alignments[0] eq 'all' or $gene_alignments[0] =~ /^max_block_[0-9]+/)) {
                        if (!$run_module{"gene_parser"}) {
                            my $dbh = PifCosm_support_subs::connect_to_database ($database);
                            my $sth = $dbh->prepare("SELECT DISTINCT(gene) FROM alignment_groups") or die;
                            $sth->execute();
                            my @genes;
                            while (my $temp = $sth->fetchrow_array()) { push(@genes,$temp); }
                            for (my $k=0; $k < scalar @gene_alignments; ++$k) {
                                my $flag='n';
                                for (my $l=0; $l < scalar @genes; ++$l) {
                                    if ($gene_alignments[$k] eq $genes[$l]) { $flag='y'; last; }
                                }
                                if ($flag eq 'n') { die "The gene $gene_alignments[$k] is not present in $database and can not be used, remove or give other gene to use (-u/--use_genes).\n"; }
                            } 
                            $dbh->disconnect();
                        }
                    }
                }
            }
        }
        if (!$run_module{"gb_parser"} and !(-e $database)) {
            die "$database does not exist, provide other database or run module gb_parser.\n";
        }
    }
    else { die "No modules specified (--modules or -o and comma separated string with module names). Nothing to do so quitting.\n"; }
    # test of variables
    if (%backup) {
        print "Will backup after these modules:\n";
        foreach (keys %backup) {
            print "$_ ";
        }
        print "\n";
    }
    if (length($path) > 0) { print "Will use the path $path for program depandencies.\n"; }
    else { print "Will use the shell's path for dependencies.\n"; }
}
my $clustered = 'n'; # flag if some sort of clustering has been done
{   # Test if any table has been clustered
    my $dbh = PifCosm_support_subs::connect_to_database($database);
    my @tables = PifCosm_support_subs::get_gene_tables ($dbh);
    foreach (@tables) {
        my $sth = $dbh->prepare("SELECT COUNT(accno) FROM $_ WHERE cluster!='empty'");
        $sth->execute();
        my $number = $sth->fetchrow_array();
        if ($number and $number > 0) { $clustered = 'y'; }
        $sth->finish();
    }
}

if ($print_batch_file eq 'y') {
    my $filename = "$database\_batchfile.txt";
    my $i = 0;
    while (-e $filename ) { $filename = "$database\_batchfile_" . (++$i) . ".txt"; }
    &print_batch_file_sub( $filename );
}

### Run modules ###
{
for (my $i=0; $i < scalar @modules; ++$i) {
    if ($modules[$i] eq "gb_parser") { 
        &print_module_title ("Starting parsing GenBank sequences from $gb_data_file to $database.");
        Gb_parser::gb_parser($database,$gb_data_file,$max_seq_length);
        print "Parsed GenBank sequences.\n";
        if ($backup{"gb_parser"} or $backup{"all"}) {
            my $backup_name = $database . ".pars";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "check_taxonomy") {
        &print_module_title ("Checking if taxonomy in field taxon_string in gb_data is consistent.");
	PifCosm_support_subs::check_taxon_string_consistency($database, 'environmental samples;mycorrhizal samples', 'r');
        if ($backup{"check_taxonomy"} or $backup{"all"}) {
            my $backup_name = $database . ".check";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "change_entries") {
        &print_module_title ("Changing entries in table gb_data.");
        PifCosm_support_subs::change_entries ($database, $change_option, @changes_files );
        if ($backup{"change_entries"} or $backup{"all"}) {
            my $backup_name = $database . ".change";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "define_individuals") {
        &print_module_title ("Starting to identify individuals.");
        Define_individuals::define_individuals( $database, scalar @individual_criteria, @individual_criteria, $individual_condition );
        if ($backup{"define_individuals"} or $backup{"all"}) {
            my $backup_name = $database . ".ind";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "gene_parser") {
        if (!(-e $database)) { die "Could not find database file: $database. Will not be able to parse genes. Quitting.\n" }
        &print_module_title ("Starting to sort sequences in gb_data table of $database into separate tables for each gene.");
        Hammer_sort::hammer_sort ($database,$path,$hmmdatabase,$min_sequence_length_out,$e_value_cut_off,$n_cut_off,$print_non_match,%intergenes);
        print "Finished sorting sequences. " . PifCosm_support_subs::date() . ".\n";
        my $dbh = PifCosm_support_subs::connect_to_database ($database);
        my @tables = PifCosm_support_subs::get_gene_tables($dbh);
        foreach (@tables) {
            my $sth = $dbh->prepare("SELECT COUNT(accno) FROM $_" ) or die "Couldn't prepare statement: " . $dbh->errstr;
            $sth->execute();
            print $sth->fetchrow_array() . " sequences belong to $_.\n";
            $sth->finish();
        }
        if ($backup{"gene_parser"} or $backup{"all"}) {
            my $backup_name = $database . ".gene";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "cluster_individuals") {
        &print_module_title ("Clustering sequences from the same individual.");
        Cluster_indv_sequences::cluster_indv_sequences( $database, $min_sequence_length_out );
        $clustered = 'y';
        if ($backup{"cluster_individuals"} or $backup{"all"}) {
            my $backup_name = $database . ".ind_clust";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "tree_clustering") {
        &print_module_title ("Clustering sequences for monophyletic units.");
        my $dbh = PifCosm_support_subs::connect_to_database($database);
        my $cut_off = PifCosm_support_subs::get_cut_off_string( $sim_cut_off, PifCosm_support_subs::get_gene_tables( $dbh ) );
        $dbh->disconnect();
        #print $cut_off, "\n";
        Tree_cluster::tree_cluster ( $database,$path,$species_tree_clust,$sim_tree_clust,$cut_off,$min_sequence_length_out, $n_threads );
        $clustered = 'y';
        if ($backup{"tree_clustering"} or $backup{"all"}) {
            my $backup_name = $database . ".tree_clust";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "CD_hit_cluster") {
        &print_module_title ("Clustering sequences based on similarity using CD-hit.");
        my $dbh = PifCosm_support_subs::connect_to_database($database);
        my $cut_off = PifCosm_support_subs::get_cut_off_string( $sim_cut_off, PifCosm_support_subs::get_gene_tables( $dbh ) );
        $dbh->disconnect();
        CD_hit_clust::CD_hit_clust ( $database,$path,$cut_off,$min_sequence_length_out, $n_threads );
        $clustered = 'y';
        if ($backup{"CD_hit_cluster"} or $backup{"all"}) {
            my $backup_name = $database . ".CD_clust";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "similarity_cluster") {
        &print_module_title ("Clustering sequences based on similarity and define groups to align based on median absolute deviation.");
        my $dbh = PifCosm_support_subs::connect_to_database($database);
        my $cut_off = PifCosm_support_subs::get_cut_off_string( $sim_cut_off, PifCosm_support_subs::get_gene_tables( $dbh ) );
        $dbh->disconnect();
        my @cut_off_array = split /,/, $cut_off;
        my $switches = "--format sqlite --group both:cut-off=$cut_off_array[1]:min_length=$min_sequence_length_out";
        if ($clustered eq 'y') { $switches .= ":only_lead"; }
        if ($n_threads > 1) { $switches .= " -T $n_threads"; }
        system "${path}$External_program::pairalign -v $switches $database";
        $clustered = 'y';
        if ($backup{"similarity_cluster"} or $backup{"all"}) {
            my $backup_name = $database . ".sim_clust";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "alignment_groups") {
        &print_module_title ("Creating or remaking alignment groups table.");
        Create_alignment_groups::create_alignment_groups($database, $min_group);
        if ($backup{"alignment_groups"} or $backup{"all"}) {
            my $backup_name = $database . ".align_gr";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "link_genes") {
        &print_module_title ("Linking genes and creating table for alignments.");
        Gene_linker::gene_linker ( $database, $min_sequence_length_out, $maximum_linking );
        if ($backup{"link_genes"} or $backup{"all"}) {
            my $backup_name = $database . ".link";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "multimotu") {
        &print_module_title ("Linking genes using graph theory and creating table for alignments.");
        CV_gene_linker::CV_gene_linker ( $database, $comp_columns, $comp_type, $comp_score );
        if ($backup{"multimotu"} or $backup{"all"}) {
            my $backup_name = $database . ".mmotu";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "latin_names") {
        &print_module_title ("Getting 'latin binominals' for the taxa in the alignment table.");
        Latin_binomials::create_species_names_column ($database);
        if ($backup{"latin_names"} or $backup{"all"}) {
            my $backup_name = $database . ".l_names";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "align_sequences") {
        &print_module_title ("Starting aligning and creating trees.");
        Align_sequences::align_sequences($database, $path, $n_threads, $create_gene_tree, $tree_method, $store_boot_trees, $rm_outliers, @linked_genes);
        if ($backup{"align_sequences"} or $backup{"all"}) {
            my $backup_name = $database . ".align";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "refine_alignments") {
        &print_module_title ("Starting to refine alignments.");
        my $dbh = PifCosm_support_subs::connect_to_database ($database);
        my $flag = PifCosm_support_subs::table_present ( $dbh, 'alignment_groups');
        if ($flag eq 'n') {
            print "Did not find table alignment_groups, creating it.\n";
            $dbh->do("CREATE TABLE alignment_groups (gene TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', tree TEXT DEFAULT 'empty', tree_method TEXT DEFAULT 'empty',alignable INTEGER, PRIMARY KEY (gene, taxon))") or die "Could not create table alignment_groups: " . $dbh->errstr;
        }
        $dbh->disconnect();
        my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
        Refine_alignment::refine_alignment ( $database, $path, $tree_method, $store_boot_trees, $stop_criterion, $max_size, $use_guide_tree, $n_threads, join (',',@genes), join (',',@linked_genes));
        if ($backup{"refine_alignments"} or $backup{"all"}) {
            my $backup_name = $database . ".refine";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "gblocks") {
        &print_module_title ("Starting to remove ambious alignment blocks.");
        my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
        Gblocks::Gblocks ($database,$path,$b0,$b1,$b2,$b3,$b4,$b5,$min_sequence_length_out,@genes);
        if ($backup{"gblocks"} or $backup{"all"}) {
            my $backup_name = $database . ".gblocks";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "concatenated_tree") {
        &print_module_title ("Starting phylogenetic analysis on concatenated gene alignments.");
        my $dbh = PifCosm_support_subs::connect_to_database ($database);
        my $flag = PifCosm_support_subs::table_present ( $dbh, 'alignment_groups');
        if ($flag eq 'n') {
            print "Did not find table alignment_groups, creating it.\n";
            $dbh->do("CREATE TABLE alignment_groups (gene TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', tree TEXT DEFAULT 'empty', tree_method TEXT DEFAULT 'empty',alignable INTEGER, PRIMARY KEY (gene, taxon))")
                or die "Could not create table alignment_groups: " . $dbh->errstr;
        }
        $dbh->disconnect();
        my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
        Concatenated_tree::concatenated_tree ($database, $path, $final_tree_method, $n_threads, $store_boot_trees, $anchor_gene, @genes); # make a tree from concatenated alignment
        if ($backup{"concatenated_tree"} or $backup{"all"}) {
            my $backup_name = $database . ".concat_tree";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq "exclude_rogues") {
        &print_module_title ("Starting search for rogue taxa.");
        Exclude_rogues::exclude_rogues( $database, $path, $print_taxa, $output_tree, $remove_from_trees, $remove_from_alignment, $out_file_stem);
        if ($backup{"exclude_rogues"} or $backup{"all"}) {
            my $backup_name = $database . ".rogues";
            my $counter = 0;
            while (-e $backup_name) {
                $backup_name =~ s/_[0-9]*$//;
                $backup_name .= '_' . (++$counter);
            }
            &backup($database,$backup_name);
            print "Database backed up in file: $backup_name\n"
        }
    }
    elsif ($modules[$i] eq 'get_distinct_entries') {
        &print_module_title ("Outputting distinct entries from table gb_data.");
        PifCosm_support_subs::get_distinct_entries ($database, $out_file_stem, @column_name );
    }
    elsif ($modules[$i] eq "get_shared_values") {
	&print_module_title ("Comparing columns for shared values.");
	print "Checking @column_name.\n";
	PifCosm_support_subs::get_shared_values ($database,$out_file_stem,@column_name);
    }
    elsif ($modules[$i] eq 'get_alignments') {
        &print_module_title ("Outputting alignment.");
        my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
        PifCosm_support_subs::print_alignment($database,$alignment_format,$get_partitions,$out_file_stem,$align_seq_name_col,$interleaved,$bases_per_row,$anchor_gene,@genes);
    }
    elsif ($modules[$i] eq 'get_otu_taxonomy') {
	&print_module_title ("Getting OTU taxonomy");
	my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
	PifCosm_support_subs::getOTUtaxonomy(*STDOUT,$database,$align_seq_name_col,$anchor_gene,@genes);
    }
    elsif ($modules[$i] eq 'get_trees') {
        &print_module_title ("Outputting trees.");
        my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
        if ($only_show eq 'y') {
            my $dbh = PifCosm_support_subs::connect_to_database ($database);
            my $sth = $dbh->prepare("SELECT gene,taxon FROM alignment_groups WHERE tree!='empty'") or die;
            $sth->execute();
            print "Trees are available for the following groups:\n";
            while (my @row = $sth->fetchrow_array()) {
                print "Gene $row[0] and taxon $row[1].\n";
            }
            $sth->finish();
            $dbh->disconnect();
        }
        else {
            PifCosm_support_subs::print_trees($database,$out_file_stem,$taxon,@genes);
        }
    }
    elsif ($modules[$i] eq 'switch_tip_names') {
        &print_module_title ("Outputting tree with new tip names.");
        my $file_name = PifCosm_support_subs::change_names ( $database, $path, $tree_file, $out_file_stem, $input_tree_file );
        print "Tree with new names printed to $file_name.\n";
    }
    elsif ($modules[$i] eq 'stats') {
        &print_module_title ("Outputting alignment statistics.");
        my @genes;
        if ( $gene_alignments[0] =~ /max_block_([0-9]+)/ ) {
            @genes = PifCosm_support_subs::find_lardgest_block($database, $1);
        }
        else { @genes = @gene_alignments; }
        PifCosm_support_subs::stats($database,$stats_type,'y',$path,$anchor_gene,@genes);
    }
    else { print STDERR "WARNING!!! $modules[$i] is not recognized. Proceeding reluctantly.\n"; }
}
print "Finished with all. " . PifCosm_support_subs::date() . ".\n"
} exit; # end of scope for main

sub print_module_title {
    print "##########################\n##########################\n";
    print $_[0];
    print " " . PifCosm_support_subs::date() . "\n";
    print "##########################\n";
}

#*******************#
#*#################*#
#*## SUBROUTINES ##*#
#*#################*#
#*******************#

sub pars_batch_file {
    my $file = shift @_;
    open BATCHFILE, "<$file" or die "Could not open batchfile ($file): $!.\n";
    while (my $infile = <BATCHFILE> ) {
        chomp ($infile);
        $infile =~ s/([^\\]|^)#.+$/$1/;
        $infile =~ s/\\#/#/g;
        if (!($infile =~ /[a-zA-Z]/)) { next; }
        if ($infile =~ /^\s*database\s+([\w\/\-\.\\,]+)/i) {
            $database = $1;
        }
        elsif ($infile =~ /^\s*gb_file\s+([\w\/\-\.\\,]+)/i) {
            $gb_data_file = $1;
	    $gb_data_file =~ s/^\s+|\s+$//g;
        }
        elsif ($infile =~ /^\s*modules\s+([\w,]+)/i) {
            @modules = split /,/, $1;
        }
        elsif ($infile =~ /^\s*hmm_database\s+([\w\.\/\-\.\\,]+)/i) {
            $hmmdatabase = $1;
        }
        elsif ($infile =~ /^\s*print_non_match\s+(\w+)/i) {
            if ($1 =~ /^(y|Y)/) { $print_non_match = 'y'; }
            else { $print_non_match = 'n'; }
        }
        elsif ($infile =~ /^\s*individual_columns\s+([\w,]+)/i) {
            @individual_criteria = split /,/, $1;
        }
        elsif ($infile =~ /^\s*individual_condition\s+(\w+)/i) {
            $individual_condition = $1;
        }
        elsif ($infile =~ /^\s*inter_genes\s+([\w,;]+)/i) {
            %intergenes = split /,/, $1;
        }
        elsif ($infile =~ /^\s*backup\s*([\w,]*)/i) {
            if ($1) {
                my @temp = split /,/,$1;
                foreach (@temp) { $backup{$_} = 'y'; }
            }
            else { $backup{"all"} = 'y'; }
        }
        elsif ($infile =~ /^\s*e_value_cut_off\s+([0-9\.]+)/i) {
            $e_value_cut_off = $1;
        }
        elsif ($infile =~ /^\s*n_cut_off\s+([0-9\.]+)/i) {
            $n_cut_off = $1;
        }
        elsif ($infile =~ /^\s*sim_cut_off\s+([\w,\.]+)/i) {
            $sim_cut_off = $1;
        }
        elsif ($infile =~ /^\s*max_seq_length\s+([0-9]+)/i) {
            $max_seq_length = $1;
        }
        elsif ($infile =~ /^\s*min_seq_length_gene\s+([0-9]+)/i) {
            $min_sequence_length_out = $1;
        }
        elsif ($infile =~ /^\s*tree_clust_type\s+(\w+)/i) {
            if ($1 =~ /^b/i) {
                $species_tree_clust = 'y';
                $sim_tree_clust = 'y';
            }
            elsif ($1 =~ /^sp/i) {
                $species_tree_clust = 'y';
                $sim_tree_clust = 'n';
            }
            elsif ($1 =~ /^si/i) {
                $species_tree_clust = 'n';
                $sim_tree_clust = 'y';
            }
        }
        elsif ($infile =~ /^\s*group_on_genera\s+(\w+)/i) {
            if ($1 =~ /^y/i) { $min_group = 'y'; }
            else { $min_group = 'n'; }
        }
        elsif ($infile =~ /^\s*comp_columns\s+([0-9A-Za-z_\.,]+)/i) {
            $comp_columns = $1;
        }
        elsif ($infile =~ /^\s*comp_kind\s+([A-Za-z,]+)/i) {
            $comp_type = $1;
        }
        elsif ($infile =~ /^\s*comp_score\s+([0-9\.,]+)/i) {
            $comp_score = $1;
        }
        elsif ($infile =~ /^\s*phylo_method\s+([\w,-]+)/i) {
            $tree_method = $1;
        }
        elsif ($infile =~ /^\s*concat_phylo_method\s+([\w,-]+)/i) {
            $final_tree_method = $1;
        }
        elsif ($infile =~ /^\s*max_alignment_group_size\s+(\w+)/i) {
            $max_size = $1;
            if ($max_size =~ /[^0-9]/) { die "max_alignment_group_size must followed by a integer value.\n"; }
        }
        elsif ($infile =~ /^\s*use_guide_tree\s+(\w+)/i) {
            if ($1 =~ /^n/i) { $use_guide_tree = 'n'; }
            elsif ($1 =~ /^y/i) { $use_guide_tree = 'y'; }
            else { die "use_guide_tree must be followed by y/n (yes or no).\n"; }
        }
        elsif ($infile =~ /^\s*stop_criterion\s+([\w,]+)/i) {
            $stop_criterion = $1;
        }
        elsif ($infile =~ /^\s*alignment_format\s+([\w,]+)/i) {
            $alignment_format = $1;
        }
        elsif ($infile =~ /^\s*interleaved\s+([\w,]*)/i) {
            my $argument = $1;
            if ($argument and $argument =~ /\w+/) {
                if ( $argument =~ /^(y|Y)/ ) { $interleaved = 'y'; }
                elsif ( $argument =~ /^(n|N)/ ) { $interleaved = 'n'; }
                else { die "Do not recognize argument given after interleaved in batch file, it should be yes or no.\n"; }
            }
            else { $interleaved = 'y'; }
        }
        elsif ($infile =~ /^\s*bases_per_row\s+([\w,]+)/i) {
            $bases_per_row = $1;
            if ($bases_per_row =~ /^gene/i) {
                $bases_per_row = 'genes';
            }
            elsif ($bases_per_row =~ /[^0-9]/) { die "bases_per_row may only be followed by the key word genes or a integer number in batch file. Try again.\n"; }
        }
        elsif ($infile =~ /^\s*partition_file\s+([\w,]+)/i) {
            if ($1 =~ /^y/i) { $get_partitions = 'y'; }
            else { $get_partitions = 'n'; }
        }
        elsif ($infile =~ /^\s*use_genes\s+([\w,]+)/i) {
            @gene_alignments = split /,/,$1;
        }
        elsif ($infile =~ /^\s*anchor_gene\s+([\w, ]+)/i) {
            $anchor_gene = lc($1);
            $anchor_gene =~ s/ //g;
        }
        elsif ($infile =~ /^\s*maximize_linking\s+([\w]+)/i) {
            if ( $1 =~ /^y/i ) { $maximum_linking = 'y'; }
            elsif ( $1 =~ /^n/i ) { $maximum_linking = 'n'; }
            else { die "Do not recognize argument given after maximize_linking, it should be yes or no.\n"; }
        }
        elsif ($infile =~ /^\s*b0\s+([0-9]+)/i) {
            $b0 = $1;
        }
        elsif ($infile =~ /^\s*b1\s+([0-9]+)/i) {
            $b1 = $1;
        }
        elsif ($infile =~ /^\s*b2\s+([0-9]+)/i) {
            $b2 = $1;
        }
        elsif ($infile =~ /^\s*b3\s+([0-9]+)/i) {
            $b3 = $1;
        }
        elsif ($infile =~ /^\s*b4\s+([0-9]+)/i) {
            $b4 = $1;
        }
        elsif ($infile =~ /^\s*b5\s+([nha])/i) {
            $b5 = lc($1);
        }
        elsif ($infile =~ /^\s*get_taxon\s+([\w,]+)/i) {
            $taxon = $1;
        }
        elsif ($infile =~ /^\s*sequence_name_column\s+([\w,]+)/i) {
            $align_seq_name_col = $1;
        }
        elsif ($infile =~ /^\s*out_file_name\s+([\w,\/\-\.\\]+)/i) {
            $out_file_stem = $1;
        }
        elsif ($infile =~ /^\s*changes_files\s+([\w,\.\/\-\.\\]+)/i) {
            @changes_files = split /,/,$1;
        }
        elsif ($infile =~ /^\s*tree_file\s+([\w,\/\-\.\\]+)/i or $infile =~ /^\s*gene_tree\s+([\w,\/\-\.\\]+)/i) {
            $tree_file = $1;
        }
        elsif ($infile =~ /^\s*tree_infile/i) {
            $input_tree_file = 'y';
        }
        elsif ($infile =~ /^\s*linked_genes\s+([\w,;]+)/i) {
            @linked_genes = split /,/,$1;
        }
        elsif ($infile =~ /^\s*store_boot_trees/i) {
            $store_boot_trees='y';
        }
        elsif ($infile =~ /^\s*path\s+([\w\/~]+)/i) {
            $path = $1;
        }
        elsif ($infile =~ /^\s*column_name\s+([\w,]+)/i) {
            @column_name = split /,/,$1;
        }
        elsif ($infile =~ /^\s*only_show\s+(y|n)/i) { $only_show = lc($1); }
        elsif ($infile =~ /^\s*only_show/i) { $only_show = 'y'; }
        elsif ($infile =~ /^\s*threads\s+([0-9]+)/i) {
            $n_threads = $1;
        }
        elsif ($infile =~ /^\s*create_gene_trees\s+(\w+)/i) {
            my $argument = $1;
            if ($argument =~ /^y/i) { $create_gene_tree = 'y'; }
            elsif ($argument =~ /^n/i) { $create_gene_tree = 'n'; }
        }
        elsif ($infile =~ /^\s*create_gene_trees/i) { $create_gene_tree = 'y'; }
        elsif ($infile =~ /^\s*remove_outliers\s+(\w+)/i) {
            my $argument = $1;
            if ($argument =~ /^y/i) { $rm_outliers = 'y'; }
            elsif ($argument =~ /^n/i) { $rm_outliers = 'n'; }
        }
        elsif ($infile =~ /^\s*remove_outliers/i) { $rm_outliers = 'y'; }
        elsif ($infile =~ /^\s*print_rogues\s+(\S+)/i) {
            my $argument = $1;
            if ($argument =~ /^y/i) { $print_taxa = 'y'; }
            elsif ($argument =~ /^n/i) { $print_taxa = 'n'; }
            else { die "Invalid argument given for print_rogues, the alternatives are y/yes/n/no.\n"; }
        }
        elsif ($infile =~ /^\s*change_option\s+(\S+)/i) {
	    if ($1 =~ /^remove$/i || $1 =~ /^rm$/i) { $change_option = 'y'; }
	    elsif ($1 =~ /^change$/i) { $change_option = 'n'; }
    	    elsif ($1 =~ /^update$/i) { $change_option = 'u'; }
            else { die "Invalid argument given for change_option, the alternatives are y/yes/n/no.\n"; }
        }

        elsif ($infile =~ /^\s*print_rogues/i) { $print_taxa = 'y'; }
        elsif ($infile =~ /^\s*output_no_rogues_tree\s+(\S+)/i) {
            my $argument = $1;
            if ($argument =~ /^y/i) { $output_tree = 'y'; }
            elsif ($argument =~ /^n/i) { $output_tree = 'n'; }
            else { die "Invalid argument given for output_no_rogues_tree, the alternatives are y/yes/n/no.\n"; }
        }
        elsif ($infile =~ /^\s*output_no_rogues_tree/i) { $output_tree = 'y'; }
        elsif ($infile =~ /^\s*rm_rogues_from_ml_tree\s+(\S+)/i) {
            my $argument = $1;
            if ($argument =~ /^y/i) { $remove_from_trees = 'y'; }
            elsif ($argument =~ /^n/i) { $remove_from_trees = 'n'; }
            else { die "Invalid argument given for rm_rogues_from_ml_tree, the alternatives are y/yes/n/no.\n"; }
        }
        elsif ($infile =~ /^\s*rm_rogues_from_ml_tree/i) { $remove_from_trees = 'y'; }
        elsif ($infile =~ /^\s*rm_rogues_from_alignment\s+(\S+)/i) {
            my $argument = $1;
            if ($argument =~ /^y/i) { $remove_from_alignment = 'y'; }
            elsif ($argument =~ /^n/i) { $remove_from_alignment = 'n'; }
            else { die "Invalid argument given for rm_rogues_from_alignment, the alternatives are y/yes/n/no.\n"; }
        }
        elsif ($infile =~ /^\s*rm_rogues_from_alignment/i) { $remove_from_alignment = 'y'; }
        elsif ($infile =~ /^\s*print_batch_file\s+([a-zA-Z]+)/i) {
            if ($1 =~ /^y/i) { $print_batch_file = 'y'; }
            elsif ($1 =~ /^n/i) { $print_batch_file = 'n'; }
        }
        elsif ($infile =~ /^\s*stats_type\s+(\w+)/i) {
            $stats_type=$1;
            if ($stats_type ne 'alignment' and $stats_type ne 'cluster') { die "Illeagal value for stats_type.\n"; }
        }
        else { die "Could not understand command $infile.\n" }
    }
    close BATCHFILE or die;
}

sub print_batch_file_sub {
    my $file = shift @_;
    open BATCHFILE, ">$file" or die "Could not open batchfile ($file) for writing: $!.\n";
    if ($database) { print BATCHFILE "database $database\n"; }
    if ($gb_data_file) { print BATCHFILE "gb_file $gb_data_file\n" }
    if ( @modules ) { my $temp = join ",", @modules; print BATCHFILE "modules $temp\n"; }
    if ($hmmdatabase) { print BATCHFILE "hmm_database $hmmdatabase\n"; }
    if ($print_non_match) { print BATCHFILE "print_non_match $print_non_match\n"; }
    if (@individual_criteria) {
            my $temp = join ",", @individual_criteria;
            print "individual_columns $temp\n";
    }
    if ($individual_condition) { print BATCHFILE "individual_condition $individual_condition\n"; }
    if (%intergenes) {
        my $temp;
        foreach (keys %intergenes) { $temp .= "$_,$intergenes{$_},"; }
        $temp =~ s/,$//;
         print BATCHFILE "inter_genes $temp\n";
    }
    if (%backup) {
        my $temp = join ",", (keys %backup);
        print "backup $temp\n";
    }
    if ($e_value_cut_off) { print BATCHFILE "e_value_cut_off $e_value_cut_off\n"; }
    if ($n_cut_off) { print BATCHFILE "n_cut_off $n_cut_off\n"; }
    if ($sim_cut_off) { print BATCHFILE "sim_cut_off $sim_cut_off\n"; }
    if ($max_seq_length) { print BATCHFILE "max_seq_length $max_seq_length\n"; }
    if ($min_sequence_length_out) { print BATCHFILE "min_seq_length_gene $min_sequence_length_out\n"; }
    if ($species_tree_clust eq 'y' and $sim_tree_clust eq 'y') { print BATCHFILE "tree_clust_type both\n"; }
    elsif ($species_tree_clust eq 'y') { print BATCHFILE "tree_clust_type species\n"; }
    elsif ($sim_tree_clust eq 'y') { print BATCHFILE "tree_clust_type similarity\n"; }
    if ($min_group) { print BATCHFILE "group_on_genera $min_group\n"; }
    if ($create_gene_tree) { print BATCHFILE "create_gene_trees $create_gene_tree\n"; }
    if ($tree_method) { print BATCHFILE "phylo_method $tree_method\n"; }
    if ($rm_outliers) { print BATCHFILE "remove_outliers $rm_outliers\n"; }
    if ($print_taxa) { print BATCHFILE "print_rogues $print_taxa\n"; }
    if ($output_tree) { print BATCHFILE "output_no_rogues_tree $output_tree\n"; }
    if ($remove_from_trees) { print BATCHFILE "rm_rogues_from_ml_tree $remove_from_trees\n"; }
    if ($remove_from_alignment) { print BATCHFILE "rm_rogues_from_alignment $remove_from_alignment\n"; }
    if ($change_option) {
	print BATCHFILE "change_option ";
	if ($change_option eq 'y') { print BATCHFILE 'remove'; }
	elsif ($change_option eq 'n') { print BATCHFILE 'change'; }
	elsif ($change_option eq 'u') { print BATCHFILE 'update'; }
    }
    if ($final_tree_method) { print BATCHFILE "concat_phylo_method $final_tree_method\n"; }
    if ($max_size) { print BATCHFILE "max_alignment_group_size $max_size\n"; }
    if ($use_guide_tree) { print BATCHFILE "use_guide_tree $use_guide_tree\n"; }
    if ($stop_criterion) { print BATCHFILE "stop_criterion $stop_criterion\n"; }
    if ($alignment_format) { print BATCHFILE "alignment_format $alignment_format\n"; }
    if ($interleaved) { print BATCHFILE "interleaved $interleaved\n"; }
    if ($bases_per_row) { print BATCHFILE "bases_per_row $bases_per_row\n" }
    if ($get_partitions) { print BATCHFILE "partition_file $get_partitions\n"; }
    if (@gene_alignments) {
        my $temp = join ",", @gene_alignments;
        print BATCHFILE "use_genes $temp\n";
    }
    if ($maximum_linking) { print BATCHFILE "maximize_linking $maximum_linking\n"; }
    if ($comp_columns) { print BATCHFILE "comp_columns $comp_columns\n"; }
    if ($comp_type) { print BATCHFILE "comp_kind $comp_type\n"; }
    if ($comp_score) { print BATCHFILE "comp_score $comp_score\n"; }
    if ($b0 > 2) { print BATCHFILE "b0 $b0\n"; }
    if ($b1 > 0) { print BATCHFILE "b1 $b1\n"; }
    if ($b2 > 0) { print BATCHFILE "b2 $b2\n"; }
    if ($b3 > 0) { print BATCHFILE "b3 $b3\n"; }
    if ($b4 > 1) { print BATCHFILE "b4 $b4\n"; }
    if ($b5 eq 'n' or $b5 eq 'h' or $b5 eq 'a') { print BATCHFILE "b5 $b5\n"; }
    if ($taxon) { print BATCHFILE "get_taxon $taxon\n"; }
    if ($anchor_gene) { print BATCHFILE "anchor_gene $anchor_gene\n"; }
    if ($align_seq_name_col) { print BATCHFILE "sequence_name_column $align_seq_name_col\n"; }
    if ($out_file_stem) { print BATCHFILE "out_file_name $out_file_stem\n"; }
    if (@changes_files) {
        my $temp = join ",", @changes_files;
        print BATCHFILE "changes_files $temp\n";
    }
    if ($tree_file) { print BATCHFILE "tree_file $tree_file\n"; }
    if ($input_tree_file eq 'y') { print BATCHFILE "tree_infile\n"; }
    if (@linked_genes) {
        my $temp = join ",", @linked_genes;
        print BATCHFILE "linked_genes $temp\n";
    }
    if ($store_boot_trees eq 'y') { print BATCHFILE "store_boot_trees\n"; }
    if ($path) { print BATCHFILE "path $path\n"; }
    if (@column_name) {
        my $temp = join ",", @column_name;
        print BATCHFILE "column_name $temp\n";
    }
    if ($only_show) { print BATCHFILE "only_show $only_show\n"; }
    if ($stats_type) { print BATCHFILE "stats_type $stats_type\n"; }
    if ($n_threads) { print BATCHFILE "threads $n_threads\n"; }
    if ($print_batch_file) { print BATCHFILE "print_batch_file $print_batch_file\n"; }
    close BATCHFILE or die;
}

sub backup {
    my $file = shift @_;
    my $backup_file = shift @_;
    open FILE, "<$file" or die "Could not open $file.\n";
    open BACKUP,">$backup_file" or die "Could not open $backup_file.\n";
    while (my $infile= <FILE>) {
        print BACKUP $infile;
    }
    close FILE or die;
    close BACKUP or die;
}
