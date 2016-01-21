#!/gsc/bin/perl

#  CLASSIFIES SPLICE-SITE MUTATIONS BASED ON EXPRESSION LEVELS AND READ COUNTS
#
#  The statistical test is designed to recognize exon skipping and intron
#  retention
#
#  VERSIONING:
#
#  0.1: original version designed for exon-skipping and intron retention
#
#  NOTES
#
#  * FDR
#
#    We've used FDR many times and code is widely available, see for example
#
#    research/completed/2011_pathscan/tsp_analysis/process_pathscan_file.pl
#    research/in_progress/blood_clonal_expansion/statistics_blood_normal_test/complete_statistical_pipeline_matched.pl
#
#####################NEWEST
#
# Date: Thu, 02 Apr 2015 07:22:07 -0500
# From: Song Cao <scao@genome.wustl.edu>
# To: Mike Wendl <mwendl@genome.wustl.edu>,
#        Reyka Jayasinghe <rjayasin@genome.wustl.edu>,
#        Li Ding <lding@genome.wustl.edu>
#Subject: rpkm for new ucec maf
#User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.9; rv:31.0) Gecko/20100101 Thunderbird/31.5.0
#
#Hi all,
#
#The rpkm values for new ucec maf are available in the following directory:
#
#1) full length
#
#case
#/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_v2/new_maf_ucec/ucec_splice_cancergene_rpkm_case.tsv
#control
#/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_v2/new_maf_ucec/ucec_splice_cancergene_rpkm_control_stringent.tsv
#
#2) 50bp
#
#case
#/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_50_v2/new_maf_ucec/ucec_splice_cancergene_rpkm_case.tsv
#control
#/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_50_v2/new_maf_ucec/ucec_splice_cancergene_rpkm_control_stringent.tsv
#
#3) 2bp
#
#case
#/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_2_v2/new_maf_ucec/ucec_splice_cancergene_rpkm_case.tsv
#control
#/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_2_v2/new_maf_ucec/ucec_splice_cancergene_rpkm_control_stringent.tsv


############
#  SET UP  #
############

#__STANDARD PACKAGES
   use strict;
   use warnings;

#__CPAN PACKAGES
   use Statistics::Descriptive;
   use Statistics::RankCorrelation;
#  use Math::Combinatorics;
#  use List::Compare;
#  use Data::Section::Simple qw/get_data_section/;

#__SPECIAL (NON-STANDARD) TGI PACKAGES
#  use lib "../integrated_analysis/scripts";
 use lib "./";
 use SpliceInator ':all';
 use lib "/gscuser/mwendl/submissions/in_progress/splice_site_reyka";  
 use Statistics::CombinePvals;

###################
#  CONFIGURATION  #
###################

#__APPLY SPOT-CHECKING TO INPUT FILES: ONLY NEEDED ON 1ST RUN WITH NEW FILES
   my $spot_check = 0;

#__THRESHOLD P-VALUE CUTOFF FOR ALL GROUPS OF CLASSIFICATION #EXCEPT GENE EXPRESSION - HARD CODED INTO GENE EXPRESSION TEST
   my $threshold_classification = 0.05;

#__SPEARMAN RANK CORRELATION COEFFICIENT THRESHOLD FOR EXCLUSION BY 3' BIAS
#  my $spearman_threshold = 0.5;
   my $spearman_threshold = 0.3;

#__THRESHOLD BELOW WHICH AN INTRON OR EXON IS PRESUMED *NOT* TO BE EXPRESSED
#  my $expression_threshold_exon = 0;
   my $expression_threshold_exon = 3;

#  my $expression_threshold_intron = 0.25;
#  my $expression_threshold_intron = 3;
#  my $expression_threshold_intron = 0.2;
   my $expression_threshold_intron = 0.5;

#__INPUT: FULL EXON RPKM EXPRESSION OF CASES AND CONTROLS
   my $path = ".";
##################
##################
##################
   my $file_cases_full_expr = "/gscuser/scao/gc2523/dinglab/splice/pipeline/combinedfile/full/all_splice_cancergene_rpkm_case.tsv";
my $file_controls_full_expr = "/gscuser/scao/gc2523/dinglab/splice/pipeline/combinedfile/full/all_splice_cancergene_control_stringent.unique.tsv";
   my $file_cases_boundary_2bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/combinedfile/2bps/all_splice_cancergene_rpkm_case.tsv";
   my $file_controls_boundary_2bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/combinedfile/2bps/all_splice_cancergene_control_stringent.unique.tsv";
   my $file_cases_boundary_50bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/combinedfile/50bps/all_splice_cancergene_rpkm_case.tsv";
   my $file_controls_boundary_50bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/combinedfile/50bps/all_splice_cancergene_control_stringent.unique.tsv";
#  my $file_cases_full_expr = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_v2/stringent/all_splice_cancergene_rpkm_case.tsv";
#  my $file_controls_full_expr = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_v2/stringent/all_splice_cancergene_control_stringent.unique.tsv";
#  my $file_cases_boundary_2bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_2_v2/stringent/all_splice_cancergene_rpkm_case.tsv";
#  my $file_controls_boundary_2bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_2_v2/stringent/all_splice_cancergene_control_stringent.unique.tsv";
#  my $file_cases_boundary_50bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_50_v2/stringent/all_splice_cancergene_rpkm_case.tsv";
#  my $file_controls_boundary_50bp = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rpkm_exon_intron_50_v2/stringent/all_splice_cancergene_control_stringent.unique.tsv";

#__INPUT: WHOLE-GENE RSEM EXPRESSION OF CASES AND CONTROLS
   my $whole_gene_rsem_file = "/gscuser/scao/gc2523/dinglab/splice/pipeline/rsem/cancer_genes_splice_rsem.tsv";

#__SET MINIMUM NUMBER OF ELEMENTS TO ACTUALLY DO PERMUTATION TEST
#  my $min_elements = 8;
#  my $min_elements = 20;
#  my $min_elements_controls = 200;
#  my $min_elements_controls = 100;
   my $min_elements_controls = 50;
   my $min_elements_samples = 10;

####################
#  PRE PROCESSING  #
####################

# hnsc      LRP1B    TCGA-CN-5366-01A-01D-1434-08       NO

#__SONG CAO LIST OF EVENTS THAT ARE NON-EXPRESSED --- THESE SHOULD BE SKIPPED
   my $group_i_hits = {};
#  while (<DATA>) {
#     chomp;
#     my ($cancer, $mutation) = split /\s+/;
#     last if $cancer eq "######" and $mutation eq "END_SONG_DATA";
#     $group_i_hits->{$cancer}->{$mutation} = [];
######print "SONG DATA: $cancer  $mutation\n";
#  }

#__READ REYKAS LIST OF UNEXPRESSED CASES (SONG'S FALSE NEGATIVES)
   my $overlooked_instance_of_no_expression = {};
#  while (<DATA>) {
#     chomp;
#     my ($cancer, $gene, $sample) = split /\s+/;
#     next if $cancer eq "######" && $gene eq "OVERLOOKED_INSTANCES_OF_NO_EXPRESSION" && $sample eq "REYKA_10_FEB_2015";
#     $overlooked_instance_of_no_expression->{$cancer}->{$gene}->{$sample} = 1;
######print "REYKA DATA: $cancer  $gene   $sample\n";
#  }

#__SET ABSOLUTE MINIMUM FOR PERMUTATION TEST: GET DIVIDE-BY-ZERO IF VAL = 1
#
#  hnsc  TCGA-CN-5360-01A-01D-1434-08  FOXA1  ENST00000250448 - 14_38061917_C_G
#      38064105-38064238,38061916-38064104,38059188-38061915
#      7.66474445246594,1.32499044006461,9.23470958537392

   $min_elements_samples = 3 unless $min_elements_samples > 3;

#__READ FULL EXPRESSION FILES
   my $cases = read_full_expression (
      $file_cases_full_expr, $spot_check, "perform_indexing", $group_i_hits
   );
   my $controls = read_full_expression (
      $file_controls_full_expr, $spot_check
   );

#####################
#  MAIN PROCESSING  #
#####################

#__DATA STRUCTURE FOR THE CATALOG OF P-VALUE RESULTS
   my $catalog_pvals = {};

#######################################################
#  GROUP I: NO OR LOW-LEVEL EXPRESSION OF WHOLE GENE  #
#######################################################

#__DIAGNOSTIC VALUES FOR REYKA --- EXPRESSION LEVELS OF CASE FOR EACH CATEGORY
   my $diagnostics = {}; 

#
  print "# <GROUP I> NO OR LOW-LEVEL EXPRESSION OF WHOLE GENE\n";
my ($rsem_data);
($rsem_data) = parse_gene_rsem($whole_gene_rsem_file);
#my ($cancertype,$genes,$samples,$caseval,$pval);
#($cancertype,$genes,$samples,$caseval,$pval) = gene_expression($rsem_data);
my $group_i_count;
   ($group_i_hits, $group_i_count,$diagnostics) = gene_expression ($cases,$controls,$catalog_pvals,$threshold_classification,$diagnostics,$rsem_data);
###Bottom part initially in there?  
 $group_i_count = 0;
   foreach my $cancer (keys %{$group_i_hits}) {
      foreach my $gene (keys %{$group_i_hits->{$cancer}}) {
         foreach my $sample (keys %{$group_i_hits->{$cancer}->{$gene}}) {
	    #__RETRIEVE AND OUTPUT
            my ($pval, $test) = @{$group_i_hits->{$cancer}->{$gene}->{$sample}};
            print "$cancer   $gene   $sample   $pval ($test)\n";
           # my ($sample, $gene) = @{$pair};
            $group_i_count++;
           #print "$cancer   $gene   $sample\n";
	 }
      }
   }
  print "# </GROUP I> TOTAL = $group_i_count\n";

#######################################
#  GROUP II: SIMPLE INTRON RETENTION  #
#######################################

#__DIAGNOSTIC VALUES FOR REYKA --- EXPRESSION LEVELS OF CASE FOR EACH CATEGORY
#   my $diagnostics = {};

#__INTRONS: TEST EACH CASE (EACH SAMPLE HAVING A SPLICE SITE MUTATION)
   print "# <GROUP II> SIMPLE INTRON RETENTION (C=CONTROL ONLY PVAL, B=CONTROL+IN_SAMPLE)\n";
   print "# cancer   gene   sample_name   P-val (basis)\n";
   print "# basis:\n";
   print "#  C = P-val calculated only against control introns\n";
   print "#  B = calculated against both controls and in-sample introns\n";
   print "# parameters:\n";
   print "#  P-value threshold for classification: $threshold_classification\n";
   print "#  Spearman threshold for 3' bias: $spearman_threshold\n";
   print "#  minimum threshold of expression: $expression_threshold_intron\n";
   my ($group_ii_hits, $group_ii_count);
   ($group_ii_hits, $group_ii_count, $diagnostics) = simple_intron_retention (
      $cases, $controls, $catalog_pvals, $min_elements_controls,
      $expression_threshold_intron, $spearman_threshold,
      $threshold_classification, $diagnostics
   );

#__PROCESS THE TEST RESULTS
   foreach my $cancer (keys %{$group_ii_hits}) {
      foreach my $gene (keys %{$group_ii_hits->{$cancer}}) {
         foreach my $sample (keys %{$group_ii_hits->{$cancer}->{$gene}}) {

         #__RETRIEVE AND OUTPUT
            my ($pval, $test) = @{$group_ii_hits->{$cancer}->{$gene}->{$sample}};
            print "$cancer   $gene   $sample   $pval ($test)\n";

#        #__REMOVE GROUP II HITS FROM THE CONTROL SET
#           delete $cases->{$cancer}->{$gene}->{$sample};
         }
      }
   }
   print "# </GROUP II> TOTAL = $group_ii_count\n";

####################################
#  GROUP III: SIMPLE EXON SKIPPING #
####################################

#__EXONS: TEST EACH CASE (EACH SAMPLE HAVING A SPLICE SITE MUTATION)
   print "# <GROUP III> SIMPLE EXON SKIPPING (C=CONTROL ONLY PVAL, B=CONTROL+IN_SAMPLE)\n";
   print "# cancer   gene   sample_name   P-val (basis)\n";
   print "# basis:\n";
   print "#  C = P-val calculated only against control exons\n";
   print "#  B = calculated against both controls and in-sample exons\n";
   print "# parameters:\n";
   print "#  P-value threshold for classification: $threshold_classification\n";
   print "#  Spearman threshold for 3' bias: $spearman_threshold\n";
   print "#  minimum threshold of expression: $expression_threshold_exon\n";
   my ($group_iii_hits, $group_iii_count);
   ($group_iii_hits, $group_iii_count, $diagnostics) = simple_exon_skipping (
      $cases, $controls, $catalog_pvals, $min_elements_controls,
      $expression_threshold_exon, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__PROCESS THE TEST RESULTS
   my $fps = 0;
   foreach my $cancer (keys %{$group_iii_hits}) {
      foreach my $gene (keys %{$group_iii_hits->{$cancer}}) {
         foreach my $sample (keys %{$group_iii_hits->{$cancer}->{$gene}}) {

         #__RETRIEVE AND OUTPUT
            my ($pval, $test) = @{$group_iii_hits->{$cancer}->{$gene}->{$sample}};
            print "$cancer   $gene   $sample   $pval ($test)";
            if (defined $overlooked_instance_of_no_expression->{$cancer}->{$gene}->{$sample}) {
               print "    FALSE POS";
               $fps++;
               $group_iii_count--;
            }
            print "\n";

#        #__REMOVE GROUP III HITS FROM THE CONTROL SET
#           delete $cases->{$cancer}->{$gene}->{$sample};
         }
      }
   }
#  print "# </GROUP III> TOTAL = $group_iii_count\n";
   print "# </GROUP III> TOTAL REAL = $group_iii_count   FPS = $fps\n";

########################################
#  GROUP IV: POSTERIOR EXON SHRINKAGE  #
########################################

#__READ BOUNDARY EXPRESSION FILES FOR THE 2BP BOUNDARY DATA
   my $cases_2bp = read_boundary_expression (
      $file_cases_boundary_2bp, 1, $spot_check
   );
   my $controls_2bp = read_boundary_expression (
      $file_controls_boundary_2bp, 0, $spot_check
   );

##_SAVE SOME MEMORY BECAUSE WE ONLY NEED INDECES NOW
## %{$controls} = ();
## undef $cases->{$cancer}->{$gene}->{$sample}->{'values'};

#die "END OF TEST";

#__PRIOR EXON SHRINKAGE EVALUATION (2BP)
   print "prior exon shrinkage evaluation (2bp)\n";
   my ($group_iv_z_hits, $group_iv_z_count);
   ($group_iv_z_hits, $group_iv_z_count, $diagnostics) = prior_exon_shrinkage (
      $cases, $cases_2bp, $controls_2bp, "prior 3' exon shrinkage (2bp)",
      $catalog_pvals, $min_elements_controls,
      $expression_threshold_exon, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__POSTERIOR EXON SHRINKAGE EVALUATION (2BP)
   print "posterior exon shrinkage evaluation (2bp)\n";
   my ($group_iv_hits, $group_iv_count);
   ($group_iv_hits, $group_iv_count, $diagnostics) = posterior_exon_shrinkage (
#     $cases, $cases_2bp, $controls_2bp, "posterior exon shrinkage (2bp)",
      $cases, $cases_2bp, $controls_2bp, "posterior 5' exon shrinkage (2bp)",
      $catalog_pvals, $min_elements_controls,
      $expression_threshold_exon, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__POSTERIOR EXON EXTENSION EVALUATION (2BP)
   print "posterior exon extension evaluation (2bp)\n";
   my ($group_v_hits, $group_v_count);
   ($group_v_hits, $group_v_count, $diagnostics) = posterior_exon_extension (
#     $cases, $cases_2bp, $controls_2bp, "posterior exon extension (2bp)",
      $cases, $cases_2bp, $controls_2bp, "posterior 5' exon extension (2bp)",
      $catalog_pvals, $min_elements_controls,
#     $expression_threshold_exon, $spearman_threshold, $threshold_classification
      0, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__PRIOR EXON EXTENSION EVALUATION (2BP)
   print "prior exon extension evaluation (2bp)\n";
   my ($group_vi_hits, $group_vi_count);
   ($group_vi_hits, $group_vi_count, $diagnostics) = prior_exon_extension (
#     $cases, $cases_2bp, $controls_2bp, "prior exon extension (2bp)",
      $cases, $cases_2bp, $controls_2bp, "prior 3' exon extension (2bp)",
      $catalog_pvals, $min_elements_controls,
#     $expression_threshold_exon, $spearman_threshold, $threshold_classification
      0, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__DUMP THE 2BP BOUNDARY EXPRESSION DATA SINCE WE'RE DONE WITH IT
#  %{$cases} = ();
#  print "UNDEFF FULL EXPR\n";
# print "READ BOUNDARY EXPR\n";
    %{$cases_2bp} = ();
    %{$controls_2bp} = ();

#__READ BOUNDARY EXPRESSION FILES FOR THE 50BP BOUNDARY DATA
   my $cases_50bp = read_boundary_expression (
      $file_cases_boundary_50bp, 1, $spot_check
   );
   my $controls_50bp = read_boundary_expression (
      $file_controls_boundary_50bp, 0, $spot_check
   );

#__PRIOR EXON SHRINKAGE EVALUATION (50BP)
   print "prior exon shrinkage evaluation (50bp)\n";
   my ($group_iv_z2_hits, $group_iv_z2_count);
   ($group_iv_z2_hits, $group_iv_z2_count, $diagnostics) = prior_exon_shrinkage (
      $cases, $cases_50bp, $controls_50bp, "prior 3' exon shrinkage (50bp)",
      $catalog_pvals, $min_elements_controls,
      $expression_threshold_exon, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__POSTERIOR EXON SHRINKAGE EVALUATION (50BP)
   print "posterior exon shrinkage evaluation (50bp)\n";
   my ($group_iv5_hits, $group_iv5_count);
   ($group_iv5_hits, $group_iv5_count, $diagnostics) = posterior_exon_shrinkage (
#     $cases, $cases_50bp, $controls_50bp, "posterior exon shrinkage (50bp)",
      $cases, $cases_50bp, $controls_50bp, "posterior 5' exon shrinkage (50bp)",
      $catalog_pvals, $min_elements_controls,
      $expression_threshold_exon, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__POSTERIOR EXON EXTENSION EVALUATION (50BP)
   print "posterior exon extension evaluation (50bp)\n";
   my ($group_v5_hits, $group_v5_count);
   ($group_v5_hits, $group_v5_count, $diagnostics) = posterior_exon_extension (
#     $cases, $cases_50bp, $controls_50bp, "posterior exon extension (50bp)",
      $cases, $cases_50bp, $controls_50bp, "posterior 5' exon extension (50bp)",
      $catalog_pvals, $min_elements_controls,
#     $expression_threshold_exon, $spearman_threshold, $threshold_classification
      0, $spearman_threshold, $threshold_classification, $diagnostics
   );

#__PRIOR EXON EXTENSION EVALUATION (50BP)
   print "prior exon extension evaluation (50bp)\n";
   my ($group_vi5_hits, $group_vi5_count);
   ($group_vi5_hits, $group_vi5_count, $diagnostics) = prior_exon_extension (
#     $cases, $cases_50bp, $controls_50bp, "prior exon extension (50bp)",
      $cases, $cases_50bp, $controls_50bp, "prior 3' exon extension (50bp)",
      $catalog_pvals, $min_elements_controls,
#     $expression_threshold_exon, $spearman_threshold, $threshold_classification
      0, $spearman_threshold, $threshold_classification, $diagnostics
   );

#####################
#  POST PROCESSING  #
#####################

#__REPORT THE CATALOG
   print "# THE CATALOG OF P-VALUES\n";
   foreach my $cancer (keys %{$catalog_pvals}) {
#     print "CANCER: $cancer\n";
      foreach my $gene (keys %{$catalog_pvals->{$cancer}}) {
#        print "  GENE: $gene\n";
         foreach my $sample (keys %{$catalog_pvals->{$cancer}->{$gene}}) {
#           print "    SAMPLE: $sample";
            print "$cancer   $gene   $sample   $cases->{$cancer}->{$gene}->{$sample}->{'mut_position'}\n";
#           print "    $cases->{$cancer}->{$gene}->{$sample}->{'strand'}";
#           print "\n";
            foreach my $mode (sort keys %{$catalog_pvals->{$cancer}->{$gene}->{$sample}}) {
#              print "      $mode : $catalog_pvals->{$cancer}->{$gene}->{$sample}->{$mode}\n" unless $catalog_pvals->{$cancer}->{$gene}->{$sample}->{$mode} =~ /no\stest/;
               print "      $mode : $catalog_pvals->{$cancer}->{$gene}->{$sample}->{$mode} (expr: $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode})\n";
            }
         }
      }
   }

###################################
#  DONE WITH SIMPLE ALTERNATIVES  #
###################################
#  my $total_remaining = 0;
#  foreach my $cancer (keys %{$cases}) {
#     foreach my $gene (keys %{$cases->{$cancer}}) {
#        foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {
#           $total_remaining++;
#        }
#     }
#  }
#  print "TOTAL REMAINING AFTER THESE CLASSIFICATIONS: $total_remaining\n";

# die "GROUPS I AND II AND III: end of test";

################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################

#  EXPLANATION OF DATA
#
#  THESE ARE THE SITES WHICH SONG CAO AD-HOC FILTERED AS BEING FROM GENES
#  THAT ARE NOT EXPRESSED (26-JAN2015) --- THIS WILL BE REPLACED LATER
#  WITH A MORE SOPHISTICATED METHOD BASED ON FISHER LDA OF ABSOLUTE RSEM
#  EXPRESSION AND PERMUTATION P-VALUE

__DATA__
lusc	2_212566892_C_A
lusc	8_88365862_A_T
lusc	7_87178835_C_T
lusc	5_21842439_T_A
lusc	2_170127413_C_A
lusc	4_66280001_C_G
lusc	1_156845871_G_T
lusc	4_66197847_C_T
lusc	12_78534061_G_T
lusc	11_106681206_T_C
laml	X_133547852_G_T
laml	11_32413611_C_T
kirp	12_78591182_G_C
blca	6_160663435_C_G
blca	9_120470999_TTTATCCAGGTAATGAATC_-
blca	11_108196036_G_A
cesc	2_141625833_C_-
cesc	5_112136975_G_A
cesc	13_28631600_C_T
cesc	5_21975495_C_T
cesc	13_32937315_G_A
lgg	10_104353823_G_T
lgg	9_289509_G_A
lgg	22_22322990_C_A
hnsc	2_212488770_C_A
hnsc	12_40631906_G_T
hnsc	20_9525142_C_A
hnsc	9_98238315_C_A
hnsc	2_212286830_C_A
hnsc	2_169791937_T_C
hnsc	2_141457818_C_A
hnsc	2_141083449_T_A
hnsc	3_89156987_G_C
hnsc	9_37014993_C_T
hnsc	2_141459289_C_-
hnsc	10_101591889_G_A
hnsc	2_141201897_A_C
hnsc	17_59770874_C_G
brca	17_41223256_C_T
brca	13_32932067_G_A
brca	2_135966549_C_A
brca	3_185823093_C_T
brca	4_106190766_G_C
brca	2_212615430_C_A
brca	19_17945532_T_C
brca	3_185165736_G_A
brca	15_88671941_C_A
brca	11_108186737_G_A
brca	19_2211749_G_T
ucs	9_8521547_C_T
ucs	13_28602426_C_A
luad	2_141135856_C_G
luad	11_92590378_A_T
luad	2_141474386_C_A
luad	4_66509147_T_A
luad	9_8465465_C_A
luad	4_66286284_C_G
luad	3_89499326_G_T
luad	16_10032409_C_A
luad	10_104593791_A_G
luad	2_140992354_C_-
luad	13_28602427_T_A
luad	2_29917881_C_T
luad	3_121230733_C_G
luad	14_45633769_G_T
luad	5_44310633_C_T
luad	9_8526626_C_A
luad	9_8492981_T_C
luad	7_87160609_C_G
luad	2_142012211_C_T
luad	2_141625833_C_-
luad	13_28644627_C_A
luad	2_141812830_T_A
luad	5_19473826_C_G
luad	2_141200193_C_A
luad	18_42449194_G_T
luad	4_66213922_C_A
luad	8_93023228_A_G
luad	5_112111325_G_T
luad	1_156845459_G_T
luad	2_141707803_C_T
luad	2_141747099_A_T
luad	17_59820497_T_A
luad	4_71066320_T_C
luad	2_140995720_C_G
luad	15_88476242_C_A
lihc	6_117658467_ATGAAGTTTTAACATGGTAAAACTCATTGTATTTCCACTAGAAAAAGAAGTCTCGATTAATATTTTTGTTTCTCTAAGAAAATATTCTTAAAACACAAAAATGTATAGGCTTATAGCAAGTGGTTAAAATTTAAGCATTCTTTTCTATTAAAG_-
lihc	2_170033099_CTAGTGGAAAAG_-
lihc	8_88364024_T_C
lihc	17_7983975_C_A
lihc	X_14876077_C_T
lihc	13_28601380_T_C
lihc	4_55604594_G_T
lihc	7_86493597_G_T
lihc	16_9934503_C_A
coadread	5_112170863_G_A
coadread	5_112162804_G_A
coadread	2_61144153_G_A
coadread	11_108114847_T_C
coadread	7_50444492_G_T
coadread	17_7978912_C_A
coadread	12_78444543_G_A
coadread	10_104591369_C_T
coadread	6_152163731_G_T
coadread	2_141356198_A_G
coadread	12_78579379_A_G
coadread	7_81335745_T_A
coadread	5_180052868_C_T
coadread	7_81388122_T_C
coadread	11_32414302_C_T
coadread	2_170092579_C_A
coadread	5_112155042_G_C
coadread	11_108178712_G_T
skcm	12_18573874_G_A
skcm	2_141128412_C_T
skcm	2_212293132_C_T
skcm	7_99247856_C_T
skcm	16_9892322_C_T
skcm	10_43595906_G_A
skcm	5_19520765_C_T
skcm	11_92258116_T_A
skcm	12_21051369_G_A
skcm	9_8470993_A_G
skcm	20_9525015_C_T
skcm	9_8341692_C_T
skcm	2_141215031_C_T
skcm	10_96826964_C_T
skcm	9_8523525_C_T
skcm	2_141122354_T_C
skcm	16_10032409_C_T
skcm	12_21036537_G_A
skcm	7_99272209_C_T
skcm	7_99270202_C_T
prad	6_117663708_C_A
prad	13_28608023_C_G
kirc	2_61147174_A_C
kirc	17_41610308_C_T
kirc	6_117674152_C_A
kirc	X_19389075_A_C
kirc	3_48609479_C_T
gbm	2_170113633_C_T
gbm	2_141092130_T_C
######     END_SONG_DATA
######   OVERLOOKED_INSTANCES_OF_NO_EXPRESSION   REYKA_10_FEB_2015
lusc   NTRK2   TCGA-18-3409-01A-01D-0983-08
lusc   TOP3B   TCGA-21-5782-01A-01D-1632-08
lusc   LRRK2   TCGA-33-4566-01A-01D-1441-08
lusc   NOTCH1   TCGA-33-4583-01A-01D-1441-08
lusc   COL7A1   TCGA-37-5819-01A-01D-1632-08
luad   ROS1   TCGA-44-8119-01A-11D-2238-08
luad   CARD11   TCGA-55-5899-01A-11D-1625-08
lusc   TSC2   TCGA-63-5131-01A-01D-1441-08
brca   TP53   TCGA-A2-A1G1-01A-21D-A13L-09
brca   TP53   TCGA-A8-A075-01A-11D-A099-09
brca   TP53   TCGA-A8-A07O-01A-11W-A019-09
kirc   KDM5C   TCGA-AK-3458-01A-01D-1501-10
brca   NCOR1   TCGA-AN-A046-01A-21W-A050-09
brca   JAK1   TCGA-AN-A0AK-01A-21W-A019-09
brca   TP53   TCGA-AN-A0FL-01A-11W-A050-09
brca   KDM6A   TCGA-AR-A24Q-01A-12D-A167-09
brca   CDH1   TCGA-AR-A24X-01A-11D-A167-09
coadread   IPO7   TCGA-AZ-4615-01A-01D-1408-10
coadread   BRWD3   TCGA-AZ-6605-01A-11D-1835-10
kirc   VHL   TCGA-B0-4842-01A-02D-1421-08
kirc   PBRM1   TCGA-B0-4847-01A-01D-1361-10
brca   GPS2   TCGA-B6-A0IM-01A-11W-A050-09
brca   MAP2K4   TCGA-B6-A0RH-01A-21D-A10Y-09
brca   TP53   TCGA-BH-A1FC-01A-11D-A13L-09
kirc   SETD2   TCGA-BP-4343-01A-02D-1366-10
kirc   VHL   TCGA-BP-5170-01A-01D-1429-08
coadread   EPHB2   TCGA-CM-5861-01A-01D-1650-10
