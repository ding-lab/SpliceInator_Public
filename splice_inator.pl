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
   my $file_cases_full_expr = $ARGV[0];
my $file_controls_full_expr = $ARGV[1];
   my $file_cases_boundary_2bp = $ARGV[2]; 
   my $file_controls_boundary_2bp = $ARGV[3]; 
   my $file_cases_boundary_50bp = $ARGV[4]; 
   my $file_controls_boundary_50bp = $ARGV[5]; 
#__INPUT: WHOLE-GENE RSEM EXPRESSION OF CASES AND CONTROLS
   my $whole_gene_rsem_file = $ARGV[6];

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

   my $group_i_hits = {};
   my $overlooked_instance_of_no_expression = {};
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

#__DIAGNOSTIC VALUES --- EXPRESSION LEVELS OF CASE FOR EACH CATEGORY
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

#__DIAGNOSTIC VALUES --- EXPRESSION LEVELS OF CASE FOR EACH CATEGORY
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

