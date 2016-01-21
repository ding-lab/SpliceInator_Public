package SpliceInator;

use Carp;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

SpliceInator - collection of handy tools for splice site mutation classification

=head1 SYNOPSIS

	use SpliceInator;
	use SpliceInator ':all';
	use SpliceInator 'your_desired_routine';

=head1 DESCRIPTION

This is simply a collection of useful tools for the Ding Lab splice
project.
It is not a
class.
Individual routines are described
below.

=head1 METHODS

The available methods are as follows.

=cut

################################################################################
##                                                                            ##
##                      P R O G R A M M E R   N O T E S                       ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
##                      P R E L I M I N A R Y   C O D E                       ##
##                                                                            ##
################################################################################

############
#  SET-UP  #
############

#__STANDARD INCLUDES
   use strict;
   use warnings;
   use Carp;
   use constant PI => 4*atan2 1, 1;

#__SET UP EXPORTING
   use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
   require Exporter;
   @ISA = qw(Exporter);
   @EXPORT = qw();
   @EXPORT_OK = qw/
      parse_gene_rsem
      simple_intron_retention
      simple_exon_skipping
      read_full_expression
      gene_expression
      read_boundary_expression
      posterior_exon_shrinkage
      prior_exon_shrinkage
      posterior_exon_extension
      prior_exon_extension
   /;
   %EXPORT_TAGS = ( all => [
      qw/
         parse_gene_rsem
         simple_intron_retention
         simple_exon_skipping
	 	 gene_expression
         read_full_expression
         read_boundary_expression
         posterior_exon_shrinkage
         prior_exon_shrinkage
         posterior_exon_extension
         prior_exon_extension
      /
   ]);
   my $pkg = 'SpliceInator';

#############################
#  CONFIGURATION VARIABLES  #
#############################

#__VARIABLES
   our ($VERSION, $errtxt);
   $VERSION = '1.0';

####################
#  PRE-PROCESSING  #
####################

################################################################################
##                                                                            ##
##              P U B L I C   R O U T I N E S   A R E   H E R E               ##
##                                                                            ##
################################################################################

######################
#  PARSING ROUTINES  #
######################

#  ========================
#  READ BOUNDARY EXPRESSION    used for both cases and controls
#  ========================
#
#  PREVIOUS PROCESSING HAS ALREADY ESTABLISHED THE INDEXING, I.E. ACCORDING TO
#
#     il =      0       1       2       3                       intron list
#     el =  0       1       2       3       4                   exon list
#          E0, I0, E1, I1, E2, I2, E3, I3, E4
#     i = | 0 | 1   2 | 3   4 | 5   6 | 7   8 |.....|2N-1 2N |  list index
#         |___|_______|_______|_______|_______|.....|________|
#       k = 0     1       2       3       4             N       physical number
#
#       
#       i = 2 * il + 1       il = (i - 1) / 2      el = il + 1
#
#       e.g. i = 7: il = (7-1)/2 = 3  index of intron in intron list is 3
#                   el = 3+1          index of following exon in exon list is 4
#
#      $element_check->{$gene}  = 2N+1  IS ODD
#      $#expr_vals              = 2N    IS EVEN
#      $#expr_vals/2            = N     NUMBER OF *COMPLETE* INTRON->EXON PAIRS
#
#  WHERE THE DATA ARE GIVEN AS A STRING
#
#     E0,I0,E1,I1,E2,I2,E3,I3,....
#
#  WHILE HERE WE HAVE THE EXPRESSION DATA IN PAIRS (FOR LEFT 2BP/50BP AND THE
#  RIGHT 2BP/50BP OF EACH ELEMENT), I.E. THE DATA ARE FURNISHED AS THE STRING:
#
#     E0L,E0R,I0L,I0R,E1L,E1R,I1L,I1R,E2L,E2R,I2L,I2R,E3L,E3R,I3L,I3R,....
#
#  SO THE INDEXING WILL BE PRESERVED IF WE SIMPLY READ (AND STORE) THE
#  INFORMATION NOT AS A LIST, BUT RATHER AS A LIST OF LIST REFERENCES, WHERE
#  EACH SUB-LIST HAS PRECISELY 2 ELEMENTS: THE EXPRESSIONS AT THE LEFT AND
#  RIGHT BOUNDARIES --- THE 2 DATA STRUCTURES ARE THEREFORE
#
#     [   E0,        I0,        E1,        I1,        E2,        I2,....]
#
#     [[E0L,E0R], [I0L,I0R], [E1L,E1R], [I1L,I1R], [E2L,E2R], [I2L,I2R],...]

sub read_boundary_expression {
   my ($file, $save_positions, $spot_check) = @_;# SPOT_CHECK NOT YET IMPLEMENTD
   my ($data, $element_check) = ({}, {}); # ELEMENT_CHECK NOT YET IMPLEMENTED

#__READ CASE DATA
   my $i_line = 0;
   warn "reading '$file'\n";
   open (F, $file) || croak "cant open file '$file'";
   while (<F>) {

   #__FILTERING
      next if /^#/;
      chomp;
      $i_line++;
      warn "   reading line $i_line\n" unless $i_line % 50000;

   #__PARSE
      my ($cancer, $sample, $gene, $transcript,
          $strand, $mutation, $coord_string, $expression_string) = split /\t/;

   #__GET LIST OF EXPRESSION VALUES AND MAKE SURE NUMBER IS EVEN
   #  (EVEN IF THE NUMBER OF INTRONS+EXONS IS ODD, THE FACT THAT THERE ARE
   #  2 POINTS OF DATA FOR EACH MEANS THE RESULT MUST BE EVEN)
      my @expr_vals = split ',', $expression_string;
      croak "$cancer $sample $gene: detected odd number of expressions"
            if scalar @expr_vals % 2;

   #__SAVE POSITIONS IF DIRECTED TO
      my $coord_pairs;
      if ($save_positions) {
         @{$coord_pairs} = split ',', $coord_string;
         croak "$cancer $sample $gene: detected odd number of position pairs"
               if scalar @{$coord_pairs} % 2;
      }

   #__REVERSE LIST(S) IF STRAND IS NEGATIVE BECAUSE THE COMPUTATIONAL
   #  REPRESENTATION VIEWS EVERYTHING IN THE CONVENTIONAL LEFT-TO-RIGHT WAY
      if ($strand eq "-") {
         @expr_vals = reverse (@expr_vals);
         @{$coord_pairs} = reverse (@{$coord_pairs}) if $save_positions;
      }

   #__DISTILL EXPRESSION VALUES IN PAIRS (DOCUMENTED ABOVE)
      my $expression_mates = [];
      while (my ($lft, $rit) = splice (@expr_vals, 0, 2)) {
         push @{$expression_mates}, [$lft, $rit];
      }

   #__DISTILL POSITIONS IN PAIRS (DOCUMENTED ABOVE)
      my $coord_mates = [];
      if ($save_positions) {
         while (my ($lft, $rit) = splice (@{$coord_pairs}, 0, 2)) {
            push @{$coord_mates}, [$lft, $rit];
         }
      }

   #__ASSIGNMENT
      $data->{$cancer}->{$gene}->{$sample}->{'values'} = $expression_mates;
      if ($save_positions) {
         $data->{$cancer}->{$gene}->{$sample}->{'coord_pairs'} = $coord_mates;
      }
   }
   close (F);
   return $data;
}

#  ====================
#  READ FULL EXPRESSION    used for both cases and controls
#  ====================

sub read_full_expression {
   my ($file, $spot_check, $perform_indexing, $group_i_hits) = @_;
   my ($data, $element_check) = ({}, {});

#__READ CASE DATA
   my $i_line = 0;
   warn "reading '$file'\n";
   open (F, $file) || croak "cant open file '$file'";
   while (<F>) {

   #__FILTERING
      next if /^#/;
      chomp;
      $i_line++;
      warn "   reading line $i_line\n" unless $i_line % 50000;

   #__PARSE
      my ($cancer, $sample, $gene, $transcript,
          $strand, $mutation, $coord_string, $expression_string) = split /\t/;

   #__SKIP IF THIS IS IN AN UNEXPRESSED GENE
   #
   #    THIS WILL BE REMOVED ONCE WE IMPLEMENT A DEDICATED METHOD FOR
   #    IDENTIFYING CASES WHERE THE EXPRESSION OF THE WHOLE GENE HAS BEEN
   #    KNOCKED-OUT

      if (defined $group_i_hits->{$cancer}->{$mutation}) {
         push @{$group_i_hits->{$cancer}->{$mutation}}, [$sample, $gene];
#warn "skipping $cancer  $mutation --- not expressed\n";
#         croak "duplicates in unexpression: $cancer  $mutation" if
#            $group_i_hits->{$cancer}->{$mutation} > 1;
         next;
      }

   #__DISTILL LISTS OF POSITIONS AND EXPRESSION VALUES
      my @expr_vals = split ',', $expression_string;
      my $coord_pairs;
      @{$coord_pairs} = split ',', $coord_string;

   #__REVERSE THESE LISTS IF STRAND IS NEGATIVE BECAUSE THE COMPUTATIONAL
   #  REPRESENTATION VIEWS EVERYTHING IN THE CONVENTIONAL LEFT-TO-RIGHT WAY
      if ($strand eq "-") {
         @expr_vals = reverse (@expr_vals);
         @{$coord_pairs} = reverse (@{$coord_pairs});
      }

   #__DETERMINE LIST INDEX OF INTRON CONTAINING THE SPLICE-SITE MUTATION
   #   il =      0       1       2       3                       intron list
   #   el =  0       1       2       3       4                   exon list
   #        E0, I0, E1, I1, E2, I2, E3, I3, E4
   #   i = | 0 | 1   2 | 3   4 | 5   6 | 7   8 |.....|2N-1 2N |  list index
   #       |___|_______|_______|_______|_______|.....|________|
   #     k = 0     1       2       3       4             N       physical number
   #
   #       
   #     i = 2 * il + 1       il = (i - 1) / 2      el = il + 1
   #
   #     e.g. i = 7: il = (7-1)/2 = 3  index of intron in intron list is 3
   #                 el = 3+1          index of following exon in exon list is 4
   #
   #    $element_check->{$gene}  = 2N+1  IS ODD
   #    $#expr_vals              = 2N    IS EVEN
   #    $#expr_vals/2            = N     NUMBER OF *COMPLETE* INTRON->EXON PAIRS

      my ($i, $exon_flag, $errtxt, $il, $el);
      if (defined $perform_indexing) {
         ($i, $exon_flag, $errtxt) = _calc_list_index_ (
            $mutation, $coord_pairs, $cancer, $sample, $gene
         );
         if (defined $errtxt) {
            print "$errtxt --- $mutation --- skipping this one\n";
            next;
         }
         $il= ($i - 1)/2;
         $el = $il + 1;
      }

#  #__DISTILL THE EXPRESSION VALUES
#     my @expr_vals = split ',', $expression_string;

   #__SPOT CHECKING
      if ($spot_check) {

      #__MAKE SURE NUMBERS OF EXONS AND INTRONS ARE CONSISTENT ACROSS SAMPLES
         if (exists $element_check->{$gene}) {
            croak "$sample does not have $element_check->{$gene} " .
                  "exons/introns for gene $gene like the others"
            unless $element_check->{$gene} == scalar @expr_vals;
         } else {
            $element_check->{$gene} = scalar @expr_vals;
         }

      #__MAKE SURE ALL EXPRESSION VALUES ARE NUMERICAL
         foreach my $expr (@expr_vals) {
            croak "bad number '$expr'" unless $expr =~ /^[-+]?[0-9]*\.?[0-9]+$/
                       || $expr =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/;
         }
      }

   #__ASSIGNMENT
      $data->{$cancer}->{$gene}->{$sample}->{'values'} = [@expr_vals];
      if (defined $perform_indexing) {
         $data->{$cancer}->{$gene}->{$sample}->{'intron_list_index'} = $il;
         $data->{$cancer}->{$gene}->{$sample}->{'exon_list_index'} = $el;
         $data->{$cancer}->{$gene}->{$sample}->{'combined_list_index'} = $i;
         $data->{$cancer}->{$gene}->{$sample}->{'coord_pairs'} = $coord_pairs;
         $data->{$cancer}->{$gene}->{$sample}->{'exon_mutation'} = $exon_flag;
         $data->{$cancer}->{$gene}->{$sample}->{'strand'} = $strand;
         $data->{$cancer}->{$gene}->{$sample}->{'mut_position'} = $mutation;
      }
   }
   close (F);
   return $data;
}

=head2 parse_gene_rsem

Parses the RSEM values for a gene from a case and its corresponding
controls.
There is no companion accessor method - data
must be retrieved directly from the data
structure.

	Cancer_Type   Sample_ID   Gene_Name   Variant   RSEM (case)   RSEM (control)
	blca   TCGA-E5-A2PC-01A-11D-A202-08   STAG2   X_123182853_A_C   2392.9020   2284.3111,2521.0872,...,539.2120

	$rsem_data = parse_gene_rsem ($file);

	$rsem_data = {
	   blca => {
	      STAG2 => {
	         TCGA-E5-A2PC-01A-11D-A202-08 => {
	            case => 2392.9020
	            controls => [2284.3111, 2521.0872, ..., 539.2120]
                 },
              },
           },
	   :
	};

=cut

sub parse_gene_rsem {
   my ($file) = @_;
   my $rsem_data = {};

#__OPEN FILE AND READ
   open (F, $file) || croak "cant open file '$file'";
   while (<F>) {

   #__FILTERING
      next if /^#/ || /^Cancer\_Type/;
      next if /^s+/;
      chomp;

   #__PARSE
      my ($cancer, $sample, $gene, $variant, $rsem_case, $rsem_control) = split /\t/;
      my @control_values = split /,/, $rsem_control;

   #__STORE
      $rsem_data->{$cancer}->{$gene}->{$sample} = {
         'case' => $rsem_case,
         'controls' => [@control_values],
      };
   }
   close (F);

#__RETURN DATA
   return $rsem_data;
}

##########################
#  CALCULATION ROUTINES  #
##########################

#  =======================
#  GENE EXPRESSION (RSEM)
#  =======================

sub gene_expression {
  # my ($rsem_data) = @_;
   #my ($cases, $controls, $catalog,$threshold_classification, $diagnostics);
   my ($cases,$controls,$catalog,$threshold_classification,$diagnostics,$rsem_data)=@_;	 
   my ($group_i_count, $group_i_hits, $mode) = (0, {}, "gene expression");
      foreach my $cancer (keys %{$rsem_data}){
         foreach my $gene (keys %{$rsem_data->{$cancer}}){
            foreach my $sample (keys %{$rsem_data->{$cancer}->{$gene}}){
               #__PLACEHOLDER AND RELEVANT INDECES
                  my $local_hash = $cases->{$cancer}->{$gene}->{$sample};
	       #__EXTRACT CONTROL RSEM AND CASE RSEM VALUES  
		  my @controls = @{$rsem_data->{$cancer}->{$gene}->{$sample}->{'controls'}};
	          my $caseval = $rsem_data->{$cancer}->{$gene}->{$sample}->{'case'};
	       #__DIAGNOSTICS
	          $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $caseval;
	       #__PERMUTATION TEST - PREP DATA
                  my $pooled_list;
                  @{$pooled_list} = @{$rsem_data->{$cancer}->{$gene}->{$sample}->{'controls'}};
                  unshift @{$pooled_list}, $caseval;
               #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" GENE EXPRESSION LOWER-THAN-AVG
                  my ($pval_controls, $pval_control_2) =_permutation_test_(0, $pooled_list, "lower_than_avg");  
                  my $pval = $pval_controls;
					#my @POOL=@{$pooled_list};
                                	#print "$cancertype\t$genes\t$samples\t$caseval\t$pval\t@POOL\n";
	       ####################
               #  CLASSIFICATION  #
               ####################

               #__CATALOG
                  $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;

#              #__SKIP OLD-STYLE EVAL IF THIS ONE HAS ALREADY BEEN EVAL'D PREVIOUSLY
#              next if _already_have_pval_ ($catalog->{$cancer}->{$gene}->{$sample});

               #__EVALUATE CLASSIFICATION BASED ON THIS PVAL
                  my $test = "C";
		  if ($pval <= $threshold_classification) {
                     $group_i_hits->{$cancer}->{$gene}->{$sample} = [$pval, $test];
                     $group_i_count++;
	       }
            }
         }
      }
   return ($group_i_hits, $group_i_count,$diagnostics);
}

#  =======================
#  SIMPLE INTRON RETENTION
#  =======================

sub simple_intron_retention {
   my ($cases, $controls, $catalog, $min_elements,
       $expression_threshold,
       $spearman_threshold, $threshold_classification, $diagnostics) = @_;
   my ($group_ii_count, $group_ii_hits, $mode) = (0, {}, "intron retention");
   foreach my $cancer (keys %{$cases}) {
      foreach my $gene (keys %{$cases->{$cancer}}) {
         foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {

         #__PLACEHOLDER AND RELEVANT INDECES
            my $local_hash = $cases->{$cancer}->{$gene}->{$sample};
            my $i  = $local_hash->{'combined_list_index'};
            my $il = $local_hash->{'intron_list_index'};
            my $el = $local_hash->{'exon_list_index'};

         ##########################
         #  ACCESS RELEVANT DATA  #
         ##########################

         #__GET EXPRESSION LEVELS OF ALL INTRONS IN THIS SAMPLE
            my $expressions_sample = _get_expressions_in_sample_ (
               1, $cases->{$cancer}->{$gene}->{$sample}->{'values'}
            );

         #__GET EXPRESSION LEVELS IN ALL CORRESPONDING INTRONS IN CONTROLS
            my $expressions_controls =
                  _get_corresponding_expressions_in_controls_ (
                        $i, $controls->{$cancer}->{$gene}
            );

         ################
         #  EXCLUSIONS  #
         ################

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $expressions_sample->[$il];

         #__SKIP IF THIS IS AN EXON MUTATION (THESE ARE NOT SIMPLE RETENTIONS)
          #  if ($local_hash->{'exon_mutation'}) {
          #     $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
          #        = "no test --- exon mutation";
          #     next;
          #  }

         #__SKIP UNLESS INTRON EXPRESSION > AVG OF CONTROLS
            my $controls_avg = _average_ ($expressions_controls);
            unless ($expressions_sample->[$il] > $controls_avg) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- expression not > controls";
               next;
            }

         #__SKIP IF NOT ENOUGH CONTROLS
            if (scalar @{$expressions_controls} < $min_elements) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- not enough controls i.e. not > $min_elements";
               next;
            }

         #__SKIP UNLESS INTRON SEEMS ACTUALLY EXPRESSED
            if ($expressions_sample->[$il] < $expression_threshold) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- intron seems not to be expressed";
               next;
            }

         ############################
         #  CONTROLS-BASED P-VALUE  #
         ############################

         #__PREPEND SUBJECT INTRON EXPRESSION TO CONTROLS LIST: REQ'D FOR METHOD
            my $pooled_list;
            @{$pooled_list} = @{$expressions_controls};
            unshift @{$pooled_list}, $expressions_sample->[$il];

         #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" INTRON HIGHER-THAN-AVG
            my ($pval_controls, $pval_control_2) = _permutation_test_ (
                0, $pooled_list, "higher_than_avg"
            );
            my $pval = $pval_controls;
            my $test = "C";

         ###################################
         #  IN-SAMPLE P-VALUE IF IT HELPS  #
         ###################################

         #__GET MIDPOINTS OF EACH INTRON IN THIS SAMPLE
            my $midpoints = _get_midpoints_ (
               1, $cases->{$cancer}->{$gene}->{$sample}->{'coord_pairs'}
            );

         #__SPEARMAN'S RANK CORRELATION COEFFICIENT FOR SAMPLE INTRON EXPRESSION
            my $rc_obj = Statistics::RankCorrelation->new(
               $midpoints, $expressions_sample, sorted => 1
            );
            my $spearman_rho = $rc_obj->spearman;

         #__CONTROLS-ONLY P-VAL IF TOO MUCH 3' BIAS
            if (abs($spearman_rho) <= $spearman_threshold) {

            #__PERMUTATION TAILED P-VAL FOR IN-SAMPLE INTRON HIGHER-THAN-AVERAGE
               @{$pooled_list} = @{$expressions_sample};
               my ($pval_sample, $pval_sample_2) = _permutation_test_ (
                   $il, $pooled_list, "higher_than_avg"
               );

            #__INITIALIZATION FOR COMBINING SAMPLE AND CONTROL PVALS FOR INTRONS
               my $combine_obj = Statistics::CombinePvals->new (
                  [$pval_sample, $pval_controls]
               );

            #__COMBINE PVALS USING LANCASTER CORRECTION OF FISHER TRANSFORM METHOD
               my $pval_combined = $combine_obj->lancaster_mixed_corrected_transform (
                  [$pval_control_2, $pval_controls],
                  [$pval_sample_2, $pval_sample],
               );

            #__USE IT IF IT IMPROVES MATTERS
               if ($pval_combined < $pval) {
                  $pval = $pval_combined;
                  $test = "B";
               }
            }

         ####################
         #  CLASSIFICATION  #
         ####################

         #__CATALOG
            $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;

#        #__SKIP OLD-STYLE EVAL IF THIS ONE HAS ALREADY BEEN EVAL'D PREVIOUSLY
#           next if _already_have_pval_ ($catalog->{$cancer}->{$gene}->{$sample});

         #__EVALUATE CLASSIFICATION BASED ON THIS PVAL
            if ($pval <= $threshold_classification) {
               $group_ii_hits->{$cancer}->{$gene}->{$sample} = [$pval, $test];
               $group_ii_count++;
            }
         }
      }
   }
   return ($group_ii_hits, $group_ii_count, $diagnostics);
}

#  ====================
#  SIMPLE EXON SKIPPING
#  ====================

sub simple_exon_skipping {
   my ($cases, $controls, $catalog, $min_elements, $expression_threshold,
       $spearman_threshold, $threshold_classification, $diagnostics) = @_;
   my ($group_iii_count, $group_iii_hits, $mode) = (0, {}, "exon skipping");
   foreach my $cancer (keys %{$cases}) {
      foreach my $gene (keys %{$cases->{$cancer}}) {
         foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {

         #__PLACEHOLDER AND RELEVANT INDECES
            my $local_hash = $cases->{$cancer}->{$gene}->{$sample};
            my $i  = $local_hash->{'combined_list_index'};
            my $il = $local_hash->{'intron_list_index'};
            my $el = $local_hash->{'exon_list_index'};

         ##########################
         #  ACCESS RELEVANT DATA  #
         ##########################

         #__GET EXPRESSION LEVELS OF ALL EXONS IN THIS SAMPLE
            my $expressions_sample = _get_expressions_in_sample_ (
               0, $cases->{$cancer}->{$gene}->{$sample}->{'values'}
            );

         #__GET EXPRESSION LEVELS IN ALL CORRESPONDING EXONS IN CONTROLS
            my $expressions_controls =
                  _get_corresponding_expressions_in_controls_ (
                        $i+1, $controls->{$cancer}->{$gene}
            );

         ################
         #  EXCLUSIONS  #
         ################

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $expressions_sample->[$el];

         #__SKIP IF THIS IS AN EXON MUTATION (THESE ARE NOT SIMPLE SKIPPING)
           # if ($local_hash->{'exon_mutation'}) {
           #    $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
           #       = "no test --- exon mutation";
           #    next;
           # }

         #__SKIP UNLESS EXON EXPRESSION < AVG OF CONTROLS
            my $controls_avg = _average_ ($expressions_controls);
            unless ($expressions_sample->[$el] < $controls_avg) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- expression not < controls";
               next;
            }

         #__SKIP IF NOT ENOUGH CONTROLS
            if (scalar @{$expressions_controls} < $min_elements) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- not enough controls i.e. not > $min_elements";
               next;
            }

         #__SKIP UNLESS EXON SEEMS ACTUALLY EXPRESSED: THIS GUARDS AGAINST
         #  SAMPLES THAT WERE MISSED BY SONG, I.E. SAMPLES THAT REALLY ARE
         #  NOT EXPRESSED BUT WHICH SONG'S PROCESSING FAILED TO REMOVE
     #      my $sample_avg = _average_ ($expressions_sample);
     #      next if $expressions_sample->[$el] < $expression_threshold;
     #      next if $sample_avg < $expression_threshold;
     #
     #   TEST                       TP   FP
     #--------------------------------------
     #    no test                   60   27
     #    min 2 on test exon        28   13
     #    min 2 on exon AVG         45   17
     #    min 3 on test exon        20   11
     #    min 3 on exon AVG         32   12

         ############################
         #  CONTROLS-BASED P-VALUE  #
         ############################

         #__PREPEND SUBJECT EXON EXPRESSION TO CONTROLS LIST: REQ'D FOR METHOD
            my $pooled_list;
            @{$pooled_list} = @{$expressions_controls};
            unshift @{$pooled_list}, $expressions_sample->[$el];

         #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" EXON LOWER-THAN-AVG
            my ($pval_controls, $pval_control_2) = _permutation_test_ (
                0, $pooled_list, "lower_than_avg"
            );
            my $pval = $pval_controls;
            my $test = "C";

         ###################################
         #  IN-SAMPLE P-VALUE IF IT HELPS  #
         ###################################

         #__GET MIDPOINTS OF EACH EXON IN THIS SAMPLE
            my $midpoints = _get_midpoints_ (
               0, $cases->{$cancer}->{$gene}->{$sample}->{'coord_pairs'}
            );

         #__SPEARMAN'S RANK CORRELATION COEFFICIENT FOR SAMPLE EXON EXPRESSION
            my $rc_obj = Statistics::RankCorrelation->new(
               $midpoints, $expressions_sample, sorted => 1
            );
            my $spearman_rho = $rc_obj->spearman;

         #__CONTROLS-ONLY P-VAL IF TOO MUCH 3' BIAS
            if (abs($spearman_rho) <= $spearman_threshold) {

            #__PERMUTATION TAILED P-VAL FOR IN-SAMPLE EXON LOWER-THAN-AVERAGE
               @{$pooled_list} = @{$expressions_sample};
               my ($pval_sample, $pval_sample_2) = _permutation_test_ (
                   $el, $pooled_list, "lower_than_avg"
               );

            #__INITIALIZATION FOR COMBINING SAMPLE AND CONTROL PVALS FOR EXONS
               my $combine_obj = Statistics::CombinePvals->new (
                  [$pval_sample, $pval_controls]
               );

            #__COMBINE PVALS USING LANCASTER CORRECTION OF FISHER TRANSFORM METHOD
               my $pval_combined = $combine_obj->lancaster_mixed_corrected_transform (
                  [$pval_control_2, $pval_controls],
                  [$pval_sample_2, $pval_sample],
               );

            #__USE IT IF IT IMPROVES MATTERS
               if ($pval_combined < $pval) {
                  $pval = $pval_combined;
                  $test = "B";
               }
            }

         ####################
         #  CLASSIFICATION  #
         ####################

         #__CATALOG
            $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;

#        #__SKIP OLD-STYLE EVAL IF THIS ONE HAS ALREADY BEEN EVAL'D PREVIOUSLY
#           next if _already_have_pval_ ($catalog->{$cancer}->{$gene}->{$sample});

         #__EVALUATE CLASSIFICATION BASED ON THIS PVAL
            if ($pval <= $threshold_classification) {
               $group_iii_hits->{$cancer}->{$gene}->{$sample} = [$pval, $test];
               $group_iii_count++;
            }
         }
      }
   }
   return ($group_iii_hits, $group_iii_count, $diagnostics);
}

#  ========================
#  PRIOR EXON EXTENSION
#  ========================

sub prior_exon_extension {
   my ($index_hash, $cases, $controls, $mode, $catalog, $min_elements,
       $expression_threshold,
       $spearman_threshold, $threshold_classification, $diagnostics) = @_;
   foreach my $cancer (keys %{$cases}) {
      foreach my $gene (keys %{$cases->{$cancer}}) {
         foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {

         #__PLACEHOLDER AND RELEVANT INDECES
            my $local_hash = $index_hash->{$cancer}->{$gene}->{$sample};
            my $i  = $local_hash->{'combined_list_index'};
            my $il = $local_hash->{'intron_list_index'};
            my $el = $local_hash->{'exon_list_index'};

         #__SKIP IF THIS ONE WAS READ IN THE BOUNDARY DATA BUT WAS DISCARDED
         #  PREVIOUSLY BY SONG'S CUTOFF
            next unless defined $i && defined $il && defined $el;

         ##########################
         #  ACCESS RELEVANT DATA  #
         ##########################

         #__GET BOUNDARY PAIR EXPRESSION LEVELS OF ALL INTRONS IN THIS SAMPLE
            my $expressions_sample_bpairs = _get_expressions_in_sample_ (
               1, $cases->{$cancer}->{$gene}->{$sample}->{'values'}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_sample_lft, $expressions_sample_rit) = ([], []);
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_sample_lft}, $pair->[0];
               push @{$expressions_sample_rit}, $pair->[1];
            }

         #__GET BOUNDARY PAIR EXPR LEVELS IN ALL CORRESP INTRONS IN CONTROLS
            my $expressions_controls_bpairs =
                  _get_corresponding_expressions_in_controls_ (
                        $i, $controls->{$cancer}->{$gene}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_controls_lft, $expressions_controls_rit) = ([],[]);
            foreach my $pair (@{$expressions_controls_bpairs}) {
               push @{$expressions_controls_lft}, $pair->[0];
               push @{$expressions_controls_rit}, $pair->[1];
            }

         ################
         #  EXCLUSIONS  #
         ################

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $expressions_sample_lft->[$il];

         #__SKIP IF THIS IS AN EXON MUTATION
         #   if ($local_hash->{'exon_mutation'}) {
         #      $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
         #         = "no test --- exon mutation";
         #      next;
         #   }

         #__SKIP UNLESS INTRON EXPRESSION > AVG OF CONTROLS
            my $controls_avg = _average_ ($expressions_controls_lft);
            unless ($expressions_sample_lft->[$il] > $controls_avg) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- expression not > controls";
               next;
            }

         #__SKIP IF NOT ENOUGH CONTROLS
            if (scalar @{$expressions_controls_lft} < $min_elements) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- not enough controls i.e. not > $min_elements";
               next;
            }

         #__SKIP UNLESS INTRON BOUNDARY SEEMS ACTUALLY EXPRESSED
            if ($expressions_sample_lft->[$il] < $expression_threshold) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- intron boundary seems not to be expressed";
               next;
            }

         ############################
         #  CONTROLS-BASED P-VALUE  #
         ############################

         #__PREPEND SUBJECT INTRON EXPRESSION TO CONTROLS LIST: REQ'D FOR METHOD
            my $pooled_list;
            @{$pooled_list} = @{$expressions_controls_lft};
            unshift @{$pooled_list}, $expressions_sample_lft->[$il];

         #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" INTRON HIGHER-THAN-AVG
            my ($pval_controls, $pval_control_2) = _permutation_test_ (
                0, $pooled_list, "higher_than_avg"
            );
            my $pval = $pval_controls;
            my $test = "C";

         ###################################
         #  IN-SAMPLE P-VALUE IF IT HELPS  #
         ###################################
         #
         # DO A SPEARMAN TEST HERE BECAUSE IT MAY NOT HAVE BEEN RUN FOR
         # THE FULL EXPRESSION CALCULATION IF THAT CASE WAS EXCLUDED

         #__GET MIDPOINTS OF EACH LEFT AND EACH RIGHT BOUNDARY ON EACH INTRON IN
         #  THIS SAMPLE (IN ORDER)
            my $midpoints = _get_midpoints_boundary_regions_ (
               1, $cases->{$cancer}->{$gene}->{$sample}->{'coord_pairs'}
            );

         #__RE-DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
         #  IN ORDER TO PREPARE FOR MATCHED SPEARMAN CORRELATION
            my $expressions_ordered = [];
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_ordered}, @{$pair};
            }

         #__SPEARMAN'S RANK CORRELATION COEFFICIENT FOR SAMPLE INTRON EXPRESSION
            my $rc_obj = Statistics::RankCorrelation->new(
               $midpoints, $expressions_ordered, sorted => 1
            );
            my $spearman_rho = $rc_obj->spearman;

         #__CONTROLS-ONLY P-VAL IF TOO MUCH 3' BIAS
            if (abs($spearman_rho) <= $spearman_threshold) {

            #__PERMUTATION TAILED P-VAL FOR IN-SAMPLE INTRON HIGHER-THAN-AVERAGE
            #  USE DOUBLE-SIZE LIST OF ALL LEFT AND ALL RIGHT BOUNDARY VALUES
            #  WITH LEFT BEING THE FIRST LIST IN WHICH THE INDEX OF THE
            #  SUBJECT EXPRESSION VALUE IS "IL" (THEN ADD THE RIGHT BOUND LIST)
               @{$pooled_list} = @{$expressions_sample_lft};
               push @{$pooled_list}, @{$expressions_sample_rit};
               my ($pval_sample, $pval_sample_2) = _permutation_test_ (
                   $il, $pooled_list, "higher_than_avg"
               );

            #__INITIALIZATION FOR COMBINING SAMPLE AND CONTROL PVALS
               my $combine_obj = Statistics::CombinePvals->new (
                  [$pval_sample, $pval_controls]
               );

            #__COMBINE PVALS USING LANCASTER CORRECTION OF FISHER TRANSFORM METHOD
               my $pval_combined = $combine_obj->lancaster_mixed_corrected_transform (
                  [$pval_control_2, $pval_controls],
                  [$pval_sample_2, $pval_sample],
               );

            #__USE IT IF IT IMPROVES MATTERS
               if ($pval_combined < $pval) {
                  $pval = $pval_combined;
                  $test = "B";
               }
            }

         ####################
         #  CLASSIFICATION  #
         ####################

         #__CATALOG
            $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;
         }
      }
   }
   return ({}, 0, $diagnostics);
}

#  ========================
#  POSTERIOR EXON EXTENSION
#  ========================

sub posterior_exon_extension {
   my ($index_hash, $cases, $controls, $mode, $catalog, $min_elements,
       $expression_threshold,
       $spearman_threshold, $threshold_classification, $diagnostics) = @_;
   foreach my $cancer (keys %{$cases}) {
      foreach my $gene (keys %{$cases->{$cancer}}) {
         foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {

         #__PLACEHOLDER AND RELEVANT INDECES
            my $local_hash = $index_hash->{$cancer}->{$gene}->{$sample};
            my $i  = $local_hash->{'combined_list_index'};
            my $il = $local_hash->{'intron_list_index'};
            my $el = $local_hash->{'exon_list_index'};

         #__SKIP IF THIS ONE WAS READ IN THE BOUNDARY DATA BUT WAS DISCARDED
         #  PREVIOUSLY BY SONG'S CUTOFF
            next unless defined $i && defined $il && defined $el;

         ##########################
         #  ACCESS RELEVANT DATA  #
         ##########################

         #__GET BOUNDARY PAIR EXPRESSION LEVELS OF ALL INTRONS IN THIS SAMPLE
            my $expressions_sample_bpairs = _get_expressions_in_sample_ (
               1, $cases->{$cancer}->{$gene}->{$sample}->{'values'}
            );
#if ($il == 3) {
#   print "$cancer  $gene   $sample: LIST OF EXPRESSIONS\n";
#   foreach my $pair (@{$expressions_sample_bpairs}) {
#      print "$pair->[0]   $pair->[1]\n";
#   }
#}

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_sample_lft, $expressions_sample_rit) = ([], []);
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_sample_lft}, $pair->[0];
               push @{$expressions_sample_rit}, $pair->[1];
            }

         #__GET BOUNDARY PAIR EXPR LEVELS IN ALL CORRESP INTRONS IN CONTROLS
            my $expressions_controls_bpairs =
                  _get_corresponding_expressions_in_controls_ (
                        $i, $controls->{$cancer}->{$gene}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_controls_lft, $expressions_controls_rit) = ([],[]);
            foreach my $pair (@{$expressions_controls_bpairs}) {
               push @{$expressions_controls_lft}, $pair->[0];
               push @{$expressions_controls_rit}, $pair->[1];
            }

         ################
         #  EXCLUSIONS  #
         ################

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $expressions_sample_rit->[$il];

         #__SKIP IF THIS IS AN EXON MUTATION
         #   if ($local_hash->{'exon_mutation'}) {
         #      $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
         #         = "no test --- exon mutation";
         #      next;
         #   }

         #__SKIP UNLESS INTRON EXPRESSION > AVG OF CONTROLS
#if ($il == 3) {
#   print "$cancer  $gene   $sample: PRIOR EXON EXTENSION    IL = 3\n";
#   print "expression value in sample = $expressions_sample_rit->[$il]\n";
#   print "corresponding expression values in controls = ",
#          join (', ', @{$expressions_controls_rit}), "\n";
#   die "END OF CURRENT TEST";
#}
            my $controls_avg = _average_ ($expressions_controls_rit);
            unless ($expressions_sample_rit->[$il] > $controls_avg) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- expression not > controls";
               next;
            }

         #__SKIP IF NOT ENOUGH CONTROLS
            if (scalar @{$expressions_controls_rit} < $min_elements) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- not enough controls i.e. not > $min_elements";
               next;
            }

         #__SKIP UNLESS INTRON BOUNDARY SEEMS ACTUALLY EXPRESSED
            if ($expressions_sample_rit->[$il] < $expression_threshold) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- intron boundary seems not to be expressed";
               next;
            }

         ############################
         #  CONTROLS-BASED P-VALUE  #
         ############################

         #__PREPEND SUBJECT INTRON EXPRESSION TO CONTROLS LIST: REQ'D FOR METHOD
            my $pooled_list;
            @{$pooled_list} = @{$expressions_controls_rit};
            unshift @{$pooled_list}, $expressions_sample_rit->[$il];

         #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" INTRON HIGHER-THAN-AVG
            my ($pval_controls, $pval_control_2) = _permutation_test_ (
                0, $pooled_list, "higher_than_avg"
            );
            my $pval = $pval_controls;
            my $test = "C";

         ###################################
         #  IN-SAMPLE P-VALUE IF IT HELPS  #
         ###################################
         #
         # DO A SPEARMAN TEST HERE BECAUSE IT MAY NOT HAVE BEEN RUN FOR
         # THE FULL EXPRESSION CALCULATION IF THAT CASE WAS EXCLUDED

         #__GET MIDPOINTS OF EACH LEFT AND EACH RIGHT BOUNDARY ON EACH INTRON IN
         #  THIS SAMPLE (IN ORDER)
            my $midpoints = _get_midpoints_boundary_regions_ (
               1, $cases->{$cancer}->{$gene}->{$sample}->{'coord_pairs'}
            );

         #__RE-DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
         #  IN ORDER TO PREPARE FOR MATCHED SPEARMAN CORRELATION
            my $expressions_ordered = [];
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_ordered}, @{$pair};
            }

# print "$cancer  $gene  $sample: MIDPOINTS OF INTRON BOUNDARY REGIONS ",
#       join (', ', @{$midpoints}), "\n";
# print "CORRESPONDING EXPRESSIONS OF INTRON BOUNDARY REGIONS ",
#       join (', ', @{$expressions_ordered}), "\n";
# die "END OF INTRON TEST";

         #__SPEARMAN'S RANK CORRELATION COEFFICIENT FOR SAMPLE INTRON EXPRESSION
            my $rc_obj = Statistics::RankCorrelation->new(
               $midpoints, $expressions_ordered, sorted => 1
            );
            my $spearman_rho = $rc_obj->spearman;

         #__CONTROLS-ONLY P-VAL IF TOO MUCH 3' BIAS
            if (abs($spearman_rho) <= $spearman_threshold) {

            #__PERMUTATION TAILED P-VAL FOR IN-SAMPLE INTRON HIGHER-THAN-AVERAGE
            #  USE DOUBLE-SIZE LIST OF ALL LEFT AND ALL RIGHT BOUNDARY VALUES
            #  WITH RIGHT BEING THE FIRST LIST IN WHICH THE INDEX OF THE
            #  SUBJECT EXPRESSION VALUE IS "IL" (THEN ADD THE LEFT BOUND LIST)
               @{$pooled_list} = @{$expressions_sample_rit};
               push @{$pooled_list}, @{$expressions_sample_lft};
               my ($pval_sample, $pval_sample_2) = _permutation_test_ (
                   $il, $pooled_list, "higher_than_avg"
               );

            #__INITIALIZATION FOR COMBINING SAMPLE AND CONTROL PVALS
               my $combine_obj = Statistics::CombinePvals->new (
                  [$pval_sample, $pval_controls]
               );

            #__COMBINE PVALS USING LANCASTER CORRECTION OF FISHER TRANSFORM METHOD
               my $pval_combined = $combine_obj->lancaster_mixed_corrected_transform (
                  [$pval_control_2, $pval_controls],
                  [$pval_sample_2, $pval_sample],
               );

            #__USE IT IF IT IMPROVES MATTERS
               if ($pval_combined < $pval) {
                  $pval = $pval_combined;
                  $test = "B";
               }
            }

         ####################
         #  CLASSIFICATION  #
         ####################

         #__CATALOG
            $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;
         }
      }
   }
   return ({}, 0, $diagnostics);
}

#  ====================
#  PRIOR EXON SHRINKAGE
#  ====================
#
#  THE INDEXING IS DIFFERENT AS COMPARED TO THE OTHER ROUTINES BECAUSE HERE THE
#  MUTATION IS IN THE EXON ITSELF: ESSENTIALLY INDECES HAVE TO BE LEFT-SHIFTED,
#  AS ILLUSTRATED IN THIS EXAMPLE   (=== is exon, --- is intron)
#
#
#    ========--------=====X=---------======---------=========------=======---->
#  i    0        1      2        3     4       5       6       7      8     9
#  il            0               1             2               3            4
#  el   0               1              2               3              4
#
#    true index: i=2
#    BUT index as set in _calc_list_index_: i=3 (so that il and el are integers)
#
#                i - 1       3 - 1
#        il  =  -------  =  -------  = 1        el  =  il + 1  =  1 + 1  =  2
#                  2           2
#    so, we have to use i-1 in the "i list" and "el - 1" in the el list to
#    refer to the data in this case

sub prior_exon_shrinkage {
   my ($index_hash, $cases, $controls, $mode, $catalog, $min_elements,
       $expression_threshold,
       $spearman_threshold, $threshold_classification, $diagnostics) = @_;
   foreach my $cancer (keys %{$cases}) {
      foreach my $gene (keys %{$cases->{$cancer}}) {
         foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {

         #__PLACEHOLDER AND RELEVANT INDECES
            my $local_hash = $index_hash->{$cancer}->{$gene}->{$sample};
            my $i  = $local_hash->{'combined_list_index'};
            my $il = $local_hash->{'intron_list_index'};
            my $el = $local_hash->{'exon_list_index'};

         #__SKIP IF THIS ONE WAS READ IN THE BOUNDARY DATA BUT WAS DISCARDED
         #  PREVIOUSLY BY SONG'S CUTOFF
            next unless defined $i && defined $il && defined $el;

         ##########################
         #  ACCESS RELEVANT DATA  #
         ##########################

         #__GET BOUNDARY PAIR EXPRESSION LEVELS OF ALL EXONS IN THIS SAMPLE
            my $expressions_sample_bpairs = _get_expressions_in_sample_ (
               0, $cases->{$cancer}->{$gene}->{$sample}->{'values'}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_sample_lft, $expressions_sample_rit) = ([], []);
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_sample_lft}, $pair->[0];
               push @{$expressions_sample_rit}, $pair->[1];
            }

         #__GET BOUNDARY PAIR EXPR LEVELS IN ALL CORRESPONDING EXONS IN CONTROLS
         #  HERE I-1, NOT I+1 AS IN OTHER SIMPLE EXON OR SHRINKAGE ROUTINES
         #  (SEE NOTE ABOVE)
            my $expressions_controls_bpairs =
                  _get_corresponding_expressions_in_controls_ (
                        $i-1, $controls->{$cancer}->{$gene}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_controls_lft, $expressions_controls_rit) = ([],[]);
            foreach my $pair (@{$expressions_controls_bpairs}) {
               push @{$expressions_controls_lft}, $pair->[0];
               push @{$expressions_controls_rit}, $pair->[1];
            }

         ################
         #  EXCLUSIONS  #
         ################

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $expressions_sample_rit->[$el-1];
###COMMENTED OUT
         #__SKIP UNLESS THIS IS AN EXON MUTATION
         #       ^^^^^^
            unless ($local_hash->{'exon_mutation'}) {
               #$catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                #  = "no test --- not an exon mutation";
               #next;
            }

         #__SKIP UNLESS EXON EXPRESSION < AVG OF CONTROLS
         #  INDEX IS "EL-1", NOT THE USUAL "EL", AS NOTED ABOVE
            my $controls_avg = _average_ ($expressions_controls_rit);
            unless ($expressions_sample_rit->[$el-1] < $controls_avg) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- expression not < controls";
               next;
            }

         #__SKIP IF NOT ENOUGH CONTROLS
            if (scalar @{$expressions_controls_rit} < $min_elements) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- not enough controls i.e. not > $min_elements";
               next;
            }

         #__SKIP UNLESS EXON SEEMS ACTUALLY EXPRESSED: THIS GUARDS AGAINST
         #  SAMPLES THAT WERE MISSED BY SONG, I.E. SAMPLES THAT REALLY ARE
         #  NOT EXPRESSED BUT WHICH SONG'S PROCESSING FAILED TO REMOVE
     #      my $sample_avg = _average_ ($expressions_sample);
     #      next if $expressions_sample->[$el] < $expression_threshold;
     #      next if $sample_avg < $expression_threshold;
     #
     #   TEST                       TP   FP
     #--------------------------------------
     #    no test                   60   27
     #    min 2 on test exon        28   13
     #    min 2 on exon AVG         45   17
     #    min 3 on test exon        20   11
     #    min 3 on exon AVG         32   12

         ############################
         #  CONTROLS-BASED P-VALUE  #
         ############################

         #__PREPEND SUBJECT EXON EXPRESSION TO CONTROLS LIST: REQ'D FOR METHOD
         #  INDEX IS "EL-1", NOT THE USUAL "EL", AS NOTED ABOVE
            my $pooled_list;
            @{$pooled_list} = @{$expressions_controls_rit};
            unshift @{$pooled_list}, $expressions_sample_rit->[$el-1];

         #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" EXON LOWER-THAN-AVG
            my ($pval_controls, $pval_control_2) = _permutation_test_ (
                0, $pooled_list, "lower_than_avg"
            );
            my $pval = $pval_controls;
            my $test = "C";

         ###################################
         #  IN-SAMPLE P-VALUE IF IT HELPS  #
         ###################################
         #
         # DO A SPEARMAN TEST HERE BECAUSE IT MAY NOT HAVE BEEN RUN FOR
         # THE FULL EXPRESSION CALCULATION IF THAT CASE WAS EXCLUDED

         #__GET MIDPOINTS OF EACH LEFT AND EACH RIGHT BOUNDARY ON EACH EXON IN
         #  THIS SAMPLE (IN ORDER)
            my $midpoints = _get_midpoints_boundary_regions_ (
               0, $cases->{$cancer}->{$gene}->{$sample}->{'coord_pairs'}
            );

         #__RE-DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
         #  IN ORDER TO PREPARE FOR MATCHED SPEARMAN CORRELATION
            my $expressions_ordered = [];
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_ordered}, @{$pair};
            }

         #__SPEARMAN'S RANK CORRELATION COEFFICIENT FOR SAMPLE EXON EXPRESSION
            my $rc_obj = Statistics::RankCorrelation->new(
               $midpoints, $expressions_ordered, sorted => 1
            );
            my $spearman_rho = $rc_obj->spearman;

         #__CONTROLS-ONLY P-VAL IF TOO MUCH 3' BIAS
            if (abs($spearman_rho) <= $spearman_threshold) {

            #__PERMUTATION TAILED P-VAL FOR IN-SAMPLE EXON LOWER-THAN-AVERAGE
            #  USE DOUBLE-SIZE LIST OF ALL LEFT AND ALL RIGHT BOUNDARY VALUES
            #  WITH RIGHT BEING THE FIRST LIST IN WHICH THE INDEX OF THE
            #  SUBJECT EXPRESSION VALUE IS "EL-1" (THEN ADD THE LEFT BOUND LIST)
            #
            #  REMINDER: INDEX IS "EL-1", NOT THE USUAL "EL", AS NOTED ABOVE
               @{$pooled_list} = @{$expressions_sample_rit};
               push @{$pooled_list}, @{$expressions_sample_lft};
               my ($pval_sample, $pval_sample_2) = _permutation_test_ (
                   $el-1, $pooled_list, "lower_than_avg"
               );

            #__INITIALIZATION FOR COMBINING SAMPLE AND CONTROL PVALS FOR EXONS
               my $combine_obj = Statistics::CombinePvals->new (
                  [$pval_sample, $pval_controls]
               );

            #__COMBINE PVALS USING LANCASTER CORRECTION OF FISHER TRANSFORM METHOD
               my $pval_combined = $combine_obj->lancaster_mixed_corrected_transform (
                  [$pval_control_2, $pval_controls],
                  [$pval_sample_2, $pval_sample],
               );

            #__USE IT IF IT IMPROVES MATTERS
               if ($pval_combined < $pval) {
                  $pval = $pval_combined;
                  $test = "B";
               }
            }

         ####################
         #  CLASSIFICATION  #
         ####################

         #__CATALOG
            $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;
         }
      }
   }
   return ({}, 0, $diagnostics);
}

#  ========================
#  POSTERIOR EXON SHRINKAGE
#  ========================

sub posterior_exon_shrinkage {
   my ($index_hash, $cases, $controls, $mode, $catalog, $min_elements,
       $expression_threshold,
       $spearman_threshold, $threshold_classification, $diagnostics) = @_;
   foreach my $cancer (keys %{$cases}) {
      foreach my $gene (keys %{$cases->{$cancer}}) {
         foreach my $sample (keys %{$cases->{$cancer}->{$gene}}) {

         #__PLACEHOLDER AND RELEVANT INDECES
            my $local_hash = $index_hash->{$cancer}->{$gene}->{$sample};
            my $i  = $local_hash->{'combined_list_index'};
            my $il = $local_hash->{'intron_list_index'};
            my $el = $local_hash->{'exon_list_index'};

         #__SKIP IF THIS ONE WAS READ IN THE BOUNDARY DATA BUT WAS DISCARDED
         #  PREVIOUSLY BY SONG'S CUTOFF
            next unless defined $i && defined $il && defined $el;

         ##########################
         #  ACCESS RELEVANT DATA  #
         ##########################

         #__GET BOUNDARY PAIR EXPRESSION LEVELS OF ALL EXONS IN THIS SAMPLE
            my $expressions_sample_bpairs = _get_expressions_in_sample_ (
               0, $cases->{$cancer}->{$gene}->{$sample}->{'values'}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_sample_lft, $expressions_sample_rit) = ([], []);
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_sample_lft}, $pair->[0];
               push @{$expressions_sample_rit}, $pair->[1];
            }

         #__GET BOUNDARY PAIR EXPR LEVELS IN ALL CORRESPONDING EXONS IN CONTROLS
            my $expressions_controls_bpairs =
                  _get_corresponding_expressions_in_controls_ (
                        $i+1, $controls->{$cancer}->{$gene}
            );

         #__DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
            my ($expressions_controls_lft, $expressions_controls_rit) = ([],[]);
            foreach my $pair (@{$expressions_controls_bpairs}) {
               push @{$expressions_controls_lft}, $pair->[0];
               push @{$expressions_controls_rit}, $pair->[1];
            }

         ################
         #  EXCLUSIONS  #
         ################

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $expressions_sample_lft->[$el];
######COMMENTED OUT
         #__SKIP IF THIS IS AN EXON MUTATION (THESE ARE NOT SIMPLE SKIPPING)
     #       if ($local_hash->{'exon_mutation'}) {
     #          $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
     #             = "no test --- exon mutation";
     #          next;
     #       }

         #__SKIP UNLESS EXON EXPRESSION < AVG OF CONTROLS
            my $controls_avg = _average_ ($expressions_controls_lft);
            unless ($expressions_sample_lft->[$el] < $controls_avg) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- expression not < controls";
               next;
            }

         #__SKIP IF NOT ENOUGH CONTROLS
            if (scalar @{$expressions_controls_lft} < $min_elements) {
               $catalog->{$cancer}->{$gene}->{$sample}->{$mode}
                  = "no test --- not enough controls i.e. not > $min_elements";
               next;
            }

         #__SKIP UNLESS EXON SEEMS ACTUALLY EXPRESSED: THIS GUARDS AGAINST
         #  SAMPLES THAT WERE MISSED BY SONG, I.E. SAMPLES THAT REALLY ARE
         #  NOT EXPRESSED BUT WHICH SONG'S PROCESSING FAILED TO REMOVE
     #      my $sample_avg = _average_ ($expressions_sample);
     #      next if $expressions_sample->[$el] < $expression_threshold;
     #      next if $sample_avg < $expression_threshold;
     #
     #   TEST                       TP   FP
     #--------------------------------------
     #    no test                   60   27
     #    min 2 on test exon        28   13
     #    min 2 on exon AVG         45   17
     #    min 3 on test exon        20   11
     #    min 3 on exon AVG         32   12

         ############################
         #  CONTROLS-BASED P-VALUE  #
         ############################

         #__PREPEND SUBJECT EXON EXPRESSION TO CONTROLS LIST: REQ'D FOR METHOD
            my $pooled_list;
            @{$pooled_list} = @{$expressions_controls_lft};
            unshift @{$pooled_list}, $expressions_sample_lft->[$el];

         #__PERMUTATION TAILED P-VAL FOR "IN CONTROLS" EXON LOWER-THAN-AVG
            my ($pval_controls, $pval_control_2) = _permutation_test_ (
                0, $pooled_list, "lower_than_avg"
            );
            my $pval = $pval_controls;
            my $test = "C";

         ###################################
         #  IN-SAMPLE P-VALUE IF IT HELPS  #
         ###################################
         #
         # DO A SPEARMAN TEST HERE BECAUSE IT MAY NOT HAVE BEEN RUN FOR
         # THE FULL EXPRESSION CALCULATION IF THAT CASE WAS EXCLUDED

         #__GET MIDPOINTS OF EACH LEFT AND EACH RIGHT BOUNDARY ON EACH EXON IN
         #  THIS SAMPLE (IN ORDER)
            my $midpoints = _get_midpoints_boundary_regions_ (
               0, $cases->{$cancer}->{$gene}->{$sample}->{'coord_pairs'}
            );

         #__RE-DISTILL ALL LEFT BOUNDARY AND RIGHT BOUNDARY EXPRESSION LEVELS
         #  IN ORDER TO PREPARE FOR MATCHED SPEARMAN CORRELATION
            my $expressions_ordered = [];
            foreach my $pair (@{$expressions_sample_bpairs}) {
               push @{$expressions_ordered}, @{$pair};
            }

# print "$cancer  $gene  $sample: MIDPOINTS OF EXON BOUNDARY REGIONS ",
#       join (', ', @{$midpoints}), "\n";
# print "CORRESPONDING EXPRESSIONS OF EXON BOUNDARY REGIONS ",
#       join (', ', @{$expressions_ordered}), "\n";
# die "END OF TEST";

         #__SPEARMAN'S RANK CORRELATION COEFFICIENT FOR SAMPLE EXON EXPRESSION
            my $rc_obj = Statistics::RankCorrelation->new(
               $midpoints, $expressions_ordered, sorted => 1
            );
            my $spearman_rho = $rc_obj->spearman;

         #__CONTROLS-ONLY P-VAL IF TOO MUCH 3' BIAS
            if (abs($spearman_rho) <= $spearman_threshold) {

            #__PERMUTATION TAILED P-VAL FOR IN-SAMPLE EXON LOWER-THAN-AVERAGE
            #  USE DOUBLE-SIZE LIST OF ALL LEFT AND ALL RIGHT BOUNDARY VALUES
            #  WITH LEFT BEING THE FIRST LIST IN WHICH THE INDEX OF THE
            #  SUBJECT EXPRESSION VALUE IS "EL" (THEN ADD THE RIGHT BOUND LIST)
               @{$pooled_list} = @{$expressions_sample_lft};
               push @{$pooled_list}, @{$expressions_sample_rit};
               my ($pval_sample, $pval_sample_2) = _permutation_test_ (
                   $el, $pooled_list, "lower_than_avg"
               );

            #__INITIALIZATION FOR COMBINING SAMPLE AND CONTROL PVALS FOR EXONS
               my $combine_obj = Statistics::CombinePvals->new (
                  [$pval_sample, $pval_controls]
               );

            #__COMBINE PVALS USING LANCASTER CORRECTION OF FISHER TRANSFORM METHOD
               my $pval_combined = $combine_obj->lancaster_mixed_corrected_transform (
                  [$pval_control_2, $pval_controls],
                  [$pval_sample_2, $pval_sample],
               );

            #__USE IT IF IT IMPROVES MATTERS
               if ($pval_combined < $pval) {
                  $pval = $pval_combined;
                  $test = "B";
               }
            }

         ####################
         #  CLASSIFICATION  #
         ####################

         #__CATALOG
            $catalog->{$cancer}->{$gene}->{$sample}->{$mode} = $pval;
         }
      }
   }
   return ({}, 0, $diagnostics);
}

################################################################################
##                                                                            ##
##             P R I V A T E   R O U T I N E S   A R E   H E R E              ##
##                                                                            ##
################################################################################

#  ================
#  PERMUTATION_TEST   CERTIFIED 19-DEC-2014
#  ================
#
#
#   ******************************************************************
#   *                                                                *
#   *  NOTE THAT THIS ROUTINE IS LIABLE TO *CHANGE* THE INPUT LIST,  *
#   *  SPECIFICALLY THE ORDER OF ITS ELEMENTS, WHICH CAN HAVE BAD    *
#   *  SIDE-EFFECTS IF THAT LIST IS TO BE USED ELSEWHERE             *
#   *                                                                *
#   ******************************************************************
#
#  
#  BASIC PERMUTATION TEST FOR A LIST WHERE ONE ELEMENT IS TESTED AGAINST
#  THE OTHERS AS A COLLECTIVE
#
#  NOTE THE EXISTENCE OF Math::Combinatorics FOR EXAMPLE AS USED IN
#  ../../in_progress/pancan_2/pvalue_mann_whitney_association.pl
#  BUT I DON'T THINK WE NEED THAT LEVEL OF SOPHISTICATION HERE FOR TESTING
#  SINGLE SAMPLES AGAINST A COLLECTIVE.
#
#        use Math::Combinatorics;
#
#     #__GET ALL COMBINATIONS OF THIS MANY CANCERS
#        my $cancer_combinations = Math::Combinatorics->new (
#           count => $num_cancers_in_test,
#           data => [@all_cancers],
#        );
#
#  HOWEVER, WE MAY FIND IT USEFUL WHEN WE TEST A SPECIFIC SPLICE SITE MUTATION
#  ACROSS MULTIPLE CASES

sub _permutation_test_ {
   my ($i_observation, $expressions, $test_type) = @_;

#__EXTRACT THE TEST VALUE FROM THE LIST (NOT ITS INDEX, BUT THE VALUE ITSELF)
   my $test_value = splice @{$expressions}, $i_observation, 1;

#__COMPUTE AVERAGE VALUE OF BACKGROUND DISTRIBUTION
   my $avg = _average_ ($expressions);

#__TEST STATISTIC IS THE DIFFERENCE OF THE AVERAGES (TEST VALUE ITSELF *IS* THE
#  AVERAGE FOR 1 SAMPLE). NOTE: THE 'DIRECTION' (ONE-TAIL HIGHER, ONE-TAIL
#  LOWER, OR TWO-TAIL) IS TAKEN INTO CONSIDERATION IN THE TALLYING PHASE BELOW
#  AND WE MAKE NO ACCOUNT HERE OF ITS SIGN

   my $statistic = $test_value - $avg;

#__POOL ALL OBSERVATIONS TO PREPARE FOR PERMUTATION: ORDER IS NOW UN-IMPORTANT
   push @{$expressions}, $test_value;

#__NULL DISTRIBUTION OF THE STATISTIC: EVERY COMBINATION OF 1 VAL VS BACKGROUND
#
#  there are a total of N elements in the pooled data, from which we
#  we are testing 1 vs the other N-1, whereby there are exactly
#  C_{N,1} = N combinations that can be tested. There are N! permutations
#  but is that just artificially scaling the null? Yes. You can see this very
#  easily graphically (not shown here).  18-DEC-2014
#
   my $h0_null_distrib = {};
   for (my $j = 0; $j <= $#{$expressions}; $j++) {

   #__MAKE THE TEST VALUE AND THE NULL LIST FOR THIS "PERMUTATION"
      my $permut_single_val = $expressions->[$j];
      my @permut_list = ();
      if ($j == 0) {
         @permut_list = @{$expressions}[1 .. $#{$expressions}];
      } elsif ($j == $#{$expressions}) {
         @permut_list = @{$expressions}[0 .. $j-1];
      } else {
         @permut_list = @{$expressions}[0 .. $j-1];
         push @permut_list, @{$expressions}[$j+1 .. $#{$expressions}];
      }

   #__GET THE AVERAGE VALUE OF THIS LIST
      my $permut_avg = _average_ ([@permut_list]);

   #__CALCULATE THE NULL TEST STATISTIC FOR THIS PERMUTATION
      my $permut_statistic = $permut_single_val - $permut_avg;

   #__TALLY THIS TO THE NULL DISTRIBUTION
      $h0_null_distrib->{$permut_statistic}++;
   }

#__PVALUE SET-UP
   my $total_permuts = scalar @{$expressions};
   my ($permuts_eq_or_more_extreme, $permuts_more_extreme) = (0, 0);

#__ONE-SIDE TEST: "OBSERVATION OR ANYTHING HIGHER" PROBABILITY - STRAIGHTFORWARD
   if ($test_type eq "higher_than_avg") {
      ($permuts_eq_or_more_extreme, $permuts_more_extreme) = _tail_right_ (
         $statistic, $h0_null_distrib
      );

#__ONE-SIDE TEST: "OBSERVATION OR ANYTHING LOWER" PROBABILITY - STRAIGHTFORWARD
   } elsif ($test_type eq "lower_than_avg") {
      ($permuts_eq_or_more_extreme, $permuts_more_extreme) = _tail_left_ (
         $statistic, $h0_null_distrib
      );

#__TWO-SIDED REQUIRES A LITTLE MORE CALCULATION TO PROPERLY DO "OTHER SIDE"
#
#  Note that the mean of the null distribution will ALWAYS (i.e. necessarily)
#  be 0 because the sum of all the elements will be 0. This is easy to
#  illustrate by an example pooled set of 4 members: a, b, c, d, the sum of the
#  null statistic elements being
#
#     a+b+c          b+c+d           c+d+a           d+a+b
#     -----  -  d  + -----  -  a  +  -----  -  b  +  -----  -  c  =  0
#       3              3               3               3
#
#  which clearly holds in general for a list of any size, N, since each element
#  will have a sum of N ratios of itself divided by N, from which will be
#  subtracted the element itself. THEREFORE, the "polarity" of the two-sided
#  calculation can be determined by just looking at the sign of the test
#  statistic, i.e. whether it is positive or negative.
 
   } elsif ($test_type eq "two_sided") {

   #__IF TEST STATISTIC IS "+", SUM THE RIGHT TAIL, NEGATE AND SUM THE LEFT TAIL
      if ($statistic > 0) {
         ($permuts_eq_or_more_extreme, $permuts_more_extreme) = _tail_right_ (
            $statistic, $h0_null_distrib
         );
         $statistic *= -1;
         my ($local_eq_or_more_extreme, $local_more_extreme) = _tail_left_ (
            $statistic, $h0_null_distrib
         );
         $permuts_eq_or_more_extreme += $local_eq_or_more_extreme;
         $permuts_more_extreme += $local_more_extreme;

   #__ELSIF STATISTIC IS "-", SUM THE LEFT TAIL, NEGATE AND SUM THE RIGHT TAIL
      } elsif ($statistic < 0) {
         ($permuts_eq_or_more_extreme, $permuts_more_extreme) = _tail_left_ (
            $statistic, $h0_null_distrib
         );
         $statistic *= -1;
         my ($local_eq_or_more_extreme, $local_more_extreme) = _tail_right_ (
            $statistic, $h0_null_distrib
         );
         $permuts_eq_or_more_extreme += $local_eq_or_more_extreme;
         $permuts_more_extreme += $local_more_extreme;

   #__ELSE STATISTIC IS EXACTLY EQUAL TO MEAN OF NULL SO PVAL=1
      } else {
         $permuts_eq_or_more_extreme += $total_permuts;
      }

#__ELSE THIS WAS A BAD ARGUMENT
   } else {
      croak "_permutation_ method does not understand the argument '$test_type'";
   }

#__PVALUE IS PROPORTION OF VALUE OF TEST STAT OR ANYTHING MORE EXTREME
   if (wantarray) {
      return ($permuts_eq_or_more_extreme / $total_permuts,
              $permuts_more_extreme / $total_permuts);
   } else {
      return $permuts_eq_or_more_extreme / $total_permuts;
   }
}

#  =====================
#  TAIL TALLYING METHODS
#  =====================
#
#  These are both configured to work in conjunction with the methods in
#  Statistics::CombinePvals, specifically with the Lancaster correction methods
#  for combining tailed P-values of discrete distributions. Basically,
#  these methods rely not only on the tail of the observed statistic, but also
#  the next most extreme P-value in the tail. This requirement is implemented
#  here in tallying for everything >= or <= the actual statistic, i.e. the
#  usual P-value tail, but also for everything > or < the actual statistic,
#  which would be next extreme P-value.
#
#  PERL ODDITY: this does not always work...value is same as $tallies
#
#       $next_extreme += $distribution->{$null_element}
#                        if $null_element < $statistic;
#
#       $tallies += $distribution->{$null_element}
#                        if $null_element <= $statistic;
#
#    this is a problem when there are values on the edge of standard
#    double-precision, e.g. 0.633943546309715, where Perl's equality and/or
#    inequality tests seem to fail. We *could* switch to arbitrary precision
#    using Math::BigFloat, but this would tend to slow processing speed
#    considerably. So, to fix this:
#
#    * in cases where we are looking for equality, we add a string comparison
#      using "eq" to the numerical comparison, so that it will catch the
#      equality in a string context
#
#    * to get the next most extreme pval (for the Lancaster correction), we
#      remove the numerical test from within the loop and simply replace that
#      with a simple posterior decrement of the equivalence value

sub _tail_left_ {
   my ($statistic, $distribution) = @_;
   my ($tallies, $next_extreme) = (0, 0);
   foreach my $null_element (sort _numerical_ keys %{$distribution}) {
      $tallies += $distribution->{$null_element} if $null_element <= $statistic
                                                 || $null_element eq $statistic;
   }
   $next_extreme = $tallies - $distribution->{$statistic};
   croak "found 0 Pvalue: impossible for permutation test" unless $tallies;
   return ($tallies, $next_extreme);
}

sub _tail_right_ {
   my ($statistic, $distribution) = @_;
   my ($tallies, $next_extreme) = (0, 0);
   foreach my $null_element (sort _numerical_ keys %{$distribution}) {
      $tallies += $distribution->{$null_element} if $null_element >= $statistic
                                                 || $null_element eq $statistic;
   }
   $next_extreme = $tallies - $distribution->{$statistic};
   croak "found 0 Pvalue: impossible for permutation test" unless $tallies;
   return ($tallies, $next_extreme);
}

######################
#  UTILITTY METHODS  #
######################

sub _numerical_ {$a <=> $b}

#  CERTIFIED 18-DEC-2014

sub _average_ {
   my ($list) = @_;
   my $sum = 0;
#  foreach my $element (sort _numerical_ @{$list}) {} # DOES NOT CURE FLOATING
#                          POINT BUG THAT GIVES DIFF ANSWER DEPENDING UPON ORDER
   foreach my $element (@{$list}) {
      $sum += $element;
   }
   return $sum / scalar @{$list};
}

sub _calc_list_index_ {
   my ($mutation, $coord_pairs, $cancer, $sample, $gene) = @_;
   my $exon_mutation_flag = 0;

#__DISTILL THE EXACT LOCATION OF THE SPLICE SITE MUTATION
   my $position;
   if ($mutation =~ /\S+\_(\d+)\_\S+\_\S+/) {
      $position = $1;
   } else {
      return (0, 0, "do not recognize locator string '$mutation' for $cancer $sample $gene");
   }

#__THE BEGIN-END COORDINATES OF EACH ELEMENT
   my @begin_end_vals = @{$coord_pairs};

#__GO THROUGH EACH BEGIN-END PAIR TO FIND THE INDEX OF THE MUTATION
#
#  THE MOST EFFICIENT METHOD OF FINDING THIS INDEX WOULD BE BY THE BISECTION
#  METHOD --- IMPLEMENT THIS AFTER THE BASIC SCRIPT IS COMPLETE
#
#  RIGHT NOW WE JUST USE NAIVE SEARCH

   my ($i, $found_it) = (0, 0);
   my ($begin_pos, $end_pos);
   foreach my $begin_end_pair (@begin_end_vals) {
      if ($begin_end_pair =~ /(\d+)\-(\d+)/) {
         ($begin_pos, $end_pos) = ($1, $2);
         croak "bad begin/end" unless $end_pos >= $begin_pos;
         if ($position >= $begin_pos && $position <= $end_pos) {
            $found_it = 1;
            last;
         }
      } else {
         return (0, 0, "do not recognize begin-end pair '$begin_end_pair' for $cancer $sample $gene");
      }
      $i++;
   }

#__CHECK THAT SOMETHING HAS BEEN FOUND
   return (0, 0, "$cancer $sample $gene cound not locate splice site mutation")
         unless $found_it;

#__WE HAVE AN EXON MUTATION UNLESS THE INDEX IS ODD (WHERE ITS AN INTRON MUT)
   unless ($i % 2) {

   #__SET INDEX TO THE NEXT INTRON
      $i++;

   #__SET FLAG THAT THIS IS AN EXON MUTATION
      $exon_mutation_flag = 1;
   }
#  return (0, 0, "$cancer $sample $gene splice does not appear to be in an intron")
#        unless $i % 2;

#__RETURN
   return ($i, $exon_mutation_flag, undef);
}

####################
#  ACCESS METHODS  #
####################
#
#  these are necessary because the input file intersperses intron and exon
#  values (expression levels and locations) and also because we need to access
#  corresponding values in other samples in the control group.

#  CERTIFIED 19-DEC-2014

sub _get_expressions_in_sample_ {
   my ($start, $mixed_expressions) = @_;
   my $needed_expressions = [];

#__RAW EXPRESSION LIST IS STAGGERED: EXON,INTRON,EXON,INTRON ETC SO SKIP BY TWO
   for (my $i = $start; $i <= $#$mixed_expressions; $i += 2) {
      push @{$needed_expressions}, $mixed_expressions->[$i];
   }
   return $needed_expressions;
}

#  CERTIFIED 20-DEC-2014

sub _get_corresponding_expressions_in_controls_ {
   my ($index, $sample_hash) = @_;
   my $expressions = [];

#__GET EXPRESSION OF DESIGNATED ELEMENT FROM EVERY SAMPLE
#  OF THIS GENE-CANCER COMBINATION
   foreach my $sample (keys %{$sample_hash}) {
      push @{$expressions}, $sample_hash->{$sample}->{'values'}->[$index];
   }
   return $expressions;
}

#  CERTIFIED 26-DEC-2014

sub _get_midpoints_ {
   my ($start, $coord_pairs) = @_;
   my $midpoints = [];

#__THE BEGIN-END COORDINATES OF EACH ELEMENT
   my @begin_end_vals = @{$coord_pairs};

#__RAW COORDINATE LIST IS STAGGERED: EXON,INTRON,EXON,INTRON ETC SO SKIP BY TWO
   my ($begin_pos, $end_pos);
   for (my $i = $start; $i <= $#begin_end_vals; $i += 2) {
      if ($begin_end_vals[$i] =~ /(\d+)\-(\d+)/) {
         ($begin_pos, $end_pos) = ($1, $2);
### print "  $begin_pos  $end_pos  -->  ", int (($begin_pos + $end_pos)/2), "\n";
#        push @{$midpoints}, int (($begin_pos + $end_pos)/2);
         push @{$midpoints}, ($begin_pos + $end_pos)/2;
      } else {
         croak "do not recognize begin-end pair --- needs diagnosis";
      }
   }
   return $midpoints;
}

#   The expression levels for boundary values are stored as
#
#      [[E0L,E0R], [I0L,I0R], [E1L,E1R], [I1L,I1R], [E2L,E2R], [I2L,I2R],...]
#
#   and the corresponding positions are stored as
#
#   [
#      [40944888-40944889, 40945361-40945362], corresponds to element [E0L,E0R]
#      [40945363-40945364, 40982722-40982723], corresponds to element [I0L,I0R]
#      [40982724-40982725, 40982976-40982977], corresponds to element [E1L,E1R]
#      [40982978-40982979, 40988251-40988252], corresponds to element [I1L,I1R]
#      [40988253-40988254, 40988397-40988398], corresponds to element [E2L,E2R]
#      [40988399-40988400, 40990708-40990709], corresponds to element [I2L,I2R]
#          :         :         :        :
#   ]
#
#   this routine returns the list of midpoints for left and right boundary
#   elements for either the exons or the introns --- for example if exons
#   (where $start = 0), it would return:
#
#      [40944888.5, 40945361.5, 40982724.5, 40982976.5, 40988253.5, 40988397.5..
#   corresp.  E0L,      E0R,        E1L,        E1R,        E2L,        E2R...

sub _get_midpoints_boundary_regions_ {
   my ($start, $coord_pairs) = @_;
   my $midpoints = [];

#__THE BEGIN-END COORDINATE MATE PAIRS FOR EACH ELEMENT
   my @begin_end_vals = @{$coord_pairs};

#__RAW LIST-OF-LISTS IS STAGGERED: EXON,INTRON,EXON,INTRON ETC SO SKIP BY TWO
   my ($begin_pos, $end_pos);
   for (my $i = $start; $i <= $#begin_end_vals; $i += 2) {
      foreach my $boundary_string (@{$begin_end_vals[$i]}) {
         if ($boundary_string =~ /(\d+)\-(\d+)/) {
            ($begin_pos, $end_pos) = ($1, $2);
            push @{$midpoints}, ($begin_pos + $end_pos)/2;
         } else {
            croak "do not recognize begin-end pair --- needs diagnosis";
         }
      }
   }
   return $midpoints;
}

# GIVEN A SUPERSET AND SUBSET, BOTH MULTISETS, RETURNS COMPLEMENT OF SUBSET

sub _multiset_complementer_ {
   my ($super_multi_set, $sub_multi_set) = @_;
   my $hash = {};
   foreach my $element (@{$super_multi_set}) {
      $hash->{$element}++;
   }
   foreach my $element (@{$sub_multi_set}) {
      $hash->{$element}--;
   }
   my @complement = ();
   foreach my $element (keys %{$hash}) {
      for (my $i = 1; $i <= $hash->{$element}; $i++) {
         push @complement, $element;
      }
   }
   return @complement;
}

sub _already_have_pval_ {
   my ($catalog_pval_hash) = @_;
   my $already_have_pval = 0;
   foreach my $mode (keys %{$catalog_pval_hash}) {
      my $val = $catalog_pval_hash->{$mode};
      if ($val =~ /^[-+]?[0-9]*\.?[0-9]+$/ ||
          $val =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
             $already_have_pval = 1;
             last;
      }
   }
   return $already_have_pval;
}


################################################################################
##                                                                            ##
##            T R A I L I N G   P O D   D O C U M E N T A T I O N             ##
##                                                                            ##
################################################################################

=head1 NOTES ON ALGORITHMS AND THEIR NUMERICAL IMPLEMENTATIONS AND PROPERTIES

The algorithms in this package are based mostly on straightforward

=cut

################################################################################
##                                                                            ##
##                                    E N D                                   ##
##                                                                            ##
################################################################################

1;

