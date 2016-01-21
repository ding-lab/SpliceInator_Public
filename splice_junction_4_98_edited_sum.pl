#!/gsc/bin/perl
use strict;
my $start = time();

my $usage =<<USAGE;
 Usage: $0 <Original MAF> <Directory of bed files>
        The output of the script will be ...
  Example: perl scripts/splice_junction_3.pl laml/laml_ss_mutation /LAML/ 
USAGE
    die $usage unless @ARGV==2;
open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
my $BEDDIR=$ARGV[1];
#bed file:
#chr	start	stop	JUNC ID	score	strand	thk_start	thk stop	rgb	block_size	overhang	?
#10	8105773	8106005	JUNC00158136	6	+	8105773	8106005	255,0,0	2	41,47	0,185
#10	8106052	8111484	JUNC00158137	3151	+	8106052	8111484	255,0,0	2	49,49	0,5383
#10	8106054	8111483	JUNC00158138	15	+	8106054	8111483	255,0,0	2	47,41	0,5388
#10	8110769	8111476	JUNC00158139	3	+	8110769	8111476	255,0,0	2	29,41	0,666
#10	8111512	8115750	JUNC00158140	3329	+	8111512	8115750	255,0,0	2	49,49	0,4189
#10	8111514	8115755	JUNC00158141	186	+	8111514	8115755	255,0,0	2	47,46	0,4195
#10	8114225	8115723	JUNC00158142	1	+	8114225	8115723	255,0,0	2	28,22	0,1476
#Load in MAF and for each splice site store 
my (@ss,$counter,@controls,%control,%sshash,$id,$cancer,$sample,$gene,$file,$linest);
$counter=0;
#########GO THROUGH MAF AND STORE EACH LINE AS VALUE IN HASH AND USE SAMPLE AND ID(GENE_CHR_START_REF_VAR) AS KEYS 
while(<$MAF>){
        next if /Hugo/;
        chomp ($linest=$_);
        @ss=split(/\t/,$linest);
		$sample=$1 if $ss[15]=~/(^TCGA-\w{2}-\w{4})/;
#Define id for MAF:GENE_17_29586049_G_A
        $id=$ss[0]."_".$ss[4]."_".$ss[5]."_".$ss[10]."_".$ss[12];
#keys SAMPLE NAME and ID
	$sshash{$sample}{$id}=$linest;
}
close $MAF;

my (%junctionthree,%junctionfive,%intron);

#FOREACH SAMPLE KEY IN THE HASH, FIND &  OPEN THE JUNCTION.BED FILE. 
#FOREACH ID IN THE HASH ....
#IF THE CHR OF THE FEATURE MATCHES THE MAF ENTRY, AND IF THE SPLICE SITE FALLS WITHIN THE START AND STOP OF THE BED FEATURE, THEN CALCULATE THE 5' AND 3' EXON INTRON JUNCTIONS USING THE OVERHANGS. 
#STORE THE 5' AND 3' CALCULATED EXON INTRON JUNCTIONS INTO A NEW HASH (USE 5' AND 3' JUNCTIONS AS KEYS) ALONG WITH ANOTHER KEY SAMPLE_ID_JUNCTID AND STORE THE BED LINE AS VALUE.

foreach my $sn (sort keys %sshash){
	$a=`ls $BEDDIR|grep $sn|head -n1`;
#if bed file exists and is a tumor bam
#TCGA-BH-A0BS-01A-11R-A12P-07.junction.bed
#TCGA-BH-A0BS-11A-11R-A12P-07.junction.bed
	if ($a=~/.bed/){
#full directory of junctions.bed file for sample
	my $file =$BEDDIR.$a;
	chomp $file;
		foreach my $id1 (keys %{$sshash{$sn}}){	
		open (my $BED,'<',$file) or die "bed file won't open!!!\n";
#find start position of splice site
		my $mafentry=$sshash{$sn}{$id1};
		my @ss1=split(/\t/,$mafentry); 
		my @con=split(/,/,$ss1[98]);
		my $splicesite=$ss1[5];
			while(my $line=<$BED>){
			next if $line=~/^track/;
			chomp $line;			
			my @bedf=split(/\t/,$line);
			$bedf[0]=~s/chr//;
			#print STDERR "$bedf[0]\t$ss1[4]\n";
#if chromosome in bed file matches chromosome in id of MAF
				if ($bedf[0] eq $ss1[4]){
#if the splice site falls within start and stop of splice junciton feature
#if splicesite>start && splicesite<stop
					if($splicesite>=$bedf[1] && $splicesite<=$bedf[2]){
					
#calculate exon intron junction using overhangs and start and stop of splice junction reads taking into account +/- strand
						chomp $bedf[5];
						my @overhang=split(/,/,$bedf[10]);
						my ($fivejxn,$threejxn);
						if($bedf[5] eq '+'){
							$fivejxn= $bedf[1]+$overhang[0]+1;
							$threejxn= $bedf[2]-$overhang[1];
						}
						if($bedf[5] eq '-'){
                                               		$fivejxn= $bedf[2]-$overhang[1];
                                              	 	$threejxn= $bedf[1]+$overhang[0]+1;
						}
#					print "$bedf[1]\t$overhang[0]\t$fivejxn\t$threejxn\t$overhang[1]\t$bedf[2]\n";	
#key = TCGA-E9-A249_GATA3_10_8111433_CA_-_JUNC00141929_100
						my $key=$sn."_".$id1."_".$bedf[3]."_".$bedf[4];
#store splice junction data based on 5' and 3' estimated exon/intron junctions
						$intron{$key}=$fivejxn."-".$threejxn."(".$bedf[4].")";
						my $controljxn=$bedf[0]."_".$fivejxn."_".$threejxn;
#For specific junction add start_stop as a key and use sample from controls as a second key that points to nothing for now but will point to coverage
#controljxn = chr_start_stop of jxn		
				       		foreach (@con){
				                	my $short=$1 if $_=~/(^TCGA-\w{2}-\w{4})/;
            						$control{$controljxn}{$short}=undef;
        					}  
						$junctionthree{$threejxn}{$key}=$line;
						$junctionfive{$fivejxn}{$key}=$line;					
					}

				}
			}
			close $BED;
		}
	}
}
############THIS IS THE SLOW PART
####GO THROUGH EACH JUNCTION START_STOP THAT IS STORED AS A KEY IN THE CONTROL HASH AND OPEN UP EACH CONTROL SAMPLE STORED AS A KEY AND ADD THE COVERAGE IF THE JUNCTION SITE EXISTS
my %con;
foreach my $jxn (sort keys %control){
my @ssinfo=split(/_/,$jxn);
my $chr=$ssinfo[0];
my $five=$ssinfo[1];
my $three=$ssinfo[2];
####GO THROUGH EACH SAMPLE KEY AND OPEN BED
	foreach my $sampl (keys %{$control{$jxn}}){
#	print "$sampl\t$jxn\t";
	if (defined $control{$jxn}{$sampl}){
#	print "EXISTS-skipping\n";
	next;
	}
#	print "Doesn't exist\n";
	$a=`ls $BEDDIR|grep $sampl|head -n1`;	
		if ($a=~/.bed/){
		my $file =$BEDDIR.$a;
		chomp $file;
	#	print STDERR "$file\n";
		open (my $BED,'<',$file) or die "bed file won't open!!!\n";
###GO THROUGH EACH LINE OF BED FILE TO FIND CHR START THAT MATCHES THE JUNCTION START 
			while(my $line=<$BED>){
			next if $line=~/^track/;		
			chomp $line;
			my @bedf=split(/\t/,$line);
			chomp $bedf[0];
			$bedf[0]=~s/chr//;
			chomp $chr;
				if ($bedf[0] eq $chr){
                                chomp $bedf[5];
###CALCULATE JUNCTION OVERHANG
                                        my @overhang=split(/,/,$bedf[10]);
                                        my ($fivejxn2,$threejxn2);
                                        if($bedf[5] eq '+'){
                                                $fivejxn2= $bedf[1]+$overhang[0]+1;
                                                $threejxn2= $bedf[2]-$overhang[1];
                                        }   
#IF THE JUNCTION IS MAPPED TO THE NEGATIVE STRAND THEN SWITCH THE 5' AND 3' JUNCTION BUT CALCULATE THE SAME WAY AS ON THE POSITIVE STRAND
                                        if($bedf[5] eq '-'){
                                               $fivejxn2= $bedf[2]-$overhang[1];
                                               $threejxn2= $bedf[1]+$overhang[0]+1;
                                        }
##STORE THE COVERAGE AS A VALUE IN CONTROL HASH WITH KEY AS CHR_START_STOP OF JUNCTION and SAMPLE and STRAND
					my $key1=$bedf[0]."_".$fivejxn2."_".$threejxn2;	
					if (exists $control{$key1}{$sampl}){
					#	if ($bedf[4]>0){
							$control{$key1}{$sampl}=$bedf[4];
					#	}
					}
				}
			}
			close $BED;
		}
	}
}

foreach my $jxn (keys %control){
       foreach my $sampl (keys %{$control{$jxn}}){
		my $val=$control{$jxn}{$sampl};
		if ($val==""){
		$control{$jxn}{$sampl}=0;
		}
              print STDERR "SAMP:$sampl\tJXN:$jxn\tCOVG:$control{$jxn}{$sampl}\n";
       }
}

#LOAD IN BED FILE TO IDENTIFY CANONICAL SITES AND ALTERNATIVE SPLICE SITES
#Example:
#chr10	8096655	8096853	GATA3:ENST00000379328:+:e1
#chr10	8096854	8097248	GATA3:ENST00000379328:+:i1
#chr10	8097249	8097858	GATA3:ENST00000379328:+:e2
#chr10	8097859	8100266	GATA3:ENST00000379328:+:i2
#chr10	8100267	8100803	GATA3:ENST00000379328:+:e3
#chr10	8100804	8105954	GATA3:ENST00000379328:+:i3
#chr10	8105955	8106100	GATA3:ENST00000379328:+:e4
#chr10	8106101	8111434	GATA3:ENST00000379328:+:i4
#chr10	8111435	8111560	GATA3:ENST00000379328:+:e5
#chr10	8111561	8115700	GATA3:ENST00000379328:+:i5
#chr10	8115701	8117160	GATA3:ENST00000379328:+:e6

open(my $BED,'<','BED');
my %bedhash;
while(my $bedline=<$BED>){
#store in hash by gene transcript start etc.
	my @B=split(/\t/,$bedline);
	my $stop=$B[3];
	my @gene=split(/:/,$B[3]);
#gene name as key and start and stop as keys	
	my $start = $B[1]+1;
	my $stop = $B[2]+1;
	my $key;
	if ($gene[2] eq '-'){
		#print "$gene[2]\n";
		$key = $gene[0]."_".$stop."_".$start;	
	
	}
	else {
		$key = $gene[0]."_".$start."_".$stop;
	}
	$bedhash{$key}=$bedline;
#	print "$key\n";
}
#GO THROUGH EACH ENTRY IN MAF AND  EXTRACT SAMPLE NAME AND INFORMATION ABOUT MUTATION
#FOR EACH KEY IN THE JUNCTION (INTRON) HASH EXTRACT THE SAMPLE NAME, JUNCTION INFORMATION AND COVERAGE
#IF THE COVERAGE IS GREATER THAN 50, THEN PRINT OUT THE 5'_3' SPLICE JUNCTION AND THE COVERAGE OF THE SAMPLE
#FURTHER SUM AND AVERAGE THE COVERAGE OF THE CONTROLS AND PRINT OUT
open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
while(my $linest=<$MAF>){
        my @ss2=split(/\t/,$linest);

#	my $gene = $ss2[0];
        my $sam=$1 if $ss2[15]=~/(^TCGA-\w{2}-\w{4})/;
        my $mafkey=$sam."_".$ss2[0]."_".$ss2[4]."_".$ss2[5]."_".$ss2[10]."_".$ss2[12];
	my @con=split(/,/,$ss2[98]);
#key = TCGA-E9-A249_GATA3_10_8111433_CA_-_JUNC00141929_100	
	print "$mafkey\t";
	foreach my $key (sort keys %intron){
		 my $samp=$1 if $key=~/(^TCGA-\w{2}-\w{4})_.*_.*_.*_.*_.*_JUNC.*_\d*$/;
		 my $intronkeymatch=$1 if $key=~/(^TCGA-\w{2}-\w{4}_.*_.*_.*_.*_.*)_JUNC.*_\d*$/;
		 my $gene=$1 if $key=~/^TCGA-\w{2}-\w{4}_(.*)_.*_.*_.*_.*_JUNC.*_\d*$/;
                 my $junc=$1 if $key=~/^TCGA-\w{2}-\w{4}_.*_.*_.*_.*_.*_(JUNC.*)_\d*$/;
                 my $covg=$1 if $key=~/^TCGA-\w{2}-\w{4}_.*_.*_.*_.*_.*_JUNC.*_(\d*$)/;
		 if ($mafkey eq $intronkeymatch){
			if ($covg>=2){
				my $J=$intron{$key};
				my $coverage=$1 if $J=~/.*-.*\((.*)\)/;
				my $jun=$1 if $J=~/(.*-.*)\(.*\)/;
				my $jun2=$jun;
				$jun2=~s/-/_/;
				my $canonical=$gene."_".$jun2;
				if (exists $bedhash{$canonical}){
					print "CANON:[$jun],MUT:$coverage(REF)";
				}
				else{
					print "ALT:[$jun],MUT:$coverage";
					my @altjun=split(/-/,$jun);
					#print "ALT1:$altjun[0]\tALT2:$altjun[1]\n";
					my @array2 = grep {$_=~/$gene\_/} keys %bedhash;	
		                     if (@array2){
#chrX	44896934	44910952	KDM6A:ENST00000377967:+:i8
#chrX	44918711	44919265	KDM6A:ENST00000377967:+:i12
#chrX	44923062	44928822	KDM6A:ENST00000377967:+:i16
#STORE ALL INFORMATION FROM BED FILE INTO HASH BASED ON INTRON # AND GENE NAME
					my (%intronbedhash);
					for my $b (0..$#array2){
                                                my $bedline2=$bedhash{$array2[$b]};
                                                my @bedarray2 = split(/\s/,$bedline2);
                                                my @introninfo2=split(/:/,$bedarray2[3]);
                                                my $keyintronbed=$gene."_".$introninfo2[3];
                                                $intronbedhash{$keyintronbed}=$bedline2;
					}	
#GO THROUGH EACH ELEMENT OF ARRAY AND DETERMINE IF THE JUNCTION IS PREDICTED TO BE AN EXON SKIPPING/SHRINKAGE EVENT
					for my $a (0..$#array2){
						my $bedline=$bedhash{$array2[$a]};
						my @bedarray = split(/\s/,$bedline);
						my @introninfo=split(/:/,$bedarray[3]);
						#POSITIVE STRAND
						if ($altjun[0]<$altjun[1]){
							#3' EXON EXTENSION
							if ($bedarray[2]>$altjun[0] && $altjun[0]>$bedarray[1]){
								if($altjun[1]==$bedarray[2]+1){	
                                                			print "(3' EXON EXTENSION +)";
								}
							}	
							#5' EXON EXTENSION
							if ($bedarray[2]>$altjun[1] && $altjun[1]>$bedarray[1]){
								if($altjun[0]==$bedarray[1]+1){
									print "(5' EXON EXTENSION +)";
								}
							}
							###3' EXON SHRINKAGE AND SKIPPING
							if ($altjun[0]<$bedarray[1]){ 
								if ($bedarray[2]+1==$altjun[1]){
									my $currentintron=$introninfo[3];
									$currentintron=~s/i//g;
									my $previousintron=$currentintron-1;
									my $keyprevious=$gene."_i".$previousintron;
									my $keycurrent=$gene."_i".$currentintron;
									my @current=split(/\s/,$intronbedhash{$keycurrent});
									 if (exists($intronbedhash{$keyprevious})){
										my @previous=split(/\s/,$intronbedhash{$keyprevious});
										my $previousstart=$previous[1];
										my $previousend=$previous[2];
										if ($altjun[0]>$previousend){
											print "(3' EXON SHRINKAGE +)";
										}
							#		if ($altjun[0]<$previousend){
							#			print "(+ EXON SKIPPING)";
							#		}
									}		
								}
							}
						 	###5' EXON SHRINKAGE AND SKIPPING
                                                        if ($altjun[1]>$bedarray[2]){ 
                                                                if ($bedarray[1]+1==$altjun[0]){
                                                                        my $currentintron=$introninfo[3];
                                                                        $currentintron=~s/i//g;
                                                                        my $nextintron=$currentintron+1;
                                                                        my $keynext=$gene."_i".$nextintron;
                                                                        my $keycurrent=$gene."_i".$currentintron;
                                                                        my @current=split(/\s/,$intronbedhash{$keycurrent});
									if (exists($intronbedhash{$keynext})){
                                                                        	my @next=split(/\s/,$intronbedhash{$keynext});
										my $nextstart=$next[1];
                                                                        	my $nextend=$next[2];
									
                                                                        	if ($altjun[1]<$nextstart){
                                                                                	print "(5' EXON SHRINKAGE +)";
                                                                        	}   
                                                                        	if ($altjun[1]>$nextstart){
                                                                               	 	print "(EXON SKIPPING +)";
                                                                        	}		
									}   
                                                                }   
                                                        }
					#		else{print "(UNKNOWN)";}
						}
						##NEGATIVE STRAND
						if ($altjun[0]>$altjun[1]){
						##5' EXON EXTENSION
							if ($bedarray[2]>$altjun[1] && $altjun[1]>$bedarray[1]){
                                                                if($altjun[0]==$bedarray[2]+1){
                                                                        print "(5' EXON EXTENSION -)";
                                                                }   
                                                        }	
						##3' EXON EXTENSION
						      if ($bedarray[2]>$altjun[0] && $altjun[0]>$bedarray[1]){
                                                                if($altjun[1]==$bedarray[1]+1){
                                                                        print "(3' EXON EXTENSION -)";
                                                                }   
                                                        }
							###5' EXON SHRINKAGE AND SKIPPING
                                                        if ($altjun[1]<$bedarray[1]){ 
                                                                if ($bedarray[2]+1==$altjun[0]){
                                                                        my $currentintron=$introninfo[3];
                                                                        $currentintron=~s/i//g;
                                                                        my $nextintron=$currentintron+1;
                                                                        my $keynext=$gene."_i".$nextintron;
                                                                        my $keycurrent=$gene."_i".$currentintron;
                                                                        my @current=split(/\s/,$intronbedhash{$keycurrent});
                                                                     	if (exists($intronbedhash{$keynext})){
                                                                        	my @next=split(/\s/,$intronbedhash{$keynext});
										my $nextstart=$next[1];
                                                                        	my $nextend=$next[2];
                                                                        	if ($altjun[1]-1>$nextstart){
										#	print "$altjun[1]>";
                                                                                	print "(5' EXON SHRINKAGE -)";
                                                                        	}   
                                                                     	   	if ($altjun[1]<$nextstart){
                                                                        	        print "(EXON SKIPPING -)";
                                                                       		}   

										if (($altjun[1]-1)==$nextstart){
                                                                                	print "(EXON SKIPPING -)";
                                                                        	}	
									}
                                                                }   
                                                        } 	
							###3' EXON SHRINKAGE AND SKIPPING
                                                        if ($altjun[0]>$bedarray[2]){ 
                                                                if ($bedarray[1]+1==$altjun[1]){
                                                                        my $currentintron=$introninfo[3];
                                                                        $currentintron=~s/i//g;
                                                                        my $previousintron=$currentintron-1;
                                                                        my $keyprevious=$gene."_i".$previousintron;
                                                                        my $keycurrent=$gene."_i".$currentintron;
                                                                        my @current=split(/\s/,$intronbedhash{$keycurrent});
                                                                        if (exists($intronbedhash{$keyprevious})){
                                                                        	my @previous=split(/\s/,$intronbedhash{$keyprevious});
										my $previousstart=$previous[1];
                                                                       	 	my $previousend=$previous[2];
                                                                        	if ($altjun[0]<$previousend){
                                                                                	print "(3' EXON SHRINKAGE -)";
                                                                        	}
                                                         #               if ($altjun[0]>$previousend){
                                                          #                      print "(- EXON SKIPPING)";
                                                           #             }
									}
                                                                }
                                                        }
					#		else{print "(UNKNOWN)";}
						}
                                	}        
				}   
		}
				my $junction=$intron{$key};
				$junction =~ tr/-/_/;
		 		my $J=$1 if $junction=~/(.*_.*)\(.*\)/;
				my $key1=$ss2[4]."_".$J;
				my $numcontrol=0;
				my $controladd=0;	
	
########################TO CHANGE FROM AVG TO SUM

			foreach (@con){
					my $S=$1 if $_=~/(^TCGA-\w{2}-\w{4})/;
					if (defined $control{$key1}{$S}){
						my $COVG=$control{$key1}{$S};
						$controladd=$controladd+$COVG;
						$numcontrol++;	
					}
				}
				my $COV;
				if ($controladd==0){
					$COV=0;
				}
				if ($controladd>0){
				
#ORGIGINAL:				$COV=$controladd/$numcontrol;
###EDITED
					$COV=$controladd;
				}
				printf ",CON:%.2f($numcontrol)\t",$COV;
			}
                 }
	}
				print "\n";
}
my $end = time();
