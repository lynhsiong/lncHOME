#!/usr/bin/perl -w

use Env qw(LD_LIBRARY_PATH);

BEGIN {
  $LD_LIBRARY_PATH = '/bin:/usr/local/gcc-4.9.2/lib64:/opt/openmpi-1.8.7/lib:/opt/intel-xe_2015u3/lib/intel64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:$HOMEZ/usr/openmpi/lib:/path/to/bcl-3.1.0-Linux-x86_64:$HOMEZ/usr/readline-6.3/lib:$HOME/tools/binutils/lib:$HOMEZ/usr/lib64:$HOMEZ/tmp/expat/lib:$HOMEZ/usr/meme/lib:$HOMEZ/shaodi/app/samtools/zlib-1.2.8/build/lib:$HOMEZ/shaodi/app/samtools/ncurses-6.0/build/lib:$HOMEZ/jellyfish/lib:$HOME/tools/stamp/gsl_dir/gsl/lib:$HOME/tools/binutils/lib:$HOME/usr/lib64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:$HOME/tmp/expat/lib:$HOME/tools/meme/lib:$HOMEZ/usr/openmpi/lib:$HOMEZ/shaodi/app/samtools/zlib-1.2.8/build/lib:$HOMEZ/shaodi/app/samtools/ncurses-6.0/build/lib:$HOMEZ/jellyfish/lib:$LD_LIBRARY_PATH';
  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
}

use Lncrnaevo::Myownperlmodule qw(DifferentialDistributionAnalysis uniq_array eucdist intersect_two_redun_array);
use List::Util qw(min max);
use experimental qw(smartmatch);

$step=shift; $wind=shift;
$fimo_redundant_motif=shift;     ### full directory
$fimo_final_filter=shift;

$j=$step;   # 2   length-stepwise: (start overlapped/) Motif Number added per step
$k=$wind;   # 120 window-size: Motif Number within a length-fixed Window (e.g. about 120bp, but end coordinates match any motif ends)

open IN, $fimo_redundant_motif;
# 0224    FIR 1   10  +   10.7241 9.93e-05        GGGCCCAGGG
while(<IN>){
     next if (/^#/);
     chomp;
     @Ln=split("\t", $_);
     push @{$hash{$Ln[1]}}, $Ln[0]."_".$Ln[4];
     push @{$hash_FL{$Ln[1]}}, join("\t", @Ln);
     push @{$hash_pos{$Ln[1]}}, $Ln[2]."_".$Ln[3];
}
close IN;

=head
foreach $nkey (sort keys %hash_pos){
  for ($i=0; $i<=$length{$nkey}; $i+=$j){
      for (0..$motifnum{$nkey}){
         last if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] > ($i+$k));
         if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] >= $i){
           push @{$pool{$nkey}{$i}}, ${$hash{$nkey}}[$_];
           push @{$pool_fimo{$nkey}{$i}{${$hash{$nkey}}[$_]}}, ${$hash_FL{$nkey}}[$_];
         }
      }
  }
}
=cut

foreach $nkey (sort keys %hash_pos){ 
	%temphash=(); @max=(); %mark=(); %m2shash=();
	%tempMotifStart=%tempMotifEnd=(); $exitsignal=0;

	@temp=@{$hash_pos{$nkey}};
	for (@temp){ push @max, ${[split("_", $_)]}[0]; } 
	$mymax=max(@max);
	for (0..$#temp){
		$temphash{${[split("_", $temp[$_])]}[0]}++;
		$tempMotifStart{$_}=${[split("_", $temp[$_])]}[0];
		$tempMotifEnd{$_}=${[split("_", $temp[$_])]}[1];
	}
	@tmpkeysort= (sort {$a<=>$b} keys %temphash);
#print join("-", @tmpkeysort),"\n";
	for (0..$#tmpkeysort){
		$mark{$_}= $temphash{$tmpkeysort[$_]};
	}
	$cummulateS=$cummulateE=0;
	foreach $ckey (sort {$a<=>$b} keys %mark){
		$cNum=$mark{$ckey}; 
		$cummulateS=($cummulateE==0)?$cummulateE:($cummulateE+1); $cummulateE=$cummulateS+$cNum-1;
		for ($cummulateS..$cummulateE){
			$m2shash{$_}=$ckey;
		}
	}		

	$motifnum=scalar(@temp)-1;

	for ($i=0; $i<=$motifnum; $i+=step($i, $j, \%mark, \%m2shash)){
#print $nkey,"... ...",$i,"... ...",step($i, $j, \%mark, \%m2shash),"\n";
		$tmpstt=${[split("_", $temp[$i])]}[0];
		$finalstt=$tmpstt+$k;
		if ($finalstt >$mymax){
			$exitsignal=1;
			$finalstt = $mymax;
		}
#print $mymax,".cmp.", $finalstt,"\n";
		for ($i..$motifnum){
		  if ($finalstt>=${[split("_", $temp[$_])]}[0]){
			#push @{$pool{$nkey}{$i}}, ${$hash{$nkey}}[$_];
			#push @{$pool_fimo{$nkey}{$i}{${$hash{$nkey}}[$_]}}, ${$hash_FL{$nkey}}[$_];
			push @{$pool{$nkey}{$tempMotifStart{$i}}}, ${$hash{$nkey}}[$_];
			push @{$pool_fimo{$nkey}{$tempMotifStart{$i}}{${$hash{$nkey}}[$_]}}, ${$hash_FL{$nkey}}[$_];
		  }
		}
		if ((defined $exitsignal)&&($exitsignal==1)){
			last;
		}
	}  

}


$mydir=basedir($fimo_redundant_motif);
open $FINALOUT, ">".$mydir."/".$fimo_final_filter;
#open $FINALOUT, ">".$fimo_final_filter;

%hash_arr_seq1=%hash_arr_seq2=();
($seq1, $seq2)=sort keys %pool;
%SEQ1=%{$pool{$seq1}}; %SEQ2=%{$pool{$seq2}};     $nnn=0;
#$seq1len=scalar(keys %SEQ1); $seq2len=scalar(keys %SEQ2); 

#print "###$seq1len\t$seq2len###\n";

for $k1 (sort {$a <=> $b} keys %SEQ1){
  $nnn++;
  @temp1=@{$SEQ1{$k1}}; $len_temp1=scalar @temp1;      $nnnid="S$k1"."S$nnn"; 
  &FileIO($mydir, $nnnid); $tempfilename=filename($tempfile);
  for $k2 (sort {$a <=> $b} keys %SEQ2){
    @temp2=@{$SEQ2{$k2}}; $len_temp2=scalar @temp2;
    ($intersect_arr, $intersect_num)=intersect_two_redun_array(\@temp1, \@temp2);
      @temp_commotif_arr=@{$intersect_arr};
      for (@temp_commotif_arr){
          push @{$hash_arr_seq1{$k1."_".$k2}}, @{$pool_fimo{$seq1}{$k1}{$_}};
          push @{$hash_arr_seq2{$k1."_".$k2}}, @{$pool_fimo{$seq2}{$k2}{$_}};
      }
      print $OUTH "$k1\t$len_temp1\t$k2\t$len_temp2\t", $intersect_num, "\n";
  }

  &FileCL($OUTH);
  `/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/ScrDB/scripts/DifferentialDistributionAnalysis_3rd.r $mydir $tempfilename $fimo_final_filter`;  

  `rm -rf $tempfile` if ( not $? );
#   next if ( not $? );
}
close $FINALOUT;

&dbscan();   ##clustering regions from $fimo_final_filter





$fullcluster=GoupCluster() if ( not $? );

&MotifRetreiving($fullcluster);



####################
### subfunctions ###
####################

#sub "DBSCAN:Density Based Clustering of Applications with Noise": redundancy_removing_from_"file fimo_final_filter", collecting pair-regions with common motifs, and clustering them
sub dbscan{
   $fimo_final_filter_dbscan=$fimo_final_filter.".dbscan";
   `/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/ScrDB/scripts/dbscan.r $mydir $fimo_final_filter $fimo_final_filter_dbscan`;
}

#sub retreive motifs from pair-regions_above without noise
sub GoupCluster{
    open INF, "$mydir/$fimo_final_filter_dbscan";
    while(<INF>){
      next if (/^#/);
      chomp;
      @NLn=split("\t", $_);
      if ($NLn[2]>0){
         push @{$cluster{$NLn[2]}}, $NLn[0]."_".$NLn[1];
      }
    }
    close INF;   
    return \%cluster;
}

#sub GrepCenterArray{
#    $tempClusterArray=shift; $threshold=shift;
#    @centerpart=();
#    $_lower=$threshold*(scalar @{$tempClusterArray});
#    $_higher=(1-$threshold)*(scalar @{$tempClusterArray});
#    for (0..$#{$tempClusterArray}){
#        if (($_>=$_lower) && ($_<=$_higher)){
#           push @centerpart, ${$tempClusterArray}[$_];
#        }
#    }
#    return \@centerpart;
#}

sub MotifRetreiving{
    $centermotifs=shift;
    foreach $_key (keys %{$centermotifs}){
       for (@{$$centermotifs{$_key}}){
           push @{$cluster_hash_FIR_common_motifs{$_key}}, @{$hash_arr_seq1{$_}};
           push @{$cluster_hash_SEC_common_motifs{$_key}}, @{$hash_arr_seq2{$_}};
           push @{$cluster_hash_FIR_ord{$_key}}, ${[split("_", $_)]}[0]; ###s=min;e=max+120
           push @{$cluster_hash_SEC_ord{$_key}}, ${[split("_", $_)]}[1];
       }       
    }

    foreach $x (sort {$a<=>$b} keys %cluster_hash_FIR_ord){

           $re= "C".$x. "_FIR_". min(@{$cluster_hash_FIR_ord{$x}}).
                            "_". min( max(@{$cluster_hash_FIR_ord{$x}})+120, $length{$seq1}+30 ).  # +120
                        "_SEC_". min(@{$cluster_hash_SEC_ord{$x}}).
                            "_". min( max(@{$cluster_hash_SEC_ord{$x}})+120, $length{$seq2}+30 );  # +120

           print "-->", 
                 $re, "\n";
           print "->",
                 join("\n->", @{[uniq_array(\@{$cluster_hash_FIR_common_motifs{$x}})]}),"\n";
           print "->",
                 join("\n->", @{[uniq_array(\@{$cluster_hash_SEC_common_motifs{$x}})]}),"\n";

           `mkdir -p $mydir/cluster`; $outf=$mydir."/cluster/".$re; $outfs=$mydir."/cluster/".$re.".sort";
           open $re, ">".$outf;
           print $re join("\n", @{[uniq_array(\@{$cluster_hash_FIR_common_motifs{$x}})]}),"\n",
                     join("\n", @{[uniq_array(\@{$cluster_hash_SEC_common_motifs{$x}})]}),"\n";
           close $re;

           `sort -k2,2 -k3,3n -k4,4n $outf >$outfs; rm -rf $outf`;
           
    }
}

sub basedir{
    $basedir=shift;
    $basedir=~s/(^\/\S+)\/\S+/$1/g;
    return $basedir;
}

sub filename{
    $filename=shift;
    $filename=~s/^\/\S+\/(\S+)/$1/g;
    return $filename;
}

sub FileIO{
    $dir=shift; $nnnmmm=shift; 
    $tempfile="${dir}/$nnnmmm.tempxtl.xxx";
    #$OUTH="tempxtl"; 
    $OUTH="$nnnmmm";
    open $OUTH, ">".$tempfile;
}

#sub FileOP{
#    $FH=shift; $MyTempFile=shift;
#    open $FH, ">".$MyTempFile;
#}

sub FileCL{
    $handle=shift;
    close $handle;
}

sub step{
    $inputpoint=shift; $mstep=shift;
    $mNumHash=shift; %mNumHash=%{$mNumHash};
    $mMarHash=shift; %mMarHash=%{$mMarHash};
    $currentpoint=$mMarHash{$inputpoint};
    $nextpoint=$currentpoint+$mstep-1;
 
    $myreturn=0;
    for($currentpoint..$nextpoint){
        $myreturn+=$mNumHash{$_};
    }
  
    return $myreturn;
   
}
