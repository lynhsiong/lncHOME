#!/usr/bin/perl -w

use Env qw(LD_LIBRARY_PATH);

BEGIN {
 $LD_LIBRARY_PATH = '/bin:/usr/local/gcc-4.9.2/lib64:/opt/openmpi-1.8.7/lib:/opt/intel-xe_2015u3/lib/intel64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:/Share/home/zhangqf2/usr/openmpi/lib:/path/to/bcl-3.1.0-Linux-x86_64:/Share/home/zhangqf2/usr/readline-6.3/lib:/Share/home/zhangqf2/xiongtl/tools/binutils/lib:/Share/home/zhangqf2/usr/lib64:/Share/home/zhangqf2/tmp/expat/lib:/Share/home/zhangqf2/usr/meme/lib:/Share/home/zhangqf2/shaodi/app/samtools/zlib-1.2.8/build/lib:/Share/home/zhangqf2/shaodi/app/samtools/ncurses-6.0/build/lib:/Share/home/zhangqf2/jellyfish/lib:/Share/home/zhangqf2/xiongtl/tools/stamp/gsl_dir/gsl/lib:/Share/home/zhangqf2/xiongtl/tools/binutils/lib:/Share/home/zhangqf2/usr/lib64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:/Share/home/zhangqf2/tmp/expat/lib:/Share/home/zhangqf2/usr/meme/lib:/Share/home/zhangqf2/usr/openmpi/lib:/Share/home/zhangqf2/shaodi/app/samtools/zlib-1.2.8/build/lib:/Share/home/zhangqf2/shaodi/app/samtools/ncurses-6.0/build/lib:/Share/home/zhangqf2/jellyfish/lib:'.$LD_LIBRARY_PATH;
  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
}

require '/Share/home/zhangqf/usr/perl/lib/site_perl/5.22.0/Parallel/ForkManager.pm';
use Lncrnaevo::Myownperlmodule qw(ModuleFinder_L2S_preprocess ModuleFinder_L2S DifferentialDistributionAnalysis uniq_array eucdist intersect_two_redun_array);
use List::Util qw(min max);
use experimental qw(smartmatch);

$step=shift; $wind=shift;
$FIR_cluster_fimo=shift;
$fimo_redundant_motif_3rd=shift;
$fimo_final_filter=shift;

$thl=1000;
($seq2_mark=$fimo_final_filter)=~s/\.\S+//;
$mydir=basedir($fimo_redundant_motif_3rd);
`mkdir -p $mydir/cluster`;

## read cluster fimo (FIR), and 3rd fimo (TRD, mutiple TrX candidates)
open INC, $FIR_cluster_fimo;
while(<INC>){
     next if (/^#/);
     chomp;
     @pLn=split("\t", $_);
     push @{$phash{$pLn[1]}}, $pLn[0]."_".$pLn[4];
     #push @{$phash_FL{$pLn[1]}}, join("\t", @pLn);
     #push @{$phash_pos{$pLn[1]}}, $pLn[2]."_".$pLn[3];
}
close INC;

open IN, $fimo_redundant_motif_3rd;
while(<IN>){
     next if (/^#/);
     chomp;
     @Ln=split("\t", $_);
     push @{$hash{$Ln[1]}}, $Ln[0]."_".$Ln[4];
     push @{$hash_FL{$Ln[1]}}, join("\t", @Ln);
     push @{$hash_pos{$Ln[1]}}, $Ln[2]."_".$Ln[3];
}
close IN;
## 

( $len_hash, $ppool_hash, $pool_hash, $poolfimo_hash )=ModuleFinder_L2S_preprocess( $thl, $step, $wind, \%phash, \%hash, \%hash_pos, \%hash_FL );

%mylength=%{$len_hash}; 
%ppool=%{$ppool_hash};
%pool=%{$pool_hash};
%poolfimo=%{$poolfimo_hash};

$mp = Parallel::ForkManager->new(5);

for $ele ( keys %pool ) {
$mp->start and next;

$prefix_temp=$fimo_final_filter.".vs.".$ele;

if($mylength{$ele}<$thl){

$re= $prefix_temp."."."C0". "_3rd";
$outf=$mydir."/cluster/".$re; $outfs=$mydir."/cluster/".$re.".sort";
if(defined $poolfimo{$ele}{"0"}){
@mytmpfimo=@{$poolfimo{$ele}{"0"}};
@tmp_pool=@{[map {s/^(\S+)\t\S+\t(.+)/$1\t$seq2_mark\t$2/; $_} @mytmpfimo]};
open RE, ">".$outf;
print RE join("\n", @{[uniq_array(\@tmp_pool)]}),"\n";
close RE;
$outf_size= -s "$outf";
if ( $outf_size<1000 ){ 
`rm -rf $outf`; 
}else{
`sort -k2,2 -k3,3n -k4,4n $outf >$outfs; rm -rf $outf`;
`perl /Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts/fimo_motif_coverage_filter.pl $outfs` if ( not $? );
}
}
}else{

$_seq1=${[keys %ppool]}[0];
#print $_seq1,"	haha\n";
($outfinalfile, $_hash_arr_seq)=ModuleFinder_L2S($mydir, $prefix_temp, $_seq1, $ele, $seq2_mark, \%ppool, \%pool, \%poolfimo);
%my_hash_arr_seq=%{$_hash_arr_seq};

#=head

$mydbscanRes=dbscan($outfinalfile, 5);

if ( not $? ){
 $fullcluster=GoupCluster($mydbscanRes) if ( -s $mydbscanRes );
}elsif(( $? >> 8 ) != 0){
 exit;
}

&MotifRetreiving($prefix_temp, $fullcluster, \%my_hash_arr_seq) if (defined $fullcluster);
#=cut

}

$mp->finish;
}

$mp->wait_all_children;



####################
### subfunctions ###
####################

#sub "DBSCAN:Density Based Clustering of Applications with Noise": redundancy_removing_from_"file fimo_final_filter", collecting pair-regions with common motifs, and clustering them
sub dbscan{
   $outfinal=shift; $outfilesize= -s "$outfinal";
   $parameter=shift;
   $fimo_final_filter_dbscan=$outfinal.".dbscan";
   if ( $outfilesize>=1500 ){    ## "250" means at least four lines (~1.5X motif length span, ); REF: 8 lines [4 X 8 = 32 base, plus " one unit average motif ~30", 62 base span], file size 504 bit 
   `python /Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts/exhdbscan.py $outfinal $parameter $fimo_final_filter_dbscan`;
     if ( not $? ){
     `rm -rf $outfinal`;
     return  $fimo_final_filter_dbscan;
     }
   }else{
     `rm -rf $outfinal`;
   }
}

#sub retreive motifs from pair-regions_above without noise
sub GoupCluster{
    $_inputF=shift;
    if ( -s $_inputF ){
    open INF, $_inputF;
    while(<INF>){
      next if (/^#/);
      chomp;
      @NLn=split("\t", $_);
      if ($NLn[2]>0){
         push @{$cluster{$NLn[2]}}, $NLn[0]."_".$NLn[1];
      }
    }
    close INF;
    `rm -rf $_inputF`;
    return \%cluster;
    }else{   
    return 0;
    }
}

sub MotifRetreiving{
    $sub_pre_tem=shift;
    $centermotifs=shift;
    #$seq1=shift;  #HUM/MOU
    #$seq2=shift;  #3rd
    $t_hash_arr_seq=shift; %hash_arr_seq=%{$t_hash_arr_seq}; 
    #$t_length=shift; %length=%{$t_length};
    foreach $_key (keys %{$centermotifs}){
       for (@{$$centermotifs{$_key}}){
           push @{$cluster_hash_SEC_common_motifs{$_key}}, @{$hash_arr_seq{$_}};
           push @{$cluster_hash_SEC_ord{$_key}}, $_;
       }       
    }

    foreach $x (sort {$a<=>$b} keys %cluster_hash_SEC_ord){

           $re= $sub_pre_tem."."."C".$x. "_3rd";

           `mkdir -p $mydir/cluster`; $outf=$mydir."/cluster/".$re; $outfs=$mydir."/cluster/".$re.".sort";
           open $re, ">".$outf;
           print $re join("\n", @{[uniq_array(\@{$cluster_hash_SEC_common_motifs{$x}})]}),"\n";
           close $re;

           $outf_size= -s "$outf";
           if ( $outf_size<1000 ){ 
               `rm -rf $outf`; 
           } 

           if( -s $outf ){
           `sort -k2,2 -k3,3n -k4,4n $outf >$outfs; rm -rf $outf`;
           `perl /Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts/fimo_motif_coverage_filter.pl $outfs` if ( not $? );
           }           
           
    }
}

sub basedir{
    $basedir=shift;
    if($basedir=~/^\S+\//){
    $basedir=~s/(^\S+)\/\S+/$1/g;
    }else{
    $basedir='./';
    }
    return $basedir;
}
