#!/usr/bin/perl -w

###################################################
#  perl fa2svg.pl all_lnc_sp.lst config.tbl out.svg  
###################################################

$inFAlst = shift;     ### include a list of lncRNA sequences and their marks (tab separated) from all species (if available) to be compared
$config = shift;      ### including settings of directories and run mode

($scrdir, $motifDB, $syngeneDB, $lncPoolDB, $basedir, $fastaDB, $svgdir, $cand2NONCODE, $run_mode)=config($config);     
$cand2NONCODE=`find $basedir -maxdepth 1 -type f -name "$cand2NONCODE"`; chomp $cand2NONCODE;

#=head
print "script_dir=>\t ", $scrdir,"\n";
print "motif_db=>\t ",$motifDB, "\n";
print "macroH_dir=>\t ",$syngeneDB, "\n";
print "allLnc_lib=>\t ",$lncPoolDB,"\n"; 
print "cand_dir=>\t ",$basedir,"\n"; 
print "cand_fas=>\t ",$fastaDB,"\n"; 
print "SVG_dir=>\t ",$svgdir,"\n"; 
print "cand2NONCODE=>\t ",$cand2NONCODE,"\n"; 
print "run_mode=>\t ",$run_mode,"\n";
#=cut

BEGIN {
  #unshift @INC,"$scrdir";
  unshift @INC, "/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/ScrDB/scripts";
}
use SVG;
use List::Util qw(min max);
use Lncrnaevo::Myownperlmodule qw(uniq_array);
#require '/Share/home/zhangqf/usr/perl/lib/site_perl/5.22.0/Parallel/ForkManager.pm';

#$procs=1;
@code_arr = ("FIR","SEC","TRD","FOR","FIF","SIX","SEV","EIG","NIN","TEN");
@sp_arr = ("hum", "mou", "smk", "tsr", "cow", "dog", "ops", "pla", "chc", "zbf");
@fsp_arr = ("humNon", "mouNon", "saiBol", "tarSyr", "bosTau", "canFam", "momDom", "ornAna", "galGal", "danRer");

for (0..$#code_arr){ $A{$code_arr[$_]}=$sp_arr[$_]; }
for (0..$#code_arr){ $B{$code_arr[$_]}=$fsp_arr[$_]; }

for (1..$#sp_arr){
    $mysp=$sp_arr[$_]; $FH=$mysp."FH";
    $synF=`find $syngeneDB -maxdepth 1 -type f -name "hum2$mysp.PostSetOperation.out"`; chomp $synF;
    open $FH, "<$synF";
    while (<$FH>){
        chomp;
        @synFilds=split("\t", $_);
        $joint_key=$mysp."_".$synFilds[0];
        $rev_joint_key=$mysp."_".$synFilds[3];
        push @{$SYN{$joint_key}}, $synFilds[3];          #### mou_NONHSAG*** (NONMMUG*** NONMMUG*** NONMMUG*** ...)
        push @{$RVSYN{$rev_joint_key}}, $synFilds[0];        #### mou_NONMMUG*** (NONHSAG*** NONHSAG*** NONHSAG*** ...)
    }
    close $FH;
}

for (0..$#fsp_arr){
    $myfsp=$fsp_arr[$_]; $FFH=$myfsp."FH"; $ssp=$sp_arr[$_];
    $fsynF=`find $lncPoolDB -maxdepth 1 -type f -name "$myfsp.gene2trx.lst"`; chomp $fsynF; 
    open $FFH, "<$fsynF";
    while (<$FFH>){
        chomp;
        @fsynFilds=split("\t", $_);
        $joint_fkey=$ssp."_".$fsynFilds[1];
        $rev_joint_fkey=$ssp."_".$fsynFilds[0];
        push @{$POOL{$joint_fkey}}, $fsynFilds[0];    #### hum_NONHSAG*** (NONHSAT*** NONHSAT*** NONHSAT*** ...) 
        $RVPOOL{$rev_joint_fkey}=$fsynFilds[1];      #### hum_NONHSAT*** NONHSAG***
    }
    close $FFH;
}

=head
@aRVPOOL=keys %RVPOOL;
foreach (0..19){
    print $aRVPOOL[$_], "\ta\t", $RVPOOL{$aRVPOOL[$_]}, "\n";
}
=cut

%inhash=(); %rev_inhash=();
open INF, "<$inFAlst";
while ($inf=<INF>){
    if ($inf!~/^#/){
    chomp $inf;
    @F=split("\t", $inf);
    $inhash{$F[0]}=$F[1];        #### TCONS_*** FIF
    $rev_inhash{$F[1]}=$F[0];     #### FIF TCONS_***
    }
}

@pair=(); @id_pair=();
@other=(); @id_other=();
for (0..$#code_arr){
    if ( defined $rev_inhash{$code_arr[$_]} ){
       if ( $_<=1 ){
          push @pair, $code_arr[$_];
          push @id_pair, $rev_inhash{$code_arr[$_]};
       }else{
          push @other, $code_arr[$_];
          push @id_other, $rev_inhash{$code_arr[$_]};
       }
    }
}

if (scalar(@pair)==0){
  die "Error in $inFAlst!\n";
}elsif ((scalar(@pair)==1) && (scalar(@other)==0)){
  $mode=1;
}elsif ((scalar(@pair)==2) && (scalar(@other)==0)){
  $mode=2;
}elsif ((scalar(@pair)==2) && (scalar(@other)>0)){
  $mode=3;
}

%trxhash=();
if ($mode==1){

	($NONCODEgene, $NONCODEtrx)=read_cand2NONCODE($cand2NONCODE);

	if ($NONCODEgene=~/NONHSAG/){ 
       @mou_arr_trx=@temp_arr_gene=();
       if (defined $SYN{"mou_".$NONCODEgene}){
          @mou_arr_gene=@{$SYN{"mou_".$NONCODEgene}};
	      for (@mou_arr_gene){ push @mou_arr_trx, @{$POOL{"mou_".$_}}; }
       }
       for (@mou_arr_trx){
           $temp_mou_trx=$_;
        for (2..$#code_arr){
           $temp_code=$code_arr[$_];
           $temp_sp=$A{$code_arr[$_]};   
           if (defined $SYN{$temp_sp."_".$NONCODEgene}){
            @temp_arr_gene=@{$SYN{$temp_sp."_".$NONCODEgene}};
            for (@temp_arr_gene){ 
                 push @{$trxhash{$NONCODEtrx."-linker-".$temp_mou_trx}{$temp_code}}, @{$POOL{$temp_sp."_".$_}}; 
            }
           }else {               
                 $trxhash{$NONCODEtrx."-linker-".$temp_mou_trx}{$temp_code}="NOTRD";
           }
        }
       }                 
       $rmode="a";     #hum to mou1 mou2 mou3......
	}elsif ($NONCODEgene=~/NONMMUG/){
       @hum_arr_trx=@temp_arr_gene=();
       if (defined $RVSYN{"mou_".$NONCODEgene}){
          @hum_arr_gene=@{$RVSYN{"mou_".$NONCODEgene}};
          for (@hum_arr_gene){ push @hum_arr_trx, @{$POOL{"hum_".$_}}; }
       }
       for (@hum_arr_trx){
          $temp_hum_trx=$_;
          $temp_hum_gen=@{$RVPOOL{"hum_".$temp_hum_trx}};
          for (2..$#code_arr){
              $temp_code=$code_arr[$_]; 
              $temp_sp=$A{$code_arr[$_]};
              if (defined $SYN{$temp_sp."_".$temp_hum_gen}){
               @temp_arr_gene=@{$SYN{$temp_sp."_".$temp_hum_gen}};
               for (@temp_arr_gene){                      
                 push @{$trxhash{$temp_hum_trx."-linker-".$NONCODEtrx}{$temp_code}}, @{$POOL{$temp_sp."_".$_}}; 
               }             
              }else {
                 $trxhash{$temp_hum_trx."-linker-".$NONCODEtrx}{$temp_code}="NOTRD";
              }              
          }
       }
       $rmode="b";     #mou to humG1t1 humG1t2 || humG2t1 humG2t2 humG2t3 .......
    }

}elsif ($mode==2){

	($NONCODEgene, $NONCODEhtrx, $NONCODEmtrx)=read_cand2NONCODE($cand2NONCODE);

    if ($NONCODEgene=~/NONHSAG/){
       @temp_arr_gene=();
       for (2..$#code_arr){
           $temp_code=$code_arr[$_];
           $temp_sp=$A{$code_arr[$_]};     
           if (defined $SYN{$temp_sp."_".$NONCODEgene}){
            @temp_arr_gene=@{$SYN{$temp_sp."_".$NONCODEgene}};
            for (@temp_arr_gene){ 
                 push @{$trxhash{$NONCODEhtrx."-linker-".$NONCODEmtrx}{$temp_code}}, @{$POOL{$temp_sp."_".$_}};
            }
           }else {
                 $trxhash{$NONCODEhtrx."-linker-".$NONCODEmtrx}{$temp_code}="NOTRD";
           }
       }
    }     
    $rmode="c";        # hum to mou

}elsif ($mode==3){

    ($NONCODEgene, $NONCODEhtrx, $NONCODEmtrx)=read_cand2NONCODE($cand2NONCODE);

	for (@other){
        $temp_code=$_;
        #$temp_sp=$A{$_}; 
        push @{$trxhash{$NONCODEhtrx."-linker-".$NONCODEmtrx}{$temp_code}}, $rev_inhash{$temp_code};
    }
    $rmode="c";

}

=head
print $rmode, "    rmode\n";
foreach $kkk (keys %trxhash){
    %ttt=%{$trxhash{$kkk}};
    for $ikkk (keys %ttt){
        print $kkk, "\tb\t", $ikkk, "\tc\t" , join(" ", @{$ttt{$ikkk}}),"\n";
    }
}
=cut

#$mp = Parallel::ForkManager->new($procs);
@pair_keys = keys %trxhash;
for  (@pair_keys){ 
#  $mp->start and next;

  $my_pair_key=$_; print "PairID\t ","$my_pair_key\n";
  %other_support=%{$trxhash{$my_pair_key}};

  &runpip($rmode, $my_pair_key, \%other_support);

  if ($run_mode <= 1){
      &pair_meme_fimo_filter($pair_dirname, $trd_sp_dirname, \@pair);
  }

  if ($run_mode >= 1){
      #$randNUM=1000000000000000*rand; $randID="RAND"."$randNUM";
      #$outsvg = $randID.".svg";
      `mkdir -p $svgdir`;
      $outsvg = $svgdir."/".$my_pair_key.".svg";
      &runSVG;
      print "\n-----------------------\n++SVG finished!++\n-----------------------\n\n";
  }  

#  $mp->finish;
}



sub runpip {

$myrunmode=shift;
$pair_two=shift;    ### "-linker-" join
@pair=split("-linker-", $pair_two);
$othersupport_hash=shift; %othersupport_hash=%{$othersupport_hash};
#$randNUM=1000000000000000*rand; $randID="RAND"."$randNUM";

$pair_dirname="${basedir}/DIR_".$pair_two;
$other_dirname=$pair_dirname."/cluster";

$trd_sp_dirname=$other_dirname."/input.dir_3rd";
$res_diag=$other_dirname."/Selected_Align";
$res_3rdfimo=$other_dirname."/compare_3rd/cluster";
$res_dir=$other_dirname."/DP_EUC_3rd";

if ( $run_mode <= 1 ){
`rm -rf $pair_dirname; mkdir -p $trd_sp_dirname`;
#`rm -rf $other_dirname|mkdir -p $trd_sp_dirname`;
}

%faHash=(); %faLength=(); 
if ($myrunmode=~/a/){
  ($faHash, $faLength)=&get_fa_len($fastaDB, "Hash", \%inhash);
  %faHash=%{$faHash}; %faLength=%{$faLength};
  $mouTrxName=$pair[1];
  $tempMoulncPool=`find $lncPoolDB -maxdepth 1 -type f -name "$B{'SEC'}.TrX.fa"`; chomp $tempMoulncPool; 

  %tempfa=();
  open FA, $tempMoulncPool;
  while(<FA>){
   chomp;  
   if(/^>(\S+)/){
     $id=$1;
   }else{  
     $tempfa{$id}.=uc($_);
   }       
  }         
  close FA; 
 
  $faHash{$mouTrxName}=$tempfa{$mouTrxName};
  $faLength{$mouTrxName}=length($tempfa{$mouTrxName});
}elsif ($myrunmode=~/b/){
  ($faHash, $faLength)=&get_fa_len($fastaDB, "Hash", \%inhash); 
   %faHash=%{$faHash}; %faLength=%{$faLength};
  $humTrxName=$pair[0];
  $tempHumlncPool=`find $lncPoolDB -maxdepth 1 -type f -name "$B{'FIR'}.TrX.fa"`; chomp $tempHumlncPool; 

  %tempfa=();             
  open FA, $tempHumlncPool;                           
  while(<FA>){           
   chomp;                
   if(/^>(\S+)/){        
     $id=$1;             
   }else{                
     $tempfa{$id}.=uc($_);                            
   }                     
  }                      
  close FA;              
 
  $faHash{$humTrxName}=$tempfa{$humTrxName};          
  $faLength{$humTrxName}=length($tempfa{$humTrxName});
}elsif ($myrunmode=~/c/){
  ($faHash, $faLength)=&get_fa_len($fastaDB, "Hash", \%inhash);
  %faHash=%{$faHash};              #### FIF       ATGCGTCGAT...
  %faLength=%{$faLength};          #### FIF       560  
}


if ( $run_mode <= 1 ){
foreach $mytempSpCode (keys %othersupport_hash){
   $mytempSpCodeLen = $mytempSpCode."LEN";
   $trd_sp_file=$trd_sp_dirname."/".$mytempSpCode.".fa"; 
   $trd_sp_lenfile=$trd_sp_dirname."/".$mytempSpCode.".len"; 
   open $mytempSpCode, ">>$trd_sp_file";
   open $mytempSpCodeLen, ">>$trd_sp_lenfile";

   $tempsplncPool=`find $lncPoolDB -maxdepth 1 -type f -name "$B{$mytempSpCode}.TrX.fa"`; chomp $tempsplncPool;
   @tempspcodeTrx= @{$othersupport_hash{$mytempSpCode}};
   ($tempfaHash, $tempfaLength)=&get_fa_len($tempsplncPool, "Array", \@tempspcodeTrx);
   foreach $templenkey (keys %{$tempfaLength}){
     print $mytempSpCode ">",$templenkey,"\n",$tempfaHash->{$templenkey},"\n";
     print $mytempSpCodeLen $templenkey,"\t",$tempfaLength->{$templenkey},"\n";

     #$faLengthTrd{$mytempSpCode}{$templenkey}=$tempfaLength->{$templenkey};

   }

   close $mytempSpCode; close $mytempSpCodeLen;
}
}

}


sub runSVG{

#if ( not $? ){

@temp_sp_dir=`find $res_dir -maxdepth 1 -type d -name "*.0001"`; chomp @temp_sp_dir; 
if (scalar(@temp_sp_dir) >=1){
for (@temp_sp_dir){
    $tempspkeydir=$_; $tempspkeyID=$tempspkeydir; $tempspkeyID=~s/\/\S+DP_EUC_3rd\///;
    $aln_cluster=$tempspkeyID; $aln_cluster=~s/\w+\.(\w+)\.0001/$1/;
    @test_dp_file=`find $tempspkeydir -maxdepth 1 -type f -name "*_3rd.dp"`; chomp @test_dp_file;
    if (scalar(@test_dp_file) >=1){
       %hashHigh=(); $hashHigh=%hashHigh;
       %hashHighscore=(); $hashHighscore=%hashHighscore;
       for (@test_dp_file){
            $tempdpfile=$_; $tempdpID=$tempdpfile; $tempdpID=~s/\/\S+0001\/\w+\.(\S+\.C\d+_3rd)\.dp/$1/; 
            ($hashHighscore, $hashHigh)=&trd_dp_Highest($tempdpfile, $tempdpID, \%hashHighscore, \%hashHigh);
       }
       %hashHighscore=%{$hashHighscore};
       %hashHigh=%{$hashHigh}; 
       @scoresort=sort {$hashHighscore{$b}<=>$hashHighscore{$a}} keys %hashHighscore;
       $tempoutHdp=$res_dir."/".$tempspkeyID.".vs.".$scoresort[0];
       open $tempspkeyID, ">$tempoutHdp" or die "Can't write to file '$tempoutHdp' [$!]\n";;
       print $tempspkeyID $hashHigh{$scoresort[0]};
       close $tempspkeyID;

       $tempoutHdp_3rdfimo=$tempoutHdp; $tempoutHdp_3rdfimo=~s/\S+\///;
       $tempoutHdp_3rdfimo_selected=$tempoutHdp_3rdfimo; $tempoutHdp_3rdfimo_selected=$tempoutHdp_3rdfimo_selected.".fimo";
       `perl $scrdir/parsing_DP_FIMO_SVG.pl $res_diag/${aln_cluster}_*.align $res_3rdfimo/$tempoutHdp_3rdfimo.sort $tempoutHdp 0001 $res_dir/$tempoutHdp_3rdfimo_selected`;
       if ( not $? ){ print "scanning dir $tempspkeyID finished\n"; }

    }


$tempoutHdp_3rdfimo_selected=`find $res_dir -maxdepth 1 -type f -name "$tempspkeyID.*.fimo.s"`; chomp $tempoutHdp_3rdfimo_selected;

if ( -s $tempoutHdp_3rdfimo_selected ){
print $tempoutHdp_3rdfimo_selected,"\n";

if ($tempspkeyID=~/(\w+)\.(\w+)\.0001/){
  $sp=$1;
  $D=$2;
}

open $tempspkeyID, "<$tempoutHdp_3rdfimo_selected";
  while (<$tempspkeyID>){
    next if (/^#/);
    chomp;
    @A=split("\t", $_);
    if ($A[1]=~/$sp/){
      push @{$domain2motif{$D}{$sp}}, $_;
      push @{$sp2domain_s{$sp}{$D}}, $A[2];
      push @{$sp2domain_e{$sp}{$D}}, $A[3];
    }
  }
close $tempspkeyID;
}

}
}


@temp_Ali_dir=`find $res_diag -maxdepth 1 -type f -name "*.fimo.s"`; chomp @temp_Ali_dir;
for (@temp_Ali_dir){

$tempAlikeydir=$_;  $tempAlikeyID=$tempAlikeydir; $tempAlikeyID=~s/\/\S+Selected_Align\///;

if ($tempAlikeyID=~/(\w+)_FIR\S+/){
  $D=$1;
}

open $tempAlikeyID, "<$tempAlikeydir";
  while (<$tempAlikeyID>){
    next if (/^#/);
    chomp;
    @B=split("\t", $_);
      push @{$domain2motif{$D}{$B[1]}}, $_;
      push @{$sp2domain_s{$B[1]}{$D}}, $B[2];
      push @{$sp2domain_e{$B[1]}{$D}}, $B[3];
      push @redun_color, $B[0];
  }
close $tempAlikeyID;

}

@color=@{[uniq_array(\@redun_color)]};
foreach $k (keys %sp2domain_s){
  %in_s=%{$sp2domain_s{$k}};
  foreach $kk (keys %in_s){
   @s=@{$in_s{$kk}};
   @e=@{$sp2domain_e{$k}{$kk}};  
   $sp2domain{$k}{$kk}=min(@s)."_".max(@e);
  }
}

#### %faLength #### $faLengthTrd{$mytempSpCode}{$templenkey}####

 %domain2block=(); # null
# %sp2domain	FIR=>{D1=>3_100; D2=>6_400; }
# %domain2motif	D1=>{FIR=>@{fimo lines}; SEC=>@{fimo lines}; }
# @color	uniq kinds of motifs
# %faLength	ok

&domainSVG( 2479, 10000, \@color, \%sp2domain, \%domain2motif, \%domain2block, \%faLength, "MYSVG" );

#1400      -> width of SVG picture
#3000      -> height of SVG picture
#%sp2domain	-> eg. {FIR=>{D1=>40_100, D2=>200_350, }, SEC=>{}, }
#%domain2motif	-> eg. D1=>FIR=>@{fimo lines}  (PS: motif coordinates are relative to "new length" of specific domain)
#%domain2block	-> eg. D1=>FIR=>block1=>3_20   (PS: motif coordinates are relative to "new length" of specific domain)
##@color    -> how many different kinds of motifs
#SVG      -> fileHandle of SVG output

#}


}


sub config{
        $conf = shift;
        open CONF, "<$conf";
        while (<CONF>){
		chomp;
                next if (/^#/);
                if (/^(\w+)\s+(\S+)/){
                        $config{$1} = $2;
                }
        }
        close CONF;

        $ENV{'ScrDBdir'}=$config{"ScrDBdir"};
        $ENV{'basedir'}=$config{"basedir"};

		$myScr = $ENV{'ScrDBdir'}.$config{"scrdir"}; 
		$myMdb = $ENV{'ScrDBdir'}.$config{"motifDB"}; 
		$mySyn = $ENV{'ScrDBdir'}.$config{"syngeneDB"}; 
		$myLnc = $ENV{'ScrDBdir'}.$config{"lncPoolDB"}; 

		$myBse = $ENV{'basedir'}; 
		$myFas = $ENV{'basedir'}.$config{"fastaDBdir"}; 
		$mySvg = $ENV{'basedir'}.$config{"svgdir"}; 

		$myCan = $config{"cand2NONCODE"}; 
		$myMod = $config{"runmode"};

		return ( $myScr, $myMdb, $mySyn, $myLnc, $myBse, $myFas, $mySvg, $myCan, $myMod );

        1;
}


sub get_fa_len{

 %fa=%_faHash=%_faLength=();
 $FA=shift; $ArrOrHash=shift;
 $temp_inhash=shift; 
 if ($ArrOrHash=~/Hash/){ 
     %temp_inhash=%{$temp_inhash}; 
 }elsif ($ArrOrHash=~/Array/){ 
     @temp_inhash=@{$temp_inhash}; 
 }

print $FA,"\n";
 open MYFA, $FA;
 while(<MYFA>){
   chomp;
   if(/^>(\S+)/){
     $id=$1;
   }else{
     $fa{$id}.=uc($_);
   }
 }
 close MYFA;

 if ($ArrOrHash=~/Hash/){ 
 foreach $fa_key (keys %fa){
   if (defined $temp_inhash{$fa_key}){
      #$_faHash{$temp_inhash{$fa_key}}=$fa{$fa_key};
      #$_faLength{$temp_inhash{$fa_key}}=length($fa{$fa_key});
      $_faHash{$fa_key}=$fa{$fa_key};
      $_faLength{$fa_key}=length($fa{$fa_key});
   }
 }
 return (\%_faHash, \%_faLength);
 }

 if ($ArrOrHash=~/Array/){
 foreach $fa_key (keys %fa){
   if ( $fa_key~~@temp_inhash ){
      $_faHash{$fa_key}=$fa{$fa_key};
      $_faLength{$fa_key}=length($fa{$fa_key});
   }
 }
 return (\%_faHash, \%_faLength);
 }

 1;

}


sub read_cand2NONCODE{
  $mycand2NONCODE = shift;
  $hExist=$mExist="no";
  open C2N, "<$mycand2NONCODE";
  while (<C2N>){
     chomp;
     next if (/^#/); 
     if(/^(\S+)\s+(NONHSAG\d+)/){
        $hExist="yes";
        $hNONintrx=$1;
        $hNONgene=$2;
     }elsif(/^(\S+)\s+(NONMMUG\d+)/){
        $mExist="yes";
        $mNONintrx=$1;
        $mNONgene=$2;
     }
  }
  close C2N;

  if (($hExist eq "yes") && ($mExist eq "yes")){
     return ($hNONgene, $hNONintrx, $mNONintrx);
  }elsif (($hExist eq "yes") && ($mExist eq "no")){
     return ($hNONgene, $hNONintrx);
  }elsif (($hExist eq "no") && ($mExist eq "yes")){
     return ($mNONgene, $mNONintrx);
  }else {
     die "Error in $mycand2NONCODE!\n";
  }

  1;

}


sub pair_meme_fimo_filter{

$dir=shift; $trddir=shift; 
$pair=shift; @_pair=@{$pair};
$pid="_".$_pair[0]."_".$_pair[1];
$motif_db=$motifDB;

$PairLncFa=$dir."/"."lnc".$pid.".fa";	$PLF="lnc".$pid;
$bash=$dir."/"."lnc".$pid.".sh";	$BASH="bash".$pid;
$meme=$dir."/meme".$pid.".out";
$fimo=$dir."/fimo".$pid.".out";
#$DP=$dir."/DP".$pid.".out";

open $PLF, ">$PairLncFa";
print $PLF ">","FIR","\n", $faHash{$_pair[0]},"\n",">","SEC","\n",$faHash{$_pair[1]},"\n";

open $BASH, ">>$bash";
print $BASH "#!/bin/bash","\n";
print $BASH "cd ",$dir,"\n";
print $BASH 'source /Share/home/zhangqf2/.bash_profile',"\n";
print $BASH "perl2\n";
print $BASH 'export LD_LIBRARY_PATH=/bin:/usr/local/gcc-4.9.2/lib64:/opt/openmpi-1.8.7/lib:/opt/intel-xe_2015u3/lib/intel64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:$HOMEZ/usr/openmpi/lib:/path/to/bcl-3.1.0-Linux-x86_64:$HOMEZ/usr/readline-6.3/lib:$HOME/tools/binutils/lib:$HOMEZ/usr/lib64:$HOMEZ/tmp/expat/lib:$HOMEZ/usr/meme/lib:$HOMEZ/shaodi/app/samtools/zlib-1.2.8/build/lib:$HOMEZ/shaodi/app/samtools/ncurses-6.0/build/lib:$HOMEZ/jellyfish/lib:$HOME/tools/stamp/gsl_dir/gsl/lib:$HOME/tools/binutils/lib:$HOME/usr/lib64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:$HOME/tmp/expat/lib:$HOME/tools/meme/lib:$HOMEZ/usr/openmpi/lib:$HOMEZ/shaodi/app/samtools/zlib-1.2.8/build/lib:$HOMEZ/shaodi/app/samtools/ncurses-6.0/build/lib:$HOMEZ/jellyfish/lib:$LD_LIBRARY_PATH',"\n";
print $BASH "\n\n";
print $BASH "/Share/home/zhangqf2/tools/meme/bin/meme ",$PairLncFa," -p 12 -nostatus -time 36000 -dna -revcomp -text -mod anr -nmotifs 20 -minw 5 -maxw 30 -maxsites 600 -maxsize 1000000 >",$meme,"\n";       ### ?? meme-chip not considering revcomp (lots of RNA segments containing 1 or 2 motifs), but meme maybe should consider revcomp (only several RNA sequences); purpose of the former is identifying significant (top) motifs, but the latter is obtaining the full landscape of de novo motifs and their (maybe also conserved) revcomp position; PS: fimo should considers both strands.      e.g.   ....atg....atg....cat....cat....   in meme mode, 1 motif 2 diretions (revcomp); 2 motifs 1 direction (no revcomp)
print $BASH "/Share/home/zhangqf2/tools/meme/bin/fimo --verbosity 1 --text ",$meme," ",$PairLncFa," >",$fimo,".denovo","\n"; 
print $BASH "/Share/home/zhangqf2/tools/meme/bin/fimo --verbosity 1 --text ",$motif_db," ",$PairLncFa," >",$fimo,".clipNdb","\n";
print $BASH "cat $fimo",".denovo ",$fimo,".clipNdb ", ">$fimo.combined","\n";  
print $BASH "cat ",$fimo,".combined","|awk -F \'\\t\' \'\$7<0.001{print \$0}\' >",$fimo,".filter","\n";
print $BASH "sort -k2,2 -k3,3n -k4,4n ",$fimo,".filter >",$fimo,".sort","\n";
print $BASH "perl $scrdir/fimo_pair_motif_filtering.pl ", $fimo,".sort >", $fimo,".sort.pair-filter-pre","\n";

print $BASH "select=\$(cat $fimo.sort.pair-filter-pre|grep -v \'#\'|awk -F \'\\t\' '!dup[\$1]++'|awk -F \'\\t\' -v ORS=\" \" \'{print \$1}\')","\n";
print $BASH "perl $scrdir/grep_meme_recognized_by_mcast_v1.0.pl ","$motif_db ",'$select >', $meme,".reselect","\n";
print $BASH "perl $scrdir/grep_meme_recognized_by_mcast_v1.0.pl $meme \$select >> $meme.reselect","\n";
  
print $BASH "rm -rf $dir/meme_temp|mkdir -p $dir/meme_temp","\n";
print $BASH "perl $scrdir/memefimo2fa.pl ", $meme,".reselect ",$fimo,".sort.pair-filter-pre ",$dir,"/meme_temp","\n";
print $BASH "/Share/home/zhangqf2/tools/meme/bin/sites2meme_new $dir/meme_temp >$meme.pre-update", "\n";
print $BASH "newselect=\$(perl $scrdir/cmp_motif_inf.pl $meme.reselect $meme.pre-update)", "\n";
print $BASH "perl $scrdir/grep_meme_recognized_by_mcast_v1.0.pl $meme.pre-update \$newselect > $meme.update","\n";

print $BASH "/Share/home/zhangqf2/tools/meme/bin/fimo --verbosity 1 --text $meme.update ",$PairLncFa," >",$fimo,".pair-filter","\n";
print $BASH "sort -k2,2 -k3,3n -k4,4n ",$fimo,".pair-filter >",$fimo,".sort.pair-filter","\n";
print $BASH "perl $scrdir/motif_block_merge_edgeR.pl 2 120 ",$fimo,".sort.pair-filter ","fimo",$pid,".out.sort.p-f.sigwindows >",$fimo,".sort.p-f.s-w.domain","\n";
print $BASH "for cl in $dir/cluster/C*.sort; do if [ -f \"\$cl\" ]; then perl $scrdir/fimo_motif_coverage_filter.pl \$cl; fi; done\n";
print $BASH "for cl in $dir/cluster/C*.sort; do if [ -f \"\$cl\" ]; then perl $scrdir/FIMO_DP_v1.0.2_MotifPvalue.pl $PairLncFa \$cl >\$cl.dp; fi; done\n";
print $BASH "for parsing in $dir/cluster/C*.sort; do if [ -f \"\$cl\" ]; then perl $scrdir/parsing_DP_FIMO_v1.0.pl \$parsing \$parsing.dp 0001 \$parsing.fimo >\$parsing.align; fi; done\n";
print $BASH "\n\n";

print $BASH "DIRdir=\"$dir\"","\n";
print $BASH "SCRdir=\"$scrdir\"","\n";
print $BASH "TRDdir=\"$trddir\"","\n";
print $BASH 'dbscan_domain=`find ${DIRdir} -maxdepth 1 -type f -name \'*.domain\'`',"\n";
print $BASH 'dbscan_dbscan=`find ${DIRdir} -maxdepth 1 -type f -name \'*.dbscan\'`',"\n";
print $BASH 'while [ ! -f $dbscan_domain ]||[ ! -s $dbscan_dbscan ]',"\n";
print $BASH 'do',"\n";
print $BASH '  dbscan_domain=`find ${DIRdir} -maxdepth 1 -type f -name \'*.domain\'`',"\n";
print $BASH '  dbscan_cluster=`find ${DIRdir} -maxdepth 1 -type f -name \'*.dbscan\'`',"\n";
print $BASH '  sleep 5',"\n";
print $BASH 'done',"\n";
print $BASH ' ',"\n";
print $BASH 'if [[ -s $dbscan_domain ]] ; then',"\n";
print $BASH '   dbscan=`find ${DIRdir} -maxdepth 1 -type f -name \'*.dbscan\'`',"\n";
#print $BASH '   while [ ! -f $dbscan ]||[ $(ls ${DIRdir}/cluster/*.align|wc -l) -ne $(cat $dbscan|tail -n1|awk -F \'\t\' \'{print $3}\') ]',"\n";
print $BASH '   while [ ! -f $dbscan ]||[ $(ls ${DIRdir}/cluster/*.align|wc -l) -ne $(ls ${DIRdir}/cluster/C*.sort|wc -l) ]',"\n";
print $BASH '   do',"\n";
print $BASH '     dbscan=`find ${DIRdir} -maxdepth 1 -type f -name \'*.dbscan\'`',"\n";
print $BASH '     echo $(cat $dbscan|tail -n1|awk -F \'\t\' \'{print $3}\')',"\n";
print $BASH '     sleep 5',"\n";
print $BASH '   done',"\n";
print $BASH ' ',"\n";
print $BASH '   Cluster_align=`find ${DIRdir}/cluster -maxdepth 1 -type f -name \'C*.align\'`',"\n";
print $BASH '   if [ -n "$Cluster_align" ]; then cat $Cluster_align >${DIRdir}/cluster/total.cluster_align; else exit 0; fi',"\n";

#print $BASH '   mkdir -p ${DIRdir}/cluster/Diagonal_Align',"\n";
#print $BASH '   perl $SCRdir/DiagonalScore_DP.pl ${DIRdir}/cluster/total.cluster_align ${DIRdir}/cluster/Diagonal_Align/diagonal_cluster.out',"\n";
#print $BASH '   for cl in `grep "^C" ${DIRdir}/cluster/Diagonal_Align/diagonal_cluster.out`; do cp ${DIRdir}/cluster/${cl}_* ${DIRdir}/cluster/Diagonal_Align; done',"\n";
#print $BASH '   cat ${DIRdir}/cluster/Diagonal_Align/C*.align >${DIRdir}/cluster/Diagonal_Align/diagonal.cluster_align',"\n";
print $BASH '   mkdir -p ${DIRdir}/cluster/Selected_Align',"\n";
print $BASH '   perl $SCRdir/ClusterSelection.pl ${DIRdir}/cluster/total.cluster_align 2 ${DIRdir}/cluster/Selected_Align/selected_cluster.out',"\n";
print $BASH '   for cl in `grep "^C" ${DIRdir}/cluster/Selected_Align/selected_cluster.out`; do cp ${DIRdir}/cluster/${cl}_* ${DIRdir}/cluster/Selected_Align; done',"\n";
print $BASH '   Selected_cluster_align=`find ${DIRdir}/cluster/Selected_Align -maxdepth 1 -type f -name \'C*.align\'`',"\n";
print $BASH '   if [ -n "$Selected_cluster_align" ]; then cat $Selected_cluster_align >${DIRdir}/cluster/Selected_Align/selected.cluster_align; else exit 0; fi',"\n"; 

print $BASH ' ',"\n\n";
print $BASH '   SP=(TRD FOR FIF SIX SEV EIG NIN TEN ELE TWE)',"\n";
print $BASH '   let spnum=${#SP[*]}-1',"\n";
print $BASH ' ',"\n";
print $BASH '   for sp in `seq 0 $spnum`',"\n";
print $BASH '   do',"\n";
print $BASH '     if [[ -f $TRDdir/${SP[$sp]}.fa ]] ; then',"\n";
print $BASH '       seq_3rd=$TRDdir/${SP[$sp]}.fa',"\n";
print $BASH '       seq_3rd_name=${seq_3rd##*/}; seq_3rd_name=${seq_3rd_name%.fa}',"\n";
print $BASH '       perl $SCRdir/compare_3rd_v1.0.pl ${DIRdir}/cluster/Selected_Align/selected.cluster_align $seq_3rd ${DIRdir}/cluster',"\n";
print $BASH ' ',"\n";
print $BASH '       Cluster_3rd_fimo=(`find ${DIRdir}/cluster/compare_3rd -maxdepth 1 -type f -name "${seq_3rd_name}*.fimo"`)',"\n";
print $BASH '       while [ ${#Cluster_3rd_fimo[*]} -ne $(find ${DIRdir}/cluster/Selected_Align -size +0c -maxdepth 1 -type f -name \'C*.align\'|wc -l) ]',"\n";
print $BASH '       do',"\n";
print $BASH '         Cluster_3rd_fimo=(`find ${DIRdir}/cluster/compare_3rd -maxdepth 1 -type f -name "${seq_3rd_name}*.fimo"`)',"\n";
print $BASH '         sleep 5',"\n";
print $BASH '       done',"\n";
print $BASH ' ',"\n";
print $BASH '       perl $SCRdir/FIMO_DP_EUC_ThirdPlus_SVG.pl ${DIRdir}/cluster/Selected_Align/selected.cluster_align $seq_3rd ${DIRdir}/cluster',"\n";  
print $BASH '     fi',"\n";
print $BASH '   done',"\n";
print $BASH ' ',"\n";
print $BASH 'fi',"\n";

`chmod 777 $bash`;

close $PLF;
close $BASH;

#`ssh bnode02 $bash`;
#`ssh node544 $bash`;
#`/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub -K -q Z -n 20 <$bash`; print "\n-----------------------\n++bsub jobs finished!++\n-----------------------\n\n";
print `$bash`,"\n"; if ( not $? ) { print "\n-----------------------\n++pair-jobs finished!++\n-----------------------\n\n"; }

}


sub domainSVG{

 $width = shift; $height = shift;
 $color = shift; @color = @{$color};
 $sp2domain = shift; %sp2domain = %{$sp2domain};
 $domain2motif = shift; %domain2motif = %{$domain2motif};
 $domain2block = shift; %domain2block = %{$domain2block};
 $faLength = shift; %faLength = %{$faLength};
 $FH = shift;

 open $FH, ">$outsvg";

 %seen=(); @unique=();
 foreach $value (@color) {
        if (!$seen{$value}) {
           push @unique, $value;
           $seen{$value} = 1;
        }
 }
 @sort = sort (@unique);
 @sort = reverse (@sort);
 #<> $motifnumber=scalar (@sort);

 foreach $sp2domain_key (keys %sp2domain){
        %sp2domain_inner_hash=%{$sp2domain{$sp2domain_key}};
        $len_sp2domain{$sp2domain_key}=scalar( keys %sp2domain_inner_hash );
 }
 $domain_maxlen_key= shift ( @{[ sort {$len_sp2domain{$b}<=>$len_sp2domain{$a}} keys %len_sp2domain ]} );
 @domain_maxlen_key_arr=sort keys %{$sp2domain{$domain_maxlen_key}};
 #<> $domain_maxlen_num=scalar @domain_maxlen_key_arr;

 %color=(); %domaincolor=();
 #@RAWcolor = ("red","yellow","blue","lightgreen","green","lightskyblue","brown","lightpink","purple","orange","black","navajowhite","orchid","rosybrown","springgreen","tomato","wheat","sienna","violet","darkgray","cyan","darkred","darkseagreen","darkblue","aquamarine","lime","PeachPuff","Peru","Linen","olive","gold","DarkKhaki","Ivory","greenyellow","beige","lime","honeydew","teal","azure","indigo","plum","navy","pink","crimson","Turquoise","MintCream","LimeGreen","OliveDrab","LemonChiffon","PapayaWhip","tan","DarkOrange","Chocolate","Coral","FireBrick","DimGray","Gainsboro","MistyRose","Tomato","SaddleBrown","Moccasin","Cornislk");
 @RAWcolor = ("black", "navy", "darkblue", "mediumblue", "blue", "darkgreen", "green", "teal", "darkcyan", "deepskyblue", "darkturquoise", "mediumspringgreen", "lime", "springgreen", "aqua", "midnightblue", "dodgerblue", "lightseagreen", "forestgreen", "seagreen", "darkslategray", "limegreen", "mediumseagreen", "turquoise", "royalblue", "steelblue", "darkslateblue", "mediumturquoise", "indigo", "darkolivegreen", "cadetblue", "cornflowerblue", "mediumaquamarine", "dimgray", "slateblue", "olivedrab", "slategray", "lightslategray", "mediumslateblue", "lawngreen", "chartreuse", "aquamarine", "maroon", "purple", "olive", "gray", "skyblue", "lightskyblue", "blueviolet", "darkred", "darkmagenta", "saddlebrown", "darkseagreen", "lightgreen", "mediumpurple", "darkviolet", "palegreen", "darkorchid", "sienna", "brown", "darkgray", "lightblue", "greenyellow", "paleturquoise", "lightsteelblue", "powderblue", "firebrick", "darkgoldenrod", "mediumorchid", "rosybrown", "darkkhaki", "silver", "mediumvioletred", "indianred", "peru", "chocolate", "tan", "lightgray", "thistle", "orchid", "goldenrod", "palevioletred", "crimson", "gainsboro", "plum", "burlywood", "lightcyan", "lavender", "darksalmon", "violet", "palegoldenrod", "lightcoral", "khaki", "aliceblue", "honeydew", "azure", "sandybrown", "wheat", "beige", "whitesmoke", "mintcream", "ghostwhite", "salmon", "antiquewhite", "linen", "lightgoldenrodyellow", "oldlace", "red", "fuchsia", "deeppink", "orangered", "tomato", "hotpink", "coral", "darkorange", "lightsalmon", "orange", "lightpink", "pink", "gold", "peachpuff", "navajowhite", "moccasin", "bisque", "mistyrose", "blanchedalmond", "papayawhip", "lavenderblush", "seashell", "cornsilk", "lemonchiffon", "floralwhite", "snow", "yellow", "lightyellow", "ivory", "white"); 

=head
 for ($i=0; $i <=@sort-1; $i++) {
     if ($i<=scalar(@RAWcolor)-1){
        $color{$sort[$i]}=$RAWcolor[$i];
     }elsif ( ($i>scalar(@RAWcolor)-1) && ($i<=2*(scalar(@RAWcolor)-1)) ) {
        $color{$sort[$i]}=$RAWcolor[$i-scalar(@RAWcolor)];
     }elsif ( ($i>2*(scalar(@RAWcolor)-1)) && ($i<=3*(scalar(@RAWcolor)-1)) ){
        $color{$sort[$i]}=$RAWcolor[$i-2*scalar(@RAWcolor)];
     }

 }
=cut

 $loop=0;
 for ($i=0; $i <=@sort-1; $i++) {
     if ( $i>($loop+1)*(scalar(@RAWcolor)-1) ){ $loop++; }
     if ( ($i>$loop*(scalar(@RAWcolor)-1)) && ($i<=($loop+1)*(scalar(@RAWcolor)-1)) ) {
        $color{$sort[$i]}=$RAWcolor[$i-$loop*scalar(@RAWcolor)];
     }elsif ($i==0){
        $color{$sort[$i]}=$RAWcolor[$i];
     }
 }

 for ($i=0; $i <=@domain_maxlen_key_arr-1; $i++) {
        $domaincolor{$domain_maxlen_key_arr[$i]}=$RAWcolor[$#RAWcolor-$i];
 }


 @tempARR=();
 foreach $_key (sort keys %faLength){
         push @tempARR, $faLength{$_key};
 }
 @tempARR = (sort{$b<=>$a} @tempARR);
 $standardLen=shift @tempARR;
 $ratio=($width - 200)/$standardLen;


$svg = SVG->new( 'width', $width, 'height', $height );
$x = $svg->group( id=>'group_x', style=>{'stroke'=>'black','stroke-width'=>'1','stroke-opacity'=>'1','fill-opacity'=>'1'} );
$y = $svg->group( id=>'group_y', style=>{'stroke'=>'black','stroke-width'=>'0.5','stroke-opacity'=>'0.5','fill-opacity'=>'0.5'} );
$z = $svg->group( id=>'group_z', style=>{'stroke'=>'black','stroke-width'=>'2','stroke-opacity'=>'1','fill-opacity'=>'1'} );

$x->line( x1=>150, y1=>50, x2=>($width-50), y2=>50 );
$sw = ($width - 200) / 10;
$value = 0.1;
for($k=0; $k <=10; $k++){
  $x->line( x1=>(150+$sw*$k), y1=>45, x2=>(150+$sw*$k), y2=>50 );
  $x->text(x=>(150+$sw*$k), y=>25, 'font-size'=>10, 'stroke', 'black', '-cdata', $value*$k );
}


$Vgap=50; $Vnick=2; $uHigh=10; $angle=0.5;
$Hnick=10; $max_total_len=0;

$i=1;


$loc1=$loc2=0;
@code_array=("FIR","SEC","TRD","FOR","FIF","SIX","SEV","EIG","NIN","TEN","ELE","TWE");
for $mykey (@code_array){

if (defined $rev_inhash{$mykey}){

  if ($loc2 > 0){
    $loc1=$loc2+$Vgap;  
  }else{
    $loc1=50+$Vgap;   
  }
  $mylen=$faLength{$mykey};
  $head=$mykey." ".$mylen; $head_add="(".$rev_inhash{$mykey}.")";
  $x->text(x=>50, y=>$loc1, 'font-size'=>15, 'stroke', 'black', '-cdata', $head );
  $x->text(x=>($width-225), y=>$loc1, 'font-size'=>15, 'stroke', 'black', '-cdata', $head_add );

  $loc2=$loc1+($Vnick+$uHigh);
  $y->line(x1=>150, y1=>$loc2, x2=>150+$mylen*$ratio, y2=>$loc2 );

  %temp_domain=%{$sp2domain{$mykey}};   #C1 3_10; C2 8_20; 
  for (0..$#domain_maxlen_key_arr){
   $up_or_down=$_;
   $temp_domain_ele=$domain_maxlen_key_arr[$_];    #C1; C2;
   if (defined $temp_domain{$temp_domain_ele}){    #
  
      @split=split("_", $temp_domain{$temp_domain_ele});
      if (($up_or_down%2)==0){
#        $svg->rect(x => 150+$split[0]*$ratio, y => $loc1+0.5*$Vnick, width => ($split[1]-$split[0]+1)*$ratio, height => $uHigh, style=>{'fill'=>$domaincolor{$temp_domain_ele},'stroke'=>'grey','stroke-width'=>'0.5','stroke-opacity'=>'1','fill-opacity'=>'1'} );
#        $svg->text(x=>150+$split[0]*$ratio, y=>$loc1+5*$Vnick, 'font-size'=>10, 'stroke', 'grey', 'stroke-width', '0.3', '-cdata', $temp_domain_ele);
              $_xv = [ 150+$split[0]*$ratio, 150+$split[0]*$ratio, 150+$split[0]*$ratio+($split[1]-$split[0]+1)*$ratio ];
              $_yv = [ $loc2, $loc2-$uHigh, $loc2-0.5*$uHigh ];
              $_points = $svg->get_path( x=>$_xv, y=>$_yv, -type=>'polygon' );
              $svg->polygon(%$_points,style=>{'fill'=>$domaincolor{$temp_domain_ele}, 'fill-opacity'=>'0.5', 'stroke'=>'grey', 'stroke-width'=>'0.5', 'stroke-opacity'=>'0.5'});
        $z->line(x1=>150+$split[0]*$ratio, y1=>$loc2, x2=>150+$split[0]*$ratio, y2=>$loc2-$uHigh );
        $svg->text(x=>150+$split[0]*$ratio, y=>$loc1, 'font-size'=>10, 'stroke', 'grey', 'stroke-width', '0.3', '-cdata', $temp_domain_ele);
      }elsif (($up_or_down%2)!=0){
#        $svg->rect(x => 150+$split[0]*$ratio, y => $loc2+0.5*$Vnick, width => ($split[1]-$split[0]+1)*$ratio, height => $uHigh, style=>{'fill'=>$domaincolor{$temp_domain_ele},'stroke'=>'grey','stroke-width'=>'0.5','stroke-opacity'=>'1','fill-opacity'=>'1'} );
#        $svg->text(x=>150+$split[0]*$ratio, y=>$loc1+5*$Vnick, 'font-size'=>10, 'stroke', 'grey', 'stroke-width', '0.3', '-cdata', $temp_domain_ele);
              $_xv = [ 150+$split[0]*$ratio, 150+$split[0]*$ratio, 150+$split[0]*$ratio+($split[1]-$split[0]+1)*$ratio ];
              $_yv = [ $loc2, $loc2+$uHigh, $loc2+0.5*$uHigh ];
              $_points = $svg->get_path( x=>$_xv, y=>$_yv, -type=>'polygon' );
              $svg->polygon(%$_points,style=>{'fill'=>$domaincolor{$temp_domain_ele}, 'fill-opacity'=>'0.5', 'stroke'=>'grey', 'stroke-width'=>'0.5', 'stroke-opacity'=>'0.5'});
        $z->line(x1=>150+$split[0]*$ratio, y1=>$loc2, x2=>150+$split[0]*$ratio, y2=>$loc2+$uHigh );
        $svg->text(x=>150+$split[0]*$ratio, y=>$loc2+2*$uHigh, 'font-size'=>10, 'stroke', 'grey', 'stroke-width', '0.3', '-cdata', $temp_domain_ele);
      }

   }
  }

} #if

} #@code_array


### motif panel data preprocessing
foreach $Dkey (keys %domain2motif){    ### D1=>{FIR=>@{fimo lines}; SEC=>@{fimo lines}; }
  %D_domains=%{$domain2motif{$Dkey}};
  foreach $sp_key (keys %D_domains){  ### FIR=>@{fimo lines}; SEC=>@{fimo lines}; 
    @sp_motifs=@{$D_domains{$sp_key}};    
    foreach $motif_ele (@sp_motifs){  ### @{fimo lines}
      @motif_ele=split("\t", $motif_ele);
      push @{$domain_motif_kinds{$Dkey}}, $motif_ele[0];
      push @{$sp_domain_motif_s{$sp_key}{$Dkey}}, $motif_ele[2];
      push @{$sp_domain_motif_e{$sp_key}{$Dkey}}, $motif_ele[3];
    }
  }
}

foreach $temp_domain_key (keys %domain_motif_kinds){
  @temp_domain_motif=@{$domain_motif_kinds{$temp_domain_key}};
  $domain_motif_kinds_uniq{$temp_domain_key}=\@{[uniq_array(\@temp_domain_motif)]}; ### ??? ref of array
}

foreach $temp_sp (keys %sp_domain_motif_s){
  %temp_domain_s=%{$sp_domain_motif_s{$temp_sp}};
  %temp_domain_e=%{$sp_domain_motif_e{$temp_sp}};
  foreach $temp_domain (keys %temp_domain_s){
    @domain_s=@{$temp_domain_s{$temp_domain}};
    @domain_e=@{$temp_domain_e{$temp_domain}};
    push @{$sp_max_domain_len{$temp_sp}}, max(@domain_e)-min(@domain_s)+1;
    $sp_ind_domain_len{$temp_sp}{$temp_domain}=max(@domain_e)-min(@domain_s)+1;
  }
}
foreach $sp2maxlen_k (keys %sp_max_domain_len){
  $sp2maxlen{$sp2maxlen_k}=max (@{$sp_max_domain_len{$sp2maxlen_k}});
  $max_total_len+=$Hnick+max(@{$sp_max_domain_len{$sp2maxlen_k}});
}


$loc1=$loc2+$Vgap; $_start=150; $m=0;
$loc3=$loc1; 
$_ratio=($width - 200)/$max_total_len;
for (@code_array){                     ### sort species within down panels
if (defined $sp2maxlen{$_}){
  $_sp=$_; $_sp_maxD_len=$sp2maxlen{$_}; $_header=$_sp." ".$_sp_maxD_len;           ### display species domain length
  $m++;
  $loc2=$loc3;
  $sp_bench_s=$_start;
  $sp_len=$sp2maxlen{$_sp};
  $R_sp_len=$sp_len*$_ratio;
  $x->line( x1=>$sp_bench_s, y1=>$loc3, x2=>$sp_bench_s+$R_sp_len, y2=>$loc3 );
  $svg->text( x=>$sp_bench_s, y=>$loc3+25, 'font-size'=>15, 'stroke', 'black', '-cdata', $_header );
  $sw = $R_sp_len / 10;
  for($k=0; $k <=10; $k++){
    $x->line( x1=>($sp_bench_s+$sw*$k), y1=>$loc3-2.5, x2=>($sp_bench_s+$sw*$k), y2=>$loc3 );
  }

  for (0..$#domain_maxlen_key_arr){
    $loc1=$loc2+$Vgap;
    $temp_domain_ele=$domain_maxlen_key_arr[$_];
    @motifs=@{$domain2motif{$temp_domain_ele}{$_sp}};
    @uniq=@{$domain_motif_kinds_uniq{$temp_domain_ele}};
    $domain_len=$sp_ind_domain_len{$_sp}{$temp_domain_ele};
    $k=0;
    for $col (@uniq){
      $_col=$color{$col};
      $loc2=$loc1+$k*($Vnick+2*$uHigh)+($Vnick+$uHigh);
      if ($m==1){
      $svg->text(x=>100, y=>$loc2+5, 'font-size'=>10, 'stroke', 'black', '-cdata', $col );
      $svg->text(x=>50, y=>$loc2+5, 'font-size'=>10, 'stroke', 'black', '-cdata', $temp_domain_ele );
      }
      $y->line(x1=>$sp_bench_s, y1=>$loc2, x2=>$sp_bench_s+$domain_len*$_ratio, y2=>$loc2 );

      foreach (0..$#motifs){    #motifs  coordinates transformation
      $ele=$motifs[$_];
      @Fs=split(/\t/, $ele);
      if ($_==0){ $origin_s=$Fs[2]; }
      if ($Fs[0] eq $col){
        $st=($Fs[2]-$origin_s)*$_ratio;
        $ed=($Fs[3]-$origin_s)*$_ratio;
        $recLen=$ed-$st+1;
        if ($Fs[4]=~/\+/){
           $xv = [ $sp_bench_s+$st, $sp_bench_s+$st+$recLen*$angle, $sp_bench_s+$ed, $sp_bench_s+$st+$recLen*$angle, $sp_bench_s+$st];
           $yv = [ $loc2-$uHigh, $loc2-$uHigh, $loc2, $loc2+$uHigh, $loc2+$uHigh];
           $points = $svg->get_path( x=>$xv, y=>$yv, -type=>'polygon' );
           $svg->polygon(%$points,style=>{'fill'=>$_col, 'stroke'=>'grey', 'stroke-width'=>'0.5'});
        }elsif ($Fs[4]=~/-/){
           $xv = [ $sp_bench_s+$st, $sp_bench_s+$st+$recLen*$angle, $sp_bench_s+$ed, $sp_bench_s+$ed, $sp_bench_s+$st+$recLen*$angle];
           $yv = [ $loc2, $loc2-$uHigh, $loc2-$uHigh, $loc2+$uHigh, $loc2+$uHigh];
           $points = $svg->get_path( x=>$xv, y=>$yv, -type=>'polygon' );
           $svg->polygon(%$points,style=>{'fill'=>$_col, 'stroke'=>'grey', 'stroke-width'=>'0.5'});
        }
      }
      }
    $k=$k+1;
    }

  }

  $_start+=($sp_len+$Hnick)*$_ratio;

}
}

print $FH $svg->xmlify;

}


sub trd_dp_Highest{

$trd_dp=shift;
$_tempdpID=shift;
$_hashHighscore=shift; %_hashHighscore=%{$_hashHighscore};
$_hashHigh=shift; %_hashHigh=%{$_hashHigh};
$num=0; $templines="";
open DPH, "<$trd_dp";
while (<DPH>){
if ($num<=1){
if(/>>>/){
  $num++;
#  $templines.=$_;
}elsif(/>>Score(\d+\.\d+)/){
  $tempkey=$1;
  $templines.=$_;
}else{
  $templines.=$_;
}
}
}

$_hashHighscore{$_tempdpID}=$tempkey;
$_hashHigh{$_tempdpID}=$templines;

return (\%_hashHighscore, \%_hashHigh);

}
