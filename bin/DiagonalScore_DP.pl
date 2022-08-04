#!/usr/bin/perl -w

use experimental qw(smartmatch);

$total_align=shift;
$outfile=shift;

open IN, $total_align;
while(<IN>){
     next if (/^#/);
     chomp;
     if (/^\//){
        $line=$_;
        @ln=split("\t", $line);
        $cluster=$ln[0]; $cluster=~s/\S+\/(C\d+)_\S+/$1/;
        #$dp_score=$ln[2];
        #$euc_score=$ln[3];
        #$fir_sted=$ln[4];
        #$sec_sted=$ln[5];
        $fir_clu_id=$cluster."FIR"; 
        while (length($fir_clu_id)!=8){
            $fir_clu_id="0".$fir_clu_id;
        }
        $sec_clu_id=$cluster."SEC";
        while (length($sec_clu_id)!=8){
            $sec_clu_id="0".$sec_clu_id;
        }
        $match_score{$fir_clu_id."_".$sec_clu_id}=$ln[2]; 
        $fir_sted{$fir_clu_id}=${[split("_",$ln[4])]}[0];
        $sec_sted{$sec_clu_id}=${[split("_",$ln[5])]}[0];
     }
}
close IN;

%micro_matrix=(); $micro_matrix=\%micro_matrix;
foreach $k1 ( sort { $fir_sted{$a} <=> $fir_sted{$b} } keys %fir_sted ){
     foreach $k2 ( sort { $sec_sted{$a} <=> $sec_sted{$b} } keys %sec_sted ){
        if (defined $match_score{$k1."_".$k2}){
           $micro_matrix->{$k1}{$k2}=$match_score{$k1."_".$k2};
        }else{
           $micro_matrix->{$k1}{$k2}=-1000;
        }
     }
}

@fir_seq=( sort { $fir_sted{$a} <=> $fir_sted{$b} } keys %fir_sted ); $fir_seq=join("", @fir_seq);
@sec_seq=( sort { $sec_sted{$a} <=> $sec_sted{$b} } keys %sec_sted ); $sec_seq=join("", @sec_seq);

open OUT, ">$outfile";

#=head
@axis=();
for $key1 (@fir_seq){
  %temp=%{$micro_matrix{$key1}};
  print OUT $key1,"\t";
  for $key2 (@sec_seq){
    print OUT $temp{$key2},"\t";
    push @axis, $key2 if (!($key2~~@axis));
  }
  print OUT "\n";
}
print OUT "\t",join("\t",@axis),"\n\n";
#=cut

$out_result = DP($fir_seq, $sec_seq, 1);

@out_result=split("\n", $out_result);

#print $out_result;
#for (0..$#out_result){
#    print $_," : ",$out_result[$_],"\n";
#}

#=head
print OUT $out_result[0],"\n",$out_result[1],"\n",$out_result[2],"\n\n";
@fir=$out_result[1]=~/.{8}/g;
@sec=$out_result[2]=~/.{8}/g;
for (0..$#fir){
    if (($fir[$_]!~/-/)&&($sec[$_]!~/-/)){
       $fc=$fir[$_]; $fc=~s/.+(C\d+)FIR/$1/;
       $sc=$sec[$_]; $sc=~s/.+(C\d+)SEC/$1/;
       print OUT "$fc\n" if ($fc eq $sc);
    }
}
#=cut



sub DP{

  $seq1=shift; $seq2=shift; $sigth=shift; $sig = 1;

  $GAP = 0; $MATCH = 0; @matrix = ();%cmp_hash = ();
  $max_i = $max_j = 0; @pool = @align1 = @align2 = ();
  $result = "";


  $matrix[0][0]{score}   = 0;
  $matrix[0][0]{pointer} = "none";
  for($j = 1; $j <= length($seq1)/8; $j++) {
    $matrix[0][$j]{score}   = 0;
    $matrix[0][$j]{pointer} = "none";
  }
  for($i = 1; $i <= length($seq2)/8; $i++) {
    $matrix[$i][0]{score}   = 0;
    $matrix[$i][0]{pointer} = "none";
  }

  for($i = 1; $i <= length($seq2)/8; $i++) {
    for($j = 1; $j <= length($seq1)/8; $j++) {
       ($diagonal_score, $left_score, $up_score)=(0,0,0);

        $letter1 = substr($seq1, ($j-1)*8, 8);
        $letter2 = substr($seq2, ($i-1)*8, 8);
           $MATCH=$micro_matrix->{$letter1}{$letter2};
           $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;

        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;

        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next;
        }

        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                if ($diagonal_score >0){
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
                }
            }
            else {
                if ($left_score >0){
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
                }
            }
        } else {
            if ($up_score >= $left_score) {
                if ($up_score >0){
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
                }
            }
            else {
                if ($left_score >0){
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
                }
            }
        }

        $cmp_hash{$i."_".$j} = $matrix[$i][$j]{score};

    }
  }



  foreach $_key_ (sort { $cmp_hash{$b} <=> $cmp_hash{$a} } keys %cmp_hash){
  
  #if ($sig<=$sigth){
    unless ($_key_~~ @pool){

      ($max_i, $max_j)=split("_", $_key_);
      if($matrix[$max_i][$max_j]{pointer} eq "diagonal"){
        $j = $max_j; $i = $max_i;

        while (1) {
          last if $matrix[$i][$j]{pointer} eq "none";
          push @pool, $i."_".$j;
          if ($matrix[$i][$j]{pointer} eq "diagonal") {
            push @align1, substr($seq1, ($j-1)*8, 8);
            push @align2, substr($seq2, ($i-1)*8, 8);
            $i--; $j--;
          }elsif ($matrix[$i][$j]{pointer} eq "left") {
            push @align1, substr($seq1, ($j-1)*8, 8);
            push @align2, "--------";
            $j--;
          }elsif ($matrix[$i][$j]{pointer} eq "up") {
            push @align1, "--------";
            push @align2, substr($seq2, ($i-1)*8, 8);
            $i--;
          }   
        }

        @align1 = reverse @align1;
        @align2 = reverse @align2;

        $result .= ">>>Score".$cmp_hash{$_key_}."\n".join("",@align1)."\n".join("",@align2)."\n\n";

        undef @align1;
        undef @align2;
      }
    }

  #}
  #$sig+=1;

  }

return $result;

}
