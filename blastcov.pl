#!/usr/bin/perl -w

use strict;
use Bio::SearchIO;

my $in = new Bio::SearchIO(-format => 'blast',
                           -file => "$ARGV[0]");
while( my $result = $in->next_result ) {
  while( my $hit = $result->next_hit ) {
    while( my $hsp = $hit->next_hsp ) {
      my $MismatchesNum    = $hsp->hsp_length-($hsp->gaps+$hsp->num_identical);
      my $QCoveragePercent = sprintf("%.2f",(($hsp->hsp_length-$MismatchesNum-$hsp->gaps)/$result->query_length)*100);
      my $SCoveragePercent = sprintf("%.2f",(($hsp->hsp_length-$MismatchesNum-$hsp->gaps)/$hit->length)*100);
      my $Percent_Identity = sprintf("%.2f",$hsp->percent_identity);

      print   $result->query_name,"\t",
              $hit->name,"\t",
              $Percent_Identity,"\t",
              $hsp->hsp_length,"\t",
              $MismatchesNum,"\t",
              $hsp->gaps,"\t",
              $hsp->start('query'),"\t",
              $hsp->end('query'),"\t",
              $hsp->start('hit'),"\t",
              $hsp->end('hit'),"\t",
              $hsp->evalue,"\t",
              $hsp->bits,"\t",
              "$QCoveragePercent\t",
              "$SCoveragePercent\n";
    }
  }
}

exit;
