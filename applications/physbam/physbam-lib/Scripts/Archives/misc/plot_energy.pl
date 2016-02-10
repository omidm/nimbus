#!/usr/bin/perl -w

use Getopt::Std;
require strict;

my %h=();
my $ok=getopts('l',\%h);

my $data_tmp="plot_data.txt";

if(@ARGV<=2){print "Usage: $0 <log-file> <out.eps> [<step>] [opts]\n";exit;}

my ($in, $out, $step) = @ARGV;

open O, ">$data_tmp";
open F, $in;
while(<F>)
{
    /total energy = ([-+.0-9eE]+).*Step $step(?:[^0-9]|$)/ && print O "$1\n";
}
close F;
close O;

my $logplot='';
if($h{'l'}){$logplot=" u (log10(\$1))"}

`echo 'set terminal postscript enhanced eps color ; set output "$out" ; plot "$data_tmp" $logplot ' | gnuplot`;

`evince $out`;
