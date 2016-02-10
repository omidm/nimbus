#!/usr/bin/perl

# USAGE:
# gen-image.pl [-bcekpsv] [-t temp-prefix] [-h header-LaTeX-code] [-g latex-flags] [-S scale] output-filename LaTeX to render here
#
# -b       The resulting image will be big.  Try to make it work.
# -c       Print the commands run.
# -e       Do not specify an empty page style.  This will leave the page number on the document.
# -g flag  Pass these commandline flags to latex.
# -h code  Add code to the LaTeX document before \begin{document}; useful to include packages.
# -k       Keep intermediate files.
# -p       Pretend.  Useful with -v to see what will happen.
# -s       Use the font for the SIGGRAPH format.  Default is LaTeX's default Computer Modern.
# -S num   Scale the image by a factor of num.  May need to specify -b as well.
# -t tmp   Use tmp as a prefix for temporary files.  Default is a random filename in the current directory.
# -v       Be verbose.  This prints out the output of commands that are run.

use Getopt::Std;
getopts("bceg:h:kpsS:t:v");

if(@ARGV < 2) {die "USAGE:\n    gen-image.pl [-bcekpsv] [-t temp-prefix] [-h header-LaTeX-code] [-g latex-flags] [-S scale] output-filename LaTeX to render here\n";}

my $tmp_prefix;
if($opt_t) {$tmp_prefix=$opt_t;}
else {$tmp_prefix="tmp-filename-" . int(rand(1000000));}

my $output_file = shift @ARGV;
my $body=join "\n", @ARGV;

my $dir_option='';
if($tmp_prefix=~/(.*)\//) {$dir_option="-output-directory=\"$1\"";}

my $format;
if($output_file=~/.*\.svg$/i) {$format='svg';}
elsif($output_file=~/.*\.sk$/i) {$format='sk';}
elsif($output_file=~/.*\.eps$/i) {$format='eps';}
elsif($output_file=~/.*\.tex$/i) {$format='tex';}
elsif($output_file=~/.*\.pdf$/i) {$format='pdf';}
elsif($output_file=~/.*\.png$/i) {$format='png';}
else {die "filename format:   *.{svg,sk,eps,tex,pdf,png}\n";}

open F, ">$tmp_prefix.tex" || die "Cannot open file \"$tmp_prefix.tex\"";
print F "\\documentclass{article}\n";
if($opt_S) {print F "\\usepackage{epsfig}\n";}
if($opt_b) {print F "\\usepackage[papersize={16383pt,16383pt},total={16383pt,16383pt}]{geometry}\n";}
if($opt_s) {print F "\\usepackage{mathptm}\n";}
if($opt_h) {print F "$opt_h\n";}
print F "\\begin{document}\n";
if(!$opt_e) {print F "\\pagestyle{empty}\n";}
if($opt_S) {print F "\\scalebox{$opt_S}{\n";}
print F "$body\n";
if($opt_S) {print F "}\n";}
print F "\\end{document}\n";

my $cmd;
my $rmlist = "\"$tmp_prefix.tex\"";
if($format eq 'pdf' || $format eq 'svg' || $format eq 'sk' || $format eq 'eps' || $format eq 'png')
{
    $cmd = "pdflatex -halt-on-error $dir_option $opt_g \"$tmp_prefix.tex\"";
    $rmlist .= " \"$tmp_prefix.pdf\" \"$tmp_prefix.aux\" \"$tmp_prefix.log\"";
}
if($format eq 'svg' || $format eq 'sk' || $format eq 'eps' || $format eq 'png')
{
    $cmd .= " && pdf2ps \"$tmp_prefix.pdf\" \"$tmp_prefix.eps\" 2>&1";
    $rmlist .= " \"$tmp_prefix.eps\"";
}
if($format eq 'svg' || $format eq 'sk')
{
    $cmd .= " && pstoedit -page 1 -dt -psarg '-r9600x9600' -f sk \"$tmp_prefix.eps\" \"$tmp_prefix.sk\" 2>&1";
    $rmlist .= " \"$tmp_prefix.sk\"";
}
if($format eq 'svg')
{
    $cmd .= " && skconvert \"$tmp_prefix.sk\" \"$tmp_prefix.svg\" 2>&1";
    $rmlist .= " \"$tmp_prefix.svg\"";
}
if($format eq 'png')
{
    $cmd .= " && convert \"$tmp_prefix.eps\" \"$tmp_prefix.png\" 2>&1";
    $rmlist .= " \"$tmp_prefix.png\"";
}

my $out;
if($opt_c) {print "$cmd\n";}
if(!$opt_p)
{
    $out = `$cmd`;
    if($opt_v) {print $out;}
}
if($opt_c) {print "mv \"$tmp_prefix.$format\" \"$output_file\" 2>&1\n";}
if(!$opt_p)
{
    $out = `mv \"$tmp_prefix.$format\" \"$output_file\" 2>&1`;
    if($opt_v) {print $out;}
}
if(!$opt_k)
{
    if($opt_c) {print "rm $rmlist 2>&1\n";}
    if(!$opt_p)
    {
        $out = `rm $rmlist 2>&1`;
        if($opt_v) {print $out;}
    }
}



