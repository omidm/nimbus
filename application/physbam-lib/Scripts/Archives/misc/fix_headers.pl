#!/usr/bin/perl -w

require strict;

my %lookup;

open F, shift @ARGV;

while(<F>){
    chomp;
    /.*\/(.*)/ or die "Failed to parse '$_' as a valid header.\n";
    (defined $lookup{$1}) and die "Duplicate file '$1' encountered\n.";
    $lookup{$1} = $_;
}
close F;

my @tools = ();
my @geometry = ();
my @solids = ();
my @fluids = ();
my @opengl = ();
my @physbam = ();
my @system = ();
my @other = ();
my %system = ("algorithm" => 1, "streambuf" => 1, "iomanip" => 1, "list" => 1, "ostream" => 1,
              "bitset" => 1, "ios" => 1, "locale" => 1, "queue" => 1, "string" => 1, 
              "complex" => 1, "typeinfo" => 1, "iosfwd" => 1, "map" => 1, "set" => 1, 
              "deque" => 1, "iostream" => 1, "memory" => 1, "sstream" => 1, "utility" => 1, 
              "exception" => 1, "valarray" => 1, "istream" => 1, "new" => 1, "stack" => 1, 
              "fstream" => 1, "iterator" => 1, "numeric" => 1, "stdexcept" => 1, "vector" => 1, 
              "functional" => 1, "limits" => 1, "cassert" => 1, "ciso646" => 1, "csetjmp" => 1, 
              "cstdio" => 1, "ctime" => 1, "cctype" => 1, "climits" => 1, "csignal" => 1, 
              "cstdlib" => 1, "cwchar" => 1, "cerrno" => 1, "clocale" => 1, "cstdarg" => 1, 
              "cstring" => 1, "cwctype" => 1, "cfloat" => 1, "cmath" => 1, "cstddef" => 1);

my $cf = "";

my $new;
my $filename;

sub dump_headers
{
    if(!@tools && !@geometry && !@solids && !@fluids && !@opengl && !@physbam && !@system && !@other){return;}
    $filename eq $cf or die "Trying to insert includes from $cf into $filename.\n";
    my $last = "";
    for $L (\@tools, \@geometry, \@solids, \@fluids, \@opengl, \@physbam, \@system, \@other)
    {
        for(sort {$x=lc $$a[0];$y=lc $$b[0];$x=~s/_/0/g;$y=~s/_/0/g;$x cmp $y;} @$L)
        {
            if($$_[0] eq $last){print STDERR "DUPLICATE HEADER: $last in $filename\n";next;}
            $last = $$_[0];
            my $end_comment=$$_[2]?" $$_[2]":"";
            $new.="$$_[1]#include $$_[0]$end_comment\n";
        }
    }
    @tools = ();
    @geometry = ();
    @solids = ();
    @fluids = ();
    @opengl = ();
    @physbam = ();
    @system = ();
    @other = ();
    $cf = "";
}

for $filenamex (@ARGV)
{
    $filename=$filenamex;
    my $orig='';
    $new='';

    open F, $filename;
    while(<F>)
    {
        $orig.=$_;
        if(!m@^\s*(//|)\s*\#\s*include\s*(([<\"])(.*/)?(.*)([>\"]))\s*(//.*\S|)@){&dump_headers();$new.=$_;next;}
        if($cf eq ""){$cf = $filename;}
        $filename eq $cf or die "Trying to mix includes from $cf into $filename.\n";
        my $ldelim = $3;
        my $rdelim = $6;
        my $full = $2;
        my $path = $4;
        my $file = $5;
        my $pre = $1;
        my $post = $7;
        if(!$path && defined $system{$file}){push @system, [$full,$pre,$post];next;}
        if(!$path){push @other, [$full,$pre,$post];next;}
        if(defined $system{$file}){print STDERR "system file specified with path: '$path$file' in '$filename'\n";push @system, [$full,$pre,$post];next;}
        if(!defined $lookup{$file}){print STDERR "include not recognized: '$file' in '$filename'\n";push @other, [$full,$pre,$post];next;}
        my $header = "<$lookup{$file}>";
        if($header =~ /PhysBAM_Tools/){push @tools, [$header,$pre,$post];next;}
        if($header =~ /PhysBAM_Geometry/){push @geometry, [$header,$pre,$post];next;}
        if($header =~ /PhysBAM_Solids/){push @solids, [$header,$pre,$post];next;}
        if($header =~ /PhysBAM_Fluids/){push @fluids, [$header,$pre,$post];next;}
        if($header =~ /PhysBAM_OpenGL/){push @opengl, [$header,$pre,$post];next;}
        push @physbam, [$header,$pre,$post];
    }
    &dump_headers();
    close F;

    if($new ne $orig)
    {
        open F, ">$filename";
        print F $new;
        close F;
    }
}

