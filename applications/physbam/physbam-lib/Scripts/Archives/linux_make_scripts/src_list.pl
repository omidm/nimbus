#!/usr/bin/perl

use Cwd 'getcwd';
use Cwd 'abs_path';

my @src_dirs = ( LOCAL_SRC, PUBLIC_LIBRARY_SRC, PERSONAL_LIBRARIES_SRC, ERAN_LIBRARY_SRC );
my %objects = {};

sub relative_directory
{
	my ($dir, $file) = @_;

	$dir = abs_path($dir);
	$file = abs_path($file);

	@dirparts = split('/', $dir);
	@fileparts = split('/', $file);

	while ((defined $dirparts[0]) && (defined $fileparts[0]) && ($dirparts[0] eq $fileparts[0]))
	{
		shift @dirparts;
		shift @fileparts;
	}

	$reldir = "";
	for (my $i = 0; $i <= $#dirparts; $i++) {
		$reldir = $reldir . "../";
	}

	$reldir = $reldir . join('/', @fileparts);

	if ($debug) {
		print "Matched at\n";
		print join('/', @dirparts),"\n";
		print join('/', @fileparts),"\n";
		print "RELDIR = $reldir\n";
	}

	return $reldir;
}

if ($ARGV[0])
	{ open (INPUT, "< $ARGV[0]"); }
else
	{ open (INPUT, "-"); } 

while ($line = <INPUT>) {
	if ($line=~/Public_Library\/(\S*).cpp/) {
		push (@{$objects{PUBLIC_LIBRARY_SRC}}, $1.".cpp");
	}
	elsif ($line=~/Personal_Libraries\/Eran_Library\/(\S*).cpp/) {
		push (@{$objects{ERAN_LIBRARY_SRC}}, $1.".cpp");
	}
	elsif ($line=~/Personal_Libraries\/(\S*).cpp/) {
		push (@{$objects{PERSONAL_LIBRARIES_SRC}}, $1.".cpp");
	}
	elsif ($line=~/(\S*).cpp/) {
		$relpath = relative_directory(getcwd(), $1);
		push (@{$objects{LOCAL_SRC}}, $relpath.".cpp");
	}
}

my $list = $objects{LOCAL_SRC};
if ($#{$list} >= 0) {
    print $src_dirs[$i]," = \\\n\t";
    my @sorted_list = sort (@{$list});
    print join(" \\\n\t", @sorted_list);
    print "\n\n";
}

my $libs = "";
if ($#{$objects{PUBLIC_LIBRARY_SRC}} >= 0) {
    $libs = "$libs PhysBAM"
}
if ($#{$objects{ERAN_LIBRARY_SRC}} >= 0) {
    $libs = "$libs Eran_Library"
}

if ($libs ne "") {
    print "LIBRARIES =$libs\n\n";
}
