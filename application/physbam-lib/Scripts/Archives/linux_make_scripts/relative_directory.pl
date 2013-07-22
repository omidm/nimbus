#!/usr/bin/perl

use Cwd 'abs_path';

my $debug = 0;

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

if ($#ARGV < 1) {
	print "Usage: $0 <dir> <targetpath>\n";
}

print relative_directory($ARGV[0], $ARGV[1]),"\n";
