#!/usr/bin/perl

use File::Spec;

my $debug = 0;

sub collapse_path
{
	my ($path) = @_;

	chomp($path);
	@paths = split('/', $path);

	@new_paths = ();
	foreach $item (@paths)
	{
		if ($#new_paths > 0 and $item eq "..")
		{
			$last_item = pop(@new_paths);
			if ($last_item eq "..") 
			{ 
				push(@new_paths, $last_item); 
			}
		}
		else
		{
			push(@new_paths, $item);
		}
	}

	$new_path = join('/', @new_paths);
	return ($new_path);
}

if (!$ARGV[0])
{
	print ("Usage: $0 <main .cpp or executable> ...\n");
	exit;
}

my @cppfiles = ();
my %main_source = ();
foreach $file (@ARGV)
{
    if (not $file=~/\.cpp$/)
    {
	    $file = $file.".cpp";
    }

    if (not -e $file)
    {
	    print ("Error: $file does not exist\n");
	    exit;
    }

    $file = File::Spec->canonpath($file);
	push (@cppfiles, $file);
	$main_source{$file} = 1;
}

my %cppfiles_done = {};
my %hfiles_done = {};
my @maincppfiles = ();

my $cc_depend = "g++";
if ($ENV{PHYSBAM_CC_DEPEND}) { $cc_depend = $ENV{PHYSBAM_CC_DEPEND}; }

while ($#cppfiles >= 0)
{
	$cppfile = pop @cppfiles;

	if (not defined $cppfiles_done{$cppfile})
	{
		$cppfiles_done{$cppfile} = 1;
		if (not defined $main_source{$cppfile}) {
			push (@maincppfiles, $cppfile);
		}

		if ($debug) { print "Processing cpp file $cppfiles\n"; }
		open (INPUT, "$cc_depend -MM $cppfile |");

		while ($line = <INPUT>)
		{
			@hfiles = grep(/\S*\.h/, split(' ', $line));
			foreach $hfile (@hfiles)
			{
				($hfile) = collapse_path($hfile);
				if (not defined $hfile_done{$hfile})
				{
					if ($debug) { print "Adding header file $hfile\n"; }
					$hfiles_done{$hfile} = 1;

					$cppfile = $hfile;
					$cppfile =~ s/\.h/\.cpp/;
					if (-e $cppfile)
					{
						push (@cppfiles, $cppfile);
					}
				}
			}
		}
	}
}

print (join("\n", @maincppfiles));
