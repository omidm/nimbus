#!/usr/bin/perl -w

use strict;

my %adj = ();

for(`find \$PUBLIC -type f -name '*.h' -o -name '*.cpp'`)
{
    chomp;
    open F, $_;
    /.*\/(.*)/ || die "Failed to parse filename $_\n";
    my $file = $1;
    my @inc = ();

    while(<F>)
    {
        if(/\#\s*include\s*[<\"][^<>\"]*\/([^\/<>\"]*?)[\">]/)
        {
            push @inc, $1;
        }
    }
    $adj{$file} = [@inc];

    close F;
}

my $i = 0;
my @S = ();
my %vi = ();
my %vlow = ();
for my $v (keys %adj)
{
    if(!defined $vi{$v})
    {
        &tarjan($v);
    }
}

my @comp;
my %hs = ();
sub tarjan
{
    my $v = $_[0];
    $vi{$v} = $i;
    $vlow{$v} = $i++;
    push @S, $v;
    $hs{$v}=1;
    for my $n (@{$adj{$v}})
    {
        if(!defined $vi{$n})
        {
            &tarjan($n);
            if($vlow{$n} < $vlow{$v})
            {
                $vlow{$v} = $vlow{$n};
            }
        }
        elsif(defined $hs{$n})
        {
            if($vi{$n} < $vlow{$v})
            {
                $vlow{$v} = $vi{$n};
            }
        }
    }
    if($vlow{$v} == $vi{$v})
    {
        my @c = ();
        my $n;
        do
        {
            $n = pop @S;
            delete $hs{$n};
            push @c, $n;
        }
        while($n ne $v);
        push @comp, [@c];
    }
}

my $id = 0;
my $ic = 0;
for(@comp)
{
    my @c = @$_;
    if(@c <= 1)
    {
        next;
    }
    if(@c == 2 && $c[0] eq "READ_WRITE_$c[1]")
    {
        print "rw: $c[1]\n";
        next;
    }
    if(@c == 2 && $c[1] eq "READ_WRITE_$c[0]")
    {
        print "rw: $c[0]\n";
        next;
    }

    my %h = ();
    open F, ">component-$ic\n";
    print F "digraph dg {\n";
    for(@c)
    {
        $h{$_}=$id++;
        print F "$h{$_} [label=\"$_\"];\n";
    }
    for $a (@c)
    {
        for $b (@{$adj{$a}})
        {
            if(defined $h{$b})
            {
                print F "$h{$a} -> $h{$b};\n";
            }
        }
    }
    print F "}\n";
    close F;
    `dot component-$ic -Tps > component-$ic.ps`;
    $ic++;
}
