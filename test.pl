#!/usr/bin/env perl

use strict;
use warnings;

my $setPath = "PATH=/bin:/usr/bin; ";
my $cmd = "ls";

my $which = qx{$setPath type -ap $cmd 2> /dev/null | head -1};
chomp $which;

printf STDERR "# which: '%s'\n", $which;
