#!/usr/bin/env perl

package setenv;

use strict;
use setenv_local;

# simstudy file format is:
# <model condition, with all params separated by _>/<R? run number>/<any method, alignest or treeest or align+treeest, with all params separated by _>
# model condition run dir has STATS file and so does method dir
# true alignment file stored at rundir

use constant BASH_COMMAND => setenv_local::BASH_COMMAND;
use constant PYTHON_COMMAND => setenv_local::PYTHON_COMMAND;

# perl availability
use constant PERL_HOME => setenv_local::PERL_HOME;

# only a few fixed files in structure
use constant STATS_FILENAME => "STATS";

# for now not used?
# state types
use constant AADATATYPE => "AA";
use constant DNADATATYPE => "DNA";

use constant PERL_COMMAND => setenv_local::PERL_COMMAND;

sub setenv {


}

1;
