#!/usr/bin/env perl


# nice this works
package setenv_local;

use strict;

# need this flag to more intelligently choose scripts/programs
# automatically
# for now, just 32 or 64 bit - can change this later as need be
use constant ARCHITECTURE_64_BIT_FLAG => 0;


# perl availability
# need this too - due to older perl interpreter only on TACC
use constant PERL_HOME_32_BIT => "/lusr/bin";
# crap - couldn't get perl 5.8.8 to build of TACC
# just use old perl 5.8.5 and import Time::HiRes manually
# huh - perl 5.8.5 works?? on TACC - gettimeofday seems to work?
# has Time::HiRes already?? check ~/test_perl
use constant PERL_HOME_64_BIT => "/usr/bin";
use constant PERL_HOME => PERL_HOME_32_BIT;

# perl - different for TACC
use constant PERL_COMMAND => PERL_HOME . "/perl";

# enforce BASH interpret all shell scripts
use constant BASH_COMMAND => "/bin/bash";

# for python support
use constant PYTHON_COMMAND => "/lusr/bin/python";

1;
