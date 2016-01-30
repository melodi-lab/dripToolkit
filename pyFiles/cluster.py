#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

import os

################ bash using /tmp
def write_cluster_script(args, outfn, num_pickle, 
                         tmpDir = '/tmp'):
    bash_script = os.path.join(args.output_dir, 'split%d.sh' % num_pickle)
    with open(bash_script, 'w') as out:
        print >> out, """#!/bin/bash
"""
        print >> out, "TMPDIR=$(mktemp -d %s)" % (os.path.join(tmpDir, 'drip.XXXXXXXXXX'))
        print >> out, """
cd $TMPDIR
"""
        print >> out, "python -OO $DRIPTOOLKIT/dripSearch.py --cluster-mode True --spectra %s --output split%d-ident" % (outfn,num_pickle)
        print >> out, "NAP=%d" % args.random_wait
        print >> out, """number=$RANDOM
while [ "$number" -gt $NAP -o "$number" -le 0 ]
do
number=$RANDOM
done
sleep $number"""
        print >> out, "cp %s %s"  % ('split%d-ident.txt' % num_pickle, 
                                     os.path.join(os.path.abspath(args.logDir), 'split%d-ident.txt' % num_pickle))

################ bash using scratch space /s0
def write_scratch_cluster_script(args, outfn, num_pickle):
    bash_script = os.path.join(args.output_dir, 'split%d.sh' % num_pickle)
    with open(bash_script, 'w') as out:
        print >> out, """#!/bin/bash
if [ ! -d /s0/$USER ]
then
    space_req /s0
fi
"""
        print >> out, "TMPDIR=$(mktemp -d %s)" % '/s0/$USER/drip.XXXXXXXXXX'
        print >> out, """
cd $TMPDIR
"""
        print >> out, "python -OO $DRIPTOOLKIT/dripSearch.py --cluster-mode True --spectra %s --output split%d-ident" % (outfn,num_pickle)
        print >> out, "NAP=%d" % args.random_wait
        print >> out, """number=$RANDOM
while [ "$number" -gt $NAP -o "$number" -le 0 ]
do
number=$RANDOM
done
sleep $number"""
        print >> out, "cp %s %s"  % ('split%d-ident.txt' % num_pickle, 
                                     os.path.join(os.path.abspath(args.logDir), args.output+ '-split%d.txt' % num_pickle))
