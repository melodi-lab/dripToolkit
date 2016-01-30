#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu> and Ajit Singh <ajit@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

"""Tools for creating a peptide database, usually using the output generated
by crux search-for-matches: i.e., search.{target,decoy}.txt.

"""

from __future__ import with_statement

__authors__ = [ 'John T. Halloran <halloj3@uw.edu>', 'Ajit Singh <ajit@ee.washington.edu>' ]

from bisect import bisect_left, bisect_right
from operator import itemgetter
from peptide import Peptide
import sys
import csv

class SimplePeptideDB(object):
    """A simpler version of a peptide database that takes the peptides.txt
    file produces by crux create-index --peptide-list T <fasta-file> <index-dir>

    The peptide database is a text file where each line is of the form
    <peptide> <neutral mass>

    Where peptide is a IUPAC sequence, and neutral mass it the mass
    of a peptide (c.f.,
    http://noble.gs.washington.edu/proj/crux/crux-search-for-matches.html)

    """
    def __init__(self, peptides):
        peptides = [peptides[pep] for pep in peptides]
        peptides.sort(key = lambda x: x.peptideMass)
        self.peptides = [ p.peptide.seq for p in peptides ]
        self.masses =  [ p.peptideMass for p in peptides ]

    def filter(self, precursor_mass, window = 3, ppm = False):
        """Return the sequence of the peptides that are within 'window'
        Daltons of the mass given: i.e., where the peptide's mass is within
        [precursor_mass - window, precursor_mass + window]

        Note: avgmass is what is reported on the neutral mass line, and
        that should be a more reasonable mass calculation since the precursor
        mass comes from MS1, which reflects the average mass.

        """
        if not ppm:
            l = bisect_left(self.masses, precursor_mass - window)
            h = bisect_right(self.masses, precursor_mass + window, lo = l)
        else:
            wl = window*1e-6
            l = bisect_left(self.masses, precursor_mass / (1.0 + float(window)*0.000001))
            h = bisect_right(self.masses, precursor_mass / (1.0 - float(window)*0.000001), lo = l)
        return self.peptides[l:h]
