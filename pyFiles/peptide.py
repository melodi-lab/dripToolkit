#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu> and Ajit Singh <ajit@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

__authors__ = [ 'John T. Halloran <halloj3@uw.edu>', 'Ajit Singh <ajit@ee.washington.edu>' ]

"""Peptides and operations on peptides.

A peptide is a short sequence of amino acids, usually a piece of a protein.
Amino acid tables are taken from [1].

[1] http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
"""

import random
import string

# Monoisotopic: elements are assumed to be the most common isotopes.
__amino_acid_mass = { 'A' : 71.03711,
                      'R' : 156.10111,
                      'N' : 114.04293,
                      'D' : 115.02694,
                      'C' : 103.00919 + 57.021464,
                      'E' : 129.04259,
                      'Q' : 128.05858,
                      'G' : 57.02146,
                      'H' : 137.05891,
                      'I' : 113.08406,
                      'L' : 113.08406,
                      'K' : 128.09496,
                      'M' : 131.04049,
                      'F' : 147.06841,
                      'P' : 97.05276,
                      'S' : 87.03203,
                      'T' : 101.04768,
                      'W' : 186.07931,
                      'Y' : 163.06333,
                      'V' : 99.06841,
                      'J' : float('inf'),
                      'O' : float('inf'),
                      'U' : float('inf'),
                      'B' : float('inf'),
                      'X' : float('inf'),
                      'Z' : float('inf'),
                      '*' : float('inf')}

__amino_acid_mass_pure = { 'A' : 71.03711,
                           'R' : 156.10111,
                           'N' : 114.04293,
                           'D' : 115.02694,
                           'C' : 103.00919,
                           'E' : 129.04259,
                           'Q' : 128.05858,
                           'G' : 57.02146,
                           'H' : 137.05891,
                           'I' : 113.08406,
                           'L' : 113.08406,
                           'K' : 128.09496,
                           'M' : 131.04049,
                           'F' : 147.06841,
                           'P' : 97.05276,
                           'S' : 87.03203,
                           'T' : 101.04768,
                           'W' : 186.07931,
                           'Y' : 163.06333,
                           'V' : 99.06841,
                           'J' : float('inf'),
                           'O' : float('inf'),
                           'U' : float('inf'),
                           'B' : float('inf'),
                           'X' : float('inf'),
                           'Z' : float('inf'),
                           '*' : float('inf')}

# Average: weighted average of isotopic masses, weighted by isotope abundance.
__amino_acid_average_mass = { 'A' : 71.0788,
                              'R' : 156.1875,
                              'N' : 114.1038,
                              'D' : 115.0886,
                              'C' : 103.1388 + 57.021464,
                              'E' : 129.1155,
                              'Q' : 128.1307,
                              'G' : 57.0519,
                              'H' : 137.1411,
                              'I' : 113.1594,
                              'L' : 113.1594,
                              'K' : 128.1741,
                              'M' : 131.1926,
                              'F' : 147.1766,
                              'P' : 97.1167,
                              'S' : 87.0782,
                              'T' : 101.1051,
                              'W' : 186.2132,
                              'Y' : 163.1760,
                              'V' : 99.1326,
                              'J' : float('inf'),
                              'O' : float('inf'),
                              'U' : float('inf'),
                              'B' : float('inf'),
                              'X' : float('inf'),
                              'Z' : float('inf'),
                              '*' : float('inf')}

__amino_acid_average_mass_pure = { 'A' : 71.0788,
                                   'R' : 156.1875,
                                   'N' : 114.1038,
                                   'D' : 115.0886,
                                   'C' : 103.1388,
                                   'E' : 129.1155,
                                   'Q' : 128.1307,
                                   'G' : 57.0519,
                                   'H' : 137.1411,
                                   'I' : 113.1594,
                                   'L' : 113.1594,
                                   'K' : 128.1741,
                                   'M' : 131.1926,
                                   'F' : 147.1766,
                                   'P' : 97.1167,
                                   'S' : 87.0782,
                                   'T' : 101.1051,
                                   'W' : 186.2132,
                                   'Y' : 163.1760,
                                   'V' : 99.1326,
                                   'J' : float('inf'),
                                   'O' : float('inf'),
                                   'U' : float('inf'),
                                   'B' : float('inf'),
                                   'X' : float('inf'),
                                   'Z' : float('inf'),
                                   '*' : float('inf')}

__amino_acid_letter_codes = {'A' : 'Ala',
                             'R' : 'Arg',
                             'N' : 'Asn',
                             'D' : 'Asp',
                             'C' : 'Cys',
                             'E' : 'Glu',
                             'Q' : 'Gln',
                             'G' : 'Gly',
                             'H' : 'His',
                             'I' : 'Ile',
                             'L' : 'Leu',
                             'K' : 'Lys',
                             'M' : 'Met',
                             'F' : 'Phe',
                             'P' : 'Pro',
                             'S' : 'Ser',
                             'T' : 'Thr',
                             'W' : 'Trp',
                             'Y' : 'Tyr',
                             'V' : 'Val'}
# ,
#                              'J' : 'NUL',
#                              'O' : 'NUL',
#                              'B' : 'NUL',
#                              'Z' : 'NUL',
#                              'U' : 'NUL',
#                              'X' : 'NUL

__amino_acids = __amino_acid_mass.keys()
__amino_acids.sort()

def mass_table(tbl = 'monoisotopic', pure = False):
    """Return a dictionary of IUPAC amino acid codes to masses (in Daltons).
    """
    if tbl == 'monoisotopic':
        if pure: return __amino_acid_mass_pure
        else: return __amino_acid_mass
    elif tbl == 'average':
        if pure: return __amino_acid_average_mass_pure
        else: return __amino_acid_average_mass
    else:
        raise ValueError("Argument must be one of ['monoisotopic', 'average'].")

def IUPAC_table(tbl = '1-to-3'):
    """Return a dictionary of IUPAC 1-letter to 3-letter codes, or IUPAC
    3-letter to 1-letter codes: e.g., 'M' -> 'Met' or 'Met' -> 'M'
    """
    if tbl == '1-to-3':
        return __amino_acid_letter_codes
    elif tbl == '3-to-1':
        return dict([(v, k) for k,v in list(__amino_acid_letter_codes.items())])
    else:
        raise ValueError("Argument must be one of ['1-to-3', '3-to-1'].")

def amino_acids():
    return __amino_acids

def amino_acids_to_indices():
    """Return a dictionary that maps amino acid letter codes to its index.

    The returned dictionary can be used to map letter like 'A' to their integer
    index, which is extensively used in CPTs which condition on amino acid
    identity.

    """
    dic = { }
    for index, aa in enumerate(amino_acids()):
        dic[aa] = index
    return dic

class Peptide(object):
    """Short sequence of amino acids.

    Peptides represent a sequence of amino acids, a string where each
    letter corresponds to the IUPAC code for an amino acid (e.g., A = Arginine).
    This class encapsulates the string sequence, and common operations on
    peptides. Peptides hash by their amino acid sequence; not their object id
    in memory.
    """

    def __init__(self, sequence = ''):
        """Test validity of peptide sequence string."""
        self.seq = sequence.upper() # String representation: e.g., 'KVRN'
        if not set(self.seq) <= set(amino_acids()):
            diff = list(set(self.seq) - set(amino_acids()))
            raise ValueError("Argument 'sequence' contains letters not in "
                             "IUPAC amino acid tables: " + str(diff))

    def __cmp__(self,other):
        if self.seq < other.seq:
            return -1
        elif self.seq > other.seq:
            return 1
        else:
            return 0

    def __hash__(self):
        return hash(self.seq)

    def __str__(self):
        return self.seq

    @property
    def length(self):
        """Number of amino acids in the peptide. Read-only computed property."""
        return len(self.seq)

    @property
    def mass(self):
        """Mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('monoisotopic', pure = False)) + 18.010564684

    @property
    def pure_mass(self):
        """Mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('monoisotopic', pure = True)) + 18.010564684

    @property
    def average_mass(self):
        """Isotope weighted mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('average', pure = False)) + 18.0153

    @property
    def pure_average_mass(self):
        """Isotope weighted mass of the peptide. Read-only computed property."""
        return self.__compute_mass(mass_table('average', pure = True)) + 18.0153

    def __compute_mass(self, mass):
        """Compute the mass of a peptide using a selected mass table."""
        #return sum(map(lambda aa: mass[aa], list(self.seq)))
        return sum([mass[aa] for aa in list(self.seq)])

    def reverse(self, tryptic = True):
        """Return a reverse copy of the peptide sequence, as a Peptide object.

        If tryptic is set to true, then we guarantee that the last position,
        which is 'R' or 'K' in a tryptic peptide, is not reversed.

        """
        if len(self.seq) == 0:
            return Peptide(self.seq)

        residue_list = list(self.seq)
        if tryptic:
            last_char = residue_list.pop()
            residue_list.reverse()
            residue_list.append(last_char)
        else:
            residue_list.reverse()
        return Peptide(string.join(residue_list, sep = ''))

    def ideal_fragment_masses(self, tbl = 'monoisotopic',
                              mass_op = lambda x: x):
        """Return the masses of the synthetic fragmentation products of self.

        Consider a peptide of length n as a character sequence: p[0..n-1].
        We return two arrays of length n+1, (left_mass, right_mass), where

        for each i = 0...n
        left_mass[i] = mass(p[0..i-1])  [ i.e., left_mass[0] = 0 ]
        right_mass[i] = mass(p[i..n-1])  [ i.e., right_mass[n] = 0 ]

        and where mass( ) is the mass of a peptide sequence. We do
        *not* consider the effect of mobile proton losses, or neutral
        losses.

        Args:
            tbl: string indicating the type of mass table to use. See
                protein.peptide.mass_table() for all possible values.
            mass_op: transformation of each of the mass table elements,
                i.e., a scalar function which takes a float and returns
                a numeric value.

        Returns:
            A tuple containing two vectors.
            left_mass:
            right_mass:

        """
        residue_mass = mass_table(tbl)
        mass = lambda aa: mass_op(residue_mass[aa])
        left_mass = [ 0 ]
        for i, aa in enumerate(self.seq):
            left_mass.append(left_mass[i] + mass(aa))

        right_mass = [ 0 ]
        for i, aa in enumerate(self.seq[::-1]): # reversed string
            right_mass.append(right_mass[i] + mass(aa))
        right_mass.reverse()
        return (left_mass, right_mass)

    def n_masses(self, tbl = 'monoisotopic'):
        """ Return the n-term masses at each amide bond location.
        Args:
            tbl:    (optional) string indicating the type of mass table to use.
                    Defaults to 'monoisotopic'; otherwise, use 'average.'
        Returns:
            A list of the n-term masses for the peptide.
        """
        residue_mass = mass_table(tbl)
        mass = lambda aa: residue_mass[aa]
        n_masses = []
        # omit the last amino acid -- just want n-masses at amide bonds
        for i, aa in enumerate(self.seq[:-1]):
            if i == 0:
                n_masses.append(mass(aa))
            else:
                n_masses.append(n_masses[i-1] + mass(aa))
        return n_masses

    def c_masses(self, tbl = 'monoisotopic'):
        """ Return the c-term masses at each amide bond location.
        Args:
            tbl:    (optional) string indicating the type of mass table to use.
                    Defaults to 'monoisotopic'; otherwise, use 'average.'
        Returns:
            A list of the c-term masses for the peptide.
        """
        residue_mass = mass_table(tbl)
        mass = lambda aa: residue_mass[aa]
        c_masses = []
        # omit the first amino acid -- just want c-masses at amide bonds
        for i, aa in enumerate(self.seq[:0:-1]):
            if i == 0:
                c_masses.append(mass(aa))
            else:
                c_masses.append(c_masses[i-1] + mass(aa))
        return c_masses
