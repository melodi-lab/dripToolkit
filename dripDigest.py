#!/usr/bin/env python
#
# Written by John Halloran <halloj3@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

import sys
import array
import tempfile
import heapq

import struct

import os
import math
import mmap
import argparse
import random
import itertools
import re
import csv
import collections
import cPickle as pickle
from pyFiles.peptide import Peptide
from shutil import rmtree

validPeps = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')-set('JOBZUX')

def boolListToString(lst):
    return ''.join(['1' if x else '0' for x in lst])

def write_binary_peptide_tuples(s, records, f):
    """ Write binary records to file
        inputs:
        s - compiled struct, instance of struct.Struct
        records - array of tuples, each tuple of which conforms to the format depicted in s
        f - output file stream, assumed open for binary writing
    """
    for r in records:
        f.write(s.pack(*r))

def read_binary_peptide_tuples(s, f):
    """ Read binary records from file
        inputs:
        s - compiled struct, instance of struct.Struct
        f - input file stream, assumed open for binary writing
    """
    chunks = iter(lambda: f.read(s.size), b'')
    return (s.unpack(chunk) for chunk in chunks)

def read_binary_peptide_tuples_buffer(s, f, buff = 10000):
    """ Read binary records from file with a fixed buffer size
        inputs:
        s - compiled struct, instance of struct.Struct
        f - input file stream, assumed open for binary writing

        this functions acts as a limited-memory queue
    """
    while True:
        l = f.read(s.size * buff)
        start_ind = 0
        end_ind = s.size
        a = []
        for i in range(len(l) / s.size):
            a.append(s.unpack(l[start_ind:end_ind]))
            start_ind += s.size
            end_ind += s.size
        if not a:
            break
        for x in a:
            yield x

class digestedPep(object):
    """ Digested peptide class containing parent protein and flanking amino acid information

    """
    def __init__(self, peptide = '', proteinName = '', 
                 ntermFlanking = '', ctermFlanking = '',
                 peptideMass = -1, 
                 missed_cleavages = 0):

        self.peptide = peptide
        self.peptideMass = peptideMass
        self.proteinName = proteinName
        self.ntermFlanking = ntermFlanking
        self.ctermFlanking = ctermFlanking
        self.missed_cleavages = missed_cleavages

        if peptideMass < 0.0:
            raise ValueError("Peptide %s mass is %f, must be non-negative" % (peptide, peptideMass))

    def __cmp__(self,other):
        if self.peptide < other.peptide:
            return -1
        elif self.peptide > other.peptide:
            return 1
        else:
            return 0

    def __hash__(self):
        # return hash((self.peptide, self.proteinName))
        return hash(self.peptide)

    def __str__(self):
        return self.peptide

    @property
    def length(self):
        """Number of amino acids in the peptide. Read-only computed property."""
        return len(self.seq)

def adapt_digestedPep(dp):
    return "%d;%s;%f;%s;%s;%s" % (int(dp.peptideMass), dp.peptide, dp.peptideMass, 
                                  dp.proteinName, dp.ntermFlanking,
                                  dp.ctermFlanking)

def convert_digestedPep(s):
    l = s.split(";")
    return digestedPep(int(l[0]), l[1], float(l[2]), l[3], l[4], l[5])

def store_db_pickle(filename, targets, decoys = None):
    if decoys:
        data = {'targets' : targets,
                'decoys' : decoys}
    else:
        data = {'targets' : targets}

    with open(filename, 'w') as out:
        pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)

def load_db_pickle(filename):
    return pickle.load(open(filename))

def check_arg_trueFalse(a):
    """ Check boolean command line options 
    
    inputs:
    a - command line option

    outputs:
    True for '1' and all lower/upper case versions of 't' or 'true'
    False for '0' and all lower/upper case versions of 'f' or 'false'

    exception is thrown if a python string is not input
    """
    try:
        a = a.lower()
    except:
        print "Supplied argument %s must be a string indicating true or false, exitting" % a
        exit(-1)
    if a == 't' or a == 'true' or a == '1':
        return True
    if a== 'f' or a == 'false' or a == '0':
        return False

def create_digest_re(enzyme, digest_rules):
    """ Create regular expression for digestion rules
    iputs (note that digest_rules takes precedence over enzyme):
    enzyme - Digesting enzyme.  From enzyme option help:
      <no-enzyme|trypsin|trypsin/p|chymotrypsin|elastase|clostripain|cyanogen-
      bromide|iodosobenzoate|proline-endopeptidase|staph-protease|asp-n|lys-c|
      lys-n|arg-c|glu-c|pepsin-a|elastase-trypsin-chymotrypsin|
      custom-enzyme> - Specify the enzyme used to digest the proteins in silico. 
      Available enzymes (with the corresponding digestion rules indicated in 
      parentheses) include no-enzyme ([X]|[X]), trypsin ([RK]|{P}), trypsin/p 
      ([RK]|[]), chymotrypsin ([FWYL]|{P}), elastase ([ALIV]|{P}), 
      clostripain ([R]|[]), cyanogen-bromide ([M]|[]), iodosobenzoate 
      ([W]|[]), proline-endopeptidase ([P]|[]), staph-protease ([E]|[]), 
      asp-n ([]|[D]), lys-c ([K]|{P}), lys-n ([]|[K]), arg-c ([R]|{P}), 
      glu-c ([DE]|{P}), pepsin-a ([FL]|{P}), elastase-trypsin-chymotrypsin 
      ([ALIVKRWFY]|{P}). Specifying --enzyme no-enzyme yields a non-enzymatic digest. 
      Warning: the resulting peptide database may be quite large. Default = trypsin.

    digest_rules - Rules detailing digest.  From custom-enzyme option help:
      <string> - Specify rules for in silico digestion of protein sequences. 
      Overrides the enzyme option. Two lists of residues are given enclosed 
      in square brackets or curly braces and separated by a |. The first 
      list contains residues required/prohibited before the cleavage site 
      and the second list is residues after the cleavage site. If the 
      residues are required for digestion, they are in square brackets, 
      '[' and ']'. If the residues prevent digestion, then they are 
      enclosed in curly braces, '{' and '}'. Use X to indicate all residues. 
      For example, trypsin cuts after R or K but not before P which is 
      represented as [RK]|{P}. AspN cuts after any residue but only before 
      D which is represented as [X]|[D]. Default = <empty>.

    outpus:
    digest_re - python regular expression to digest proteins
    """
    enzyme = enzyme.lower()

    if not digest_rules:
        if enzyme == 'no-enzyme':
            digest_rules = '[X]|[X]'
        elif enzyme == 'trypsin':
            digest_rules = '[RK]|{P}'
        elif enzyme == 'trypsin/p':
            digest_rules = '[RK]|[]'
        elif enzyme == 'chymotrypsin':
            digest_rules = '[FWYL]|{P}'
        elif enzyme == 'elastase':
            digest_rules = '[ALIV]|{P}'
        elif enzyme == 'clostripain':
            digest_rules = '[R]|[]'
        elif enzyme == 'cyanogen-bromide':
            digest_rules = '[M]|[]'
        elif enzyme == 'iodosobenzoate':
            digest_rules = '[W]|[]'
        elif enzyme == 'proline-endopeptidase':
            digest_rules = '[P]|[]'
        elif enzyme == 'staph-protease':
            digest_rules = '[E]|[]'
        elif enzyme == 'asp-n':
            digest_rules = '[]|[D]'
        elif enzyme == 'lys-c':
            digest_rules = '[K]|{P}'
        elif enzyme == 'lys-n':
            digest_rules = '[]|[K]'
        elif enzyme == 'arg-c':
            digest_rules = '[R]|{P}'
        elif enzyme == 'glu-c':
            digest_rules = '[DE]|{P}'
        elif enzyme == 'pepsin-a':
            digest_rules = '[FL]|{P}'
        elif enzyme == 'elastase-trypsin-chymotrypsin':
            digest_rules = '[ALIVKRWFY]|{P}'
        else:
            print "%s invalid digestion enzyme.  For option enzyme:"
            print enzyme_help
            exit(1)

    # parse digestion rule
    h = digest_rules.split('|')
    if len(h) != 2:
        print "Invalid custom enzyme specified.  For option custom-enzyme:"
        print custom_enzyme_help
        exit(2)

    nterm = h[0]
    if nterm[0]=='[':
        if nterm[-1] != ']':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(3)

        include_nterm = True
        if nterm[1:-1]:
            nterm_aa = nterm[1:-1]
        else:
            nterm_aa = ''
    elif nterm[0]=='{':
        if nterm[-1] != '}':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(4)

        include_nterm = False
        if nterm[1:-1]:
            nterm_aa = nterm[1:-1]
        else:
            nterm_aa = ''
            
    cterm = h[1]
    if cterm[0]=='[':
        if cterm[-1] != ']':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(5)

        include_cterm = True
        if cterm[1:-1]:
            cterm_aa = cterm[1:-1]
        else:
            cterm_aa = ''
    elif cterm[0]=='{':
        if cterm[-1] != '}':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(6)

        include_cterm = False
        if cterm[1:-1]:
            cterm_aa = cterm[1:-1]
        else:
            cterm_aa = ''


    if nterm_aa == 'X':
        nterm_aa = ''
    if cterm_aa == 'X':
        cterm_aa = ''
    # form digestion regular expression
    # tryptic digest with lysine, arginine secondary digestion suppression ".(?:(?<![KR](?!P|K|R)).)*"
    if include_nterm:
        if not include_cterm:
            digest_re = ".(?:(?<!" 
            if nterm_aa:
                digest_re += "[" + nterm_aa + "]" 
            if cterm_aa:
                digest_re += "(?!"
                for aa in cterm_aa[:-1]:
                    digest_re += aa
                    digest_re += "|"
                digest_re += cterm_aa[-1]
                digest_re += ")"
            digest_re += ").)*"
        else:
            digest_re = ".(?:(?<!" 
            if nterm_aa:
                digest_re += "[" + nterm_aa + "]" 
            if cterm_aa:
                digest_re += "["
                digest_re += cterm_aa
                digest_re += "]"
            digest_re += ").)*"
    else:
        if not include_cterm:
            digest_re = ".(?:(?<!" 
            if nterm_aa:
                digest_re += "[^" + nterm_aa + "]" 
            if cterm_aa:
                digest_re += "(?!"
                for aa in cterm_aa[:-1]:
                    digest_re += aa
                    digest_re += "|"
                digest_re += cterm_aa[-1]
                digest_re += ")"
            digest_re += ").)*"
        else:
            digest_re = ".(?:(?<!" 
            if nterm_aa:
                digest_re += "[^" + nterm_aa + "]" 
            if cterm_aa:
                digest_re += "["
                digest_re += cterm_aa
                digest_re += "]"
            digest_re += ").)*"

    return digest_re

def enzyme_enzSet(enzyme, digest_rules):
    """ Create regular expression for digestion rules
    iputs (note that digest_rules takes precedence over enzyme):
    enzyme - Digesting enzyme.  From enzyme option help:
      <no-enzyme|trypsin|trypsin/p|chymotrypsin|elastase|clostripain|cyanogen-
      bromide|iodosobenzoate|proline-endopeptidase|staph-protease|asp-n|lys-c|
      lys-n|arg-c|glu-c|pepsin-a|elastase-trypsin-chymotrypsin|
      custom-enzyme> - Specify the enzyme used to digest the proteins in silico. 
      Available enzymes (with the corresponding digestion rules indicated in 
      parentheses) include no-enzyme ([X]|[X]), trypsin ([RK]|{P}), trypsin/p 
      ([RK]|[]), chymotrypsin ([FWYL]|{P}), elastase ([ALIV]|{P}), 
      clostripain ([R]|[]), cyanogen-bromide ([M]|[]), iodosobenzoate 
      ([W]|[]), proline-endopeptidase ([P]|[]), staph-protease ([E]|[]), 
      asp-n ([]|[D]), lys-c ([K]|{P}), lys-n ([]|[K]), arg-c ([R]|{P}), 
      glu-c ([DE]|{P}), pepsin-a ([FL]|{P}), elastase-trypsin-chymotrypsin 
      ([ALIVKRWFY]|{P}). Specifying --enzyme no-enzyme yields a non-enzymatic digest. 
      Warning: the resulting peptide database may be quite large. Default = trypsin.

    digest_rules - Rules detailing digest.  From custom-enzyme option help:
      <string> - Specify rules for in silico digestion of protein sequences. 
      Overrides the enzyme option. Two lists of residues are given enclosed 
      in square brackets or curly braces and separated by a |. The first 
      list contains residues required/prohibited before the cleavage site 
      and the second list is residues after the cleavage site. If the 
      residues are required for digestion, they are in square brackets, 
      '[' and ']'. If the residues prevent digestion, then they are 
      enclosed in curly braces, '{' and '}'. Use X to indicate all residues. 
      For example, trypsin cuts after R or K but not before P which is 
      represented as [RK]|{P}. AspN cuts after any residue but only before 
      D which is represented as [X]|[D]. Default = <empty>.

    outpus:
    digest_re - python regular expression to digest proteins
    """
    enzyme = enzyme.lower()

    if not digest_rules:
        if enzyme == 'no-enzyme':
            digest_rules = '[X]|[X]'
        elif enzyme == 'trypsin':
            digest_rules = '[RK]|{P}'
        elif enzyme == 'trypsin/p':
            digest_rules = '[RK]|[]'
        elif enzyme == 'chymotrypsin':
            digest_rules = '[FWYL]|{P}'
        elif enzyme == 'elastase':
            digest_rules = '[ALIV]|{P}'
        elif enzyme == 'clostripain':
            digest_rules = '[R]|[]'
        elif enzyme == 'cyanogen-bromide':
            digest_rules = '[M]|[]'
        elif enzyme == 'iodosobenzoate':
            digest_rules = '[W]|[]'
        elif enzyme == 'proline-endopeptidase':
            digest_rules = '[P]|[]'
        elif enzyme == 'staph-protease':
            digest_rules = '[E]|[]'
        elif enzyme == 'asp-n':
            digest_rules = '[]|[D]'
        elif enzyme == 'lys-c':
            digest_rules = '[K]|{P}'
        elif enzyme == 'lys-n':
            digest_rules = '[]|[K]'
        elif enzyme == 'arg-c':
            digest_rules = '[R]|{P}'
        elif enzyme == 'glu-c':
            digest_rules = '[DE]|{P}'
        elif enzyme == 'pepsin-a':
            digest_rules = '[FL]|{P}'
        elif enzyme == 'elastase-trypsin-chymotrypsin':
            digest_rules = '[ALIVKRWFY]|{P}'
        else:
            print "%s invalid digestion enzyme.  For option enzyme:"
            print enzyme_help
            exit(1)

    # parse digestion rule
    h = digest_rules.split('|')
    if len(h) != 2:
        print "Invalid custom enzyme specified.  For option custom-enzyme:"
        print custom_enzyme_help
        exit(2)

    nterm = h[0]
    if nterm[0]=='[':
        if nterm[-1] != ']':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(3)

        include_nterm = True
        if nterm[1:-1]:
            nterm_aa = nterm[1:-1]
        else:
            nterm_aa = ''
    elif nterm[0]=='{':
        if nterm[-1] != '}':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(4)

        include_nterm = False
        if nterm[1:-1]:
            nterm_aa = nterm[1:-1]
        else:
            nterm_aa = ''
            
    cterm = h[1]
    if cterm[0]=='[':
        if cterm[-1] != ']':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(5)

        include_cterm = True
        if cterm[1:-1]:
            cterm_aa = cterm[1:-1]
        else:
            cterm_aa = ''
    elif cterm[0]=='{':
        if cterm[-1] != '}':
            print "Error reading digesting enzyme, missing closing bracket ]"
            exit(6)

        include_cterm = False
        if cterm[1:-1]:
            cterm_aa = cterm[1:-1]
        else:
            cterm_aa = ''


    if nterm_aa == 'X' or nterm_aa == '':
        left_cleavage_site = validPeps
    else:
        left_cleavage_site = set(nterm_aa)
    if cterm_aa == 'X' or cterm_aa == '':
        right_cleavage_site = validPeps
    else:
        right_cleavage_site = set(cterm_aa)

    if not include_nterm:
        left_cleavage_site = validPeps - left_cleavage_site
    if not include_cterm:
        right_cleavage_site = validPeps - right_cleavage_site

    return [left_cleavage_site,right_cleavage_site]

def check_peptide(p, min_length, max_length,
                  protein_name, isMonoMass,
                  min_mass, max_mass, 
                  nterm_mods, cterm_mods, mods,
                  nterm_flanking, cterm_flanking):
    # check peptide mass
    if len(p) < min_length or len(p) > max_length:
        return ''
    pep = Peptide(p)
    try:
        pep = Peptide(p)
    except:
        print p
        print protein
        exit(7)

    if isMonoMass:
        pm = pep.mass
    else:
        pm = pep.average_mass
    if pm < min_mass or pm > max_mass:
        return ''

    # modifications
    if p[0] in nterm_mods:
        pm += nterm_mods[p[0]]
    elif p[0] in mods:
        pm += mods[p[0]]
    if p[-1] in cterm_mods:
        pm += cterm_mods[p[-1]]
    elif p[-1] in mods:
        pm += mods[p[-1]]
    if mods:
        for aa in p[1:-1]:
            if aa in mods:
                pm += mods[aa]

    return(digestedPep(p, protein_name,
                       nterm_flanking, 
                       cterm_flanking, 
                       pm))

def check_peptide_dbEntry(p, min_length, max_length,
                          protein_name, isMonoMass,
                          min_mass, max_mass, 
                          nterm_mods, cterm_mods, mods,
                          nterm_flanking, cterm_flanking):
    # check peptide mass
    if len(p) < min_length or len(p) > max_length:
        return ''
    try:
        pep = Peptide(p)
    except:
        print p
        print protein
        exit(7)

    if isMonoMass:
        pm = pep.mass
    else:
        pm = pep.average_mass
    if pm < min_mass or pm > max_mass:
        return ''

    # modifications
    if p[0] in nterm_mods:
        pm += nterm_mods[p[0]]
    elif p[0] in mods:
        pm += mods[p[0]]
    if p[-1] in cterm_mods:
        pm += cterm_mods[p[-1]]
    elif p[-1] in mods:
        pm += mods[p[-1]]
    if mods:
        for aa in p[1:-1]:
            if aa in mods:
                pm += mods[aa]

    return(int(pm), p, pm, protein_name, nterm_flanking, cterm_flanking)

def check_peptide_tuple(p, min_length, max_length,
                        protein_name, isMonoMass,
                        min_mass, max_mass, 
                        nterm_mods, cterm_mods, mods,
                        nterm_flanking, cterm_flanking):
    # check peptide mass
    if len(p) < min_length or len(p) > max_length:
        return ''
    pep = Peptide(p)
    # try:
    #     pep = Peptide(p)
    # except:
    #     print p
    #     print protein_name
    #     exit(7)

    if isMonoMass:
        pm = pep.mass
    else:
        pm = pep.average_mass

    # modifications
    if p[0] in nterm_mods:
        pm += nterm_mods[p[0]]
    elif p[0] in mods:
        pm += mods[p[0]]
    if p[-1] in cterm_mods:
        pm += cterm_mods[p[-1]]
    elif p[-1] in mods:
        pm += mods[p[-1]]
    if mods:
        for aa in p[1:-1]:
            if aa in mods:
                pm += mods[aa]

    if pm < min_mass or pm > max_mass:
        return ''

    return (pm, p, protein_name, nterm_flanking, cterm_flanking)

def check_peptide_tuple_var_mods(p, min_length, max_length,
                                 isMonoMass, min_mass, max_mass, 
                                 nterm_mods, cterm_mods, mods,
                                 nterm_var_mods, cterm_var_mods, var_mods,
                                 max_mods, min_mods):
    # check peptide mass
    if len(p) < min_length or len(p) > max_length:
        return ''
    pep = Peptide(p)

    num_stat_mods = 0

    if isMonoMass:
        pm = pep.mass
    else:
        pm = pep.average_mass

    if pm == float('inf'): # string contains character outside of 20 amino acids
        return ''

    # modifications
    if p[0] in nterm_mods:
        pm += nterm_mods[p[0]]
        num_stat_mods += 1
    elif p[0] in mods:
        pm += mods[p[0]]
        num_stat_mods += 1
    if p[-1] in cterm_mods:
        pm += cterm_mods[p[-1]]
        num_stat_mods += 1
    elif p[-1] in mods:
        pm += mods[p[-1]]
        num_stat_mods += 1
    if mods:
        for aa in p[1:-1]:
            if aa in mods:
                pm += mods[aa]
                num_stat_mods += 1

    # identify indices of possible variable mods
    var_aas = []
    for i,aa in enumerate(p):
        if aa in var_mods:
            var_aas.append(i)

    valid_peps = []

    if not var_aas and p[0] not in nterm_var_mods and p[-1] not in cterm_var_mods:
        if pm >= min_mass and pm <= max_mass:
            return [((), pm)]
        else:
            return valid_peps

    if p[0] not in nterm_var_mods and p[-1] not in cterm_var_mods:
        check_peptide_var_mod_subsets(valid_peps, p, var_aas, pm,
                                      max_mods, min_mods,
                                      var_mods, 
                                      min_mass, max_mass,
                                      num_stat_mods)

    elif p[0] in nterm_var_mods and p[-1] in cterm_var_mods:
        check_peptide_var_mod_subsets_nterm_cterm(valid_peps, p, var_aas, pm,
                                                  max_mods, min_mods,
                                                  var_mods, nterm_var_mods, cterm_var_mods, 
                                                  min_mass, max_mass,
                                                  num_stat_mods)
    elif p[0] in nterm_var_mods:
        check_peptide_var_mod_subsets_nterm(valid_peps, p, var_aas, pm,
                                            max_mods, min_mods,
                                            var_mods, nterm_var_mods, 
                                            min_mass, max_mass,
                                            num_stat_mods)
    else:
        check_peptide_var_mod_subsets_cterm(valid_peps, p, var_aas, pm,
                                            max_mods, min_mods,
                                            var_mods, cterm_var_mods, 
                                            min_mass, max_mass,
                                            num_stat_mods)
    return (valid_peps)

def check_peptide_var_mod_subsets(valid_peps, p, var_aas, pm,
                                  max_mods, min_mods,
                                  var_mods, 
                                  min_mass, max_mass,
                                  num_stat_mods = 0):
    # consider all possible variable mods, which is the powerset
    # of var_inds restricted to sets of at least size min_mods and
    # at most size max_mods
    for var_subset in itertools.chain.from_iterable(itertools.combinations(var_aas, r) for r in range(len(var_aas)+1)):
        if len(var_subset) > (max_mods + num_stat_mods): 
            break
        elif len(var_subset) >= min_mods:
            aa_tally = {}
            pm_var_mod_offset = pm
            valid_aas = True
            for ind in var_subset:
                aa = p[ind]
                if aa not in aa_tally:
                    aa_tally[aa] = 1
                else:
                    aa_tally[aa] += 1

                if aa_tally[aa] > var_mods[aa][0]: # aa is modified more than desired number of times
                    valid_aas = False
                    break

                pm_var_mod_offset += var_mods[aa][1]

            if pm_var_mod_offset >= min_mass or pm_var_mod_offset <= max_mass:
                if valid_aas:
                    valid_peps.append((var_subset, pm_var_mod_offset))

def check_peptide_var_mod_subsets_nterm(valid_peps, p, var_aas, pm,
                                        max_mods, min_mods,
                                        var_mods, nterm_var_mods, 
                                        min_mass, max_mass,
                                        num_stat_mods = 0):
    # consider all possible variable mods, which is the powerset
    # of var_inds restricted to sets of at least size min_mods and
    # at most size max_mods

    # pre: p[0] is in specified n-terminal mods, 
    for var_subset in itertools.chain.from_iterable(itertools.combinations(var_aas, r) for r in range(len(var_aas)+1)):
        if len(var_subset) > (max_mods + num_stat_mods): 
            break
        elif len(var_subset) >= min_mods:
            aa_tally = {}
            pm_var_mod_offset = pm
            valid_aas = True
            for ind in var_subset:
                aa = p[ind]
                if aa not in aa_tally:
                    aa_tally[aa] = 1
                else:
                    aa_tally[aa] += 1

                if aa_tally[aa] > var_mods[aa][0]: # aa is modified more than desired number of times
                    valid_aas = False
                    break

                pm_var_mod_offset += var_mods[aa][1]

            if pm_var_mod_offset < min_mass or pm_var_mod_offset > max_mass:
                continue

            if valid_aas:
                valid_peps.append((var_subset, pm_var_mod_offset))

                # add n-term var mod
                var_subset = list(var_subset)
                if len(var_subset) < max_mods and (0 not in var_subset):
                    if p[0] in nterm_var_mods:
                        pm_var_mod_offset += nterm_var_mods[p[0]][1]
                        var_subset.append(-1)
                        if pm_var_mod_offset >= min_mass and pm_var_mod_offset <= max_mass:
                            # we have a valid variable modified peptide
                            valid_peps.append((var_subset, pm_var_mod_offset))

def check_peptide_var_mod_subsets_cterm(valid_peps, p, var_aas, pm,
                                        max_mods, min_mods,
                                        var_mods, cterm_var_mods, 
                                        min_mass, max_mass,
                                        num_stat_mods = 0):
    # consider all possible variable mods, which is the powerset
    # of var_inds restricted to sets of at least size min_mods and
    # at most size max_mods

    # pre: p[0] is in specified n-terminal mods, 
    for var_subset in itertools.chain.from_iterable(itertools.combinations(var_aas, r) for r in range(len(var_aas)+1)):
        if len(var_subset) > (max_mods + num_stat_mods): 
            break
        elif len(var_subset) >= min_mods:
            aa_tally = {}
            pm_var_mod_offset = pm
            valid_aas = True
            for ind in var_subset:
                aa = p[ind]
                if aa not in aa_tally:
                    aa_tally[aa] = 1
                else:
                    aa_tally[aa] += 1

                if aa_tally[aa] > var_mods[aa][0]: # aa is modified more than desired number of times
                    valid_aas = False
                    break

                pm_var_mod_offset += var_mods[aa][1]

            if pm_var_mod_offset < min_mass or pm_var_mod_offset > max_mass:
                continue

            if valid_aas:
                valid_peps.append((var_subset, pm_var_mod_offset))

                var_subset = list(var_subset)
                # add n-term var mod
                if len(var_subset) < max_mods and (len(p)-1 not in var_subset):
                    if p[-1] in cterm_var_mods:
                        pm_var_mod_offset += cterm_var_mods[p[-1]][1]
                        var_subset.append(-2)
                        if pm_var_mod_offset >= min_mass and pm_var_mod_offset <= max_mass:
                            valid_peps.append((var_subset, pm_var_mod_offset))

def check_peptide_var_mod_subsets_nterm_cterm(valid_peps, p, var_aas, pm,
                                              max_mods, min_mods,
                                              var_mods, 
                                              nterm_var_mods, cterm_var_mods, 
                                              min_mass, max_mass,
                                              num_stat_mods = 0):
    # consider all possible variable mods, which is the powerset
    # of var_inds restricted to sets of at least size min_mods and
    # at most size max_mods

    # pre: p[0] is in specified n-terminal mods, 
    for var_subset in itertools.chain.from_iterable(itertools.combinations(var_aas, r) for r in range(len(var_aas)+1)):
        if len(var_subset) > (max_mods + num_stat_mods): 
            break
        elif len(var_subset) >= min_mods:
            aa_tally = {}
            pm_var_mod_offset = pm
            valid_aas = True
            for ind in var_subset:
                aa = p[ind]
                if aa not in aa_tally:
                    aa_tally[aa] = 1
                else:
                    aa_tally[aa] += 1

                if aa_tally[aa] > var_mods[aa][0]: # aa is modified more than desired number of times
                    valid_aas = False
                    break

                pm_var_mod_offset += var_mods[aa][1]

            if pm_var_mod_offset < min_mass or pm_var_mod_offset > max_mass:
                continue

            if valid_aas:
                valid_peps.append((var_subset, pm_var_mod_offset))

                var_subset2 = list(var_subset)
                pmo = pm_var_mod_offset

                # add n-term var mod
                if len(var_subset2) < max_mods and (0 not in var_subset2):
                    if p[0] in nterm_var_mods:
                        pmo += nterm_var_mods[p[0]][1]
                        if pmo >= min_mass or pmo <= max_mass:
                            var_subset2.append(-1)
                            valid_peps.append((var_subset2, pmo))
                            # add c-term variable mod in tandem with n-term variable mode
                            if len(var_subset2) < max_mods-1 and (len(p)-1 not in var_subset2):
                                if p[-1] in cterm_var_mods:
                                    pmo += cterm_var_mods[p[-1]][1]
                                    if pmo >= min_mass and pmo <= max_mass:
                                        var_subset2.append(-2)
                                        valid_peps.append((var_subset2, pmo))

                var_subset2 = list(var_subset)
                pmo = pm_var_mod_offset

                # add n-term var mod
                if len(var_subset2) < max_mods and (len(p)-1 not in var_subset2):
                    if p[-1] in cterm_var_mods:
                        pmo += cterm_var_mods[p[-1]][1]
                        if pmo >= min_mass and pmo <= max_mass:
                            var_subset2.append(-2)
                            valid_peps.append((var_subset2, pmo))

#################### peptides have been digested and written to a binary file, 
#################### read from binary file and evaluate variable mods
def check_digested_peptide_tuple_var_mods(p, pm, 
                                          min_length, max_length,
                                          min_mass, max_mass, 
                                          nterm_var_mods, cterm_var_mods, var_mods,
                                          max_mods, min_mods):
    # identify indices of possible variable mods
    var_aas = []
    for i,aa in enumerate(p):
        if aa in var_mods:
            var_aas.append(i)

    if not var_aas and p[0] not in nterm_var_mods and p[-1] not in cterm_var_mods:
        if pm >= min_mass and pm <= max_mass:
            return [((), pm)]
        else:
            return valid_peps

    valid_peps = []

    if p[0] not in nterm_var_mods and p[-1] not in cterm_var_mods:
        check_peptide_var_mod_subsets(valid_peps, p, var_aas, pm,
                                      max_mods, min_mods,
                                      var_mods, 
                                      min_mass, max_mass)
        
    elif p[0] in nterm_var_mods and p[-1] in cterm_var_mods:
        check_peptide_var_mod_subsets_nterm_cterm(valid_peps, p, var_aas, pm,
                                                  max_mods, min_mods,
                                                  var_mods, nterm_var_mods, cterm_var_mods, 
                                                  min_mass, max_mass)
    elif p[0] in nterm_var_mods:
        check_peptide_var_mod_subsets_nterm(valid_peps, p, var_aas, pm,
                                            max_mods, min_mods,
                                            var_mods, nterm_var_mods, 
                                            min_mass, max_mass)
    else:
        check_peptide_var_mod_subsets_cterm(valid_peps, p, var_aas, pm,
                                            max_mods, min_mods,
                                            var_mods, cterm_var_mods, 
                                            min_mass, max_mass)
    return (valid_peps)

def load_proteins(fastaFile, digest_regexp,
                  min_length, max_length,
                  min_mass, max_mass,
                  mods, nterm_mods, cterm_mods,
                  max_mods, min_mods, missed_cleavages = 0,
                  isMonoMass = True, 
                  digestion = 'full-digest'):
    """ Parse fasta file and digest proteins into peptides
    inputs:
    fastaFile - FASTA file containing proteins
    digest_regexp - re compiled regular expression dictating digestion rules
    min_length - minimum residue length
    max_length - maximum residue length

    output:
    peptides - dict whose keys are the digested peptide string, elements 
               are instances of digestedPep class with fields:
               peptide - instance of pyFiles.peptide's Peptide
               proteinName - name of parent protein from Fasta file
               ntermFlanking - nterm (left flanking) amino acid
               ctermFlanking - cterm (right flanking) amino acid
               peptideMass - calculated peptide mass

    todo: add missed cleavages, static/variable modification masses
    """
    peptides = {}
    
    with open(fastaFile) as f:
        for isname, group in itertools.groupby(f, lambda x: x[0]=='>'):
            if isname: # descriptor line
                l = group.next()
                protein_name = l[1:].split()[0]
            else: # protein sequence lines
                protein = ''.join(l.strip() for l in group)
                # protein = protein.replace('*', '')
                # h = digest_regexp.findall(protein)
                h = []
                for cutProt in protein.split('*'):
                    h += digest_regexp.findall(cutProt)

                hl = len(h)-1
                if len(h)==1:
                    p = h[0]
                    nterm_flanking = '-'
                    cterm_flanking = '-'
                    valid_p = check_peptide(p, min_length, max_length,
                                            protein_name, isMonoMass,
                                            min_mass, max_mass, 
                                            nterm_mods, cterm_mods, mods,
                                            nterm_flanking, cterm_flanking)
                    if valid_p:
                        peptides[p] = valid_p
                    if digestion == 'partial-digest':
                        # check all partial digests
                        
                        # left to right string flanking amino acid information
                        ltr_cterm_flanking = cterm_flanking
                        # right to left string flanking amino acid information
                        rtl_nterm_flanking = nterm_flanking
                        # for ind in range(min_length, len(p)+1):
                        for ind in range(len(p)+1):
                            curr_end = min(len(p)-1,ind)
                            if ind > 0:
                                # right to left string
                                rtl_nterm_flanking = p[ind-1]
                                if (rtl_nterm_flanking in c_res and p[ind] in n_res) or (p[-1] in c_res and cterm_flanking in n_res):
                                    valid_p = check_peptide(p[ind:], min_length, max_length,
                                                            protein_name, isMonoMass,
                                                            min_mass, max_mass, 
                                                            nterm_mods, cterm_mods, mods,
                                                            rtl_nterm_flanking, cterm_flanking)
                                    if valid_p:
                                        peptides[p[ind:]] = valid_p

                            if ind < len(p):
                                # left to right string
                                ltr_cterm_flanking = p[ind]
                                if (rtl_nterm_flanking in c_res and p[ind] in n_res) or (p[-1] in c_res and cterm_flanking in n_res):
                                    valid_p = check_peptide(p[:ind], min_length, max_length,
                                                            protein_name, isMonoMass,
                                                            min_mass, max_mass, 
                                                            nterm_mods, cterm_mods, mods,
                                                            nterm_flanking, ltr_cterm_flanking)
                                    if valid_p:
                                        peptides[p[:ind]] = valid_p
                else:
                    for i in range(len(h)):
                        p = ''
                        for j in range(i, min(hl, i+missed_cleavages)+1):
                            p += h[j]
                            # get flanking information
                            if i == 0:
                                nterm_flanking = '-'
                            else:
                                nterm_flanking = h[i-1][-1]
                            if j == hl:
                                cterm_flanking = '-'
                            else:
                                cterm_flanking = h[j+1][0]

                            valid_p = check_peptide(p, min_length, max_length,
                                                    protein_name, isMonoMass,
                                                    min_mass, max_mass, 
                                                    nterm_mods, cterm_mods, mods,
                                                    nterm_flanking, cterm_flanking)
                            if valid_p:
                                peptides[p] = valid_p

                        if digestion == 'partial-digest':

                            # if p == 'HDDGDIDGDDFYK':
                            #     print 'hi'
                            # p is now the longest substring including all allowed missed cleavages
                            # check all partial digests

                            # left to right string flanking amino acid information
                            ltr_cterm_flanking = cterm_flanking
                            # right to left string flanking amino acid information
                            rtl_nterm_flanking = nterm_flanking
                            for ind in range(len(p)+1):
                                if ind > 0:
                                    # right to left string
                                    rtl_nterm_flanking = p[ind-1]
                                    if (rtl_nterm_flanking in c_res and p[ind] in n_res) or (p[-1] in c_res and cterm_flanking in n_res):
                                        valid_p = check_peptide(p[ind:], min_length, max_length,
                                                                protein_name, isMonoMass,
                                                                min_mass, max_mass, 
                                                                nterm_mods, cterm_mods, mods,
                                                                rtl_nterm_flanking, cterm_flanking)
                                        if valid_p:
                                            peptides[p[ind:]] = valid_p
                                if ind < len(p):
                                    # left to right string
                                    ltr_cterm_flanking = p[ind]
                                    if (rtl_nterm_flanking in c_res and p[ind] in n_res) or (p[-1] in c_res and cterm_flanking in n_res):
                                        valid_p = check_peptide(p[:ind], min_length, max_length,
                                                                protein_name, isMonoMass,
                                                                min_mass, max_mass, 
                                                                nterm_mods, cterm_mods, mods,
                                                                nterm_flanking, ltr_cterm_flanking)
                                        if valid_p:
                                            peptides[p[:ind]] = valid_p
                        # # old code
                        # p = h[i]
                        # # check peptide length
                        # if len(p) < min_length or len(p) > max_length:
                        #     continue
                        # # check peptide mass
                        # try:
                        #     pep = Peptide(p)
                        # except:
                        #     print p
                        #     print protein
                        #     exit(7)
                        # if isMonoMass:
                        #     pm = pep.mass
                        # else:
                        #     pm = pep.average_mass
                        # if pm < min_mass or pm > max_mass:
                        #     continue

                        # # modifications
                        # if p[0] in nterm_mods:
                        #     pm += nterm_mods[p[0]]
                        # elif p[0] in mods:
                        #     pm += mods[p[0]]
                        # if p[-1] in cterm_mods:
                        #     pm += cterm_mods[p[-1]]
                        # elif p[-1] in mods:
                        #     pm += mods[p[-1]]
                        # if mods:
                        #     for aa in p[1:-1]:
                        #         if aa in mods:
                        #             pm += mods[aa]

                        # # get flanking information
                        # if i == 0:
                        #     nterm_flanking = '-'
                        # else:
                        #     nterm_flanking = h[i-1][-1]
                        # if i == hl:
                        #     cterm_flanking = '-'
                        # else:
                        #     cterm_flanking = h[i+1][0]
                        # peptides[p] = digestedPep(pep, protein_name,
                        #                           nterm_flanking, 
                        #                           cterm_flanking, 
                        #                           pm)
    return peptides

def load_proteins_transStop(fastaFile, digest_regexp,
                            n_res, c_res,
                            min_length, max_length,
                            min_mass, max_mass,
                            mods, nterm_mods, cterm_mods,
                            max_mods, min_mods, missed_cleavages = 0,
                            isMonoMass = True, 
                            digestion = 'full-digest'):
    """ Parse fasta file and digest proteins into peptides
    inputs:
    fastaFile - FASTA file containing proteins
    digest_regexp - re compiled regular expression dictating digestion rules
    min_length - minimum residue length
    max_length - maximum residue length

    output:
    peptides - dict whose keys are the digested peptide string, elements 
               are instances of digestedPep class with fields:
               peptide - instance of pyFiles.peptide's Peptide
               proteinName - name of parent protein from Fasta file
               ntermFlanking - nterm (left flanking) amino acid
               ctermFlanking - cterm (right flanking) amino acid
               peptideMass - calculated peptide mass

    todo: add missed cleavages, static/variable modification masses
    """
    peptides = {}
    
    with open(fastaFile) as f:
        for isname, group in itertools.groupby(f, lambda x: x[0]=='>'):
            if isname: # descriptor line
                l = group.next()
                protein_name = l[1:].split()[0]
            else: # protein sequence lines
                protein = ''.join(l.strip() for l in group)
                # protein = protein.replace('*', '')
                # h = digest_regexp.findall(protein)
                split_col = []
                for cutProt in protein.split('*'):
                    h = digest_regexp.findall(cutProt)
                    if h:
                        split_col.append(digest_regexp.findall(cutProt))

                numTranslations = len(split_col)
                for currSplit, h in enumerate(split_col):

                    hl = len(h)-1
                    if len(h)==1:
                        p = h[0]
                        nterm_flanking = '-'
                        cterm_flanking = '-'
                        if not currSplit:
                            valid_p = check_peptide(p, min_length, max_length,
                                                    protein_name, isMonoMass,
                                                    min_mass, max_mass, 
                                                    nterm_mods, cterm_mods, mods,
                                                    nterm_flanking, cterm_flanking)
                            if valid_p:
                                peptides[p] = valid_p
                        if digestion == 'partial-digest':
                            # prefix strings
                            prefix_nterm_flanking = nterm_flanking # n-terminus doesn't change
                            prefix_nterm_aa = p[0]
                            for ind in range(min_length, len(p)+1): # p[ind:-1]
                                if ind == len(p):
                                    prefix_cterm_flanking = cterm_flanking
                                else:
                                    prefix_cterm_flanking = p[ind]
                                prefix_cterm_aa = p[ind-1] # handle case where user sets min-length to 1, in which
                                                                  # case an AA is both it's c-terminus and n-terminus
                                if currSplit == 0\
                                        or (prefix_nterm_flanking in n_res and prefix_nterm_aa in c_res) \
                                        or (prefix_cterm_aa in n_res and prefix_cterm_flanking in c_res):
                                    valid_p = check_peptide(p[:ind], min_length, max_length,
                                                            protein_name, isMonoMass,
                                                            min_mass, max_mass, 
                                                            nterm_mods, cterm_mods, mods,
                                                            prefix_nterm_flanking, prefix_cterm_flanking)
                                    if valid_p:
                                        peptides[p[:ind]] = valid_p
                            # suffix strings
                            suffix_cterm_flanking = cterm_flanking
                            suffix_cterm_aa = p[-1]
                            for ind in range(0,len(p) - min_length+1): # p[0:ind]
                                suffix_nterm_aa = p[ind]
                                if ind > 0:
                                    suffix_nterm_flanking = p[ind-1]
                                else:
                                    suffix_nterm_flanking = nterm_flanking
                                if currSplit== numTranslations-1 \
                                        or (suffix_cterm_flanking in n_res and suffix_nterm_aa in c_res) \
                                        or (suffix_cterm_aa in n_res and suffix_cterm_flanking in c_res):
                                    valid_p = check_peptide(p[ind:], min_length, max_length,
                                                            protein_name, isMonoMass,
                                                            min_mass, max_mass, 
                                                            nterm_mods, cterm_mods, mods,
                                                            suffix_nterm_flanking, suffix_cterm_flanking)
                                    if valid_p:
                                        peptides[p[ind:]] = valid_p
                    else:
                        for i in range(len(h)):
                            p = ''
                            for j in range(i, min(hl, i+missed_cleavages)+1):
                                p += h[j]
                                # get flanking information
                                if i == 0:
                                    nterm_flanking = '-'
                                else:
                                    nterm_flanking = h[i-1][-1]
                                if j == hl:
                                    cterm_flanking = '-'
                                else:
                                    cterm_flanking = h[j+1][0]

                                if currSplit:
                                    if p[-1] not in n_res:
                                        continue
                                # if currSplit != numTranslations-1:
                                #     if p[-1] not in n_res:
                                #         continue
                                valid_p = check_peptide(p, min_length, max_length,
                                                        protein_name, isMonoMass,
                                                        min_mass, max_mass, 
                                                        nterm_mods, cterm_mods, mods,
                                                        nterm_flanking, cterm_flanking)
                                if valid_p:
                                    peptides[p] = valid_p

                            if digestion == 'partial-digest':
                                # p is now the longest substring including all allowed missed cleavages
                                # prefix strings
                                prefix_nterm_flanking = nterm_flanking # n-terminus doesn't change
                                prefix_nterm_aa = p[0]
                                for ind in range(min_length, len(p)+1): # p[0:ind]
                                    if ind == len(p):
                                        prefix_cterm_flanking = cterm_flanking
                                    else:
                                        prefix_cterm_flanking = p[ind]
                                    prefix_cterm_aa = p[ind-1] # handle case where user sets min-length to 1, in which

                                    # case an AA is both it's c-terminus and n-terminus
                                    if currSplit == 0 \
                                            or (prefix_nterm_flanking in n_res and prefix_nterm_aa in c_res) \
                                            or (prefix_cterm_aa in n_res and prefix_cterm_flanking in c_res):
                                        valid_p = check_peptide(p[:ind], min_length, max_length,
                                                                protein_name, isMonoMass,
                                                                min_mass, max_mass, 
                                                                nterm_mods, cterm_mods, mods,
                                                                prefix_nterm_flanking, prefix_cterm_flanking)
                                        if valid_p:
                                            peptides[p[:ind]] = valid_p
                                # suffix strings
                                suffix_cterm_flanking = cterm_flanking
                                suffix_cterm_aa = p[-1]
                                for ind in range(0,len(p) - min_length+1): # p[ind:]
                                    suffix_nterm_aa = p[ind]
                                    if ind > 0:
                                        suffix_nterm_flanking = p[ind-1]
                                    else:
                                        suffix_nterm_flanking = nterm_flanking

                                    # if p[ind:] == 'TFPLVHSHLK':
                                    #     print n_res
                                    #     print c_res
                                    #     print suffix_nterm_flanking
                                    #     print suffix_nterm_aa
                                    #     print suffix_cterm_aa
                                    #     print suffix_cterm_flanking
                                    #     print currSplit
                                    #     print i
                                    #     print j
                                    #     print hl
                                    #     print h[j+1]
                                    #     print h[j+1][0]

                                    if currSplit == numTranslations -1 \
                                            or (suffix_cterm_flanking in n_res and suffix_nterm_aa in c_res) \
                                            or (suffix_cterm_aa in n_res and suffix_cterm_flanking in c_res):
                                        valid_p = check_peptide(p[ind:], min_length, max_length,
                                                                protein_name, isMonoMass,
                                                                min_mass, max_mass, 
                                                                nterm_mods, cterm_mods, mods,
                                                                suffix_nterm_flanking, suffix_cterm_flanking)
                                        if valid_p:
                                            peptides[p[ind:]] = valid_p
    return peptides

def ooc_write_merge_peptides(fastaFile, digest_regexp,
                             output,
                             n_res, c_res,
                             min_length, max_length,
                             min_mass, max_mass,
                             mods, nterm_mods, cterm_mods,
                             max_mods, min_mods, missed_cleavages = 0,
                             isMonoMass = True, 
                             digestion = 'full-digest', 
                             peptide_buffer = 100000):

    curr = 0
    # todo: map protein name to integer
    s = struct.Struct('f %ds I s s' % (max_length))

    peptides = []
    iters = []
    num_buffered = 0
    pb_bound = peptide_buffer-1

    for peptide in load_proteins_generator(fastaFile, digest_regexp,
                                            n_res, c_res,
                                            min_length, max_length,
                                            min_mass, max_mass,
                                            mods, nterm_mods, cterm_mods,
                                            max_mods, min_mods, missed_cleavages,
                                            isMonoMass, 
                                            digestion):
        if num_buffered < pb_bound:
            peptides.append(peptide)
            num_buffered += 1
        else:
            peptides.append(peptide)
            # sort in memory
            # for peptide in peptides:
            #    peptide[0] = peptide mass (float)
            #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
            #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
            #    peptide[3] = nterm_flanking (character)
            #    peptide[4] = cterm_flanking (character)
            peptides.sort(key = lambda r: r[0]) # sort tuples based on peptide mass
            # write to temporary file
            f = tempfile.TemporaryFile()
            write_binary_peptide_tuples(s, peptides, f)
            f.seek(0)

            iters.append(read_binary_peptide_tuples_buffer(s, f))

            del peptides[:]
            num_buffered = 0

    if num_buffered:
        peptides.sort(key = lambda r: r[0]) # sort tuples based on peptide mass
        # write to temporary file
        f = tempfile.TemporaryFile()
        write_binary_peptide_tuples(s, peptides, f)
        f.seek(0)
        
        iters.append(read_binary_peptide_tuples_buffer(s, f))

        del peptides[:]
        

    peptides = []
    write_buffer = peptide_buffer

    with open(output, 'wb') as f:
        for p in heapq.merge(*iters):
            peptides.append(p)
            if len(peptides) >= write_buffer:
                write_binary_peptide_tuples(s, peptides, f)
                del peptides[:]
        if peptides:
            write_binary_peptide_tuples(s, peptides, f)

def ooc_write_merge_peptides_var_mods(fastaFile, digest_regexp,
                                      digest_dir,
                                      n_res, c_res,
                                      min_length, max_length,
                                      min_mass, max_mass,
                                      mods, nterm_mods, cterm_mods,
                                      var_mods, nterm_var_mods, cterm_var_mods,
                                      max_mods, min_mods, missed_cleavages = 0,
                                      isMonoMass = True, 
                                      digestion = 'full-digest', 
                                      outFileName = 'targets.bin', 
                                      peptide_buffer = 100000):

    curr = 0
    # todo: map protein name to integer
    s = struct.Struct('f %ds I s s %ds' % (max_length, max_length))

    peptides = []
    iters = []
    num_buffered = 0
    pb_bound = peptide_buffer-1

    for peptide in load_proteins_generator_var_mods(fastaFile, digest_regexp,
                                                    n_res, c_res,
                                                    min_length, max_length,
                                                    min_mass, max_mass,
                                                    mods, nterm_mods, cterm_mods,
                                                    var_mods, nterm_var_mods, cterm_var_mods,
                                                    max_mods, min_mods, missed_cleavages,
                                                    isMonoMass, digestion):
        if num_buffered < pb_bound:
            peptides.append(peptide)
            num_buffered += 1
        else:
            peptides.append(peptide)
            # sort in memory
            # for peptide in peptides:
            #    peptide[0] = peptide mass (float)
            #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
            #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
            #    peptide[3] = nterm_flanking (character)
            #    peptide[4] = cterm_flanking (character)
            #    peptide[5] = list of bools
            peptides.sort(key = lambda r: r[0]) # sort tuples based on peptide mass
            # write to temporary file
            # curr_f = os.path.join(digest_dir, 'peps-%d.bin' % curr)
            f = tempfile.TemporaryFile()
            # with open(file_list, 'wb') as f:
            write_binary_peptide_tuples(s, peptides, f)
            f.seek(0)

            iters.append(read_binary_peptide_tuples_buffer(s, f))

            del peptides[:]
            num_buffered = 0

    if num_buffered:
        peptides.sort(key = lambda r: r[0]) # sort tuples based on peptide mass
        # write to temporary file
        f = tempfile.TemporaryFile()
        write_binary_peptide_tuples(s, peptides, f)
        f.seek(0)
        
        iters.append(read_binary_peptide_tuples_buffer(s, f))

        del peptides[:]
        

    peptides = []
    write_buffer = peptide_buffer

    output = os.path.join(digest_dir, outFileName)
    with open(output, 'wb') as f:
        for p in heapq.merge(*iters):
            peptides.append(p)
            if len(peptides) >= write_buffer:
                write_binary_peptide_tuples(s, peptides, f)
                del peptides[:]
        if peptides:
            write_binary_peptide_tuples(s, peptides, f)

def digested_peptides_var_mods(targets, output,
                               min_length, max_length,
                               min_mass, max_mass,
                               var_mods, nterm_var_mods, cterm_var_mods,
                               max_mods, min_mods,
                               peptide_buffer = 100000):
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    #    peptide[5] = sequence denoting variable modifications
    s = struct.Struct('f %ds I s s %ds' % (max_length, max_length))

    curr = 0
    peptides = []
    iters = []
    num_buffered = 0
    pb_bound = peptide_buffer-1

    for peptide in load_digested_peptides_var_mods_generator(targets, 
                                                             min_length, max_length,
                                                             min_mass, max_mass,
                                                             var_mods, nterm_var_mods, cterm_var_mods,
                                                             max_mods, min_mods):
        # current fields: p[0] = peptide mass, p[1] = peptide string, p[2] = protein number, p[3] = nterm_flanking, p[4] = cterm_flanking
        if num_buffered < pb_bound:
            peptides.append(peptide)
            num_buffered += 1
        else:
            peptides.append(peptide)
            # sort in memory
            # for peptide in peptides:
            #    peptide[0] = peptide mass (float)
            #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
            #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
            #    peptide[3] = nterm_flanking (character)
            #    peptide[4] = cterm_flanking (character)
            #    peptide[5] = list of bools
            peptides.sort(key = lambda r: r[0]) # sort tuples based on peptide mass
            # write to temporary file
            # curr_f = os.path.join(digest_dir, 'peps-%d.bin' % curr)
            f = tempfile.TemporaryFile()
            # with open(file_list, 'wb') as f:
            write_binary_peptide_tuples(s, peptides, f)
            f.seek(0)

            iters.append(read_binary_peptide_tuples_buffer(s, f))

            del peptides[:]
            num_buffered = 0

    if num_buffered:
        peptides.sort(key = lambda r: r[0]) # sort tuples based on peptide mass
        # write to temporary file
        f = tempfile.TemporaryFile()
        write_binary_peptide_tuples(s, peptides, f)
        f.seek(0)
        
        iters.append(read_binary_peptide_tuples_buffer(s, f))

        del peptides[:]
        

    peptides = []
    write_buffer = peptide_buffer

    with open(output, 'wb') as f:
        for p in heapq.merge(*iters):
            peptides.append(p)
            if len(peptides) >= write_buffer:
                write_binary_peptide_tuples(s, peptides, f)
                del peptides[:]
        if peptides:
            write_binary_peptide_tuples(s, peptides, f)

def load_digested_peptides_var_mods_generator(targets, 
                                              min_length, max_length,
                                              min_mass, max_mass,
                                              var_mods, nterm_var_mods, cterm_var_mods,
                                              max_mods, min_mods):
    """
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    s = struct.Struct('f %ds I s s' % (max_length))

    with open(targets, 'rb') as fid:
        for peptide in read_binary_peptide_tuples(s, fid):
            static_mod_mass = peptide[0]
            p = peptide[1].split('\x00')[0]
            protein_num = peptide[2]
            nterm_flanking = peptide[3]
            cterm_flanking = peptide[4]
            valid_peps = check_digested_peptide_tuple_var_mods(p, static_mod_mass, 
                                                               min_length, max_length,
                                                               min_mass, max_mass, 
                                                               nterm_var_mods, cterm_var_mods, var_mods,
                                                               max_mods, min_mods)
            # (var_subset, pm)
            if valid_peps:
                for mod_pep in valid_peps:
                    modded_aas = ['0'] * max_length
                    for ma in mod_pep[0]:
                        if ma == -1: # variable n-term modification
                            modded_aas[0] = '2'
                        elif ma == -2: # variable c-term modification
                            modded_aas[len(p)-1] = '3'
                        else: # variable mod
                            modded_aas[ma] = '1'
                                        
                    yield (mod_pep[1], p, protein_num, nterm_flanking, cterm_flanking, ''.join(modded_aas))

def load_proteins_generator_var_mods(fastaFile, digest_regexp, 
                                     n_res, c_res,
                                     min_length, max_length,
                                     min_mass, max_mass,
                                     mods, nterm_mods, cterm_mods,
                                     var_mods, nterm_var_mods, cterm_var_mods,
                                     max_mods, min_mods,
                                     missed_cleavages = 0,
                                     isMonoMass = True, 
                                     digestion = 'full-digest'):
    """ Parse fasta file and digest proteins into peptides
    inputs:
    fastaFile - FASTA file containing proteins
    digest_regexp - re compiled regular expression dictating digestion rules
    min_length - minimum residue length
    max_length - maximum residue length

    output:
    peptides - dict whose keys are the digested peptide string, elements 
               are instances of digestedPep class with fields:
               peptide - instance of pyFiles.peptide's Peptide
               proteinName - name of parent protein from Fasta file
               ntermFlanking - nterm (left flanking) amino acid
               ctermFlanking - cterm (right flanking) amino acid
               peptideMass - calculated peptide mass

    todo: add missed cleavages, static/variable modification masses
    """

    protein_num = -1

    with open(fastaFile) as f:
        for isname, group in itertools.groupby(f, lambda x: x[0]=='>'):
            protein_num += 1
            if isname: # descriptor line
                l = group.next()
                protein_name = l[1:].split()[0]
            else: # protein sequence lines
                protein = ''.join(l.strip() for l in group)
                if protein[-1] == '*':
                    h = digest_regexp.findall(protein[:-1])
                else:
                    h = digest_regexp.findall(protein)

                if h:
                    hl = len(h)-1
                    if len(h)==1:
                        p = h[0]
                        nterm_flanking = '-'
                        cterm_flanking = '-'
                        valid_peps = check_peptide_tuple_var_mods(p, min_length, max_length,
                                                                  isMonoMass, min_mass, max_mass, 
                                                                  nterm_mods, cterm_mods, mods,
                                                                  nterm_var_mods, cterm_var_mods, var_mods,
                                                                  max_mods, min_mods)

                        # (var_subset, pm)
                        if valid_peps:
                            for mod_pep in valid_peps:
                                modded_aas = ['0'] * max_length
                                for ma in mod_pep[0]:
                                    if ma == -1: # variable n-term modification
                                        modded_aas[0] = '2'
                                    elif ma == -2: # variable c-term modification
                                        modded_aas[len(p)-1] = '3'
                                    else: # variable mod
                                        modded_aas[ma] = '1'
                                        
                                yield (mod_pep[1], p, protein_num, nterm_flanking, cterm_flanking, ''.join(modded_aas))

                        if digestion == 'partial-digest':
                            # prefix strings
                            prefix_nterm_flanking = nterm_flanking # n-terminus doesn't change
                            prefix_nterm_aa = p[0]
                            for ind in range(min_length, len(p)+1): # p[ind:-1]
                                if ind == len(p):
                                    prefix_cterm_flanking = cterm_flanking
                                else:
                                    prefix_cterm_flanking = p[ind]
                                prefix_cterm_aa = p[ind-1] # handle case where user sets min-length to 1, in which
                                                                  # case an AA is both it's c-terminus and n-terminus
                                pep_str = p[:ind]
                                valid_peps = check_peptide_tuple_var_mods(pep_str, min_length, max_length,
                                                                          isMonoMass, min_mass, max_mass, 
                                                                          nterm_mods, cterm_mods, mods,
                                                                          nterm_var_mods, cterm_var_mods, var_mods,
                                                                          max_mods, min_mods)
                                # (var_subset, pm)
                                if valid_peps:
                                    for mod_pep in valid_peps:
                                        modded_aas = ['0'] * max_length
                                        for ma in mod_pep[0]:
                                            if ma == -1: # variable n-term modification
                                                modded_aas[0] = '2'
                                            elif ma == -2: # variable c-term modification
                                                modded_aas[len(pep_str)-1] = '3'
                                            else:
                                                modded_aas[ma] = '1'
                                        yield (mod_pep[1], pep_str, protein_num, prefix_nterm_flanking, prefix_cterm_flanking, ''.join(modded_aas))

                            # suffix strings
                            suffix_cterm_flanking = cterm_flanking
                            suffix_cterm_aa = p[-1]
                            for ind in range(0,len(p) - min_length+1): # p[0:ind]
                                suffix_nterm_aa = p[ind]
                                if ind > 0:
                                    suffix_nterm_flanking = p[ind-1]
                                else:
                                    suffix_nterm_flanking = nterm_flanking
                                    
                                pep_str = p[ind:]
                                valid_peps = check_peptide_tuple_var_mods(pep_str, min_length, max_length,
                                                                          isMonoMass, min_mass, max_mass, 
                                                                          nterm_mods, cterm_mods, mods,
                                                                          nterm_var_mods, cterm_var_mods, var_mods,
                                                                          max_mods, min_mods)
                                # (var_subset, pm)
                                if valid_peps:
                                    for mod_pep in valid_peps:
                                        modded_aas = ['0'] * max_length
                                        for ma in mod_pep[0]:
                                            if ma == -1: # variable n-term modification
                                                modded_aas[0] = '2'
                                            elif ma == -2: # variable c-term modification
                                                modded_aas[len(pep_str)-1] = '3'
                                            else:
                                                modded_aas[ma] = '1'
                                        yield (mod_pep[1], pep_str, protein_num, suffix_nterm_flanking, suffix_cterm_flanking, ''.join(modded_aas))

                    else:
                        for i in range(len(h)):
                            p = ''
                            for j in range(i, min(hl, i+missed_cleavages)+1):
                                p += h[j]
                                # get flanking information
                                if i == 0:
                                    nterm_flanking = '-'
                                else:
                                    nterm_flanking = h[i-1][-1]
                                if j == hl:
                                    cterm_flanking = '-'
                                else:
                                    cterm_flanking = h[j+1][0]

                                valid_peps = check_peptide_tuple_var_mods(p, min_length, max_length,
                                                                          isMonoMass, min_mass, max_mass, 
                                                                          nterm_mods, cterm_mods, mods,
                                                                          nterm_var_mods, cterm_var_mods, var_mods,
                                                                          max_mods, min_mods)
                                # (var_subset, pm)
                                if valid_peps:
                                    for mod_pep in valid_peps:
                                        modded_aas = ['0'] * max_length
                                        for ma in mod_pep[0]:
                                            if ma == -1: # variable n-term modification
                                                modded_aas[0] = '2'
                                            elif ma == -2: # variable c-term modification
                                                modded_aas[len(p)-1] = '3'
                                            else:
                                                modded_aas[ma] = '1'

                                        yield (mod_pep[1], p, protein_num, nterm_flanking, cterm_flanking, ''.join(modded_aas))

                            if digestion == 'partial-digest':
                                # p is now the longest substring including all allowed missed cleavages
                                # prefix strings
                                prefix_nterm_flanking = nterm_flanking # n-terminus doesn't change
                                prefix_nterm_aa = p[0]
                                for ind in range(min_length, len(p)+1): # p[0:ind]
                                    if ind == len(p):
                                        prefix_cterm_flanking = cterm_flanking
                                    else:
                                        prefix_cterm_flanking = p[ind]
                                    prefix_cterm_aa = p[ind-1] # handle case where user sets min-length to 1, in which

                                    pep_str = p[:ind]
                                    # case an AA is both it's c-terminus and n-terminus
                                    if  i==0 \
                                            or (prefix_nterm_flanking in n_res and prefix_nterm_aa in c_res) \
                                            or (prefix_cterm_aa in n_res and prefix_cterm_flanking in c_res):
                                        valid_peps = check_peptide_tuple_var_mods(pep_str, min_length, max_length,
                                                                                  isMonoMass, min_mass, max_mass, 
                                                                                  nterm_mods, cterm_mods, mods,
                                                                                  nterm_var_mods, cterm_var_mods, var_mods,
                                                                                  max_mods, min_mods)
                                        # (var_subset, pm)
                                        if valid_peps:
                                            for mod_pep in valid_peps:
                                                modded_aas = ['0'] * max_length
                                                for ma in mod_pep[0]:
                                                    if ma == -1: # variable n-term modification
                                                        modded_aas[0] = '2'
                                                    elif ma == -2: # variable c-term modification
                                                        modded_aas[len(pep_str)-1] = '3'
                                                    else:
                                                        modded_aas[ma] = '1'
                                                yield (mod_pep[1], pep_str, protein_num, prefix_nterm_flanking, prefix_cterm_flanking, ''.join(modded_aas))

                                # suffix strings
                                suffix_cterm_flanking = cterm_flanking
                                suffix_cterm_aa = p[-1]
                                for ind in range(0,len(p) - min_length+1): # p[ind:]
                                    suffix_nterm_aa = p[ind]
                                    if ind > 0:
                                        suffix_nterm_flanking = p[ind-1]
                                    else:
                                        suffix_nterm_flanking = nterm_flanking

                                    pep_str = p[ind:]
                                    if i==hl \
                                    or (suffix_cterm_flanking in n_res and suffix_nterm_aa in c_res) \
                                            or (suffix_cterm_aa in n_res and suffix_cterm_flanking in c_res):
                                        valid_peps = check_peptide_tuple_var_mods(pep_str, min_length, max_length,
                                                                                  isMonoMass, min_mass, max_mass, 
                                                                                  nterm_mods, cterm_mods, mods,
                                                                                  nterm_var_mods, cterm_var_mods, var_mods,
                                                                                  max_mods, min_mods)
                                        # (var_subset, pm)
                                        if valid_peps:
                                            for mod_pep in valid_peps:
                                                modded_aas = ['0'] * max_length
                                                for ma in mod_pep[0]:
                                                    if ma == -1: # variable n-term modification
                                                        modded_aas[0] = '2'
                                                    elif ma == -2: # variable c-term modification
                                                        modded_aas[len(pep_str)-1] = '3'
                                                    else:
                                                        modded_aas[ma] = '1'
                                                yield (mod_pep[1], pep_str, protein_num, suffix_nterm_flanking, suffix_cterm_flanking, ''.join(modded_aas))

def load_proteins_generator(fastaFile, digest_regexp, 
                            n_res, c_res,
                            min_length, max_length,
                            min_mass, max_mass,
                            mods, nterm_mods, cterm_mods,
                            max_mods, min_mods,
                            missed_cleavages = 0,
                            isMonoMass = True, 
                            digestion = 'full-digest'):
    """ Parse fasta file and digest proteins into peptides
    inputs:
    fastaFile - FASTA file containing proteins
    digest_regexp - re compiled regular expression dictating digestion rules
    min_length - minimum residue length
    max_length - maximum residue length

    output:
    peptides - dict whose keys are the digested peptide string, elements 
               are instances of digestedPep class with fields:
               peptide - instance of pyFiles.peptide's Peptide
               proteinName - name of parent protein from Fasta file
               ntermFlanking - nterm (left flanking) amino acid
               ctermFlanking - cterm (right flanking) amino acid
               peptideMass - calculated peptide mass

    todo: add missed cleavages, static/variable modification masses
    """

    protein_num = -1

    with open(fastaFile) as f:
        for isname, group in itertools.groupby(f, lambda x: x[0]=='>'):
            protein_num += 1
            if isname: # descriptor line
                l = group.next()
                protein_name = l[1:].split()[0]
            else: # protein sequence lines
                protein = ''.join(l.strip() for l in group)
                if protein[-1] == '*':
                    h = digest_regexp.findall(protein[:-1])
                else:
                    h = digest_regexp.findall(protein)

                if h:
                    hl = len(h)-1
                    if len(h)==1:
                        p = h[0]
                        nterm_flanking = '-'
                        cterm_flanking = '-'
                        valid_p = check_peptide_tuple(p, min_length, max_length,
                                                      protein_num, isMonoMass,
                                                      min_mass, max_mass, 
                                                      nterm_mods, cterm_mods, mods,
                                                      nterm_flanking, cterm_flanking)
                        if valid_p:
                            yield valid_p
                                
                        if digestion == 'partial-digest':
                            # prefix strings
                            prefix_nterm_flanking = nterm_flanking # n-terminus doesn't change
                            prefix_nterm_aa = p[0]
                            for ind in range(min_length, len(p)+1): # p[ind:-1]
                                if ind == len(p):
                                    prefix_cterm_flanking = cterm_flanking
                                else:
                                    prefix_cterm_flanking = p[ind]
                                prefix_cterm_aa = p[ind-1] # handle case where user sets min-length to 1, in which
                                                                  # case an AA is both it's c-terminus and n-terminus
                                valid_p = check_peptide_tuple(p[:ind], min_length, max_length,
                                                              protein_num, isMonoMass,
                                                              min_mass, max_mass, 
                                                              nterm_mods, cterm_mods, mods,
                                                              nterm_flanking, cterm_flanking)
                                if valid_p:
                                    yield valid_p
                            # suffix strings
                            suffix_cterm_flanking = cterm_flanking
                            suffix_cterm_aa = p[-1]
                            for ind in range(0,len(p) - min_length+1): # p[0:ind]
                                suffix_nterm_aa = p[ind]
                                if ind > 0:
                                    suffix_nterm_flanking = p[ind-1]
                                else:
                                    suffix_nterm_flanking = nterm_flanking
                                valid_p = check_peptide_tuple(p[ind:], min_length, max_length,
                                                              protein_num, isMonoMass,
                                                              min_mass, max_mass, 
                                                              nterm_mods, cterm_mods, mods,
                                                              suffix_nterm_flanking, suffix_cterm_flanking)
                                if valid_p:
                                    yield valid_p
                    else:
                        for i in range(len(h)):
                            p = ''
                            for j in range(i, min(hl, i+missed_cleavages)+1):
                                p += h[j]
                                # get flanking information
                                if i == 0:
                                    nterm_flanking = '-'
                                else:
                                    nterm_flanking = h[i-1][-1]
                                if j == hl:
                                    cterm_flanking = '-'
                                else:
                                    cterm_flanking = h[j+1][0]

                                valid_p = check_peptide_tuple(p, min_length, max_length,
                                                                protein_num, isMonoMass,
                                                                min_mass, max_mass, 
                                                                nterm_mods, cterm_mods, mods,
                                                                nterm_flanking, cterm_flanking)
                                if valid_p:
                                    yield valid_p

                            if digestion == 'partial-digest':
                                # p is now the longest substring including all allowed missed cleavages
                                # prefix strings
                                prefix_nterm_flanking = nterm_flanking # n-terminus doesn't change
                                prefix_nterm_aa = p[0]
                                for ind in range(min_length, len(p)+1): # p[0:ind]
                                    if ind == len(p):
                                        prefix_cterm_flanking = cterm_flanking
                                    else:
                                        prefix_cterm_flanking = p[ind]
                                    prefix_cterm_aa = p[ind-1] # handle case where user sets min-length to 1, in which

                                    # case an AA is both it's c-terminus and n-terminus
                                    if  i==0 \
                                            or (prefix_nterm_flanking in n_res and prefix_nterm_aa in c_res) \
                                            or (prefix_cterm_aa in n_res and prefix_cterm_flanking in c_res):
                                        valid_p = check_peptide_tuple(p[:ind], min_length, max_length,
                                                                      protein_num, isMonoMass,
                                                                      min_mass, max_mass, 
                                                                      nterm_mods, cterm_mods, mods,
                                                                      prefix_nterm_flanking, prefix_cterm_flanking)
                                        if valid_p:
                                            yield valid_p
                                # suffix strings
                                suffix_cterm_flanking = cterm_flanking
                                suffix_cterm_aa = p[-1]
                                for ind in range(0,len(p) - min_length+1): # p[ind:]
                                    suffix_nterm_aa = p[ind]
                                    if ind > 0:
                                        suffix_nterm_flanking = p[ind-1]
                                    else:
                                        suffix_nterm_flanking = nterm_flanking

                                    if i==hl \
                                    or (suffix_cterm_flanking in n_res and suffix_nterm_aa in c_res) \
                                            or (suffix_cterm_aa in n_res and suffix_cterm_flanking in c_res):
                                        valid_p = check_peptide_tuple(p[ind:], min_length, max_length,
                                                                      protein_num, isMonoMass,
                                                                      min_mass, max_mass, 
                                                                      nterm_mods, cterm_mods, mods,
                                                                      suffix_nterm_flanking, suffix_cterm_flanking)
                                        if valid_p:
                                            yield valid_p

def write_peptides(f, peptides):
    # first sort peptides by mass
    peps = [peptides[pep] for pep in peptides]
    peps.sort(key = lambda x: x.peptideMass)

    # write digested peptides
    fid = open(f, 'w')
    fid.write("Protein\tPeptide\tMass\tNterm_flanking\tCterm_flanking\n")

    for p in peps:
        fid.write("%s\t%s\t%f\t%c\t%c\n" % 
                  (p.proteinName, p.peptide, 
                   p.peptideMass, p.ntermFlanking, p.ctermFlanking))
    fid.close()

def filter_double_binarydb(output, targets, max_length, peptide_buffer = 100000):
    """ pre: ooc_write_merge_peptides has been called and peptides are written
             in binary format according to struct.Struct instance below.
             Note that after ooc_write_merge_peptides, peptides should be
             sorted by mass
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    s = struct.Struct('f %ds I s s' % (max_length))

    peptides = []
    num_buffered = 0

    prev_mass = -1.0

    eq_mass_set = set([])

    # not as simple as looking at two adjacent peptides, since several peptides could have the same mass
    # and lead to a sequence of peptides, each with the same mass, with redundant peptides

    with open(output, 'wb') as f:
        with open(targets, 'rb') as fid:
            for p in read_binary_peptide_tuples(s, fid):
                if p[0] == prev_mass:
                    if p[1] not in eq_mass_set:
                        peptides.append(p)
                        eq_mass_set.add(p[1])
                        num_buffered += 1
                else:
                    peptides.append(p)
                    num_buffered += 1
                    prev_mass = p[0]
                    eq_mass_set = set([p[1]])

                if num_buffered >= peptide_buffer:
                    write_binary_peptide_tuples(s, peptides, f)
                    del peptides[:]
                    num_buffered = 0

            if num_buffered:
                write_binary_peptide_tuples(s, peptides, f)

def filter_double_binarydb_var_mods(output, targets, max_length,
                                    peptide_buffer = 100000):
    """ pre: ooc_write_merge_peptides has been called and peptides are written
             in binary format according to struct.Struct instance below.
             Note that after ooc_write_merge_peptides, peptides should be
             sorted by mass
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    #    peptide[5] = boolean sequence of variable modificiations
    # s = struct.Struct('f %ds I s s' % (max_length))
    s = struct.Struct('f %ds I s s %ds' % (max_length, max_length))

    peptides = []
    num_buffered = 0

    prev_mass = -1.0

    eq_mass_set = set([])

    # not as simple as looking at two adjacent peptides, since several peptides could have the same mass
    # and lead to a sequence of peptides, each with the same mass, with redundant peptides.
    # with the advent of variable mods, looking for doubles is significantly more complicated.  now,
    # we construct a dictionary whose key's keys are pairs (peptide string, peptide mass) and whose 
    # entries are the boolean sequences of variable modifications

    with open(output, 'wb') as f:
        with open(targets, 'rb') as fid:
            for p in read_binary_peptide_tuples(s, fid):
                # boolVarMods = boolListToString(p[5])
                if p[0] == prev_mass:
                    # construct string based on binary sequence of variable mods
                    if (p[1], p[5]) not in eq_mass_set: # dictionary key is (peptide string, peptide mass)
                        peptides.append(p)
                        eq_mass_set.add((p[1], p[5]))
                        num_buffered += 1
                else:
                    peptides.append(p)
                    num_buffered += 1
                    prev_mass = p[0]
                    eq_mass_set = set([(p[1], p[5])])

                if num_buffered >= peptide_buffer:
                    write_binary_peptide_tuples(s, peptides, f)
                    del peptides[:]
                    num_buffered = 0

            if num_buffered:
                write_binary_peptide_tuples(s, peptides, f)

def create_decoys_from_binarydb(output, targets, max_length,
                                decoy_format, keep_terminal_aminos, s0, peptide_buffer = 100000):
    """ pre: ooc_write_merge_peptides has been called and peptides are written
             in binary format according to struct.Struct instance below.
             Note that after ooc_write_merge_peptides, peptides should be
             sorted by mass
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    s = struct.Struct('f %ds I s s' % (max_length))

    if s0.lower == 'time':
        s0 = None
    else:
        try:
            s0 = int(s0)
        except:
            print "Seed must be either an unsigned integer or time, exitting"
            exit(9)
    # seed random number generator
    random.seed(s0)

    # check terminal amino acids option
    kta = keep_terminal_aminos
    if kta != 'N' and kta != 'C' and kta != 'NC' and kta.lower != 'none':
        print "keep-terminal-aminos value supplied is %s, must be one of N, C, NC, or none.  Exitting" % kta
        exit(10)

    peptides = []
    num_buffered = 0

    prev_mass = 0

    # fill eq_mass_dict with targets of equal mass.  when we move to the next mass, shuffle targets in 
    # eq_mass_dict
    eq_mass_dict = {}

    # not as simple as looking at two adjacent peptides, since several peptides could have the same mass
    # and lead to a sequence of peptides, each with the same mass, with redundant peptides

    decoy_format = decoy_format.lower()

    if decoy_format != 'shuffle' and decoy_format != 'reverse':
        print "valid options for decoy-format are shuffle or reverse, exitting."
        exit(11)

    curr_visited = 0
    decoy_eq_mass_set = set([])
    with open(output, 'wb') as f:
        with open(targets, 'rb') as fid:
            for p in read_binary_peptide_tuples(s, fid):
                pep_str = p[1].split('\x00')[0] # remove null characters
                if p[0] == prev_mass:
                    # assume only unique peptides in target database (i.e., filter function has removed all duplicates)
                    eq_mass_dict[pep_str] = p
                else: # we've filled up the queue of targets to be shuffled/reversed with equal mass
                    decoy_eq_mass_set.clear()
                    for t in eq_mass_dict:
                        if decoy_format == 'shuffle':
                            d = create_shuffled_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                        elif decoy_format == 'reverse':
                            d = create_reversed_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                        if d:
                            decoy_eq_mass_set.add(d[1])
                            peptides.append(d)
                            num_buffered +=1

                            if num_buffered >= peptide_buffer:
                                write_binary_peptide_tuples(s, peptides, f)
                                del peptides[:]
                                num_buffered = 0

                    prev_mass = p[0]
                    eq_mass_dict.clear()

                    eq_mass_dict[pep_str] = p

            # we have to shuffle the last set of targets with equal mass, since we fill the buffer first than create
            # decoys
            decoy_eq_mass_set.clear()
            for t in eq_mass_dict:
                if decoy_format == 'shuffle':
                    d = create_shuffled_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                elif decoy_format == 'reverse':
                    d = create_reversed_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                if d:
                    decoy_eq_mass_set.add(d[1])
                    peptides.append(d)
                    num_buffered +=1

                    if num_buffered >= peptide_buffer:
                        write_binary_peptide_tuples(s, peptides, f)
                        del peptides[:]
                        num_buffered = 0

            if num_buffered:
                write_binary_peptide_tuples(s, peptides, f)

def create_recalibrateDecoys_from_binarydb(decoyFile, targets, recalibrateDecoyFile,
                                           max_length,
                                           decoy_format, keep_terminal_aminos, s0,
                                           peptide_buffer = 100000):
    """ pre: ooc_write_merge_peptides has been called and peptides are written
             in binary format according to struct.Struct instance below.
             Note that after ooc_write_merge_peptides, peptides should be
             sorted by mass
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    s = struct.Struct('f %ds I s s' % (max_length))

    if s0.lower == 'time':
        s0 = None
    else:
        try:
            s0 = int(s0)
        except:
            print "Seed must be either an unsigned integer or time, exitting"
            exit(9)
    # seed random number generator
    random.seed(s0)

    # check terminal amino acids option
    kta = keep_terminal_aminos
    if kta != 'N' and kta != 'C' and kta != 'NC' and kta.lower != 'none':
        print "keep-terminal-aminos value supplied is %s, must be one of N, C, NC, or none.  Exitting" % kta
        exit(10)

    peptides = []
    peptide_decoys2 = []
    num_buffered = 0

    prev_mass = 0

    # fill eq_mass_dict with targets of equal mass.  when we move to the next mass, shuffle targets in 
    # eq_mass_dict
    eq_mass_dict = {}

    # not as simple as looking at two adjacent peptides, since several peptides could have the same mass
    # and lead to a sequence of peptides, each with the same mass, with redundant peptides
    decoy_format = decoy_format.lower()

    if decoy_format != 'shuffle' and decoy_format != 'reverse':
        print "valid options for decoy-format are shuffle or reverse, exitting."
        exit(11)

    curr_visited = 0

    decoy_eq_mass_set = set([])

    f2 = open(recalibrateDecoyFile, 'wb')
    with open(decoyFile, 'wb') as f:
        with open(targets, 'rb') as fid:
            for p in read_binary_peptide_tuples(s, fid):
                pep_str = p[1].split('\x00')[0] # remove null characters
                if p[0] == prev_mass:
                    # assume only unique peptides in target database (i.e., filter function has removed all duplicates)
                    eq_mass_dict[pep_str] = p
                else: # we've filled up the queue of targets to be shuffled/reversed with equal mass, 
                    decoy_eq_mass_set.clear()
                    for t in eq_mass_dict:
                        if decoy_format == 'shuffle':
                            d = create_shuffled_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                        elif decoy_format == 'reverse':
                            d = create_reversed_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                        if d:
                            decoy_eq_mass_set.add(d[1])
                            peptides.append(d)
                            num_buffered +=1

                            # create second set of decoys for recalibration
                            d = create_reversed_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                            if d:
                                decoy_eq_mass_set.add(d[1])
                                peptide_decoys2.append(d)
                                num_buffered +=1
                                

                            if num_buffered >= peptide_buffer:
                                write_binary_peptide_tuples(s, peptides, f)
                                del peptides[:]
                                # write recalibration decoy set
                                write_binary_peptide_tuples(s, peptide_decoys2, f2)
                                del peptide_decoys2[:]

                                num_buffered = 0

                    prev_mass = p[0]
                    eq_mass_dict.clear()

                    eq_mass_dict[pep_str] = p

            # we have to shuffle the last set of targets with equal mass, since we fill the buffer first than create
            # decoys
            decoy_eq_mass_set.clear()
            for t in eq_mass_dict:
                if decoy_format == 'shuffle':
                    d = create_shuffled_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                elif decoy_format == 'reverse':
                    d = create_reversed_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                if d:
                    decoy_eq_mass_set.add(d[1])
                    peptides.append(d)
                    num_buffered +=1

                # create second set of decoys for recalibration
                d = create_reversed_decoy(t, eq_mass_dict[t], kta, eq_mass_dict, decoy_eq_mass_set)
                if d:
                    decoy_eq_mass_set.add(d[1])
                    peptide_decoys2.append(d)
                    num_buffered +=1


                    if num_buffered >= peptide_buffer:
                        write_binary_peptide_tuples(s, peptides, f)
                        del peptides[:]
                        # write recalibration decoy set
                        write_binary_peptide_tuples(s, peptide_decoys2, f2)
                        del peptide_decoys2[:]

                        num_buffered = 0

            if num_buffered:
                write_binary_peptide_tuples(s, peptides, f)
                write_binary_peptide_tuples(s, peptide_decoys2, f2)
    f2.close()

def write_binary_peptide_to_ascii(output, max_length, targets):
    """ pre: ooc_write_merge_peptides has been called and peptides are written
             in binary format according to struct.Struct instance below.
             Note that after ooc_write_merge_peptides, peptides should be
             sorted by mass
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    s = struct.Struct('f %ds I s s' % (max_length))
    
    with open(output, 'w') as f:
        with open(targets, 'rb') as fid:
            f.write("Protein\tPeptide\tMass\tNterm_flanking\tCterm_flanking\n")
            for p in read_binary_peptide_tuples(s, fid):
                f.write("%d\t%s\t%f\t%c\t%c\n" % (p[2], p[1].split('\x00')[0], p[0], p[3], p[4]))

# def write_binary_peptide_to_ascii(output, max_length, targets):
#     """ pre: ooc_write_merge_peptides has been called and peptides are written
#              in binary format according to struct.Struct instance below.
#              Note that after ooc_write_merge_peptides, peptides should be
#              sorted by mass
#     """

#     # targets = os.path.join(digest_dir, inFile)
#     #    peptide[0] = peptide mass (float)
#     #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
#     #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
#     #    peptide[3] = nterm_flanking (character)
#     #    peptide[4] = cterm_flanking (character)
#     s = struct.Struct('f %ds I s s' % (max_length))
    
#     with open(output, 'w') as f:
#         with open(targets, 'rb') as fid:
#             f.write("Protein\tPeptide\tMass\tNterm_flanking\tCterm_flanking\n")
#             for p in read_binary_peptide_tuples(s, fid):
#                 f.write("%d\t%s\t%f\t%c\t%c\n" % (p[2], p[1].split('\x00')[0], p[0], p[3], p[4]))

def write_binary_peptide_to_ascii_var_mods(output, max_length, targets,
                                           mods, nterm_mods, cterm_mods,
                                           var_mods, nterm_var_mods, cterm_var_mods):
    """ pre: ooc_write_merge_peptides has been called and peptides are written
             in binary format according to struct.Struct instance below.
             Note that after ooc_write_merge_peptides, peptides should be
             sorted by mass
    """
    #    peptide[0] = peptide mass (float)
    #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
    #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
    #    peptide[3] = nterm_flanking (character)
    #    peptide[4] = cterm_flanking (character)
    #    peptide[5] = binary string deoting variable modifications
    s = struct.Struct('f %ds I s s %ds' % (max_length, max_length))
    
    with open(output, 'w') as f:
        with open(targets, 'rb') as fid:
            f.write("Protein\tPeptide\tMass\tNterm_flanking\tCterm_flanking\n")
            for p in read_binary_peptide_tuples(s, fid):
                pep_str = []
                # pep = p[1].split('\x00')[0]
                for c, vm in zip(p[1].split('\x00')[0], p[5]):
                    pep_str += c
                    if vm == '1': # variable mods
                        pep_str += str('[%1.0e]' % var_mods[c][1])
                    elif vm == '2':
                        pep_str += str('[%1.0e]' % nterm_var_mods[c][1])
                    elif vm == '3':
                        pep_str += str('[%1.0e]' % cterm_var_mods[c][1])
                        
                f.write("%d\t%s\t%f\t%c\t%c\n" % (p[2], ''.join(pep_str), p[0], p[3], p[4]))
                
def load_peptides(f):
    peptides = {}
    fid = open(f, 'r')
    for l in csv.DictReader(fid, delimiter = '\t'):
        peptides[l["Peptide"]] = digestedPep(l["Peptide"], 
                                             l["Protein"],
                                             l["Nterm_flanking"], 
                                             l["Cterm_flanking"], 
                                             float(l["Mass"]))
    fid.close()
    return peptides

def create_decoys(targets, decoy_format, keep_terminal_aminos, s0):
    """ Create decoy peptides
    """
    if s0.lower == 'time':
        s0 = None
    else:
        try:
            s0 = int(s0)
        except:
            print "Seed must be either an unsigned integer or time, exitting"
            exit(9)
    # seed random number generator
    random.seed(s0)

    # check terminal amino acids option
    kta = keep_terminal_aminos
    if kta != 'N' and kta != 'C' and kta != 'NC' and kta.lower != 'none':
        print "keep-terminal-aminos value supplied is %s, must be one of N, C, NC, or none.  Exitting" % kta
        exit(10)

    decoys = {}
    if decoy_format.lower() == 'shuffle':
        for p in targets:
            cterm = [p[-1]]
            nterm = [p[0]]
            if kta == 'NC':
                l = list(p[1:-1])
                random.shuffle(l)
                shuffled_target = nterm + l + cterm
            elif kta == 'N':
                l = list(p[1:])
                random.shuffle(l)
                shuffled_target = nterm + l
            elif kta == 'C':
                l = list(p[:-1])
                random.shuffle(l)
                shuffled_target = l + cterm
            else:
                l = list(p)
                random.shuffle(l)
                shuffled_target = l

            d = ''.join(shuffled_target)
            perms = 1
            properly_shuffled = 1
            max_iters = 100*len(l)
            while d in targets or d in decoys:
                if perms >= max_iters:
                    properly_shuffled = 0
                    break
                if kta == 'NC':
                    l = list(p[1:-1])
                    random.shuffle(l)
                    shuffled_target = nterm + l + cterm
                elif kta == 'N':
                    l = list(p[1:])
                    random.shuffle(l)
                    shuffled_target = nterm + l
                elif kta == 'C':
                    l = list(p[:-1])
                    random.shuffle(l)
                    shuffled_target = l + cterm
                else:
                    l = list(p)
                    random.shuffle(l)
                    shuffled_target = l
                d = ''.join(shuffled_target)
                perms += 1
            if properly_shuffled:
                t = targets[p]
                protein_name = 'decoy' + t.proteinName
                decoys[d] = digestedPep(Peptide(d),
                                        protein_name,
                                        t.ntermFlanking,
                                        t.ctermFlanking,
                                        t.peptideMass,
                                        t.missed_cleavages)
    elif decoy_format.lower() == 'reverse':
        for p in targets:
            cterm = p[-1]
            nterm = p[0]
            if kta == 'NC':
                l = p[1:-1]
                d = nterm + l[::-1] + cterm
            elif kta == 'N':
                l = p[1:]
                d = nterm + l[::-1]
            elif kta == 'C':
                l = p[:-1]
                d = l[::-1] + cterm
            else:
                d = p[::-1]
            if d not in targets and d not in decoys:
                t = targets[p]
                protein_name = 'decoy' + t.proteinName
                decoys[d] = digestedPep(Peptide(d),
                                        protein_name,
                                        t.ntermFlanking,
                                        t.ctermFlanking,
                                        t.peptideMass,
                                        t.missed_cleavages)
    else:
        print "valid options for decoy-format are shuffle or reverse, exitting."
        exit(11)
    return decoys

def create_shuffled_decoy(p, pep_tuple, kta, targets, decoys):
    cterm = [p[-1]]
    nterm = [p[0]]
    if kta == 'NC':
        l = list(p[1:-1])
        random.shuffle(l)
        shuffled_target = nterm + l + cterm
    elif kta == 'N':
        l = list(p[1:])
        random.shuffle(l)
        shuffled_target = nterm + l
    elif kta == 'C':
        l = list(p[:-1])
        random.shuffle(l)
        shuffled_target = l + cterm
    else:
        l = list(p)
        random.shuffle(l)
        shuffled_target = l

    d = ''.join(shuffled_target)
    perms = 1
    properly_shuffled = 1
    max_iters = 100*len(l)
    while d in targets or d in decoys:
        if perms >= max_iters:
            properly_shuffled = 0
            break
        if kta == 'NC':
            l = list(p[1:-1])
            random.shuffle(l)
            shuffled_target = nterm + l + cterm
        elif kta == 'N':
            l = list(p[1:])
            random.shuffle(l)
            shuffled_target = nterm + l
        elif kta == 'C':
            l = list(p[:-1])
            random.shuffle(l)
            shuffled_target = l + cterm
        else:
            l = list(p)
            random.shuffle(l)
            shuffled_target = l
        d = ''.join(shuffled_target)
        perms += 1
    if properly_shuffled:
        t = pep_tuple
        return (t[0], d, t[2], t[3], t[4])
    return ''

def create_reversed_decoy(p, pep_tuple, kta, targets, decoys):
    cterm = p[-1]
    nterm = p[0]
    if kta == 'NC':
        l = p[1:-1]
        d = nterm + l[::-1] + cterm
    elif kta == 'N':
        l = p[1:]
        d = nterm + l[::-1]
    elif kta == 'C':
        l = p[:-1]
        d = l[::-1] + cterm
    else:
        d = p[::-1]
    if d not in targets and d not in decoys:
        t = pep_tuple
        return (t[0], d, t[2], t[3], t[4])
    return ''

def parse_var_mods(mod_spec, mainMods = True):
    # parse modifications
    mods = {}

    var_mods = {}

    if not mod_spec:
        return mods, var_mods

    pattern = re.compile('(?P<maxPerPeptide>[0-9]*)(?P<residues>[a-zA-z]*)(?P<massChange>\S+)')

    for m in mod_spec.split(','):
        x = pattern.match(m)
        if x.group('maxPerPeptide'):
            mpp = int(x.group('maxPerPeptide'))

            # these are the variable mods
            assert x.group('residues'), "No residues supplied in variable modification definition %s, exitting" % m
            
            for aminoAcid in x.group('residues'):
                if aminoAcid == 'X':
                    for aa in validPeps:
                        if aa in var_mods or aa in mods:
                            print "Two modifications supplied for %c, exitting" % aa
                            exit(12)
                        var_mods[aa] = (mpp, float(x.group('massChange')))
                else:
                    if aminoAcid in var_mods or aminoAcid in mods:
                        print "Two modifications supplied for %c, exitting" % m[0]
                        exit(13)
                    var_mods[aminoAcid] = (mpp, float(x.group('massChange')))
        else: # static mod
            # these are the mods
            assert x.group('residues'), "No residues supplied in modification definition %s, exitting" % m
            
            for aminoAcid in x.group('residues'):
                if aminoAcid == 'X':
                    for aa in validPeps:
                        if aa in mods or aa in var_mods:
                            print "Two modifications supplied for %c, exitting" % aa
                            exit(12)
                        mods[aa] = float(x.group('massChange'))
                else:
                    if aminoAcid in mods or aminoAcid in var_mods:
                        print "Two modifications supplied for %c, exitting" % m[0]
                        exit(13)
                    mods[aminoAcid] = float(m[1:])
    if mainMods:
        if 'C' in mods:
            if mods['C'] != 0:
                mods.pop('C')
            else:
                mods['C'] = -57.021464 # cancel the default cysteine modification
    return mods, var_mods

def parse_mods(mod_spec, mainMods = True):
    # parse modifications
    mods = {}
    
    if not mod_spec:
        return mods

    for m in mod_spec.split(','):
        # to do, add parsing for variable mods, which begin with a number
        if m[0] != '1':
            if m[0] == 'X':
                for aa in validPeps:
                    if aa in mods:
                        print "Two modifications supplied for %c, exitting" % aa
                        exit(12)
                    mods[aa] = float(m[1:])
            else:
                if m[0] in mods:
                    print "Two modifications supplied for %c, exitting" % m[0]
                    exit(13)
                mods[m[0]] = float(m[1:])

    if mainMods:
        if 'C' in mods:
            if mods['C'] != 0:
                mods.pop('C')
            else:
                mods['C'] = -57.021464 # cancel the default cysteine modification
    return mods    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Digest FASTA file, producing proteins.")
    ############## input and output options
    iFileGroup = parser.add_argument_group('iFileGroup', 'Necessary input files.')
    help_digest_dir = '<string> - Digestion output file.  Default=dripDigest-output.'
    iFileGroup.add_argument('--digest-dir', type = str, action = 'store',
                        help = help_digest_dir, default = 'dripDigest-output')
    help_pepdb = '<string> - Protein FASTA file.'
    iFileGroup.add_argument('--fasta', type = str, action = 'store',
                        help = help_pepdb)
    help_peptide_list = '<T|F> -  Create in the output directory a text file listing of all the peptides in the database, along with their neutral masses, one per line. If decoys are generated, then a second file will be created containing the decoy peptides. Decoys that also appear in the target database are marked with an asterisk in a third column. Default = False'
    iFileGroup.add_argument('--peptide-list', type = str, action = 'store',
                        help = help_peptide_list, default = 'False')
    ############## peptide properties
    peptidePropertiesGroup = parser.add_argument_group('peptidePropertiesGroup', 'Options for digested peptids.')
    max_length_help = """<integer> - The maximum length of peptides to consider. Default = 50."""
    peptidePropertiesGroup.add_argument('--max-length', type = int, action = 'store',
                                        default = 50, help = max_length_help)
    max_mass_help = """<float> - The maximum mass (in Da) of peptides to consider. Default = 7200."""
    peptidePropertiesGroup.add_argument('--max-mass', type = float, action = 'store',
                                        default = 7200.0, help = max_mass_help)
    min_length_help = """<integer> - The minimum length of peptides to consider. Default = 6."""
    peptidePropertiesGroup.add_argument('--min-length', type = int, action = 'store',
                                        default = 6, help = min_length_help)
    min_mass_help = """<float> - The minimum mass (in Da) of peptides to consider. Default = 7200."""
    peptidePropertiesGroup.add_argument('--min-mass', type = float, action = 'store',
                                        default = 200.0, help = min_mass_help)
    monoisotopic_precursor_help = """<T|F> - When computing the mass of a peptide, use monoisotopic masses rather than average masses. Default = true."""
    peptidePropertiesGroup.add_argument('--monoisotopic-precursor', type = str, action = 'store', 
                                        default = 'True', help = monoisotopic_precursor_help)
    peptide_buffer_help = """<int> - The maximum number of peptides to keep in memory."""
    peptidePropertiesGroup.add_argument('--peptide-buffer', type = int, action = 'store',
                                        default = 100000, help = peptide_buffer_help)
    ############## amino acid modifications
    aaModsGroup = parser.add_argument_group('aaModsGroup', 'Options for amino acid modifications.')
    aaModsGroup.add_argument('--mods-spec', type = str, action = 'store',
                             default = 'C+57.02146')
                        # default = '')
    aaModsGroup.add_argument('--cterm-peptide-mods-spec', type = str, action = 'store',
                             default = '')
    aaModsGroup.add_argument('--nterm-peptide-mods-spec', type = str, action = 'store',
                             default = '')
    aaModsGroup.add_argument('--max-mods', type = int, action = 'store',
                             default = 255)
    aaModsGroup.add_argument('--min-mods', type = int, action = 'store',
                             default = 0)
    ############## decoy database generation
    decoyDBGroup = parser.add_argument_group('decoyDBGroup', 'Options decoy database generation.')
    help_decoy_format = '<shuffle|peptide-reverse> - Include a decoy version of every peptide by shuffling or reversing the target sequence or protein. In shuffle or peptide-reverse mode, each peptide is either reversed or shuffled, leaving the N-terminal and C-terminal amino acids in place. Note that peptides appear multiple times in the target database are only shuffled once. In peptide-reverse mode, palindromic peptides are shuffled. Also, if a shuffled peptide produces an overlap with the target or decoy database, then the peptide is re-shuffled up to 100n times, where n is the length of the string to be shuffled. Note that, despite this repeated shuffling, homopolymers will appear in both the target and decoy database. The protein-reverse mode reverses the entire protein sequence, irrespective of the composite peptides. Default = shuffle.'
    decoyDBGroup.add_argument('--decoy-format', type = str, action = 'store', 
                              default = 'shuffle', help = help_decoy_format)
    help_keep_terminal_aminos = '<N|C|NC|none> - When creating decoy peptides using decoy-format=shuffle or decoy-format=peptide-reverse, this option specifies whether the N-terminal and C-terminal amino acids are kept in place or allowed to be shuffled or reversed. For a target peptide "EAMPK" with decoy-format=peptide-reverse, setting keep-terminal-aminos to "NC" will yield "EPMAK"; setting it to "C" will yield "PMAEK"; setting it to "N" will yield "EKPMA"; and setting it to "none" will yield "KPMAE". Default = NC.'
    decoyDBGroup.add_argument('--keep-terminal-aminos', type = str, action = 'store', 
                              default = 'NC', help = help_keep_terminal_aminos)
    help_seed = '<seed> - When given a unsigned integer value seeds the random number generator with that value. When given the string "time" seeds the random number generator with the system time. Default = 1.'
    decoyDBGroup.add_argument('--seed', type = str, action = 'store', 
                              default = '1', help = help_seed)
    help_decoys = '<T|F> - whether to create (shuffle target peptides) and search decoy peptides. Default = True'
    decoyDBGroup.add_argument('--decoys', type = str, action = 'store', 
                              default = 'True', help = help_decoys)
    help_recalibrate = '<T|F> - whether to create second set of decoys (shuffling only). Default = False'
    decoyDBGroup.add_argument('--recalibrate', type = str, action = 'store', 
                              default = 'False', help = help_recalibrate)
    ############ enzymatic digestion
    enzymeDigestGroup = parser.add_argument_group('enzymeDigestGroup', 'Options for the digesting enzyme.')
    enzyme_help = """"<no-enzyme|trypsin|trypsin/p|chymotrypsin|elastase|clostripain|cyanogen-bromide|iodosobenzoate|proline-endopeptidase|staph-protease|asp-n|lys-c|lys-n|arg-c|glu-c|pepsin-a|elastase-trypsin-chymotrypsin|custom-enzyme> - Specify the enzyme used to digest the proteins in silico. Available enzymes (with the corresponding digestion rules indicated in parentheses) include no-enzyme ([X]|[X]), trypsin ([RK]|{P}), trypsin/p ([RK]|[]), chymotrypsin ([FWYL]|{P}), elastase ([ALIV]|{P}), clostripain ([R]|[]), cyanogen-bromide ([M]|[]), iodosobenzoate ([W]|[]), proline-endopeptidase ([P]|[]), staph-protease ([E]|[]), asp-n ([]|[D]), lys-c ([K]|{P}), lys-n ([]|[K]), arg-c ([R]|{P}), glu-c ([DE]|{P}), pepsin-a ([FL]|{P}), elastase-trypsin-chymotrypsin ([ALIVKRWFY]|{P}). Specifying --enzyme no-enzyme yields a non-enzymatic digest. Warning: the resulting peptide database may be quite large. Default = trypsin."""
    enzymeDigestGroup.add_argument('--enzyme', type = str, action = 'store',
                                   default = 'trypsin', help = enzyme_help)
    custom_enzyme_help = """<string> - Specify rules for in silico digestion of protein sequences. Overrides the enzyme option. Two lists of residues are given enclosed in square brackets or curly braces and separated by a |. The first list contains residues required/prohibited before the cleavage site and the second list is residues after the cleavage site. If the residues are required for digestion, they are in square brackets, '[' and ']'. If the residues prevent digestion, then they are enclosed in curly braces, '{' and '}'. Use X to indicate all residues. For example, trypsin cuts after R or K but not before P which is represented as [RK]|{P}. AspN cuts after any residue but only before D which is represented as [X]|[D]. Default = <empty>."""
    enzymeDigestGroup.add_argument('--custom-enzyme', type = str, action = 'store',
                                   default = '', help = custom_enzyme_help)
    missed_cleavages_help = """<integer> - Maximum number of missed cleavages per peptide to allow in enzymatic digestion. Default = 0."""
    enzymeDigestGroup.add_argument('--missed-cleavages', type = int, action = 'store',
                                   default = 0, help = custom_enzyme_help)
    digestion_help = """<full-digest|partial-digest> - Specify whether every peptide in the database must have two enzymatic termini (full-digest) or if peptides with only one enzymatic terminus are also included (partial-digest). Default=full-digest."""
    enzymeDigestGroup.add_argument('--digestion', type = str, action = 'store',
                                   default = "full-digest", help = digestion_help)
    args = parser.parse_args()
    
    if not args.fasta:
        parser.print_help()
        exit(8)

    # set true or false strings to booleans
    args.decoys = check_arg_trueFalse(args.decoys)
    args.monoisotopic_precursor = check_arg_trueFalse(args.monoisotopic_precursor)
    args.recalibrate = check_arg_trueFalse(args.recalibrate)
    args.peptide_list = check_arg_trueFalse(args.peptide_list)

    # formulate digestion regular expression
    digest_re = create_digest_re(args.enzyme, args.custom_enzyme)
    print digest_re
    r = re.compile(digest_re)
    
    # parse modifications
    mods, var_mods = parse_var_mods(args.mods_spec, True)
    print "Static mods:"
    print mods
    print "Variable mods:"
    print var_mods
    nterm_mods, nterm_var_mods = parse_var_mods(args.nterm_peptide_mods_spec, False)
    print "Static n-terminal mods:"
    print nterm_mods
    print "Variable n-terminal mods:"
    print nterm_var_mods
    cterm_mods, cterm_var_mods = parse_var_mods(args.cterm_peptide_mods_spec, False)
    print "Static c-terminal mods:"
    print cterm_mods
    print "Variable c-terminal mods:"
    print cterm_var_mods

    args.digestion = args.digestion.lower()
    if args.digestion != "full-digest" and args.digestion != "partial-digest":
        print "Option --digestion specified as %s, must be either full-digest or partial-digest.  Exitting" % (args.digestion)
        exit(-1)

    print args.digestion

    try:
        base = os.path.abspath(args.digest_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
        args.digest_dir = os.path.abspath(args.digest_dir)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.digest_dir)
        exit(-1)

    # calculate enzymatic markers
    n_res, c_res = enzyme_enzSet(args.enzyme, args.custom_enzyme)    

    # pickle parameters to be used by other programs
    with open(os.path.join(args.digest_dir, 'options.pickle'), 'w') as f:
        pickle.dump(args, f, pickle.HIGHEST_PROTOCOL)

    ##################################################################
    ################## out of core sort
    ##################################################################
    initDbFile = os.path.join(args.digest_dir, 'tempTargets.bin')
    ooc_write_merge_peptides(args.fasta, r,
                             initDbFile,
                             n_res, c_res,
                             args.min_length, args.max_length,
                             args.min_mass, args.max_mass,
                             mods, nterm_mods, cterm_mods,
                             args.max_mods, args.min_mods,
                             args.missed_cleavages,
                             args.monoisotopic_precursor,
                             args.digestion,
                             args.peptide_buffer)

    # the sort does not get rid of doubles in the database, so filter these
    filteredTargets = os.path.join(args.digest_dir, 'digestedTargets.bin')
    filter_double_binarydb(filteredTargets, initDbFile, args.max_length, args.peptide_buffer)
    # delete temporary file
    os.remove(initDbFile)

    if args.decoys:
        decoyFile = os.path.join(args.digest_dir, 'digestedDecoys.bin')
        create_decoys_from_binarydb(decoyFile, filteredTargets, args.max_length,
                                    args.decoy_format, args.keep_terminal_aminos, 
                                    args.seed, args.peptide_buffer)
        if args.recalibrate:
            recalibrateDecoyFile = os.path.join(args.digest_dir, 'recalibrateDecoys.bin')
            create_recalibrateDecoys_from_binarydb(decoyFile, filteredTargets, recalibrateDecoyFile,
                                                   args.max_length, 
                                                   args.decoy_format, args.keep_terminal_aminos, args.seed, args.peptide_buffer)
            # write_binary_peptide_to_ascii('recalibrationDecoys.txt', args.max_length, args.digest_dir, 'recalibrationDecoys.bin')

    if not var_mods and not nterm_var_mods and not cterm_var_mods:
        # no variable mods supplied, stop program
        os.rename(filteredTargets, os.path.join(args.digest_dir, 'targets.bin'))
        os.rename(decoyFile, os.path.join(args.digest_dir, 'decoys.bin'))

        if args.peptide_list: # write ascii files with digested peptides
            write_binary_peptide_to_ascii(os.path.join(args.digest_dir, 'targets.txt'), 
                                          args.max_length,
                                          os.path.join(args.digest_dir, 'targets.bin'))

            write_binary_peptide_to_ascii(os.path.join(args.digest_dir, 'decoys.txt'), 
                                          args.max_length,
                                          os.path.join(args.digest_dir, 'decoys.bin'))

            if args.recalibrate:
                # os.rename(recalibrateDecoyFile, os.path.join(args.digest_dir, 'recalibrateDecoys.bin'))
                write_binary_peptide_to_ascii(os.path.join(args.digest_dir, 'recalibrateDecoys.txt'), 
                                              args.max_length,
                                              os.path.join(args.digest_dir, 'recalibrateDecoys.bin'))
        exit(0)
        
    # create sets of variable mods
    tempVarModsTargets = os.path.join(args.digest_dir, 'tempVarModsTargets.bin')
    digested_peptides_var_mods(filteredTargets, tempVarModsTargets,
                               args.min_length, args.max_length,
                               args.min_mass, args.max_mass,
                               var_mods, nterm_var_mods, cterm_var_mods,
                               args.max_mods, args.min_mods,
                               args.peptide_buffer)
    # delete original file
    os.remove(filteredTargets)
    # the sort does not get rid of doubles in the database, so filter these
    varModsTargets = os.path.join(args.digest_dir, 'targets.bin')
    filter_double_binarydb_var_mods(varModsTargets, tempVarModsTargets, args.max_length, args.peptide_buffer)
    # delete temporary file
    os.remove(tempVarModsTargets)

    if args.peptide_list:
        write_binary_peptide_to_ascii_var_mods(os.path.join(args.digest_dir, 'targets.txt'), 
                                               args.max_length, varModsTargets,
                                               mods, nterm_mods, cterm_mods,
                                               var_mods, nterm_var_mods, cterm_var_mods)

    if args.decoys:
        # create sets of variable mods
        tempVarModsDecoys = os.path.join(args.digest_dir, 'tempVarModsDecoys.bin')
        digested_peptides_var_mods(decoyFile, tempVarModsDecoys,
                                   args.min_length, args.max_length,
                                   args.min_mass, args.max_mass,
                                   var_mods, nterm_var_mods, cterm_var_mods,
                                   args.max_mods, args.min_mods,
                                   args.peptide_buffer)

        # delete original file
        os.remove(decoyFile)

        varModsDecoys = os.path.join(args.digest_dir, 'decoys.bin')
        filter_double_binarydb_var_mods(varModsDecoys, tempVarModsDecoys, args.max_length, args.peptide_buffer)
        # delete temporary file and old decoy file
        os.remove(tempVarModsDecoys)

        if args.peptide_list:
            write_binary_peptide_to_ascii_var_mods(os.path.join(args.digest_dir, 'decoys.txt'), 
                                                   args.max_length, varModsDecoys,
                                                   mods, nterm_mods, cterm_mods,
                                                   var_mods, nterm_var_mods, cterm_var_mods)

        if args.recalibrate:
        # create sets of variable mods
            tempVarModsDecoys = os.path.join(args.digest_dir, 'tempVarModsRecalibrateDecoys.bin')
            digested_peptides_var_mods(recalibrateDecoyFile, tempVarModsDecoys,
                                       args.min_length, args.max_length,
                                       args.min_mass, args.max_mass,
                                       var_mods, nterm_var_mods, cterm_var_mods,
                                       args.max_mods, args.min_mods,
                                       args.peptide_buffer)
            varModsDecoys = os.path.join(args.digest_dir, 'recalibrateDecoys.bin')
            filter_double_binarydb_var_mods(varModsDecoys, tempVarModsDecoys, args.max_length, args.peptide_buffer)
            # delete temporary file
            os.remove(tempVarModsDecoys)
            if args.peptide_list:
                write_binary_peptide_to_ascii_var_mods(os.path.join(args.digest_dir, 'recalibrateDecoys.txt'), 
                                                       args.max_length, varModsDecoys,
                                                       mods, nterm_mods, cterm_mods,
                                                       var_mods, nterm_var_mods, cterm_var_mods)
