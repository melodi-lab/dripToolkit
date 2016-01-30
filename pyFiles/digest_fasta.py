#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

"""Support for absolute ranking plots"""

import argparse
import random
import itertools
import re
import csv
import collections
import cPickle as pickle
from peptide import Peptide

validPeps = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')-set('JOBZUX')
# struct for proteins
# digestedPep = collections.namedtuple('dPep',
#                                      'peptide proteinName ntermFlanking ctermFlanking peptideMass')

class digestedPep(object):
    """ Digested peptide class containing parent protein and flanking amino acid information

    """
    def __init__(self, peptide = '', proteinName = '', 
                 ntermFlanking = '', ctermFlanking = '',
                 peptideMass = -1, 
                 missed_cleavages = 0):

        self.peptide = peptide
        self.proteinName = proteinName
        self.ntermFlanking = ntermFlanking
        self.ctermFlanking = ctermFlanking
        self.peptideMass = peptideMass 
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

    return(digestedPep(pep, protein_name,
                       nterm_flanking, 
                       cterm_flanking, 
                       pm))


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


def write_peptides(f, peptides):
    # first sort peptides by mass
    peps = [peptides[pep] for pep in peptides]
    peps.sort(key = lambda x: x.peptideMass)

    # write digested peptides
    fid = open(f, 'w')
    fid.write("Protein\tPeptide\tMass\tNterm_flanking\tCterm_flanking\n")

    for p in peps:
        fid.write("%s\t%s\t%f\t%c\t%c\n" % 
                  (p.proteinName, p.peptide.seq, 
                   p.peptideMass, p.ntermFlanking, p.ctermFlanking))
    fid.close()

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

def parse_mods(mod_spec, mainMods = True):
    # parse modifications
    mods = {}

    var_mods = {}

    if not mod_spec:
        return mods

    pattern = re.compile('(?P<maxPerPeptide>[0-9]*)(?P<residues>[a-zA-z]*)(?P<massChange>\S+)')

    for m in mod_spec.split(','):
        x = pattern.match(m)
        if x.group('maxPerPeptide'):
            
            mpp = int(x.group('maxPerPeptide'))

            # these are the variable mods
            assert x.group('residues'), "No residues supplied in variable modification definition %s, exitting" % m
            
            for aminoAcid in x.group('residues'):
                if m[0] == 'X':
                    for aa in validPeps:
                        if aa in var_mods:
                            print "Two modifications supplied for %c, exitting" % aa
                            exit(12)
                        var_mods[aa] = (mpp, float(x.group('massChange')))

        else: # static mod
            # these are the mods
            assert x.group('residues'), "No residues supplied in modification definition %s, exitting" % m
            
            for aminoAcid in x.group('residues'):
                if m[0] == 'X':
                    for aa in validPeps:
                        if aa in mods:
                            print "Two modifications supplied for %c, exitting" % aa
                            exit(12)
                        mods[aa] = float(x.group('massChange'))
    if mainMods:
        if 'C' in mods:
            if mods['C'] != 0:
                mods.pop('C')
            else:
                mods['C'] = -57.021464 # cancel the default cysteine modification
    return mods, var_mods

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Digest FASTA file, producing proteins.")
    ############## input and output options
    help_output_file = '<string> - Output file for resulting peptides.'
    parser.add_argument('--output_file', type = str, action = 'store', 
                        default = 'peptides.txt', help = help_output_file)
    help_pepdb = '<string> - Protein FASTA file.'
    parser.add_argument('--fasta', type = str, action = 'store',
                        help = help_pepdb)

    ############## peptide properties
    max_length_help = """<integer> - The maximum length of peptides to consider. Default = 50."""
    parser.add_argument('--max-length', type = int, action = 'store',
                        default = 50, help = max_length_help)
    max_mass_help = """<foat> - The maximum mass (in Da) of peptides to consider. Default = 7200."""
    parser.add_argument('--max-mass', type = float, action = 'store',
                        default = 7200.0, help = max_mass_help)
    min_length_help = """<integer> - The minimum length of peptides to consider. Default = 6."""
    parser.add_argument('--min-length', type = int, action = 'store',
                        default = 6, help = min_length_help)
    min_mass_help = """<foat> - The minimum mass (in Da) of peptides to consider. Default = 7200."""
    parser.add_argument('--min-mass', type = float, action = 'store',
                        default = 200.0, help = min_mass_help)
    monoisotopic_precursor_help = """<T|F> - When computing the mass of a peptide, use monoisotopic masses rather than average masses. Default = true."""
    parser.add_argument('--monoisotopic-precursor', type = str, action = 'store', 
                        default = 'True', help = monoisotopic_precursor_help)

    ############## amino acid modifications
    parser.add_argument('--mods-spec', type = str, action = 'store',
                        default = 'C+57.02146')
                        # default = '')
    parser.add_argument('--cterm-peptide-mods-spec', type = str, action = 'store',
                        default = '')
    parser.add_argument('--nterm-peptide-mods-spec', type = str, action = 'store',
                        default = '')
    parser.add_argument('--max-mods', type = int, action = 'store',
                        default = 255)
    parser.add_argument('--min-mods', type = int, action = 'store',
                        default = 0)
    ############## decoy database generation
    help_decoys = '<T|F> - whether to create (shuffle target peptides) and search decoy peptides. Default = True'
    parser.add_argument('--decoys', type = str, action = 'store', 
                        default = 'True', help = help_decoys)
    help_decoy_format = '<shuffle|peptide-reverse> - Include a decoy version of every peptide by shuffling or reversing the target sequence or protein. In shuffle or peptide-reverse mode, each peptide is either reversed or shuffled, leaving the N-terminal and C-terminal amino acids in place. Note that peptides appear multiple times in the target database are only shuffled once. In peptide-reverse mode, palindromic peptides are shuffled. Also, if a shuffled peptide produces an overlap with the target or decoy database, then the peptide is re-shuffled up to 100n times, where n is the length of the string to be shuffled. Note that, despite this repeated shuffling, homopolymers will appear in both the target and decoy database. The protein-reverse mode reverses the entire protein sequence, irrespective of the composite peptides. Default = shuffle.'
    parser.add_argument('--decoy-format', type = str, action = 'store', 
                        default = 'shuffle', help = help_decoy_format)
    help_keep_terminal_aminos = '<N|C|NC|none> - When creating decoy peptides using decoy-format=shuffle or decoy-format=peptide-reverse, this option specifies whether the N-terminal and C-terminal amino acids are kept in place or allowed to be shuffled or reversed. For a target peptide "EAMPK" with decoy-format=peptide-reverse, setting keep-terminal-aminos to "NC" will yield "EPMAK"; setting it to "C" will yield "PMAEK"; setting it to "N" will yield "EKPMA"; and setting it to "none" will yield "KPMAE". Default = NC.'
    parser.add_argument('--keep-terminal-aminos', type = str, action = 'store', 
                        default = 'NC', help = help_keep_terminal_aminos)
    help_seed = '<seed> - When given a unsigned integer value seeds the random number generator with that value. When given the string "time" seeds the random number generator with the system time. Default = 1.'
    parser.add_argument('--seed', type = str, action = 'store', 
                        default = '1', help = help_seed)
    ############ enzymatic digestion
    enzyme_help = """"<no-enzyme|trypsin|trypsin/p|chymotrypsin|elastase|clostripain|cyanogen-bromide|iodosobenzoate|proline-endopeptidase|staph-protease|asp-n|lys-c|lys-n|arg-c|glu-c|pepsin-a|elastase-trypsin-chymotrypsin|custom-enzyme> - Specify the enzyme used to digest the proteins in silico. Available enzymes (with the corresponding digestion rules indicated in parentheses) include no-enzyme ([X]|[X]), trypsin ([RK]|{P}), trypsin/p ([RK]|[]), chymotrypsin ([FWYL]|{P}), elastase ([ALIV]|{P}), clostripain ([R]|[]), cyanogen-bromide ([M]|[]), iodosobenzoate ([W]|[]), proline-endopeptidase ([P]|[]), staph-protease ([E]|[]), asp-n ([]|[D]), lys-c ([K]|{P}), lys-n ([]|[K]), arg-c ([R]|{P}), glu-c ([DE]|{P}), pepsin-a ([FL]|{P}), elastase-trypsin-chymotrypsin ([ALIVKRWFY]|{P}). Specifying --enzyme no-enzyme yields a non-enzymatic digest. Warning: the resulting peptide database may be quite large. Default = trypsin."""
    parser.add_argument('--enzyme', type = str, action = 'store',
                        default = 'trypsin', help = enzyme_help)
    custom_enzyme_help = """<string> - Specify rules for in silico digestion of protein sequences. Overrides the enzyme option. Two lists of residues are given enclosed in square brackets or curly braces and separated by a |. The first list contains residues required/prohibited before the cleavage site and the second list is residues after the cleavage site. If the residues are required for digestion, they are in square brackets, '[' and ']'. If the residues prevent digestion, then they are enclosed in curly braces, '{' and '}'. Use X to indicate all residues. For example, trypsin cuts after R or K but not before P which is represented as [RK]|{P}. AspN cuts after any residue but only before D which is represented as [X]|[D]. Default = <empty>."""
    parser.add_argument('--custom-enzyme', type = str, action = 'store',
                        default = '', help = custom_enzyme_help)
    missed_cleavages_help = """<integer> - Maximum number of missed cleavages per peptide to allow in enzymatic digestion. Default = 0."""
    parser.add_argument('--missed-cleavages', type = int, action = 'store',
                        default = 0, help = custom_enzyme_help)
    digestion_help = """<full-digest|partial-digest> - Specify whether every peptide in the database must have two enzymatic termini (full-digest) or if peptides with only one enzymatic terminus are also included (partial-digest). Default=full-digest."""
    parser.add_argument('--digestion', type = str, action = 'store',
                        default = "full-digest", help = digestion_help)
    args = parser.parse_args()
    
    if not args.output_file or not args.fasta:
        parser.print_help()
        exit(8)

    # set true or false strings to booleans
    args.decoys = check_arg_trueFalse(args.decoys)
    args.monoisotopic_precursor = check_arg_trueFalse(args.monoisotopic_precursor)

    # formulate digestion regular expression
    digest_re = create_digest_re(args.enzyme, args.custom_enzyme)
    print digest_re
    r = re.compile(digest_re)
    
    # parse modifications
    mods = parse_mods(args.mods_spec, True)
    print mods
    nterm_mods = parse_mods(args.nterm_peptide_mods_spec, False)
    print nterm_mods
    cterm_mods = parse_mods(args.cterm_peptide_mods_spec, False)
    print cterm_mods

    args.digestion = args.digestion.lower()
    if args.digestion != "full-digest" and args.digestion != "partial-digest":
        print "Option --digestion specified as %s, must be either full-digest or partial-digest.  Exitting" % (args.digestion)
        exit(-1)

    print args.digestion

    # calculate enzymatic markers
    n_res, c_res = enzyme_enzSet(args.enzyme, args.custom_enzyme)    

    # parse FASTA file and digest proteins
    peptides = load_proteins_transStop(args.fasta, r, 
                                       n_res, c_res,
                                       args.min_length, args.max_length,
                                       args.min_mass, args.max_mass,
                                       mods, nterm_mods, cterm_mods,
                                       args.max_mods, args.min_mods,
                                       args.missed_cleavages,
                                       args.monoisotopic_precursor,
                                       args.digestion)
    # write digested peptides
    write_peptides('targets.txt', peptides)

    # if args.decoys:
    #     decoys = create_decoys(peptides, args.decoy_format, args.keep_terminal_aminos, args.seed)
    #     write_peptides('decoys.txt', decoys)

