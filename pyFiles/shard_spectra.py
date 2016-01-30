#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>, Ajit Singh <ajit.ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

from __future__ import with_statement

__authors__ = [ 'John T. Halloran <halloj3@uw.edu>', 'Ajit Singh <ajit@ee.washington.edu>' ]

import stat
import struct
import cPickle as pickle
import collections
import csv
import itertools
import math
import time
import argparse
import os
import sys
import re

from shutil import rmtree
from cluster import write_cluster_script, write_scratch_cluster_script
from random import shuffle
from bisect import bisect_left, bisect_right
from peptide_db import SimplePeptideDB
from spectrum import MS2Iterator, SampleSpectra
from parsers import load_ident

import digest_fasta as df

def izip_longest(*args, **kwds):
    # izip_longest('ABCD', 'xy', fillvalue='-') --> Ax By C- D-
    fillvalue = kwds.get('fillvalue')
    def sentinel(counter = ([fillvalue]*(len(args)-1)).pop):
        yield counter()         # yields the fillvalue, or raises IndexError
    fillers = itertools.repeat(fillvalue)
    iters = [itertools.chain(it, sentinel(), fillers) for it in args]
    try:
        for tup in itertools.izip(*iters):
            yield tup
    except IndexError:
        pass

def grouper(n, iterable, fillvalue=None):
    """Split a list into evenly sized segments of size n.

    Examples:

    >>> list(grouper(3, 'ABCDEFG'))
    [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', None, None)]

    >>> list(grouper(2, xrange(5), 'x'))
    [(0, 1), (2, 3), (4, 'x')]

    Args:
        n: Length of each segment
        iterable: Sequence to split up into segments
        fillvalue: What to pad the last segment with, if padding is required.

    Returns:
        A list of segments.
    """
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def prune_peptide_buffer(peptides, massKey):
    """
    """
    masses = [p[0] for p in peptides]
    # h = bisect_right(masses, massKey)
    h = bisect_left(masses, massKey)
    del peptides[:h]

def mass_filter(peptides, precursor_mass, window = 3, ppm = False):
    """Return the sequence of the peptides that are within 'window'
    Daltons of the mass given: i.e., where the peptide's mass is within
    [precursor_mass - window, precursor_mass + window]

    Note: avgmass is what is reported on the neutral mass line, and
    that should be a more reasonable mass calculation since the precursor
    mass comes from MS1, which reflects the average mass.

    """
    masses = [p[0] for p in peptides]

    # peptide_strings = [p[1].split('\x00')[0] for p in peptides]

    if not ppm:
        l = bisect_left(masses, precursor_mass - window)
        h = bisect_right(masses, precursor_mass + window, lo = l)
    else:
        wl = window*1e-6
        l = bisect_left(masses, precursor_mass / (1.0 + float(window)*0.000001))
        h = bisect_right(masses, precursor_mass / (1.0 - float(window)*0.000001), lo = l)

    # return peptide_strings[l:h]
    return peptides[l:h]

def load_target_decoy_db(args, r, mods, nterm_mods, cterm_mods):
    digest = 1
    # check if we should load an existing file
    if args.load_peptide_database_file:
        data = df.load_db_pickle(args.peptide_database_file)
        try:
            targets = data['targets']
            digest = 0
        except KeyError:
            print "Target peptides not present in supplied pickle, digesting protein database"
            digest = 1
        # load decoys if necessary
        if args.decoys and digest == 0:
            try:
                decoys = data['decoys']
            except KeyError:
                print "Creating decoy peptides"
                decoys = df.create_decoys(targets, 
                                          args.decoy_format, 
                                          args.keep_terminal_aminos, 
                                          args.seed)
    if digest:
        # parse FASTA file and digest proteins
        targets = df.load_proteins(args.fasta, r, 
                                    args.min_length, args.max_length,
                                    args.min_mass, args.max_mass,
                                    mods, nterm_mods, cterm_mods,
                                    args.max_mods, args.min_mods,
                                    args.monoisotopic_precursor)

        # write targets in ASCII
        if args.target_output_file:
            df.write_peptides(args.target_output_file, targets)

        if args.decoys:
            decoys = df.create_decoys(targets, 
                                      args.decoy_format, 
                                      args.keep_terminal_aminos, 
                                      args.seed)
            # write decoys in ASCII
            if args.decoy_output_file:
                df.write_peptides(args.decoy_output_file, decoys)

            if args.peptide_database_file:
                df.store_db_pickle(args.peptide_database_file, targets, decoys)

        else:
            if args.peptide_database_file:
                df.store_db_pickle(args.peptide_database_file, targets)

    if args.decoys:
        return targets, decoys
    else:
        return targets, []

def load_target_decoy_db_no_reshuffle(args, r, mods, nterm_mods, cterm_mods):
    """ Exit if decoys are desired but not present in the datafile
        Do not rewrite database to output file
    """
    digest = 1
    # check if we should load an existing file
    if args.load_peptide_database_file:
        data = df.load_db_pickle(args.peptide_database_file)
        try:
            targets = data['targets']
            digest = 0
        except KeyError:
            print "Target peptides not present in supplied pickle, digesting protein database"
            digest = 1
        # load decoys if necessary
        if args.decoys and digest == 0:
            try:
                decoys = data['decoys']
            except KeyError:
                print "Could not load decoy database, exitting"
                exit(-1)

    if digest:
        # parse FASTA file and digest proteins
        targets = df.load_proteins(args.fasta, r, 
                                    args.min_length, args.max_length,
                                    args.min_mass, args.max_mass,
                                    mods, nterm_mods, cterm_mods,
                                    args.max_mods, args.min_mods,
                                    args.monoisotopic_precursor)

        if args.decoys:
            print "Decoys are requested but not encountered."
            print "Reshuffling on the fly does not guarantee the same set of run decoys, exitting."
            exit(-1)

    if args.decoys:
        return targets, decoys
    else:
        return targets, []


def calc_minMaxMz(ms2File, charges, ident = None):
    """
    """
    try:
        if ms2File[-3:] == '.gz':
            gzcat = True
        else:
            gzcat = False
    except:
        gzcat = False

    minMz = 2000.0
    maxMz = -1
    # Consider pre-specified charges
    validcharges = set(int(c) for c in charges.split(','))
    for s in MS2Iterator(ms2File, gzcat):
        if (set(s.charges) & validcharges):
            if s.mz[0] < minMz:
                minMz = spec.mz[0]
            if s.mz[-1] > maxMz:
                maxMz = s.mz[-1]
    return minMz, maxMz

def load_spectra(ms2File, charges, ident = None, randOrder = True):
    """
    """
    try:
        if ms2File[-3:] == '.gz':
            gzcat = True
        else:
            gzcat = False
    except:
        gzcat = False
    if charges.lower() == 'all':
        all_c = True
        validcharges = {}
        spectra = list(s for s in MS2Iterator(ms2File, gzcat)) 
    else:
        all_c = False
        validcharges = set(int(c) for c in charges.split(','))
        spectra = list(s for s in MS2Iterator(ms2File, gzcat) if
                       set(s.charges) & validcharges) 

    print >> sys.stderr, '%s has %d spectra with charges %s' % (
        ms2File, len(spectra),
        '{' + ','.join('+%d' % c for c in validcharges) + '}')

    minMz = 2000.0
    maxMz = -1
    if ident:
        print "Filtering spectra from ident %s" % ident
        targets, decoys = load_ident(ident)
        sids = [ t[0] for t in targets ]

        s = []
        for spec in spectra:
            if(spec.spectrum_id in sids):
                s.append(spec)
                if all_c:
                    for c,m in spec.charge_lines:
                        validcharges[c] = 1

            # find minimum and maximum mz
            if spec.mz[0] < minMz:
                minMz = spec.mz[0]
            if spec.mz[-1] > maxMz:
                maxMz = spec.mz[-1]
    else:
        s = []
        for spec in spectra:
            s.append(spec)
            if all_c:
                for c,m in spec.charge_lines:
                    validcharges[c] = 1
            # find minimum and maximum mz
            if spec.mz[0] < minMz:
                minMz = spec.mz[0]
            if spec.mz[-1] > maxMz:
                maxMz = spec.mz[-1]
                
    spectra = s
    # print "minMz=%f" % minMz
    # print "maxMz=%f" % maxMz

    if len(spectra) == 0:
        print >> sys.stderr, 'There are no spectra with charges %s.' % charges
        exit(-1)

    if randOrder:
        shuffle(spectra)

    return spectra, minMz, maxMz, validcharges

def load_spectra_ret_dict(ms2File, charges):
    """
    """
    try:
        if ms2File[-3:] == '.gz':
            gzcat = True
        else:
            gzcat = False
    except:
        gzcat = False
    if charges.lower() == 'all':
        all_c = True
        validcharges = {}
        spectra = list(s for s in MS2Iterator(ms2File, gzcat)) 
    else:
        all_c = False
        validcharges = set(int(c) for c in charges.split(','))
        spectra = list(s for s in MS2Iterator(ms2File, gzcat) if
                       set(s.charges) & validcharges) 

    sidChargePreMass = [] # tuple of sid, charge, precursor mass

    minMz = 2000.0
    maxMz = -1
    s = {}
    for spec in spectra:
        s[spec.spectrum_id] = spec
        if all_c:
            for c,m in spec.charge_lines:
                sidChargePreMass.append((spec.spectrum_id, c, m))
                validcharges[c] = 1
        else:
            for c,m in spec.charge_lines:
                if c in validcharges:
                    sidChargePreMass.append((spec.spectrum_id, c, m))
        # find minimum and maximum mz
        if spec.mz[0] < minMz:
            minMz = spec.mz[0]
        if spec.mz[-1] > maxMz:
            maxMz = spec.mz[-1]
                
    if len(s) == 0:
        print >> sys.stderr, 'There are no spectra with charges %s.' % charges
        exit(-1)

    sidChargePreMass.sort(key = lambda r: r[2])

    return s, minMz, maxMz, validcharges, sidChargePreMass

def pickle_candidate_spectra(args, r, 
                             mods, nterm_mods, cterm_mods):
    """ Create pickles for each spectrum, each containing the observed spectrum and
        target/decoy candidates
    """
    targets, decoys = load_target_decoy_db(args, r, mods, nterm_mods, cterm_mods)
    # create target/decoy database instances
    targets = SimplePeptideDB(targets)
    if decoys:
        decoys = SimplePeptideDB(decoys)
    # Compute the number of spectra in each slice.
    # Should be ok to load all spectra into memory so long as we 
    # don't also load all PSMs into memory
    spectra, minMz, maxMz = load_spectra(args.spectra,
                                         args.charges, 
                                         args.ident, 
                                         True)
    n = len(spectra)
    if n < args.shards:
        print('More shards than spectra: %d vs. %d, defaulting to %d shards' % (
        args.shards, n, max(int(n/4), 1)))
        args.shards = max(int(n/4), 1)
    # assert n >= args.shards, 'More shards than spectra: %d vs. %d' % (
    #     args.shards, n)
    sz = int(math.floor(float(n)/args.shards))
    print args.shards
    print >> sys.stderr, 'Each shard has at most %d spectra' % sz


    # calculate candidate peptides, create pickles
    mass_h = 1.00727646677
    validcharges = set(int(c) for c in args.charges.split(','))
    for part, group in enumerate(grouper(sz, spectra)):
        #print >> sys.stderr, '%d' % part
        spectra_app = []
        spectra_list = list(s for s in group if s)
        emptySpectra = 0
        # Use neutral mass to select peptides (Z-lines report M+H+ mass)
        t = collections.defaultdict(list)
        if args.decoys:
            d = collections.defaultdict(list)
        for s in spectra_list:
            spec_added = 0
            repSpec = 0
            for c, m in s.charge_lines:
                if c in validcharges:
                    tc = targets.filter(m - mass_h, args.precursor_window, args.ppm)
                    if args.decoys:
                        dc = decoys.filter(m - mass_h, args.precursor_window, args.ppm)
                        if len(tc) > 0 and len(dc) > 0:
                            print "sid=%d, charge=%d, num targets=%d, num decoys=%d" % (s.spectrum_id, c, len(tc), len(dc))
                            if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                                t[(s.spectrum_id,c)] += tc
                                d[(s.spectrum_id,c)] += dc
                                repSpec = 1
                            else:
                                t[(s.spectrum_id,c)] = tc
                                d[(s.spectrum_id,c)] = dc
                            if not spec_added:
                                spectra_app.append(s)
                                spec_added = 1
                        else:
                            emptySpectra += 1
                    else:
                        if len(tc) > 0:
                            print "sid=%d, charge=%d, num targets=%d" % (s.spectrum_id, c, len(tc))
                            if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                                t[(s.spectrum_id,c)] += tc
                                repSpec = 1
                            else:
                                t[(s.spectrum_id,c)] = tc
                            if not spec_added:
                                spectra_app.append(s)
                                spec_added = 1
                        else:
                            emptySpectra += 1

            if repSpec: # prune any multiply added peptide candidates
                for c in validcharges:
                    t[(s.spectrum_id, c)] = list(set(t[(s.spectrum_id, c)]))
                    if args.decoys:
                        d[(s.spectrum_id, c)] = list(set(d[(s.spectrum_id, c)]))
                

        print "%d spectra with empty candidate sets" % emptySpectra

        if args.decoys:
            data = { 'spectra' : spectra_app,
                     'target' : t,
                     'decoy' : d,
                     'shardname' : 'shard-%d' % part }
        else:
            data = { 'spectra' : spectra_app,
                     'target' : t,
                     'shardname' : 'shard-%d' % part }

        outfn = os.path.join(args.output_dir, 'shard-%d.pickle' % part)
        with open(outfn, 'w') as out:
            pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)

def candidate_spectra_generator(args, r,
                                mods, nterm_mods, cterm_mods):
    """ Load all spectra and candidate targets/decoys into memory,
        yield each spectrum and candidates
    """

    targets, decoys = load_target_decoy_db(args, r, mods, nterm_mods, cterm_mods)
    # create target/decoy database instances
    targets = SimplePeptideDB(targets)
    if decoys:
        decoys = SimplePeptideDB(decoys)
    # Compute the number of spectra in each slice.
    # Should be ok to load all spectra into memory so long as we 
    # don't also load all PSMs into memory
    spectra, minMz, maxMz, validcharges = load_spectra(args.spectra,
                                                       args.charges, 
                                                       args.ident,
                                                       True)

    # set global variable to parsed/all encountered charges
    args.charges = validcharges

    n = len(spectra)
    if n < args.shards:
        print('More shards than spectra: %d vs. %d, defaulting to %d shards' % (
        args.shards, n, max(int(n/4), 1)))
        args.shards = max(int(n/4), 1)
    sz = int(math.floor(float(n)/args.shards))
    print args.shards
    print >> sys.stderr, 'Each shard has at most %d spectra' % sz

    # calculate candidate peptides, create pickles
    mass_h = 1.00727646677
    for part, group in enumerate(grouper(sz, spectra)):
        spectra_app = []
        spectra_list = list(s for s in group if s)
        emptySpectra = 0
        # Use neutral mass to select peptides (Z-lines report M+H+ mass)
        t = collections.defaultdict(list)
        if args.decoys:
            d = collections.defaultdict(list)
        for s in spectra_list:
            spec_added = 0
            repSpec = 0
            for c, m in s.charge_lines:
                if c in validcharges:
                    tc = targets.filter(m - mass_h, args.precursor_window, args.ppm)
                    if args.decoys:
                        dc = decoys.filter(m - mass_h, args.precursor_window, args.ppm)
                        if len(tc) > 0 and len(dc) > 0:
                            print "sid=%d, charge=%d, num targets=%d, num decoys=%d" % (s.spectrum_id, c, len(tc), len(dc))
                            if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                                t[(s.spectrum_id,c)] += tc
                                d[(s.spectrum_id,c)] += dc
                                repSpec = 1
                            else:
                                t[(s.spectrum_id,c)] = tc
                                d[(s.spectrum_id,c)] = dc
                            if not spec_added:
                                spectra_app.append(s)
                                spec_added = 1
                        else:
                            emptySpectra += 1
                    else:
                        if len(tc) > 0:
                            print "sid=%d, charge=%d, num targets=%d" % (s.spectrum_id, c, len(tc))
                            if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                                t[(s.spectrum_id,c)] += tc
                                repSpec = 1
                            else:
                                t[(s.spectrum_id,c)] = tc
                            if not spec_added:
                                spectra_app.append(s)
                                spec_added = 1
                        else:
                            emptySpectra += 1

            if repSpec: # prune any multiply added peptide candidates
                for c in validcharges:
                    if (s.spectrum_id, c) in t:
                        t[(s.spectrum_id, c)] = list(set(t[(s.spectrum_id, c)]))
                        if args.decoys:
                            d[(s.spectrum_id, c)] = list(set(d[(s.spectrum_id, c)]))

        print "%d spectra with empty candidate sets" % emptySpectra
        if args.decoys:
            data = { 'spectra' : spectra_app,
                     'target' : t,
                     'decoy' : d,
                     'minMz' : minMz, 
                     'maxMz' : maxMz}
        else:
            data = { 'spectra' : spectra_app,
                     'target' : t,
                     'minMz' : minMz, 
                     'maxMz' : maxMz}
        yield data

def find_curr_candidate_peps(t, tfid, m, precursor_window, ppm,
                             lowerMassKey, upperMassKey,
                             stst, t_l_buffered, buff):
    # load in targets
    run = 1
        # see if old targets are within mass range
    if t:
        if t[-1][0] < lowerMassKey: # none are in mass range
            del t[:]
        elif t[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
            # subsequent peptides to be read are also greater
            print "Cached peptides are greater than current precursor mass, exitting"
            exit(-1)
        else: # some must be in mass range, determine which and clear some of the buffer
            # do bisection
            prune_peptide_buffer(t, lowerMassKey)

    # first check buffered peptides
    if t_l_buffered:
        start_ind = 0
        end_ind = stst.size
        for i in range(len(t_l_buffered) / stst.size):
            p = stst.unpack(t_l_buffered[start_ind:end_ind])
            start_ind += stst.size
            end_ind += stst.size
            if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                t.append(p)
            elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                # del t_l_buffered[:start_ind]
                start_ind -= stst.size
                t_l_buffered = t_l_buffered[start_ind:]
                run = 0 
                break

    while run:
        # read in buffer of peptides
        l = tfid.read(stst.size * buff)
        start_ind = 0
        end_ind = stst.size
        t_l_buffered = []
        for i in range(len(l) / stst.size):
            p = stst.unpack(l[start_ind:end_ind])
            start_ind += stst.size
            end_ind += stst.size
            if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                t.append(p)
            elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                # include current peptide in buffer
                start_ind -= stst.size
                t_l_buffered = l[start_ind:] # buffer remaining peptides for next iteration
                run = 0
                break
            # else, just keep reading in peptides from the binary file until
            # we get within the mass range

    return mass_filter(t, m, precursor_window, ppm)                

# def candidate_binarydb_spectra_generator(args,
#                                          mods, nterm_mods, cterm_mods,
#                                          var_mods, nterm_var_mods, cterm_var_mods):
#     """ Binary database assumed to be created.  Load spectra into memory, sort by 
#         precursor mass.  Simultaneously traverse binary peptide database(s) (both target
#         and decoy, if decoys have been generated) in order of succeeding mass, thus keeping
#         memory usage down for when the number of peptides is exponentially large compared
#         to the number of protein strings in the FASTA database (i.e., many missed cleavages,
#         partial-digestion, variable modifications).
#     """

#     # Compute the number of spectra in each slice.
#     # Should be ok to load all spectra into memory so long as we 
#     # don't also load all PSMs into memory

#     buff = args.peptide_buffer

#     # unfortunately, we can no longer shuffle spectra
#     spectra, minMz, maxMz, validcharges, sidChargePreMass = load_spectra_ret_dict(args.spectra,
#                                                                                   args.charges)

#     # set global variable to parsed/all encountered charges
#     args.charges = validcharges

#     t_buff = []
#     d_buff = []

#     if not var_mods and not nterm_var_mods and not cterm_var_mods:
#         stst = struct.Struct('f %ds I s s' % (args.max_length))
#     else:    
#         # we have some variable mods, so var_mods, nterm_var_mods, or cterm_var_mods are non-empty
#         # serialized data:
#         #    peptide[0] = peptide mass (float)
#         #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
#         #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
#         #    peptide[3] = nterm_flanking (character)
#         #    peptide[4] = cterm_flanking (character)
#         #    peptide[5] = binary string deoting variable modifications
#         stst = struct.Struct('f %ds I s s %ds' % (args.max_length, args.max_length))


#     t_l_buffered = []

#     d_l_buffered = []

#     t = []
#     d = []

#     t_out = {}
#     d_out = {}

#     spectra_app = []
#     spectra_buffer = 1000

#     buffered_spec = 0
    
#     tfid = open(os.path.join(args.digest_dir, 'targets.bin'), 'rb')
#     if args.decoys:
#         dfid = open(os.path.join(args.digest_dir, 'decoys.bin'), 'rb')
#         if args.recalibrate:
#             recal_dfid = open(os.path.join(args.digest_dir, 'recalibrateDecoys.bin'), 'rb')
#             recal_d_out = {}

#     # calculate candidate peptides, yield a buffer of spectra and candidates
#     mass_h = 1.00727646677

#     if not args.ppm:
#         tol_buf = int(math.ceil(args.precursor_window))
#     else:
#         tol_buf = 1

#     sid_visited = set([])

#     for sid, c, pm in sidChargePreMass:

#         if buffered_spec > spectra_buffer: # clear memory
#             d_out.clear()
#             t_out.clear()
#             del spectra_app[:]
#             buffered_spec = 0 # this will have been yielded at the end of the previous iteration

#             if args.recalibrate:
#                 recal_d_out.clear()

#         m = pm - mass_h

#         # would use read_binary_peptide_tuples_buffer function for the following, but much
#         # care is needed to carefully traverse the lists of targets and decoys

#         lowerMassKey = math.floor(m-tol_buf)
#         upperMassKey = math.ceil(m+tol_buf)

#         # print "%d, %f" % (sid, pm)
#         # if t_l_buffered:
#         #     p = stst.unpack(t_l_buffered[0:stst.size])
#         #     print "%d, %d, %f, %f, %f, %f, %d" % (sid, c, m, lowerMassKey, upperMassKey, p[0], len(t))

#         curr_peps = find_curr_candidate_peps(t, tfid, m, args.precursor_window, args.ppm,
#                                              lowerMassKey, upperMassKey,
#                                              stst, t_l_buffered,
#                                              buff)
#         if curr_peps:
#             t_out[sid, c] = curr_peps

#         curr_peps = find_curr_candidate_peps(d, dfid, m, args.precursor_window, args.ppm,
#                                              lowerMassKey, upperMassKey,
#                                              stst, d_l_buffered,
#                                              buff)
#         if curr_peps:
#             d_out[sid, c] = curr_peps

#         if sid not in sid_visited:
#             spectra_app.append(spectra[sid])
#             sid_visited.add(sid)

#         buffered_spec += 1

#         if buffered_spec == spectra_buffer:
#             data = { 'spectra' : spectra_app,
#                      'target' : t_out,
#                      'decoy' : d_out,
#                      'minMz' : minMz, 
#                      'maxMz' : maxMz}
#             yield data


#     tfid.close()
#     dfid.close()
#     if buffered_spec:
#         data = { 'spectra' : spectra_app,
#                  'target' : t_out,
#                  'decoy' : d_out,
#                  'minMz' : minMz, 
#                  'maxMz' : maxMz}
#         yield data        

def pickle_candidate_binarydb_spectra(args,
                                      mods, nterm_mods, cterm_mods,
                                      var_mods, nterm_var_mods, cterm_var_mods):
    """ Binary database assumed to be created.  Load spectra into memory, sort by 
        precursor mass.  Simultaneously traverse binary peptide database(s) (both target
        and decoy, if decoys have been generated) in order of succeeding mass, thus keeping
        memory usage down for when the number of peptides is exponentially large compared
        to the number of protein strings in the FASTA database (i.e., many missed cleavages,
        partial-digestion, variable modifications).
    """

    # Compute the number of spectra in each slice.
    # Should be ok to load all spectra into memory so long as we 
    # don't also load all PSMs into memory

    # unfortunately, we can no longer shuffle spectra
    spectra, minMz, maxMz, validcharges, sidChargePreMass = load_spectra_ret_dict(args.spectra,
                                                                                  args.charges)
    # set global variable to parsed/all encountered charges
    args.charges = validcharges

    buff = args.peptide_buffer

    if args.write_cluster_scripts:
        script_list = os.path.join(args.output_dir, 'clusterJobs.txt')
        listFid = open(script_list, 'w')

    if not var_mods and not nterm_var_mods and not cterm_var_mods:
        # only static mods
        # serialized data:
        #    peptide[0] = peptide mass (float)
        #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
        #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
        #    peptide[3] = nterm_flanking (character)
        #    peptide[4] = cterm_flanking (character)
        stst = struct.Struct('f %ds I s s' % (args.max_length))
    else:    
        # we have some variable mods, so var_mods, nterm_var_mods, or cterm_var_mods are non-empty
        # serialized data:
        #    peptide[0] = peptide mass (float)
        #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
        #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
        #    peptide[3] = nterm_flanking (character)
        #    peptide[4] = cterm_flanking (character)
        #    peptide[5] = binary string deoting variable modifications
        stst = struct.Struct('f %ds I s s %ds' % (args.max_length, args.max_length))

    t_l_buffered = ''
    d_l_buffered = ''

    t = []
    d = []

    t_out = {}
    d_out = {}

    spectra_app = []
    spectra_buffer = int(math.ceil(float(len(sidChargePreMass))/args.num_jobs))

    buffered_spec = 0
    
    tfid = open(os.path.join(args.digest_dir, 'targets.bin'), 'rb')
    if args.decoys:
        dfid = open(os.path.join(args.digest_dir, 'decoys.bin'), 'rb')
        if args.recalibrate:
            recal_dfid = open(os.path.join(args.digest_dir, 'recalibrateDecoys.bin'), 'rb')
            recal_d_l_buffered = ''
            recal_d = []
            recal_d_out = {}

    # calculate candidate peptides, yield a buffer of spectra and candidates
    mass_h = 1.00727646677

    if not args.ppm:
        tol_buf = int(math.ceil(args.precursor_window))
    else:
        tol_buf = 1

    sid_visited = set([])

    num_pickle = 0

    for sid, c, pm in sidChargePreMass:
        # if buffered_spec >= spectra_buffer: # clear memory
        #     d_out.clear()
        #     t_out.clear()
        #     del spectra_app[:]
        #     buffered_spec = 0 # this will have been yielded at the end of the previous iteration

        #     if args.recalibrate:
        #         recal_d_out.clear()

        m = pm - mass_h

        lowerMassKey = math.floor(m-tol_buf)
        upperMassKey = math.ceil(m+tol_buf)

        # load in targets
        run = 1
        # see if old targets are within mass range
        if t:
            if t[-1][0] < lowerMassKey: # none are in mass range
                del t[:]
            elif t[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
                # subsequent peptides to be read are also greater
                print "Cached peptides are greater than current precursor mass, exitting"
                exit(-1)
            else: # some must be in mass range, determine which and clear some of the buffer
                # do bisection
                prune_peptide_buffer(t, lowerMassKey)

        # first check buffered peptides
        if t_l_buffered:
            start_ind = 0
            end_ind = stst.size
            for i in range(len(t_l_buffered) / stst.size):
                p = stst.unpack(t_l_buffered[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    t.append(p)
                elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                    # del t_l_buffered[:start_ind]
                    start_ind -= stst.size
                    t_l_buffered = t_l_buffered[start_ind:]
                    run = 0 
                    break
                # else, keep checking the buffered peptides til we get into the mass range

        while run:
            # read in buffer of peptides
            l = tfid.read(stst.size * buff)
            start_ind = 0
            end_ind = stst.size
            t_l_buffered = ''
            for i in range(len(l) / stst.size):
                p = stst.unpack(l[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    t.append(p)
                elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                    # include current peptide in buffer
                    start_ind -= stst.size
                    t_l_buffered = l[start_ind:] # buffer remaining peptides for next iteration
                    run = 0
                    break
                # else, just keep reading in peptides from the binary file until
                # we get within the mass range

        curr_peps = mass_filter(t, m, args.precursor_window, args.ppm)                
        if curr_peps:
            t_out[sid, c] = curr_peps

        # load in decoys
        run = 1
        # see if old decoys are within mass range
        if d:
            if d[-1][0] < lowerMassKey: # none are in mass range
                del d[:]
            elif d[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
                # subsequent peptides to be read are also greater
                print "Cached peptides are greater than current precursor mass, exitting"
                exit(-1)
            else: # some must be in mass range, determine which and clear some of the buffer
                # do bisection
                prune_peptide_buffer(d, lowerMassKey)

        # first checked buffered peptides
        if d_l_buffered:
            start_ind = 0
            end_ind = stst.size
            for i in range(len(d_l_buffered) / stst.size):
                p = stst.unpack(d_l_buffered[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    d.append(p)
                elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                    start_ind -= stst.size
                    d_l_buffered = d_l_buffered[start_ind:]                    
                    run = 0
                    break
                # else, keep checking the buffered peptides til we get into the mass range
        while run:
            # read in buffer of peptides
            l = dfid.read(stst.size * buff)
            start_ind = 0
            end_ind = stst.size
            d_l_buffered = ''
            for i in range(len(l) / stst.size):
                p = stst.unpack(l[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    d.append(p)
                elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                    start_ind -= stst.size
                    d_l_buffered = l[start_ind:]
                    run = 0
                    break
                # else, just keep reading in peptides from the binary file until
                # we get within the mass range

        curr_peps = mass_filter(d, m, args.precursor_window, args.ppm)
        if curr_peps:
            d_out[sid, c] = curr_peps

        if args.recalibrate:
            # load in recalibration decoys
            run = 1
            # see if old decoys are within mass range
            if recal_d:
                if recal_d[-1][0] < lowerMassKey: # none are in mass range
                    del recal_d[:]
                elif recal_d[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
                    # subsequent peptides to be read are also greater
                    print "Cached peptides are greater than current precursor mass, exitting"
                    exit(-1)
                else: # some must be in mass range, determine which and clear some of the buffer
                    # do bisection
                    prune_peptide_buffer(recal_d, lowerMassKey)

            # first checked buffered peptides
            if recal_d_l_buffered:
                start_ind = 0
                end_ind = stst.size
                for i in range(len(recal_d_l_buffered) / stst.size):
                    p = stst.unpack(recal_d_l_buffered[start_ind:end_ind])
                    start_ind += stst.size
                    end_ind += stst.size
                    if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                        recal_d.append(p)
                    elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                        start_ind -= stst.size
                        recal_d_l_buffered = recal_d_l_buffered[start_ind:]                    
                        run = 0
                        break
                    # else, keep checking the buffered peptides til we get into the mass range
            while run:
                # read in buffer of peptides
                l = recal_dfid.read(stst.size * buff)
                start_ind = 0
                end_ind = stst.size
                recal_d_l_buffered = ''
                for i in range(len(l) / stst.size):
                    p = stst.unpack(l[start_ind:end_ind])
                    start_ind += stst.size
                    end_ind += stst.size
                    if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                        recal_d.append(p)
                    elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                        start_ind -= stst.size
                        recal_d_l_buffered = l[start_ind:]
                        run = 0
                        break
                    # else, just keep reading in peptides from the binary file until
                    # we get within the mass range

            curr_peps = mass_filter(recal_d, m, args.precursor_window, args.ppm)
            if curr_peps:
                recal_d_out[sid, c] = curr_peps

        spectra_app.append(spectra[sid])
        # if sid not in sid_visited:
        #     spectra_app.append(spectra[sid])
        #     sid_visited.add(sid)

        buffered_spec += 1

        if buffered_spec == spectra_buffer:
            if args.recalibrate:
                data = { 'args' : args,
                         'mods' : mods,
                         'nterm_mods' : nterm_mods,
                         'cterm_mods' : cterm_mods,
                         'var_mods' : var_mods,
                         'nterm_var_mods' : nterm_var_mods,
                         'cterm_var_mods' : cterm_var_mods,
                         'spectra' : spectra_app,
                         'target' : t_out,
                         'decoy' : d_out,
                         'recal_decoy' : recal_d_out,
                         'minMz' : minMz, 
                         'maxMz' : maxMz}
            else:
                data = { 'args' : args, 
                         'mods' : mods,
                         'nterm_mods' : nterm_mods,
                         'cterm_mods' : cterm_mods,
                         'var_mods' : var_mods,
                         'nterm_var_mods' : nterm_var_mods,
                         'cterm_var_mods' : cterm_var_mods,
                         'spectra' : spectra_app,
                         'target' : t_out,
                         'decoy' : d_out,
                         'minMz' : minMz, 
                         'maxMz' : maxMz}

            outfn = os.path.join(args.output_dir, 'split%d.pickle' % num_pickle)
            with open(outfn, 'w') as out:
                pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)

            if args.write_cluster_scripts:
                write_cluster_script(args, outfn, num_pickle, args.cluster_dir)
                # write script to master list
                print >> listFid, "%s" % os.path.join(os.path.abspath(args.output_dir), 'split%d.sh' % num_pickle)

                os.chmod(os.path.join(os.path.abspath(args.output_dir), 'split%d.sh' % num_pickle), stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR )

            num_pickle += 1
            # clear buffers
            d_out.clear()
            t_out.clear()
            del spectra_app[:]
            buffered_spec = 0 # this will have been yielded at the end of the previous iteration

            if args.recalibrate:
                recal_d_out.clear()

    tfid.close()
    dfid.close()
    if args.recalibrate:
        recal_dfid.close()

    if buffered_spec:
        if args.recalibrate:
            data = { 'args' : args,
                     'mods' : mods,
                     'nterm_mods' : nterm_mods,
                     'cterm_mods' : cterm_mods,
                     'var_mods' : var_mods,
                     'nterm_var_mods' : nterm_var_mods,
                     'cterm_var_mods' : cterm_var_mods,
                     'spectra' : spectra_app,
                     'target' : t_out,
                     'decoy' : d_out,
                     'recal_decoy' : recal_d_out,
                     'minMz' : minMz, 
                     'maxMz' : maxMz}
        else:
            data = { 'args' : args,
                     'mods' : mods,
                     'nterm_mods' : nterm_mods,
                     'cterm_mods' : cterm_mods,
                     'var_mods' : var_mods,
                     'nterm_var_mods' : nterm_var_mods,
                     'cterm_var_mods' : cterm_var_mods,
                     'spectra' : spectra_app,
                     'target' : t_out,
                     'decoy' : d_out,
                     'minMz' : minMz, 
                     'maxMz' : maxMz}
        outfn = os.path.join(args.output_dir, 'split%d.pickle' % num_pickle)
        with open(outfn, 'w') as out:
            pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)

        if args.write_cluster_scripts:
            write_cluster_script(args, outfn, num_pickle, args.cluster_dir)
            # write script to master list
            print >> listFid, "%s" % os.path.join(os.path.abspath(args.output_dir), 'split%d.sh' % num_pickle)

            os.chmod(os.path.join(os.path.abspath(args.output_dir), 'split%d.sh' % num_pickle), stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR )

            listFid.close()


def candidate_binarydb_spectra_generator(args,
                                         mods, nterm_mods, cterm_mods,
                                         var_mods, nterm_var_mods, cterm_var_mods):
    """ Binary database assumed to be created.  Load spectra into memory, sort by 
        precursor mass.  Simultaneously traverse binary peptide database(s) (both target
        and decoy, if decoys have been generated) in order of succeeding mass, thus keeping
        memory usage down for when the number of peptides is exponentially large compared
        to the number of protein strings in the FASTA database (i.e., many missed cleavages,
        partial-digestion, variable modifications).
    """

    # Compute the number of spectra in each slice.
    # Should be ok to load all spectra into memory so long as we 
    # don't also load all PSMs into memory

    # unfortunately, we can no longer shuffle spectra
    spectra, minMz, maxMz, validcharges, sidChargePreMass = load_spectra_ret_dict(args.spectra,
                                                                                  args.charges)

    # set global variable to parsed/all encountered charges
    args.charges = validcharges

    buff = args.peptide_buffer

    if not var_mods and not nterm_var_mods and not cterm_var_mods:
        stst = struct.Struct('f %ds I s s' % (args.max_length))
    else:    
        # we have some variable mods, so var_mods, nterm_var_mods, or cterm_var_mods are non-empty
        # serialized data:
        #    peptide[0] = peptide mass (float)
        #    peptide[1] = peptide string (string of max_length character, possibly many of which are null)
        #    peptide[2] = protein name (mapped to an integer for the protein value encountered in the file)
        #    peptide[3] = nterm_flanking (character)
        #    peptide[4] = cterm_flanking (character)
        #    peptide[5] = binary string deoting variable modifications
        stst = struct.Struct('f %ds I s s %ds' % (args.max_length, args.max_length))


    t_l_buffered = ''
    d_l_buffered = ''

    t = []
    d = []

    t_out = {}
    d_out = {}

    spectra_app = []
    spectra_buffer = 1000

    buffered_spec = 0
    
    tfid = open(os.path.join(args.digest_dir, 'targets.bin'), 'rb')
    if args.decoys:
        dfid = open(os.path.join(args.digest_dir, 'decoys.bin'), 'rb')
        if args.recalibrate:
            recal_dfid = open(os.path.join(args.digest_dir, 'recalibrateDecoys.bin'), 'rb')
            recal_d_l_buffered = ''
            recal_d = []
            recal_d_out = {}

    # calculate candidate peptides, yield a buffer of spectra and candidates
    mass_h = 1.00727646677

    if not args.ppm:
        tol_buf = int(math.ceil(args.precursor_window))
    else:
        tol_buf = 1

    sid_visited = set([])

    for sid, c, pm in sidChargePreMass:
        # if buffered_spec > spectra_buffer: # clear memory
        #     d_out.clear()
        #     t_out.clear()
        #     del spectra_app[:]
        #     buffered_spec = 0 # this will have been yielded at the end of the previous iteration

        #     if args.recalibrate:
        #         recal_d_out.clear()

        m = pm - mass_h

        # would use read_binary_peptide_tuples_buffer function for the following, but much
        # care is needed to carefully traverse the lists of targets and decoys

        lowerMassKey = math.floor(m-tol_buf)
        upperMassKey = math.ceil(m+tol_buf)

        # print "%d, %f" % (sid, pm)
        # if t_l_buffered:
        #     p = stst.unpack(t_l_buffered[0:stst.size])
        #     print "%d, %d, %f, %f, %f, %f, %d" % (sid, c, m, lowerMassKey, upperMassKey, p[0], len(t))

        # load in targets
        run = 1
        # see if old targets are within mass range
        if t:
            if t[-1][0] < lowerMassKey: # none are in mass range
                del t[:]
            elif t[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
                # subsequent peptides to be read are also greater
                print "Cached peptides are greater than current precursor mass, exitting"
                exit(-1)
            else: # some must be in mass range, determine which and clear some of the buffer
                # do bisection
                prune_peptide_buffer(t, lowerMassKey)

        # if sid==48 and c==2:
        #     print len(t)
        #     print m
        #     print lowerMassKey
        #     print upperMassKey
        #     print len(t_l_buffered)

        # first check buffered peptides
        if t_l_buffered:
            start_ind = 0
            end_ind = stst.size
            for i in range(len(t_l_buffered) / stst.size):
                p = stst.unpack(t_l_buffered[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    t.append(p)
                elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                    # del t_l_buffered[:start_ind]
                    start_ind -= stst.size
                    t_l_buffered = t_l_buffered[start_ind:]
                    run = 0 
                    break
                # else, keep checking the buffered peptides til we get into the mass range

        while run:
            # read in buffer of peptides
            l = tfid.read(stst.size * buff)
            start_ind = 0
            end_ind = stst.size
            t_l_buffered = ''
            for i in range(len(l) / stst.size):
                p = stst.unpack(l[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    t.append(p)
                elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                    # include current peptide in buffer
                    start_ind -= stst.size
                    t_l_buffered = l[start_ind:] # buffer remaining peptides for next iteration
                    run = 0
                    break
                # else, just keep reading in peptides from the binary file until
                # we get within the mass range

        curr_peps = mass_filter(t, m, args.precursor_window, args.ppm)                
        if curr_peps:
            t_out[sid, c] = curr_peps

        # load in decoys
        run = 1
        # see if old decoys are within mass range
        if d:
            if d[-1][0] < lowerMassKey: # none are in mass range
                del d[:]
            elif d[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
                # subsequent peptides to be read are also greater
                print "Cached peptides are greater than current precursor mass, exitting"
                exit(-1)
            else: # some must be in mass range, determine which and clear some of the buffer
                # do bisection
                prune_peptide_buffer(d, lowerMassKey)

        # first checked buffered peptides
        if d_l_buffered:
            start_ind = 0
            end_ind = stst.size
            for i in range(len(d_l_buffered) / stst.size):
                p = stst.unpack(d_l_buffered[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    d.append(p)
                elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                    start_ind -= stst.size
                    d_l_buffered = d_l_buffered[start_ind:]                    
                    run = 0
                    break
                # else, keep checking the buffered peptides til we get into the mass range
        while run:
            # read in buffer of peptides
            l = dfid.read(stst.size * buff)
            start_ind = 0
            end_ind = stst.size
            d_l_buffered = ''
            for i in range(len(l) / stst.size):
                p = stst.unpack(l[start_ind:end_ind])
                start_ind += stst.size
                end_ind += stst.size
                if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                    d.append(p)
                elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                    start_ind -= stst.size
                    d_l_buffered = l[start_ind:]
                    run = 0
                    break
                # else, just keep reading in peptides from the binary file until
                # we get within the mass range

        curr_peps = mass_filter(d, m, args.precursor_window, args.ppm)
        if curr_peps:
            d_out[sid, c] = curr_peps

        if args.recalibrate:
            # load in recalibration decoys
            run = 1
            # see if old decoys are within mass range
            if recal_d:
                if recal_d[-1][0] < lowerMassKey: # none are in mass range
                    del recal_d[:]
                elif recal_d[0][0] > upperMassKey: # this shouldn't happen since this would mean all 
                    # subsequent peptides to be read are also greater
                    print "Cached peptides are greater than current precursor mass, exitting"
                    exit(-1)
                else: # some must be in mass range, determine which and clear some of the buffer
                    # do bisection
                    prune_peptide_buffer(recal_d, lowerMassKey)

            # first checked buffered peptides
            if recal_d_l_buffered:
                start_ind = 0
                end_ind = stst.size
                for i in range(len(recal_d_l_buffered) / stst.size):
                    p = stst.unpack(recal_d_l_buffered[start_ind:end_ind])
                    start_ind += stst.size
                    end_ind += stst.size
                    if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                        recal_d.append(p)
                    elif p[0] > upperMassKey: # we're outside the mass tolerance range now
                        start_ind -= stst.size
                        recal_d_l_buffered = recal_d_l_buffered[start_ind:]                    
                        run = 0
                        break
                    # else, keep checking the buffered peptides til we get into the mass range
            while run:
                # read in buffer of peptides
                l = recal_dfid.read(stst.size * buff)
                start_ind = 0
                end_ind = stst.size
                recal_d_l_buffered = ''
                for i in range(len(l) / stst.size):
                    p = stst.unpack(l[start_ind:end_ind])
                    start_ind += stst.size
                    end_ind += stst.size
                    if p[0] >= lowerMassKey and p[0] <= upperMassKey:
                        recal_d.append(p)
                    elif p[0] > upperMassKey: # stop reading in peptides, we're outside mass tolerance range now
                        start_ind -= stst.size
                        recal_d_l_buffered = l[start_ind:]
                        run = 0
                        break
                    # else, just keep reading in peptides from the binary file until
                    # we get within the mass range

            curr_peps = mass_filter(recal_d, m, args.precursor_window, args.ppm)
            if curr_peps:
                recal_d_out[sid, c] = curr_peps

        spectra_app.append(spectra[sid])
        # if sid not in sid_visited:
        #     spectra_app.append(spectra[sid])
        #     sid_visited.add(sid)

        buffered_spec += 1

        if buffered_spec == spectra_buffer:
            if args.recalibrate:
                data = { 'spectra' : spectra_app,
                         'target' : t_out,
                         'decoy' : d_out,
                         'recal_decoy' : recal_d_out,
                         'minMz' : minMz, 
                         'maxMz' : maxMz}
            else:
                data = { 'spectra' : spectra_app,
                         'target' : t_out,
                         'decoy' : d_out,
                         'minMz' : minMz, 
                         'maxMz' : maxMz}
            yield data

            # clear buffers
            d_out.clear()
            t_out.clear()
            del spectra_app[:]
            buffered_spec = 0 # this will have been yielded at the end of the previous iteration

            if args.recalibrate:
                recal_d_out.clear()

    tfid.close()
    dfid.close()
    if args.recalibrate:
        recal_dfid.close()

    if buffered_spec:
        if args.recalibrate:
            data = { 'spectra' : spectra_app,
                     'target' : t_out,
                     'decoy' : d_out,
                     'recal_decoy' : recal_d_out,
                     'minMz' : minMz, 
                     'maxMz' : maxMz}
        else:
            data = { 'spectra' : spectra_app,
                     'target' : t_out,
                     'decoy' : d_out,
                     'minMz' : minMz, 
                     'maxMz' : maxMz}
        yield data        

def candidate_spectra_memeffic_generator(args, r,
                                         mods, nterm_mods, cterm_mods):
    """ Generator yielding spectra and target/decoy candidates per spectrum.
        Does not load all spectra into memory
    """

    # parse FASTA file and digest proteins
    peptides = df.load_proteins(args.fasta, r, 
                               args.min_length, args.max_length,
                               args.min_mass, args.max_mass,
                               mods, nterm_mods, cterm_mods,
                               args.max_mods, args.min_mods,
                               args.monoisotopic_precursor)
    df.write_peptides('targets.txt', peptides)

    targets = SimplePeptideDB(peptides)
    if args.decoys:
        decoys = df.create_decoys(peptides, 
                                  args.decoy_format, 
                                  args.keep_terminal_aminos, 
                                  args.seed)
        df.write_peptides('decoys.txt', decoys)
        decoys = SimplePeptideDB(decoys)

    # calculate candidate peptides, create pickles
    mass_h = 1.00727646677
    validcharges = set(int(c) for c in args.charges.split(','))
    for s in MS2Iterator(ms2File, gzcat):
        if set(s.charges) & validcharges:
            repSpec = 0
            t = collections.defaultdict(list)
            if args.decoys:
                d = collections.defaultdict(list)
            for c, m in s.charge_lines:
                if c in validcharges:
                    tc = targets.filter(m - mass_h, args.precursor_window, args.ppm)
                    if args.decoys:
                        dc = decoys.filter(m - mass_h, args.precursor_window, args.ppm)
                        if len(tc) > 0 and len(dc) > 0:
                            if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                                t[(s.spectrum_id,c)] += tc
                                d[(s.spectrum_id,c)] += dc
                                repSpec = 1
                            else:
                                t[(s.spectrum_id,c)] = tc
                                d[(s.spectrum_id,c)] = dc
                    else:
                        if len(tc) > 0:
                            if (s.spectrum_id,c) in t: # high res MS1 may have multiple precursor charge estimates
                                t[(s.spectrum_id,c)] += tc
                                repSpec = 1
                            else:
                                t[(s.spectrum_id,c)] = tc

            if repSpec: # prune any multiply added peptide candidates
                for c in validcharges:
                    t[(s.spectrum_id, c)] = list(set(t[(s.spectrum_id, c)]))
                    if args.decoys:
                        d[(s.spectrum_id, c)] = list(set(d[(s.spectrum_id, c)]))
                
        if args.decoys:
            yield s, t, d
        else:
            yield s, t
