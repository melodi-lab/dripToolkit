#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

from __future__ import with_statement

__authors__ = [ 'John Halloran <halloj3@ee.washington.edu>' ]

import csv
import re
import math
import optparse
import cPickle as pickle
import collections

from operator import itemgetter
from itertools import izip
from pyFiles.peptide import Peptide
from pyFiles.dripEncoding import load_drip_means

# psm = collections.namedtuple('psm',
#                              'kind scan peptide score insertions deletions sumScoredIntensities sumScoredMzDist numFrames delSequence ch')
psm = collections.namedtuple('psm',
                             'kind scan peptide score numBY fragments usedFragments insSequence ch')

def write_psm_ins_dels(psm, s, drip_means, identFid):
    numDels = psm.numBY - len(psm.usedFragments)
    sum_scored_intensities = 0.0
    sum_scored_mz_dists = 0.0
    numIns = 0
    numFrames = len(psm.insSequence)
        
    assert len(s.mz) == numFrames
    
    for intensity,mz,ins,fragment in zip(s.intensity,s.mz,psm.insSequence,psm.fragments):
        if not ins:
            sum_scored_intensities += intensity
            sum_scored_mz_dists += abs(mz - drip_means[fragment])
        else:
            numIns += 1


    # ident file header:
    # (1)Kind(2)Sid(3)Frames(4)Score(5)Peptide(6)Obs_Inserts(7)Theo_Deletes(8)Obs_peaks_scored(9)Theo_peaks_used(10)Sum_obs_intensities
    # (11)Sum_scored_mz_dist(12)Charge
    try:
        identFid.write('%c\t%d\t%d\t%f\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\n' % (psm.kind, 
                                                                             s.spectrum_id,
                                                                             numFrames,
                                                                             psm.score,
                                                                             psm.peptide,
                                                                             numIns,
                                                                             numDels,
                                                                             numFrames - numIns,
                                                                             psm.numBY - numDels,
                                                                             sum_scored_intensities,
                                                                             sum_scored_mz_dists,
                                                                             psm.ch))
    except IOError:
        print "Could not write to ident stream, exitting"
        exit(-1)

    return 1

def parse_dripExtract(filename, pepDBlist):
    """Parse the DRIP Viterbi path output by gmtkViterbi.  This is done spectrum by spectrum where each PSM
       is created as an instance of the psm collection.namedtuple.  Although some of this namedtuple's fields
       involve lists, we only retain the top topMatch per spectrum so that we shouldn't hit any memory problems
    """
    try:
        f = open(filename, 'r')
    except IOError:
        print "Could not open vitValsFile %s for reading, exitting" % filename
        exit(-1)
    log = f.read()
    f.close()

    pattern = 'Segment (?P<segment>\d+), number of frames = (?P<frame>\d+), viterbi-score = (?P<score>\S+)'
    vit_insert_pattern = '.*FRAGMENT_MASS\(\d*\)=(?P<theo_peak>\d+),.*INSERTION\(\d*\)=(?P<insert>\d+)'
    insert = [] # record all inserts
    fragment_masses = []
    for match in re.finditer(vit_insert_pattern, log, re.MULTILINE):
        insert.append(match.group('insert')=='1')
        fragment_masses.append(int(match.group('theo_peak')))

    target_psms = {}
    decoy_psms = {}

    pepLookUp = open(pepDBlist, "r")
    # header: (1)Kind (2)Peptide (3) NumBY (4) Charge
    pepDB = [pepRow for pepRow in csv.DictReader(pepLookUp, delimiter = '\t')]
    pepLookUp.close()
    # for pepRow in csv.DictReader(pepLookUp, delimiter = '\t'):
    #     pepDB.append((pepRow['Kind'], pepRow['Peptide'], pepRow['Charge']))
    
    insert_curr = 0
    count = 0
    
    for match in re.finditer(pattern, log, re.MULTILINE):
        # current GMTK segment number
        curr_segment = int(match.group('segment'))
        # look up peptide in database
        sid = int(pepDB[curr_segment]['Sid'])
        curr_kind = pepDB[curr_segment]['Kind']
        curr_pep_seq = pepDB[curr_segment]['Peptide']

        try:
            curr_numby = int(pepDB[curr_segment]['NumBY'])
        except ValueError:
            print "Trouble trying to convert numBy %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], sid, pepDB[curr_segment]['NumBY'])

        try:
            curr_charge = int(pepDB[curr_segment]['Charge'])
        except ValueError:
            print "Trouble trying to convert charge %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], sid, pepDB[curr_segment]['Peptide'])
            
        # calculate 
        curr_frames = int(match.group('frame'))
        curr_score = float(match.group('score'))/float(curr_frames)
        fragments = []
        used_peaks = {}
        ins_sequence = []

        for i, f in zip(insert[insert_curr:(insert_curr+curr_frames)],fragment_masses[insert_curr:(insert_curr+curr_frames)]):
            ins_sequence.append(i)
            fragments.append(f)
            if not i:
                used_peaks[f] = 1

        insert_curr += curr_frames

        # log all PSM info
        currPsm = psm(kind = curr_kind,
                      scan = sid,
                      ch = curr_charge,
                      peptide = curr_pep_seq,
                      score = curr_score,
                      numBY = curr_numby,
                      fragments = fragments,
                      usedFragments = used_peaks,
                      insSequence = ins_sequence)

        if(curr_kind=='t'):
            if (sid,curr_charge) not in target_psms:
                target_psms[sid,curr_charge] = []
            target_psms[sid,curr_charge].append(currPsm)
        else:
            if (sid,curr_charge) not in decoy_psms:
                decoy_psms[sid,curr_charge] = []
            decoy_psms[sid,curr_charge].append(currPsm)
                           
        count += 1

    return target_psms, decoy_psms

def parse_segments_persid(filename, dripMeans, 
                          topMatch, sid, 
                          currSpec, numPeps, 
                          pepDBlist, identFid,
                          writeDelSequences = False, writeAllDelSequences = False, 
                          delSequencesFid = None):
    """Parse the DRIP Viterbi path output by gmtkViterbi.  This is done spectrum by spectrum where each PSM
       is created as an instance of the psm collection.namedtuple.  Although some of this namedtuple's fields
       involve lists, we only retain the top topMatch per spectrum so that we shouldn't hit any memory problems
    """
    pattern = 'Segment (?P<segment>\d+), number of frames = (?P<frame>\d+), viterbi-score = (?P<score>\S+)'
    try:
        f = open(filename, 'r')
    except IOError:
        return [(None, float("-inf")), (None, float("-inf"))]
    log = f.read()
    f.close()

    vit_insert_pattern = '.*FRAGMENT_MASS\(\d*\)=(?P<theo_peak>\d+),.*INSERTION\(\d*\)=(?P<insert>\d+)'
    insert = [] # record all inserts
    fragment_masses = []
    for match in re.finditer(vit_insert_pattern, log, re.MULTILINE):
        insert.append(match.group('insert')=='1')
        fragment_masses.append(int(match.group('theo_peak')))

    target_psms = {}
    decoy_psms = {}

    pepLookUp = open(pepDBlist, "r")
    # header: (1)Kind (2)Peptide (3) NumBY (4) Charge
    pepDB = [pepRow for pepRow in csv.DictReader(pepLookUp, delimiter = '\t')]
    pepLookUp.close()
    # for pepRow in csv.DictReader(pepLookUp, delimiter = '\t'):
    #     pepDB.append((pepRow['Kind'], pepRow['Peptide'], pepRow['Charge']))
    
    insert_curr = 0
    count = 0
    
    for match in re.finditer(pattern, log, re.MULTILINE):
        # current GMTK segment number
        curr_segment = int(match.group('segment'))
        # look up peptide in database
        curr_kind = pepDB[curr_segment]['Kind']
        curr_pep_seq = pepDB[curr_segment]['Peptide']

        try:
            curr_numby = int(pepDB[curr_segment]['NumBY'])
        except ValueError:
            print "Trouble trying to convert numBy %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], sid, pepDB[curr_segment]['NumBY'])

        try:
            curr_charge = int(pepDB[curr_segment]['Charge'])
        except ValueError:
            print "Trouble trying to convert charge %s to int for sid %d, peptide %s, exitting" % (pepDB[curr_segment]['Charge'], sid, pepDB[curr_segment]['Peptide'])
            
        # calculate 
        curr_frames = int(match.group('frame'))
        curr_score = float(match.group('score'))/float(curr_frames)
        fragments = []
        used_peaks = {}
        ins_sequence = []

        if curr_frames != len(currSpec.mz):
            print "Spectrum %d number of frames does not match up with number of observations, exitting" % sid
            exit(-1)

        for i, f in zip(insert[insert_curr:(insert_curr+curr_frames)],fragment_masses[insert_curr:(insert_curr+curr_frames)]):
            ins_sequence.append(i)
            fragments.append(f)
            if not i:
                used_peaks[f] = 1

        insert_curr += curr_frames

        # log all PSM info
        currPsm = psm(kind = curr_kind,
                      scan = sid,
                      ch = curr_charge,
                      peptide = curr_pep_seq,
                      score = curr_score,
                      numBY = curr_numby,
                      fragments = fragments,
                      usedFragments = used_peaks,
                      insSequence = ins_sequence)

        if(curr_kind=='t'):
            if curr_charge not in target_psms:
                target_psms[curr_charge] = []
            target_psms[curr_charge].append(currPsm)
        else:
            if curr_charge not in decoy_psms:
                decoy_psms[curr_charge] = []
            decoy_psms[curr_charge].append(currPsm)
                           
        count += 1

    if count != numPeps: #didn't read all spectra scores, return empty scores
        td = [(None, float("-inf")), (None, float("-inf"))]

    sid_td_psms = []
    for charge in target_psms:
        # find maximum scoring target and decoy
        # todo: make the number of top matches returned a user input parameter
        target_psms[charge].sort(key = lambda x: x.score, reverse = True)
        decoy_psms[charge].sort(key = lambda x: x.score, reverse = True)
        t = target_psms[charge][0]
        d = decoy_psms[charge][0]
        sid_td_psms.append((target_psms[charge][0],decoy_psms[charge][0]))
        for i in range(min(topMatch, len(target_psms[charge]))):
            write_psm_ins_dels(target_psms[charge][i], currSpec, 
                               dripMeans, identFid)

        for i in range(min(topMatch, len(decoy_psms[charge]))):
            write_psm_ins_dels(decoy_psms[charge][i], currSpec, 
                               dripMeans, identFid)
    return sid_td_psms

def make_ident_persid(output, logDir, 
                      topMatch, highResMs2, 
                      spec_dict, meanFile,
                      pepDBperSidList, output_dir,
                      writeDelSequences = False,
                      writeAllDelSequences = False):
    """Driver function
    """

    # means are constant for low-res
    if not highResMs2:
        dripMeans = load_drip_means(meanFile)

    try:
        identFile = open(output, "w")
    except IOError:
        print "Could not open file %s for writing, exitting" % output

    sid_numPeps_lookup = open(pepDBperSidList, "r")
    reader = csv.DictReader(sid_numPeps_lookup, delimiter = '\t')
    identFile.write('Kind\tSid\tFrames\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\n')
    
    if(writeDelSequences):
        try:
            delSequencesFid=open(output[0:-4]+'_deletionSequences'+output[-4:], "w")
            delSequencesFid.write('Kind\tSid\tPeptide\tNum_deletes\tCharge\n')
        except IOError:
            print('Could not open %s for writing' % (output[0:-4]+'_deletionSequences'+output[-4:]))
            writeDelSequences=False
            delSequencesFid = None
    else:
        delSequencesFid = None
        
    target_peptide = None
    decoy_peptide = None
    for row in reader:
        sid = int(row['sid'])
        currSpec = spec_dict[sid]
        numPeps = int(row['numPeps'])
        testOutputFile = '%s/vitVals-sid%d.txt' % (logDir, sid)
        pepDBlist = '%s/sid%d-pepDB.txt' % (output_dir, sid)

        if highResMs2:
            # load this spectrum's collection of means
            mean_file = meanFile[:-4] + '-sid' + str(sid) + '.txt'
            dripMeans = load_drip_means(mean_file)

        td = parse_segments_persid(testOutputFile, dripMeans, 
                                   topMatch, sid, 
                                   currSpec, numPeps, 
                                   pepDBlist, identFile,
                                   writeDelSequences, writeAllDelSequences, 
                                   delSequencesFid)
    identFile.close()
