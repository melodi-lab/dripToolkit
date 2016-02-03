#!/usr/bin/env python
#
# Written by John Halloran <halloj3@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

from __future__ import with_statement

__authors__ = ['John Halloran <halloj3@uw.edu>' ]

import os
import re
import sys
import argparse
import pyFiles.digest_fasta as df
import shlex
import math
import string
import random
import pyFiles.psm as ppsm
import copy

from shutil import rmtree
from dripExtract import runDripExtract
from pyFiles.spectrum import MS2Spectrum
from pyFiles.peptide import Peptide, amino_acids_to_indices
from pyFiles.normalize import pipeline
from pyFiles.dripEncoding import (create_drip_structure,
                                  dripExtractParams,
                                  dripGaussianCollectionNames,
                                  create_drip_master,
                                  triangulate_drip,
                                  write_covar_file,
                                  drip_peptide_sentence,
                                  drip_spectrum_sentence,
                                  make_master_parameters,
                                  make_master_parameters_lowres,
                                  interleave_b_y_ions,
                                  interleave_b_y_ions_lowres,
                                  filter_theoretical_peaks,
                                  filter_theoretical_peaks_lowres,
                                  load_drip_means)

from pyFiles.constants import allPeps
from pyFiles.ioPsmFunctions import load_psms, load_pin_file, load_pin_return_dict
from pyFiles.shard_spectra import load_spectra_ret_dict

from pyFiles.process_vitvals import parse_dripExtract
from subprocess import call, check_output

# open GMTK file streams
stdo = open(os.devnull, "w")
stde = sys.stdout

def peptide_sentence_flatascii(pep_dt, peptide, bNy, pep_num, spectra_id, max_mass,
                               filename, istarget, max_candidate_theo_peaks):
    """Add a sentence/segment to the peptide ascii file.

    Arguments:
        peptide_pfile: The file stream pointer, already opened for writing.
        peptide: The peptide in the PSM, instance of protein.peptide.Peptide.
        tbl: The amino acid mass table to use ('average' | 'mono').
        pep_num: (i-1)th peptide in dt to be written
        spectra_id: sid, for dt naming purposes

    Effects:
        Adds a sentence to 'pep_dt'.

    """
    if(istarget):
        pep_dt.write('%d\nsid%dtargetseq%s\n1\n0 %d ' % (pep_num, spectra_id, peptide, len(bNy)))
    else:
        pep_dt.write('%d\nsid%ddecoyseq%s\n1\n0 %d ' % (pep_num, spectra_id, peptide, len(bNy)))

    for i in range(0, len(bNy)-1): 
        pep_dt.write('%d ' % i)
    pep_dt.write('default\n')

    for ion in bNy:
        pep_dt.write('\t-1 { %d }\n' % ion)
    pep_dt.write('\n')

    #write number of peaks and number of spurious peaks
    segment = 0
    frame = 0
    try:
        fid = open(filename, 'w')
    except IOError:
        print "Could not open file %s for writing, exitting" % filename
        exit(-1)

    fid.write("%d %d %d\n" % (segment, frame, max_candidate_theo_peaks))
    fid.close()

def spectrum_sentence_flatascii(filename, mz, intensity):
    """
    """
    try:
        fid = open(filename, 'w')
    except IOError:
        print "Could not open file %s for writing, exitting" % filename
        exit(-1)

    segment = 0
    for frame, (m, i) in enumerate(zip(mz, intensity)):
        fid.write("%d %d %f %f\n" % (segment, frame, m, i))
    fid.close()

def load_spectra_minMaxMz(spectra):
    """
    """
    # currently ignore ident file input for spectra filtering
    spectra, minMz, maxMz, validcharges, _ = load_spectra_ret_dict(spectra, 'all')
    return spectra, minMz, maxMz, validcharges

def load_spectra(spectra):
    """
    """
    # currently ignore ident file input for spectra filtering
    spectra, minMz, maxMz, validcharges, _ = load_spectra_ret_dict(spectra, 'all')
    return spectra

def plot_psms(psmFile, spectrumFile, plotList = 'currPsms.html',
              highResMs2 = False,
              dripLearnedMeans = 'dripLearned.means',
              dripLearnedCovars = 'dripLearned.covars',
              mods = '', ntermMods = '', ctermMods = ''):
    """
    """
    # parse modifications
    mods = df.parse_mods(mods, True)
    ntermMods = df.parse_mods(ntermMods, False)
    ctermMods = df.parse_mods(ctermMods, False)

    # initialize arguments for dripExtract
    args = dripExtractParams(psmFile, spectrumFile, 'all', 
                             mods, ntermMods, ctermMods, 
                             highResMs2, 
                             dripLearnedMeans, dripLearnedCovars)

    stde = open('gmtk_err', "w")
    # stdo = sys.stdout
    stdo = stde

    args.normalize = 'top300TightSequest'

    t, d, spectra0 = runDripExtract(args, stdo, stde)
    
    spectra, minMz, maxMz, validCharges = load_spectra_minMaxMz(spectrumFile)

    # get original intensity values to plot
    for sid in spectra:
        spectra[sid].mz = list(spectra0[sid].mz)
        mz_vals = set(spectra0[sid].mz)
        z = max(spectra0[sid].intensity)
        spectra[sid].intensity = [i/z for mz, i in zip(spectra[sid].mz, spectra[sid].intensity)
                                  if mz in mz_vals]

    for sid in spectra:
        s = spectra[sid]

    if not highResMs2:
        dripMeans = load_drip_means(dripLearnedMeans)
    else:
        dripMeans = {}
        for sid, c in t:
            p = t[sid,c]
            bNy = interleave_b_y_ions(Peptide(p.peptide), c, mods,
                                      ntermMods, ctermMods)
            filter_theoretical_peaks(bNy, minMz, maxMz, 0.1)
            for i, ion in enumerate(bNy):
                dripMeans[i] = ion
        for sid, c in d:
            p = d[sid,c]
            bNy = interleave_b_y_ions(Peptide(p.peptide), c, mods,
                                      ntermMods, ctermMods)
            filter_theoretical_peaks(bNy, minMz, maxMz, 0.1)
            for i, ion in enumerate(bNy):
                dripMeans[i] = ion

    ion_to_index_map = {} # reverse mapping, from ions to indices
    for ind in dripMeans:
        ion_to_index_map[dripMeans[ind]] = ind

    all_psms = []
    for sid, c in t:
        s = spectra[sid]
        for p in t[sid,c]:
            p.add_obs_spectrum(s)
            p.calculate_drip_features(dripMeans)
            p.calc_by_sets(c, mods,
                           ntermMods, ctermMods, highResMs2, 
                           ion_to_index_map)
        all_psms.append(p)
    for sid, c in d:
        s = spectra[sid]
        for p in d[sid,c]:
            p.add_obs_spectrum(s)
            p.calculate_drip_features(dripMeans)
            p.calc_by_sets(c, mods,
                           ntermMods, ctermMods, highResMs2, 
                           ion_to_index_map)
        all_psms.append(p)

    fid = open(plotList, "w")

    all_psms.sort(key = lambda r: r.score, reverse = True)
    # for p in all_psms.sort(key = lambda r: r.score):
    for p in all_psms:
        if p.kind == 't':
            kind = 'target'
        elif p.kind == 'd':
            kind = 'decoy'
        else:
            continue

        plotName = kind + 'Scan' + str(p.scan) + \
            'Charge' + str(p.charge) + \
            p.peptide + '.png'

        p.plot_drip_viterbi(plotName)
        fid.write("<a href=\"%s\">%s Scan %d Charge %d %s</a><br>\n" %
                  (plotName, kind, p.scan, p.charge, p.peptide))

    fid.close()

def psm(p, s0, c = 2, highResMs2 = False,
        dripLearnedMeans = 'dripLearned.means',
        dripLearnedCovars = 'dripLearned.covars',
        mods = '', ntermMods = '', ctermMods = ''):
    """ Inputs:
               p = peptide string
               s = observed spectrum, instance of class MS2Spectrum
               c = psm charge
               mods = static mods
               ntermMods = static nterm-mods
               ctermMods = static cterm-mods
    """

    s = copy.deepcopy(s0)

    args = dripGaussianCollectionNames()
    sid = s.spectrum_id

    # parse modifications
    mods = df.parse_mods(mods, True)
    ntermMods = df.parse_mods(ntermMods, False)
    ctermMods = df.parse_mods(ctermMods, False)

    # set observed spectrum preprocessing
    normalize = 'top300TightSequest'
    preprocess = pipeline(normalize)
    preprocess(s)

    # get original intensity values to plot
    s0.mz = list(s.mz)
    mz_vals = set(s.mz)
    z = max(s0.intensity)    
    s0.intensity = [i/z for mz, i in zip(s0.mz, s0.intensity)
                                  if mz in mz_vals]
    num_psms = 1

    max_obs_mass = 2001

    # dirBase = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(10))
    dirBase = 'dtk'

    # output_dir = os.path.abspath('dripEncode_' + dirBase)
    output_dir = os.path.abspath('encode')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    obs_dir = 'obs' # sub directory of output_dir
    pfile_dir = os.path.join(output_dir, obs_dir)
    if not os.path.exists(pfile_dir):
        os.mkdir(pfile_dir)

    # log_dir = os.path.abspath('dripLog_' + dirBase)
    log_dir = os.path.abspath('log')
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    if not highResMs2:
        dripMeans = load_drip_means(dripLearnedMeans)
        bNy = interleave_b_y_ions_lowres(Peptide(p), c, mods,
                                         ntermMods, ctermMods)
        l = len(bNy)
        filter_theoretical_peaks_lowres(bNy, 
                                        dripMeans, s.mz[0], s.mz[-1])
    else:
        # calculate b- and y-ions, filter peaks outside of spectrum range
        bNy = interleave_b_y_ions(Peptide(p), c, mods,
                                  ntermMods, ctermMods)
        l = len(bNy)
        filter_theoretical_peaks(bNy, s.mz[0], s.mz[-1], 0.1)
        # now construct means based on this
        dripMeans = {}
        for i, ion in enumerate(bNy):
            dripMeans[i] = ion

    ion_to_index_map = {} # reverse mapping, from ions to indices
    for ind in dripMeans:
        ion_to_index_map[dripMeans[ind]] = ind

    # make collection per spectrum
    make_master_parameters_lowres(args, dripMeans)
    peptide_obs_file = os.path.join(pfile_dir,'pep-lengths')
    spectrum_obs_file = os.path.join(pfile_dir,'spectrum')

    pep_dt = open(os.path.join(output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    # write peptide database to parse and identify GMTK segments later
    pepdb_list = open(os.path.join(output_dir, 'pepDB.txt'), "w")
    pepdb_list.write("Kind\tSid\tPeptide\tNumBY\tCharge\n")

    pep_num = 0
    # create iterable dt and peptide pfile
    peptide_sentence_flatascii(pep_dt, p, bNy, 
                               pep_num, sid, max_obs_mass,
                               peptide_obs_file, True, len(bNy))
    # create spectrum pfile
    spectrum_sentence_flatascii(spectrum_obs_file, s.mz, s.intensity)
    pepdb_list.write("t\t%d\t%s\t%d\t%d\n" % (sid, 
                                              p, l, c))
        
    # close streams for this spectrum
    pep_dt.close()
    pepdb_list.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)
         # stdout = sys.stderr, stderr = sys.stderr)

    # create structure and master files then triangulate
    try:
        create_drip_structure(highResMs2, args.structure_file, 
                              max_obs_mass)
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        create_drip_master(highResMs2, args.master_file, 
                           max_obs_mass,
                           "DRIP_MZ",
                           "drip_collection/covar.txt",
                           "DRIP_GAUSSIAN_COMPONENTS",
                           "DRIP_GAUSSIAN_MIXTURES",
                           "DRIP_MZ_GAUSSIANS")
    except:
        print "Could not create DRIP master file %s, exitting" % args.master_file
        exit(-1)

    try:
        triangulate_drip(args.structure_file, args.master_file)
    except:
        print "Could not create triangulate structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        write_covar_file(highResMs2, args.covar_file, 
                         dripLearnedCovars)
    except:
        print "Could not create covariance file %s, exitting" % args.covar_file
        exit(-1)

    # run GMTK
    dtFile = os.path.join(output_dir, 'iterable.dts')
    cppCommand = '\'-DITERABLE_DT=' + dtFile \
        + ' -DDRIP_MZ=' + args.mean_file \
        + ' -DDRIP_GAUSSIAN_COMPONENTS=' + args.gauss_file \
        + ' -DDRIP_GAUSSIAN_MIXTURES=' + args.mixture_file \
        + ' -DDRIP_MZ_GAUSSIANS=' + args.collection_file \
        + '\''

    # call gmtkViterbi
    vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
        + " -fdiffact2 rl" \
        + " -inputMasterFile " + args.master_file + " -inputTrainableParameters trained.params -failOnZeroClique F"
    # gmtkViterbi command line
    vitValsFile = os.path.join(log_dir, 'vitVals.txt')
    vitStr = vitStr0 + ' -vitValsFile ' +  vitValsFile \
        + ' -of1 ' + spectrum_obs_file \
        + ' -fmt1 flatascii ' \
        + ' -of2 ' + peptide_obs_file \
        + ' -fmt2 flatascii ' \
        + ' -cppCommand ' + cppCommand
    # call(shlex.split(vitStr), stdout = sys.stdout, stderr = sys.stdout)
    call(shlex.split(vitStr), stdout = stdo, stderr = stde)

    # parse output
    t,d = ppsm.parse_dripExtract(vitValsFile, os.path.join(output_dir, 'pepDB.txt'))

    t = t[sid,c][0]
    # calculate insertions and deletions
    t.add_obs_spectrum(s0)
    t.calculate_drip_features(dripMeans)
    t.calc_by_sets(c, mods,
                   ntermMods, ctermMods, highResMs2, 
                   ion_to_index_map)
    return t
        
if __name__ == '__main__':
    # process input arguments
    args = process_args()
    targets, decoys, spec_dict = runDripExtract(args)
