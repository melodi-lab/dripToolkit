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
import cPickle as pickle
import pyFiles.digest_fasta as df
import multiprocessing
import shlex
import math
from pyFiles.psm import dripPSM, parse_dripExtract

from pyFiles.gmtkUtils import parse_params
from random import shuffle
from shutil import rmtree
from pyFiles.spectrum import MS2Spectrum
from pyFiles.peptide import Peptide, amino_acids_to_indices
from pyFiles.normalize import pipeline
from pyFiles.pfile.wrapper import PFile
# from pyFiles.args import make_dir_callback
from pyFiles.dripEncoding import (create_drip_structure,
                                  create_drip_master,
                                  triangulate_drip,
                                  write_covar_file,
                                  create_pfile, 
                                  drip_peptide_sentence,
                                  drip_spectrum_sentence,
                                  training_spectrum_sentence,
                                  make_master_parameters,
                                  make_master_parameters_lowres,
                                  interleave_b_y_ions,
                                  interleave_b_y_ions_lowres,
                                  filter_theoretical_peaks,
                                  filter_theoretical_peaks_lowres,
                                  create_empty_master,
                                  create_params_no_train,
                                  load_drip_means,
                                  dripMean)

from pyFiles.constants import allPeps
from pyFiles.ioPsmFunctions import load_psm_library
from pyFiles.shard_spectra import (calc_minMaxMz,
                                   pickle_candidate_spectra, 
                                   candidate_spectra_generator,
                                   candidate_spectra_memeffic_generator,
                                   load_spectra_ret_dict,
                                   load_target_decoy_db_no_reshuffle)
from subprocess import call, check_output, STDOUT

debug=1

bw = 1.0005079
bo = 0.68
########## would have to rework all code to use below, but this would make working with the DRIP means
########## considerably simpler.  However, mean class is overkill as the means are scalars
def load_dripMeansClass(master_file):
    """ Load DRIP learned means from master file.
        Return dictionary of gmtkMean objects whose keys are the mean names
    """
    mean_file = open(master_file)
    log = mean_file.read()
    mean_file.close()
    mean_pattern = '[^%]mean(?P<meanInd>\d+) 1 (?P<meanVal>\S+)'

    means = {}
    for match in re.finditer(mean_pattern, log):
        name = 'mean' + str(match.group('meanInd'))
        dimension = 1
        val = float(match.group('meanVal'))
        currMean = dripMean(name, 
                            val, 
                            int(match.group('meanInd')), 
                            dimension)
        means[name] = currMean
        means[int(match.group('meanInd'))] = float(match.group('meanVal'))

    return means

def process_args(args):
    """ Check whether relevant directories already exist (remove them if so), create
        necessary directories.  Adjust parameters.

        pre:
        - args has been created by calling parse_args()
        
        post:
        - Relevant directories will be first removed then created, boolean parameters
          will be adjusted to reflect their boolean command line strings
    """
    # check input arguments, create necessary directories, 
    # set necessary global variables based on input arguments
    # constant string
    if args.max_obs_mass < 0:
        print "Supplied max-mass %d, must be greater than 0, exitting" % args.max_obs_mass
        exit(-1)
    # check arguments
    # create obervation directories
    # vit vals output directory
    try:
        base = os.path.abspath(args.logDir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.logDir)
        exit(-1)

    # top level encode directory
    try:
        base = os.path.abspath(args.output_dir)
        args.output_dir = base
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.output_dir)
        exit(-1)
    # observation directory
    try:
        base = os.path.join(os.path.abspath(args.output_dir), args.obs_dir)
        os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.join(os.path.abspath(args.output_dir), args.obs_dir)
        exit(-1)
    # Gaussian collection directory
    try:
        base = os.path.abspath(args.collection_dir)
        if not os.path.exists(base):
            os.mkdir(base)
        else:
            rmtree(base)
            os.mkdir(base)
    except:
        print "Could not create observation directory %s, exitting" % os.path.abspath(args.collection_dir)
        exit(-1)

    args.covar_file = os.path.join(os.path.abspath(args.collection_dir), args.covar_file)
    args.mean_file = os.path.join(os.path.abspath(args.collection_dir), args.mean_file)
    args.gauss_file = os.path.join(os.path.abspath(args.collection_dir), args.gauss_file)
    args.mixture_file = os.path.join(os.path.abspath(args.collection_dir), args.mixture_file)
    args.collection_file = os.path.join(os.path.abspath(args.collection_dir), args.collection_file)


    args.filt_theo_peaks = True
    # set true or false strings to booleans
    args.high_res_ms2 = df.check_arg_trueFalse(args.high_res_ms2)

    args.num_threads = 1
    # # make sure number of input threads does not exceed number of supported threads
    # if args.num_threads > multiprocessing.cpu_count():
    #     args.num_threads = max(multiprocessing.cpu_count()-1,1)

def remove_theoretical_holes(dripMeans, usedTheoPeaks):
    """ A theoretical hole is a DRIP mean not accessed by any training PSM
        (i.e., such a theoretical peak is not present in any theoretical
        spectrum amongst the training PSMs).  Such holes will, by definition,
        receive zero probability during training and could cause unexpected
        behavior due to numerical instability (the EM update for the means
        involves the reciprocal of the posterior probability) so it is 
        best to remove these and add them back in after training.

        Note: "theoretical hole" is not to be confused with the notion of a deletion.
              A deletion pertains to a peak in a theoretical spectrum where,
              given an observed spectrum, a DRIP alignment does not contain
              the scoring of any observed peaks with this theoretical peak. 
              In contrast, a hole pertains to a theoretical peak that is not
              contained in any of the union of theoretical spectra amongst
              a collection of PSMs.

        inputs: 
                dripMeans - dictionary whose keys are the index of the theoretical peak and whose values are the DRIP mean value
                usedTheoPeaks - set of theoretical peaks that are not holes

        outputs: None

        pre:
             -dripMeans is constructed before calling this function (either 
             read from a file of prior means or initialized to default 
             values) and should be a python dictionary object
             -usedTheoPeaks is constructed and should be a python set object
             
        post:
             -theoretical holes are deleted from dripMeans.  Note that
             dripMeans is passed by reference and adjusted herein
    """
    
    meanKeys = list(dripMeans.iterkeys())
    # print sorted(list(set(meanKeys).difference(set(usedTheoPeaks))))
    for m in meanKeys:
        if m not in usedTheoPeaks:
            hole = dripMeans.pop(m, None)

def remap_mean_indices(args, spectra,
                       dripMeans, usedTheoPeaks, theo_dict):
    """ Remap theoretical peaks to indices under collection of
        non-theoretical-holes

        todo: add high-res ms2 to this, even though learning
              high-res ms2 parameters hasn't lead to any 
              improvements

        pre:
            -drip means have been loaded into dripMeans and
            theoretical holes have been removed
            -non-theoretical-holes have been properly denoted
            in the set usedTheoPeaks
            -dictionary of theoretical spectra have been constructed
            but have not been remapped under the collection of 
            active theoretical peaks
        post:
             -each theoretical peak in theo_dict is an index
             into array of drip means sorted by key in dripMeans
    """
    
    new_theo_index_mapping = {}
    for new_ind, old_ind in enumerate(sorted(dripMeans.iterkeys())):
        new_theo_index_mapping[old_ind] = new_ind

    # spec_considered = set([])
    for sid,charge,pep in theo_dict:
        # s = spectra[sid]
        # if sid not in spec_considered:
        #     preprocess(s)
        #     spec_considered.add(sid)
        # calculate b- and y-ions under new mapping
        bNy = theo_dict[sid,charge,pep]
        # bNy = [new_theo_index_mapping[i] for i in bNy]
        theo_dict[sid,charge,pep] = [new_theo_index_mapping[i] for i in bNy]

def inject_mean_evidence(args, spectra,
                         dripMeans,
                         fa_psms, theo_dict, lowProbMeans,
                         observedHoles):
    """ Remap theoretical peaks to indices under collection of
        non-theoretical-holes

        todo: add high-res ms2 to this, even though learning
              high-res ms2 parameters hasn't lead to any 
              improvements

        pre:
            -drip means have been loaded into dripMeans and
            theoretical holes have been removed
            -non-theoretical-holes have been properly denoted
            in the set usedTheoPeaks
            -dictionary of theoretical spectra have been constructed
            but have not been remapped under the collection of 
            active theoretical peaks
        post:
             -each theoretical peak in theo_dict is an index
             into array of drip means sorted by key in dripMeans
    """

    new_theo_index_mapping = {}
    new_to_old_mapping = {}
    for new_ind, old_ind in enumerate(sorted(dripMeans.iterkeys())):
        new_theo_index_mapping[old_ind] = new_ind
        new_to_old_mapping[new_ind] = old_ind

    remapped_lowProbMeans = set([new_theo_index_mapping[m] for m in observedHoles])
    remapped_lowProbMeans |= set([new_theo_index_mapping[m] for m in lowProbMeans])
    # remapped_lowProbMeans = set([new_theo_index_mapping[m] for m in lowProbMeans])

    for sid,charge,pep in theo_dict:
        dripPsm = []
        ind = -1
        for psm in fa_psms[sid,charge]:
            ind += 1
            if psm.peptide == pep:
                dripPsm = psm
                break
        assert dripPsm

        bNy = set(theo_dict[sid,charge,pep])
        bNy_original = [new_to_old_mapping[theoPeak] for theoPeak in theo_dict[sid,charge,pep]]
        overlap = bNy & remapped_lowProbMeans
        if overlap:
            s = spectra[sid]
            temp_spec = [(mz, intensity) for mz, intensity in zip(s.real_mz, s.real_intensity)]
            for m in overlap:
                k = new_to_old_mapping[m]
                synthetic_mz = (bw*float(2*k+1)-2*bo)/2
                temp_spec.append((synthetic_mz,1.0))
            temp_spec.sort(key = lambda r: r[0])
            s.mz = [mz for mz, _ in temp_spec]
            s.intensity = [intensity for _, intensity in temp_spec]
            dripPsm.insertion_sequence = [0 if calc_bin(mz,bo,bw) in bNy_original else 1 for mz in s.mz]
            dripPsm.add_obs_spectrum(s)

            spectra[sid] = s
            fa_psms[sid,charge][ind] = dripPsm

    # # create dictionary of charges with sid keys
    # sidChargeDict = {}
    # pepDict = {}
    # for sid, charge, pep in theo_dict:
    #     if sid in sidChargeDict:
    #         sidChargeDict[sid].append(charge)
    #     else:
    #         sidChargeDict[sid] = [charge]

    #     if (sid,charge) not in pepDict:
    #         pepDict[sid,charge] = [pep]
    #     else:
    #         pepDict[sid,charge].append(pep)

    # # add synthetic peaks to observed spectrum holes,
    # # estimate forced alignment
    # for sid in sidChargeDict:
    #     s = spectra[sid]
    #     # temp_spec = [(mz, intensity) for mz, intensity in zip(s.mz, s.intensity)]
    #     temp_spec = [(mz, intensity) for mz, intensity in zip(s.real_mz, s.real_intensity)]
    #     synthetic_peak_locations = set([])
    #     for charge in sidChargeDict[sid]:
    #         for pep in pepDict[sid,charge]:
    #             bNy = theo_dict[sid,charge,pep]
    #             synthetic_peak_locations |= set(bNy).intersection(remapped_lowProbMeans)
    #     for mz in synthetic_peak_locations:
    #         # synthetic_mz = dripMeans[new_to_old_mapping[mz]]
    #         synthetic_mz = (bw*float(2*mz+1)-2*bo)/2

    #         temp_spec.append((synthetic_mz, 1.0))

    #     temp_spec.sort(key = lambda r: r[0])

    #     s.mz = [mz for mz, intensity in temp_spec]
    #     s.intensity = [intensity for mz, intensity in temp_spec]

    #     curr_ins_seq = [False if calc_bin(mz,bo,bw) in bNy else True for mz in s.mz]
    #     # update forced alignment estimate
    #     for charge in sidChargeDict[sid]:
    #         for i in range(len(fa_psms[sid,charge])):
    #             pep = fa_psms[sid,charge][i].peptide
    #             bNy = theo_dict[sid,charge,pep]
    #             curr_ins_seq = [False if calc_bin(mz,bo,bw) in bNy else True for mz in s.mz]
    #             fa_psms[sid,charge][i].insertion_sequence = curr_ins_seq

def remap_learned_means(dripMeans_og, learnedMeans):
    """ Remap theoretical peaks to indices under collection of
        non-theoretical-holes

        todo: add high-res ms2 to this, even though learning
              high-res ms2 parameters hasn't lead to any 
              improvements

        pre:
            -drip means have been loaded into dripMeans and
            theoretical holes have been removed
            -non-theoretical-holes have been properly denoted
            in the set usedTheoPeaks
            -dictionary of theoretical spectra have been constructed
            but have not been remapped under the collection of 
            active theoretical peaks
        post:
             -each theoretical peak in theo_dict is an index
             into array of drip means sorted by key in dripMeans
    """
    pat = re.compile('\d+')

    for m in learnedMeans:
        if str(m[0])=='intensity_mean':
            continue
        oldInd = int(pat.findall(str(m[0]))[0])
        assert m[1]==1 and len(m[2])==1, "Learned mean %s has dimension %d, DRIP means should only have dimension 1.  Exitting" % (m[0], m[1])
        dripMeans_og[oldInd] = m[2][0]
                

def write_learned_means_covars(learnedMeans, learnedCovars, meanName, covarName):
    # meanName = baseName + '.means'
    # covarName = baseName + '.covars'

    # check that means are increasing in m/z value
    prevMean = -1.0
    resortedMeans = []
    resort = 0
    for m in sorted(learnedMeans.iterkeys()):
        newMean = learnedMeans[m]
        if newMean < prevMean:
            resort = 1 # if we have to resort, the indices no longer have meaning
        resortedMeans.append(learnedMeans[m])
        prevMean = newMean

    if resort:
        print "Sorting learned means in increasing order"
        resortedMeans.sort()

    # write mean file
    meanFid = open(meanName, 'w')
    meanFid.write("%d\n" % len(resortedMeans))
    for i,m in enumerate(resortedMeans):
        meanFid.write("%d mean%d 1 %.10e\n" % (i,i,m))
    meanFid.close()
    # # write mean file
    # meanFid = open(meanName, 'w')
    # meanFid.write("%d\n" % len(learnedMeans))
    # for m in sorted(learnedMeans.iterkeys()):
    #     meanFid.write("%d mean%d 1 %f\n" % (m,m,learnedMeans[m]))
    # meanFid.close()

    # write covar file
    covarFid = open(covarName, 'w')
    covarFid.write("%d\n" % len(learnedCovars))
    for i, covar in enumerate(learnedCovars):
        covarFid.write("%d %s %d %.10e\n" % (i,covar[0],
                                             covar[1],
                                             covar[2][0]))
    covarFid.close()

def calc_bin(mz, bin_offset = 0.68, bin_width = 1.0005079, ma = False):
    if ma:
        return int(math.floor((mz+bin_offset)/bin_width))
    else:
        return int(math.floor(mz/bin_width-bin_offset+1))

def calcTheoDict_estimateFa(args, spectra):
    """Dictionary of theoretical spectra for each training PSM

       don't run this for high-res; theoretical spectrum
       preprocessing is explicitly handled in high-res data
       generation
    """

    assert not args.high_res_ms2, "calculate_theoretical_dictionary shouldn't be called in high-res ms2 mode"

    # parse modifications
    mods = df.parse_mods(args.mods_spec, True)
    print "mods:"
    print mods
    ntermMods = df.parse_mods(args.nterm_peptide_mods_spec, False)
    print "n-term mods:"
    print ntermMods
    ctermMods = df.parse_mods(args.cterm_peptide_mods_spec, False)
    print "c-term mods:"
    print ctermMods

    target,num_psms = load_psm_library(args.psm_library)
    sid_charges = list(target.iterkeys())

    theo_dict = {} # keys: sid, charge, peptide string
    usedTheoPeaks = set([])
    validcharges = args.charges # this should have been updated by loading
                                # the spectra into memory first

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    max_mass = 2001
    not_observed_holes = [0] * max_mass

    observedPeaks = set([])

    minMz = 0
    maxMz = max_mass-1

    spec_considered = set([])
    for sid, charge in sid_charges:
        s = spectra[sid]
        if sid not in spec_considered:
            preprocess(s)
            spec_considered.add(sid)
            s.real_mz = list(s.mz)
            s.real_intensity = list(s.intensity)

            spectra[sid] = s

        for p in target[sid,charge]:
            pep = p.peptide
            # calculate theoretical spectra
            bNy = interleave_b_y_ions_lowres(Peptide(pep), charge, mods,
                                             ntermMods, ctermMods)
            # calculate complement of observed holes

            observedPeaks |= set([calc_bin(mz, bo, bw) for mz in s.mz])
            if args.filt_theo_peaks:
                filter_theoretical_peaks(bNy, minMz, maxMz)
            theo_dict[sid,charge,pep] = bNy
            # update set of active theoretical peaks
            usedTheoPeaks |= set(bNy)

    observedHoles = usedTheoPeaks.difference(observedPeaks)

    # create dictionary of charges with sid keys
    sidChargeDict = {}
    for sid, charge in sid_charges:
        if sid in sidChargeDict:
            sidChargeDict[sid].append(charge)
        else:
            sidChargeDict[sid] = [charge]
    
    # add synthetic peaks to observed spectrum holes,
    # estimate forced alignment
    for sid in sidChargeDict:            
        s = spectra[sid]
        temp_spec = [(mz, intensity) for mz, intensity in zip(s.mz, s.intensity)]
        synthetic_peak_locations = set([])
        for charge in sidChargeDict[sid]:
            for p in target[sid,charge]:
                pep = p.peptide
                bNy = theo_dict[sid,charge,pep]
                synthetic_peak_locations |= set(bNy).intersection(observedHoles)
        for mz in synthetic_peak_locations:
            synthetic_mz = (bw*float(2*mz+1)-2*bo)/2
            # synthetic_mz = float(mz)/2
            # synthetic_mz = dripMeans[new_to_old_mapping[mz]]
            temp_spec.append((synthetic_mz, 1.0))

        temp_spec.sort(key = lambda r: r[0])

        s.mz = [mz for mz, _ in temp_spec]
        s.intensity = [intensity for _, intensity in temp_spec]
        spectra[sid] = s

    fa_psms = {}
    # iteratre through new m/z values, estimate forced alignment
    for sid in sidChargeDict:
        s = spectra[sid]
        for charge in sidChargeDict[sid]:
            for p in target[sid,charge]:
                pep = p.peptide
                bNy = set(theo_dict[sid,charge,pep])
                curr_ins_seq = [False if calc_bin(mz,bo,bw) in bNy else True for mz in s.mz]
                
                currPsm = dripPSM(pep, 0.0, sid, 't', 
                                  charge, len(bNy), [100] * len(curr_ins_seq),
                                  curr_ins_seq)
                if (sid, charge) in fa_psms:
                    fa_psms[sid,charge].append(currPsm)
                else:
                    fa_psms[sid,charge] = [currPsm]

    return theo_dict, usedTheoPeaks, target, num_psms, observedHoles, fa_psms

def make_fa_data_highres(args, spectra, target, num_psms, stdo, stde):
    """Generate test data .pfile. and create job scripts for cluster use (if num_jobs > 1).
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go.

       inputs:
       args - output of parsed input arguments (struct)

       outputs:
       sids - list of scan IDs for the generated data

       pre:
       - args has been created by parse_args(), directories have been created/checked for existence,
         relevant arguments have been processed (Booleans, mods, digesting enzyme, etc)
       - data has been created by candidate_spectra_generate() and contains the above mentioned fields

       post:
       - args.{mean_file, gauss_file, mixture_file, collection_file} will all be adjusted
       - args.max_mass will be updated to the size of the number of unique theoretical fragmentation locations (floating point if high-res ms2, integers if low-res ms2)
    """
    # parse modifications
    mods = df.parse_mods(args.mods_spec, True)
    print "mods:"
    print mods
    ntermMods = df.parse_mods(args.nterm_peptide_mods_spec, False)
    print "n-term mods:"
    print ntermMods
    ctermMods = df.parse_mods(args.cterm_peptide_mods_spec, False)
    print "c-term mods:"
    print ctermMods

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)
    sid_charges = list(target.iterkeys())

    # assume that we should randomize PSMs for multithreading purposes; only reason
    # why we are currently assuming this is that there is already a parameter for dripSearch
    # which signifies whether we should shuffle the data
    shuffle(sid_charges)

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    validcharges = args.charges

    ion_dict = {} # global dictionary for used fragment ions
    theo_spec_dict = {}
    numBY_dict_per_sid = {}
    # construct ion_dict
    for sid in spectra:
        s = spectra[sid]
        # preprocess(s)
        for charge in validcharges:
            if (s.spectrum_id, charge) not in target:
                continue
            # check if we're filtering theoretical peaks outside observed m/z values
            if args.filt_theo_peaks:
                if args.per_spectrum_mz_bound:
                    minMz = s.mz[0]
                    maxMz = s.mz[-1]
                else:
                    minMz = args.mz_lb
                    maxMz = args.mz_ub

            # calculate maximum decoy and target theoretical spectra cardinalities
            for p in target[s.spectrum_id, charge]:
                pep = p.peptide
                bNy = interleave_b_y_ions(Peptide(pep), charge, mods,
                                          ntermMods, ctermMods)
                numBY_dict_per_sid[sid, pep] = len(bNy)
                if args.filt_theo_peaks:
                    filter_theoretical_peaks(bNy, minMz, maxMz)
                # to be backwards compatible with make_training_data_lowres,
                # we must add a dict key charge.  This makes sense for training since
                # we a high-confidence PSM may be desired to train over different
                # charge states.  This doesn't make sense for testing (dripSearch/dripExtract)
                # since it's impossible that one peptide will be in multiple charge-candidate-peptide-sets
                theo_spec_dict[s.spectrum_id, charge, pep] = bNy

                for i in bNy:
                    ion_dict[i] = 1
            for d in decoy[s.spectrum_id, charge]:
                pep = d.peptide
                bNy = interleave_b_y_ions(Peptide(pep), charge, mods, 
                                          ntermMods, ctermMods)
                numBY_dict_per_sid[sid, pep] = len(bNy)
                if args.filt_theo_peaks:
                    filter_theoretical_peaks(bNy, minMz, maxMz)
                theo_spec_dict[s.spectrum_id, charge, pep] = bNy
                for i in bNy:
                    ion_dict[i] = 1

    ions = list(ion_dict.iterkeys())
    ions.sort()
    dripMeans = {}
    for i, ion in enumerate(ions):
        ion_dict[ion] = i
        dripMeans[i] = ion

    # make collection per spectrum
    make_master_parameters(args, ion_dict, ions)
    peptide_pfile = create_pfile(pfile_dir,
                                 'pep-lengths.pfile',
                                 0, 1)
            
    spectrum_pfile = create_pfile(pfile_dir,
                                  'spectrum.pfile',
                                  2,0)

    pep_dt = open(os.path.join(args.output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    # write peptide database to parse and identify GMTK segments later
    pepdb_list = open(os.path.join(args.output_dir, 'pepDB.txt'), "w")
    pepdb_list.write("Kind\tSid\tPeptide\tNumBY\tCharge\n")

    spec_dict = {}
    pep_num = 0
    for sid, charge in sid_charges:
        if sid not in spec_dict:
            s = spectra[sid]
            preprocess(s)
            spec_dict[sid] = s
        else:
            s = spec_dict[sid]

        for p in target[sid,charge]:
            pep = p.peptide
            bNy = theo_spec_dict[s.spectrum_id, charge, pep]
            bNy = [ion_dict[bOrY] for bOrY in bNy]
            drip_peptide_sentence(pep_dt, pep, bNy, 
                                  pep_num, s.spectrum_id, args.max_obs_mass,
                                  peptide_pfile, True, len(bNy)-1)
            drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
            pepdb_list.write("t\t%d\t%s\t%d\t%d\n" % (sid, 
                                                      pep, 
                                                      numBY_dict_per_sid[sid, pep],
                                                      charge))
            pep_num += 1

        if (sid,charge) in decoy:
            for d in decoy[sid,charge]:
                pep = d.peptide
                bNy = theo_spec_dict[s.spectrum_id, charge, pep]
                bNy = [ion_dict[bOrY] for bOrY in bNy]
                drip_peptide_sentence(pep_dt, pep, bNy, 
                                      pep_num, s.spectrum_id, args.max_obs_mass,
                                      peptide_pfile, False, len(bNy)-1)
                drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
                pepdb_list.write("d\t%d\t%s\t%d\t%d\n" % (sid, 
                                                          pep, 
                                                          numBY_dict_per_sid[sid, pep],
                                                          charge))
                pep_num += 1
                theo_spec_dict[s.spectrum_id, charge, pep] = bNy

    # close streams for this spectrum
    pep_dt.close()
    pepdb_list.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(args.output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)

    return theo_spec_dict, dripMeans

def make_fa_data_lowres(args, spectra, dripMeans, theo_dict, 
                        target, num_psms, stdo, stde):
    """Generate test data .pfile. and create job scripts for cluster use.
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go
    """
    # make master file
    make_master_parameters_lowres(args, dripMeans)

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)
    sid_charges = list(target.iterkeys())

    # shuffle(sid_charges) # shuffle if we're parallelizing GMTK jobs
    validcharges = args.charges

    # write peptide database to parse and identify GMTK segments later
    pepdb_list = open(os.path.join(args.output_dir, 'pepDB.txt'), "w")
    pepdb_list.write("Kind\tSid\tPeptide\tNumBY\tCharge\n")

    peptide_pfile = create_pfile(pfile_dir,
                                 'pep-lengths.pfile',
                                 0, 1)
            
    spectrum_pfile = create_pfile(pfile_dir,
                                  'spectrum.pfile',
                                  2,0)

    pep_dt = open(os.path.join(args.output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    pep_num = 0
    for sid, charge in sid_charges:
        s = spectra[sid]
        for p in target[sid,charge]:
            pep = p.peptide
            bNy = theo_dict[sid,charge,pep]
            pepdb_list.write("t\t%d\t%s\t%d\t%d\n" % (sid, pep, len(bNy), charge))
            drip_peptide_sentence(pep_dt, pep, bNy, 
                                  pep_num, s.spectrum_id, args.max_obs_mass,
                                  peptide_pfile, True, len(bNy)-1)
            drip_spectrum_sentence(spectrum_pfile, s.mz, s.intensity)
            pep_num += 1
    # close streams for this spectrum
    pep_dt.close()
    pepdb_list.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(args.output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)

def make_training_data_highres(args, spectra, 
                               theo_spec_dict, ion_dict,
                               target, num_psms,
                               fa_psms, stdo, stde):
    """Generate test data .pfile. and create job scripts for cluster use (if num_jobs > 1).
       Decrease number of calls to GMTK by only calling once per spectrum
       and running for all charge states in one go.

       inputs:
       args - output of parsed input arguments (struct)

       outputs:
       sids - list of scan IDs for the generated data

       pre:
       - args has been created by parse_args(), directories have been created/checked for existence,
         relevant arguments have been processed (Booleans, mods, digesting enzyme, etc)
       - data has been created by candidate_spectra_generate() and contains the above mentioned fields

       post:
       - args.{mean_file, gauss_file, mixture_file, collection_file} will all be adjusted
       - args.max_mass will be updated to the size of the number of unique theoretical fragmentation locations (floating point if high-res ms2, integers if low-res ms2)
    """
    # parse modifications
    mods = df.parse_mods(args.mods_spec, True)
    print "mods:"
    print mods
    ntermMods = df.parse_mods(args.nterm_peptide_mods_spec, False)
    print "n-term mods:"
    print ntermMods
    ctermMods = df.parse_mods(args.cterm_peptide_mods_spec, False)
    print "c-term mods:"
    print ctermMods

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)
    sid_charges = list(target.iterkeys())

    # assume that we should randomize PSMs for multithreading purposes; only reason
    # why we are currently assuming this is that there is already a parameter for dripSearch
    # which signifies whether we should shuffle the data
    shuffle(sid_charges)

    if(args.normalize != 'filter0'):
        preprocess = pipeline(args.normalize)

    validcharges = args.charges

    ions = list(ion_dict.iterkeys())
    ions.sort()
    for i, ion in enumerate(ions):
        ion_dict[ion] = i

    # make collection per spectrum
    make_master_parameters(args, ion_dict, ions)
    peptide_pfile = create_pfile(pfile_dir,
                                 'pep-lengths.pfile',
                                 0, 1)
            
    spectrum_pfile = create_pfile(pfile_dir,
                                  'spectrum.pfile',
                                  2,1)

    pep_dt = open(os.path.join(args.output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    pep_num = 0
    for sid, charge in sid_charges:
        s = spec_dict[sid]
        for p in target[sid,charge]:
            pep = p.peptide
            # look for corresponding dripPSM
            dripPsm = []
            for t_psm in training_psms[sid,charge]:
                if t_psm.peptide == pep:
                    dripPsm = t_psm
                    break
            assert dripPsm

            bNy = theo_spec_dict[s.spectrum_id, pep]
            drip_peptide_sentence(pep_dt, pep, bNy, 
                                  pep_num, s.spectrum_id, args.max_obs_mass,
                                  peptide_pfile, True, len(bNy)-1)
            training_spectrum_sentence(spectrum_pfile, s.mz, s.intensity, dripPsm)
            pep_num += 1

    # close streams for this spectrum
    pep_dt.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(args.output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)

    return pep_num

def make_training_data_lowres(args, spectra, dripMeans, theo_dict,
                              target, num_psms, fa_psms, stdo, stde):
    """Generate .pfiles for each training PSM.

       inputs:
              args - user input options
              spectra - dictionary of spectra loaded into memory
              theo_dict - dictionary of theoretical spectra per training PSM
              target - training psms, instances of class PSM
              num_psms - number of training PSMs
              fa_psms - dictionary forced alignment PSMs, each an instance of 
                        class dripPSM and containing sequences of insertions and
                        deletions per PSM (used for regularization during training)

       outputs: None

       pre:
           args - user inputs have been checked process_args
           spectra - spectra have been loaded but not preprocessed
           theo_dict - dictionary have theoretical peaks has been created
           target - training PSMs have been loaded
           fa_psms - a forced alignment has been run and the sequences of 
                     insertions and deletions has been computed 
                     (calculateForcedAlignment has been run)

       post:
            GMTK training iterable dt will be created 
            in args.output_dir and training .pfiles will be created
            in args.output_dir / args.obs_dir
    """
    # make master file
    make_master_parameters_lowres(args, dripMeans)

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)
    sid_charges = list(target.iterkeys())

    sid_charges.sort(key = lambda r: r[0])

    # # assume that we should randomize PSMs for multithreading purposes; only reason
    # # why we are currently assuming this is that there is already a parameter for dripSearch
    # # which signifies whether we should shuffle the data
    # shuffle(sid_charges)
    peptide_pfile = create_pfile(pfile_dir,
                                 'pep-lengths.pfile',
                                 0, 1)
            
    spectrum_pfile = create_pfile(pfile_dir,
                                  'spectrum.pfile',
                                  2,1)

    pep_dt = open(os.path.join(args.output_dir, 'iterable.dts'), "w")
    pep_dt.write('%d\n\n' % (num_psms))

    pep_num = 0
    for sid, charge in sid_charges:
        s = spectra[sid]
        for p in target[sid,charge]:
            pep = p.peptide
            bNy = theo_dict[sid,charge,pep]
            # look for corresponding dripPSM
            dripPsm = []
            for t_psm in fa_psms[sid,charge]:
                if t_psm.peptide == pep:
                    dripPsm = t_psm
                    break
            assert dripPsm
            drip_peptide_sentence(pep_dt, pep, bNy, 
                                  pep_num, s.spectrum_id, args.max_obs_mass,
                                  peptide_pfile, True, len(bNy)-1)
            training_spectrum_sentence(spectrum_pfile, s.mz, s.intensity,
                                       dripPsm)
            pep_num += 1
    # close streams for this spectrum
    pep_dt.close()
    # compile dt using gmtkDTIndex
    call(['gmtkDTindex', '-decisionTreeFiles', 
          os.path.join(args.output_dir,'iterable.dts')], 
         stdout = stdo, stderr = stde)

    return pep_num

def calculateForcedAlignment(args, spectra, dripMeans, theo_dict,
                             hq_psms, num_psms, stde, stdo):
    """ Run Viterbi to compute a forced alignment, i.e., most probable
        sequence of insertions which are used as observations during 
        training for boot model.  
        Boot model is directly analogous to a speech boot model, which
        acts to regularize learned parameters, greatly decreases overall 
        training time and makes sure that "junk" (too many hypotheses 
        that consist of random matches) do not influence the learned 
        parameters.  Indeed, the number of junk hypotheses is much 
        larger than the number of "useful" hypotheses (hypotheses 
        where true peak matches are considered).
    """
    # create constant gmtkViterbi command line string
    # don't need frame/segment difference actions since each PSM corresponds to a specific spectrum, 
    # so that there isn't much redudandancy to exploit
    vitStr0 = "gmtkViterbi -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 0 -nf1 2 -ni2 1 -nf2 0" \
        + " -fdiffact2 rl" \
        + " -inputMasterFile " + args.master_file + " -inputTrainableParameters trained.params -failOnZeroClique F"

    if args.high_res_ms2:
        ################# unlike low-res mode, define theo_dict and dripMeans herein
        theo_dict, dripMeans = make_fa_data_highres(args, spectra,
                                                    hq_psms, num_psms,
                                                    stde, stdo)
    else:
        make_fa_data_lowres(args, spectra, dripMeans, theo_dict,
                            hq_psms, num_psms, stde, stdo)

    # create structure file and triangulate
    # create structure and master files then triangulate
    try:
        create_drip_structure(args.high_res_ms2, args.structure_file, 
                              args.max_obs_mass)
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        create_drip_master(args.high_res_ms2, args.master_file, 
                           args.max_obs_mass,
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

    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # run GMTK
    dtFile = os.path.join(args.output_dir, 'iterable.dts')
    cppCommand = '\'-DITERABLE_DT=' + dtFile \
        + ' -DDRIP_MZ=' + args.mean_file \
        + ' -DDRIP_GAUSSIAN_COMPONENTS=' + args.gauss_file \
        + ' -DDRIP_GAUSSIAN_MIXTURES=' + args.mixture_file \
        + ' -DDRIP_MZ_GAUSSIANS=' + args.collection_file \
        + '\''

    # call gmtkViterbi
    # gmtkViterbi command line
    vitValsFile = os.path.join(args.logDir, 'vitVals.txt')
    vitStr = vitStr0 + ' -vitValsFile ' +  vitValsFile \
        + ' -of1 ' + pfile_dir + '/spectrum.pfile' \
        + ' -of2 ' + pfile_dir + '/pep-lengths.pfile' \
        + ' -cppCommand ' + cppCommand
    call(shlex.split(vitStr), stdout = stdo, stderr = stde)

    t,_ = parse_dripExtract(vitValsFile, os.path.join(args.output_dir, 'pepDB.txt'))
    return t
    # if not args.high_res_ms2:
    #     return t
    # else:
    #     return t, theo_dict, ion_dict

def runDripTrain(args, spectra, dripMeans, theo_dict, 
                 hq_psms, num_psms, fa_psms, stdo, stde,
                 maxIters = 100):
    """ Run gmtkEMtrain
    """
    make_training_data_lowres(args, spectra, dripMeans, theo_dict,
                              hq_psms, num_psms, fa_psms, stdo, stde)
    # we don't need
    # if args.high_res_ms2:
    #     make_training_data_highres(args, spectra, fa_psms, stdo, stde)
    # else:
    #     make_training_data_lowres(args, spectra, dripMeans, theo_dict,
    #                               hq_psms, num_psms, fa_psms, stdo, stde)

    # create structure and master files then triangulate
    try:
        create_drip_structure(args.high_res_ms2, args.structure_file, 
                              args.max_obs_mass, True, True)
    except:
        print "Could not create DRIP structure file %s, exitting" % args.structure_file
        exit(-1)

    try:
        create_drip_master(args.high_res_ms2, args.master_file, 
                           args.max_obs_mass,
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

    # create static gmtkEMtrain command line string
    emTrainStr0 = "gmtkEMtrain -strFile " + args.structure_file \
        + " -triFile " + args.structure_file + ".trifile -ni1 1 -nf1 2 -ni2 1 -nf2 0" \
        + " -fdiffact2 rl" \
        + " -dirichletPriors T" \
	+ " -objsNotToTrain params.notrain" \
        + " -allocateDenseCpts 1" \
        + " -lldp 1.0e-5" \
        + " -maxEmIters " + str(maxIters) \
        + " -inputMasterFile " + args.master_file + " -inputTrainableParameters trained.params -failOnZeroClique F"

    # create GMTK observation files
    pfile_dir = os.path.join(args.output_dir, args.obs_dir)

    # run GMTK
    dtFile = os.path.join(args.output_dir, 'iterable.dts')
    cppCommand = '\'-DITERABLE_DT=' + dtFile \
        + ' -DDRIP_MZ=' + args.mean_file \
        + ' -DDRIP_GAUSSIAN_COMPONENTS=' + args.gauss_file \
        + ' -DDRIP_GAUSSIAN_MIXTURES=' + args.mixture_file \
        + ' -DDRIP_MZ_GAUSSIANS=' + args.collection_file \
        + '\''

    # call gmtkEMtrain
    # gmtkEMtrain command line
    emTrainStr = emTrainStr0 \
        + ' -outputTrainableParameters ' +  args.dripTrain_file \
        + ' -of1 ' + pfile_dir + '/spectrum.pfile' \
        + ' -of2 ' + pfile_dir + '/pep-lengths.pfile' \
        + ' -cppCommand ' + cppCommand
    # call(shlex.split(emTrainStr), stdout = stdo, stderr = stde)
    return check_output(shlex.split(emTrainStr), stderr=STDOUT)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(conflict_handler='resolve')
    ############## input
    iFileGroup = parser.add_argument_group('iFileGroup', 'Necessary input files.')
    help_psm_library = """<string> - Collection of high-confidence peptide-spectrum matches (PSMs). File must be in either tab-delimited (with fields Peptide and Scan), PIN, pepXML, or mzIdentML format."""
    iFileGroup.add_argument('--psm-library', type = str, action = 'store',
                            help = help_psm_library)
    help_spectra = """<string> - The name of the file from which to parse fragmentation spectra, in ms2 format."""
    iFileGroup.add_argument('--spectra', type = str, action = 'store',
                            help = help_spectra)
    ############## training parameters
    trainingParamsGroup = parser.add_argument_group('trainingParamsGroup', 'Search parameter options.')
    help_dripTrain_file = """<string> - Name of output file for learned parameters. Default = dripLearned.params"""
    trainingParamsGroup.add_argument('--dripTrain-file', type = str, action = 'store', default = "dripLearned.params",
                            help = help_dripTrain_file)
    help_output_mean_file = """<string> - Name of output file for learned Gaussian means. Default = dripLearned.means"""
    trainingParamsGroup.add_argument('--output-mean-file', type = str, action = 'store', default = "dripLearned.means",
                            help = help_output_mean_file)
    help_output_covar_file = """<string> - Name of output file for learned Gaussian covariances. Default = dripLearned.covars"""
    trainingParamsGroup.add_argument('--output-covar-file', type = str, action = 'store', default = "dripLearned.covars",
                            help = help_output_covar_file)
    help_high_res_ms2 = """<T|F> - boolean, whether the search is over high-res ms2 (high-high) spectra. When this parameter is true, DRIP used the real valued masses of candidate peptides as its Gaussian means. For low-res ms2 (low-low or high-low), the observed m/z measures are much less accurate so these Gaussian means are learned using training data (see dripTrain). Default=False."""
    trainingParamsGroup.add_argument('--high-res-ms2', type = str, action = 'store', default = 'false', help = help_high_res_ms2)
    # help_num_threads = '<integer> - the number of threads to run on a multithreaded CPU. If supplied value is greater than number of supported threads, defaults to the maximum number of supported threads minus one. Multithreading is not suppored for cluster use as this is typically handled by the cluster job manager. Default=1.'
    # trainingParamsGroup.add_argument('--num-threads', type = int, action = 'store', 
    #                                default = 1, help = help_num_threads)
    ############## amino acid modifications
    aaModsGroup = parser.add_argument_group('aaModsGroup', 'Options for amino acid modifications.')
    aaModsGroup.add_argument('--mods-spec', type = str, action = 'store',
                             default = 'C+57.02146')
    aaModsGroup.add_argument('--cterm-peptide-mods-spec', type = str, action = 'store',
                             default = '')
    aaModsGroup.add_argument('--nterm-peptide-mods-spec', type = str, action = 'store',
                             default = '')
    ############ shard spectra options
    parser.add_argument('--output-dir', type = str, default = 'encode')
    parser.add_argument('--obs-dir', type = str, action = 'store',
                        default = 'obs')
    parser.add_argument('-z', '--use_gzcat', action = 'store_true',
                        dest = 'gzcat', default = False,
                        help = "Use gzcat to decompress .ms2.gz files, if avail.")
    ######### encoding options
    parser.add_argument('--logDir', type = str, 
                        help = 'directory with output log files',
                        default = 'log')
    parser.add_argument("--max_obs_mass", type = int, default = 2001)
    parser.add_argument('--normalize',  type = str, action = 'store',
                        help = "Name of the spectrum preprocessing pipeline.", 
                        default = 'top300Sequest')
                        # default = 'top300TightSequest')
    help_filt_theo_peaks = "Filter theoretical peaks outside of observed spectra min and max"
    parser.add_argument('--filt_theo_peaks', action = 'store_false', dest = 'filt_theo_peaks', 
                        default = True, help = help_filt_theo_peaks)
    parser.add_argument('--per_spectrum_mz_bound', action = 'store_true', 
                        default = False, help = "Calculate observed m/z bound per spectrum")
    parser.add_argument('--mz_lb', type = float, action = 'store',
                        dest = 'mz_lb', help = 'Lower m/z bound' ,default = 150.0)
    parser.add_argument('--mz_ub', type = float, action = 'store',
                        dest = 'mz_ub', help = 'Upper m/z bound' ,default = 2000.0)
    parser.add_argument('--structure-file', type = str, action = 'store',
                        default = 'drip.str',
                        help = "DRIP GMTK structure file")
    parser.add_argument('--master-file', type = str, action = 'store',
                        default = 'drip.mtr',
                        help = "DRIP GMTK master file")
    parser.add_argument('--collection-dir', action = 'store', dest = 'collection_dir',
                        type = str, help = 'Where to store Gaussian files', 
                        default = 'drip_collection')
    parser.add_argument('--learned-means', action = 'store', dest = 'learned_means',
                        type = str, help = 'Previously learned means', default = "")
    parser.add_argument('--mean-file', action = 'store', dest = 'mean_file',
                        type = str, help = 'Where to put mean file.',
                        default = 'drip_mz.txt')
    parser.add_argument('--covar-file', action = 'store',
                        type = str, help = 'Covariance file name', 
                        default = "covar.txt")
    parser.add_argument('--gauss_file', action = 'store', dest = 'gauss_file',
                        type = str, help = 'Gaussians.', 
                        default = 'drip_Gaussian_components.txt')

    parser.add_argument('--mixture_file', action = 'store', dest = 'mixture_file',
                        type = str, help = 'Gaussian mixtures.',
                        default = 'drip_Gaussian_mixtures.txt')

    parser.add_argument('--collection_file', action = 'store', dest = 'collection_file',
                        type = str, help = 'Where to put collection file.',
                        default = 'drip_mz_Gaussians.txt')

    parser.add_argument('--covar_name', action = 'store', dest = 'covar_name',
                        type = str, help = 'Variance name.',
                        default = 'covar0')

    parser.add_argument('--mixture_name', action = 'store', dest = 'mixture_name',
                        type = str, help = 'Gaussian mixture name.',
                        default = 'mixture')

    parser.add_argument('--gaussian_component_name', action = 'store', dest = 'gaussian_component_name',
                        type = str, help = 'Gaussian mixture name.',
                        default = 'gc')

    parser.add_argument('--collection_name', action = 'store', dest = 'collection_name',
                        type = str, help = 'Gaussian mixture name.', 
                        default = 'drip_mz_Gaussians')

    parser.add_argument('--dpmf_name', action = 'store', dest = 'dpmf_name',
                        type = str, help = 'DPMF name.', 
                        default = 'unityDPMF')
    args = parser.parse_args()

    # process input arguments
    process_args(args)

    stdo = sys.stdout
    stde = stdo
    # # set stdout and stderr for subprocess
    # if debug:
    #     # stdo = sys.stdout
    #     stdo = open(os.devnull, "w")
    #     stde = stdo
    # else:
    #     stdo = open(os.devnull, "w")
    #     # stdo = sys.stdout
    #     stde = stdo
    #     sys.stdout = open("dripTrain_output", "w")

    # currently ignore ident file input for spectra filtering
    spectra, minMz, maxMz, validcharges, sidChargePreMass = load_spectra_ret_dict(args.spectra, 'all')
    
    # update encountered charges
    args.charges = validcharges
    args.mz_lb = minMz
    args.mz_ub = maxMz

    print args.covar_file
    try:
        write_covar_file(args.high_res_ms2, args.covar_file, '', False)
    except:
        print "Could not create covariance file %s, exitting" % args.covar_file
        exit(-1)

    if not args.high_res_ms2:
        # compute collection of theoretical spectra
        # theo_dict, usedTheoPeaks, hq_psms, num_psms, observedHoles = calculate_theoretical_dictionary(args, spectra)
        theo_dict, usedTheoPeaks, hq_psms, num_psms, observedHoles, fa_psms = calcTheoDict_estimateFa(args, spectra)

        # load drip means
        # load prior learned means
        if args.learned_means: # prior on the means
            dripMeans = load_drip_means(args.learned_means)
        else:
            dripMeans = {}
            for x in range(args.max_obs_mass+1):
                dripMeans[x] = x+0.5

        dripMeans0 = dripMeans.copy()

        # remove theoretical holes
        remove_theoretical_holes(dripMeans, usedTheoPeaks)

        # filter theoretical peaks out of observed spectra range 
        # and remap theoretical peaks to space of active theoretical peaks
        remap_mean_indices(args, spectra,
                           dripMeans, usedTheoPeaks, theo_dict)
    else:
        dripMeans = []
        theo_dict = []
        observedHoles = []
        hq_psms,num_psms = load_psm_library(args.psm_library)
        # compute forced alignment
        stde = open('dripTrain_fa', "w")
        fa_psms = calculateForcedAlignment(args, spectra, dripMeans, theo_dict, 
                                           hq_psms, num_psms, stdo, stde)
        stde.close()

    # after call to calculateForcedAlignment, should be no difference betweent running
    # lowres and high-res (from the perspective of all means being loaded into memory
    # with the same data-type

    # create necessary gmtkEMtrain files
    create_empty_master()
    create_params_no_train()

    meanPat = 'WARNING: Mean vec \'mean(?P<meanNum>\d+)\''

    # stde = open('dripTrain', "w")

    # run one iteration of EM to see whether we must inject 
    # synthetic peaks
    k = runDripTrain(args, spectra, dripMeans, theo_dict,
                     hq_psms, num_psms, fa_psms, stdo, stde,
                     1)

    lowProbMeans = []
    for meanNum in re.findall(meanPat, k):
        lowProbMeans.append(int(meanNum))

    # if debug:
    #     new_theo_index_mapping = {}
    #     new_to_old_mapping = {}
    #     for new_ind, old_ind in enumerate(sorted(dripMeans.iterkeys())):
    #         new_theo_index_mapping[old_ind] = new_ind
    #         new_to_old_mapping[new_ind] = old_ind

    #     for m in lowProbMeans:
    #         print new_theo_index_mapping[m]

    while lowProbMeans:        
        print "%d means without enough evidence, adding synthetic peaks" % len(lowProbMeans)
        inject_mean_evidence(args, spectra,
                             dripMeans,
                             fa_psms, theo_dict, lowProbMeans,
                             observedHoles)

        k = runDripTrain(args, spectra, dripMeans, theo_dict,
                         hq_psms, num_psms, fa_psms, stdo, stde,
                         1)
        lowProbMeans = []
        for meanNum in re.findall(meanPat, k):
            lowProbMeans.append(int(meanNum))

    k = runDripTrain(args, spectra, dripMeans, theo_dict,
                     hq_psms, num_psms, fa_psms, stdo, stde,
                     100)
    if debug:
        print k

    # parse trained output
    gm, gc, dp, m, c, dc = parse_params(args.dripTrain_file)
    # parsed fields/output format:
    # mean[0]: name, mean[1]: dimension, mean[2]: value
    # covar[0]: name, covar[1]: dimension, covar[2]: value
    # dpmf[0]: name, dpmf[1]: dimension, dpmf[2]: value
    # cpt[0] - name, cpt[1] - dimension(s), cpt[2] - flattened cpt

    # gaussian_component parameters: 
    # gaussian_component[0]: name, gaussian_component[1]: dimension, 
    # gaussian_component[2]: type
    # gaussian_component[3]: mean
    # gaussian_component[4]: covariance

    # gaussian_mixture[0]: name, gaussian_mixture[1]: mixture dim, gaussian_mixture[2]: num components, gaussian_mixture[3]: dpmf name, gaussian_mixture[4]: components

    remap_learned_means(dripMeans0, m)

    baseName = args.dripTrain_file.split('.')

    assert len(baseName)==2 or len(baseName)==1, 'Output name %s has multiple periods in it\'s name, cannot distinguish file type.  Exitting' % args.dripTrain_file


    # write learned means and covariances
    write_learned_means_covars(dripMeans0, c, args.output_mean_file, args.output_covar_file)

    if stdo:
        stdo.close()
    if stde:
        stde.close()
