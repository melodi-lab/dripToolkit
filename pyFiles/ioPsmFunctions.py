#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

import sys
import csv
import itertools
import types
import numpy
import math
import scipy.interpolate
import collections

from pyFiles.psm import PSM
####### how to use linear interpolation:
####### y_interp = scipy.interpolate.interp1d(x,y)
####### new_range = range(10)
####### y_interp(new_range)

def load_peptides(filename):
    """ Return dictionary
    """
    peptides = {}
    ind = 0
    with open(filename) as inputf:
        for line in inputf:
            peptides[line.split()[0]] = ind
            ind += 1
    return peptides

def load_peptide_masses(filename):
    """ Return array of peptide strings
    """
    peptides = []
    with open(filename) as inputf:
        for line in inputf:
            peptides.append(float(line.split()[1]))
    return peptides

def load_peptide_strings(filename):
    """ Return array of peptide strings
    """
    peptides = []
    with open(filename) as inputf:
        for line in inputf:
            peptides.append(line.split()[0])
    return peptides

def write_merged_ident(output_file, targets, decoys, tc, dc, 
                       target_db, decoy_db):
    targets.sort(key = lambda r: r[0])
    decoys.sort(key = lambda r: r[0])
    # load target and decoy databases
    target_peptides = load_peptide_strings(target_db)
    decoy_peptides = load_peptide_strings(decoy_db)
    # write new ident
    print "Writing new ident file %s" % (output_file)
    outfile = open(output_file, "w")
    outfile.write("Kind\tSid\tPeptide\tScore\tCharge\n")
    for t1, d1, ct, cd in zip(targets, decoys, tc, dc):
        outfile.write("t\t%d\t%s\t%f\t%d\n" % (t1[0], target_peptides[t1[1]], t1[2], ct))
        outfile.write("d\t%d\t%s\t%f\t%d\n" % (d1[0], decoy_peptides[d1[1]], d1[2], cd))

    outfile.close()

def merge_idents(targets1, decoys1, targets2, decoys2, c1, c2):
        # (targets, decoys): each is a list of records of the form (s, p, f)
        # where s = spectrum id, an integer; p = best peptide, an instance of
        # protein.peptide.Peptide; and, f = score, a floating point score for
        # the peptide spectrum match.

    ident1_sids = [targets1[i][0] for i in range(len(targets1))]
    ident2_sids = [targets2[i][0] for i in range(len(targets2))]
    ident_sid_intersection = set(ident1_sids)&set(ident2_sids)

    merged_targets = []
    merged_decoys = []
    tcharges = []
    dcharges = []
    # append PSMs not in intersection
    # create two vectors of the intersection
    ident1_target_intersection = []
    ident1_decoy_intersection = []
    ident2_target_intersection = []
    ident2_decoy_intersection = []
    if ident_sid_intersection:
        for t,d in zip(targets1, decoys1):
            if t[0] not in ident_sid_intersection:
                merged_targets.append(t)
                merged_decoys.append(d)
                tcharges.append(c1)
                dcharges.append(c1)
            else:
                ident1_target_intersection.append(t)
                ident1_decoy_intersection.append(d)
        for t,d in zip(targets2, decoys2):
            if t[0] not in ident_sid_intersection:
                merged_targets.append(t)
                merged_decoys.append(d)
                tcharges.append(c2)
                dcharges.append(c2)                
            else:
                ident2_target_intersection.append(t)
                ident2_decoy_intersection.append(d)

    ident1_target_intersection.sort(key=lambda r: r[0])
    ident1_decoy_intersection.sort(key=lambda r: r[0])
    ident2_target_intersection.sort(key=lambda r: r[0])
    ident2_decoy_intersection.sort(key=lambda r: r[0])
    if len(ident1_target_intersection) != len(ident2_target_intersection):
        print "Error: target ident intersections not of equal length"
        exit(-11)
    if len(ident1_decoy_intersection) != len(ident2_decoy_intersection):
        print "Error: decoy ident intersections not of equal length"
        exit(-11)

    for t1,d1,t2,d2 in zip(ident1_target_intersection, ident1_decoy_intersection, 
                           ident2_target_intersection, ident2_decoy_intersection):
        if t1[2]>=t2[2]: # t1 score is greater than t2 score
            merged_targets.append(t1)
            tcharges.append(c1)
        else:
            merged_targets.append(t2)
            tcharges.append(c2)
        if d1[2]>=d2[2]: # d1 score is greater than d2 score
            merged_decoys.append(d1)
            dcharges.append(c1)
        else:
            merged_decoys.append(d2)
            dcharges.append(c2)

    # can't sort since then charges won't be correct
    # merged_targets.sort(key=lambda r: r[0])
    # merged_decoys.sort(key=lambda r: r[0])

    return merged_targets, merged_decoys, tcharges, dcharges

def write_ident(output_file, targets, decoys, target_db, decoy_db):
    # load target and decoy databases
    target_peptides = load_peptide_strings(target_db)
    decoy_peptides = load_peptide_strings(decoy_db)
    # write new ident with maximal shift
    print "Writing new ident file %s" % (output_file)
    outfile = open(output_file, "w")
    outfile.write("Kind\tSid\tPeptide\tScore\n")
    for t1, d1 in zip(targets, decoys):
        outfile.write("t\t%d\t%s\t%f\n" % (t1[0], target_peptides[t1[1]], t1[2]))
        outfile.write("d\t%d\t%s\t%f\n" % (d1[0], decoy_peptides[d1[1]], d1[2]))

    outfile.close()

def comp_psm_files(file0, file1):
    p0 = load_drip_output(file0)
    p1 = load_drip_output(file1)

    disc = 0

    print "Comparing PSMs in %s and %s" % (file0, file1)

    for sidCharge in (set(p0.iterkeys())-set(p1.iterkeys())):
        print "Scan %d, charge %d PSMs not in %s but in %s" % (sidCharge[0],
                                                               sidCharge[1],
                                                               file1, file0)
        disc += 1
    for sidCharge in (set(p1.iterkeys())-set(p0.iterkeys())):
        print "Scan %d, charge %d PSMs not in %s but in %s" % (sidCharge[0],
                                                               sidCharge[1],
                                                               file0, file1)
        disc += 1
    # set order for outputting keys in PSM().other dict
    keys = []
    for sid, charge in p0:
        for p in p0[sid,charge]:
            for k in p.other:
                keys.append(k)
            break
        break

    for sid, charge in p0:
        if (sid,charge) in p1:
            for p in (set([p for p in p0[sid,charge]]) - set([p for p in p1[sid,charge]])):
                if p.kind == 'd':
                    continue
                print "Scan %d charge %d, peptide %s not in file %s" % (sid, charge, 
                                                                        p.peptide, file1)
                disc += 1
            for p in (set([p for p in p1[sid,charge]]) - set([p for p in p0[sid,charge]])):
                if p.kind == 'd':
                    continue
                print "Scan %d charge %d, peptide %s not in file %s" % (sid, charge, 
                                                                        p.peptide, file0)
                disc += 1
            currP0Dict = {}
            for p in p0[sid,charge]:
                if p.kind == 'd':
                    continue
                currP0Dict[p] = p

            for p in p1[sid,charge]:
                if p in currP0Dict:
                    q = currP0Dict[p]
                    if p.score != q.score:
                        print "Scan %d charge %d: f1 score=%f, f2 score=%f" % (sid, charge,
                                                                               p.score, q.score)
                        disc += 1
                    po = p.other
                    qo = q.other
                    for k in keys:
                        if po[k] != qo[k]:
                            print "Scan %d charge %d: f1 %s=%s, f2 %s=%s" % (sid, charge,
                                                                             po[k], qo[k])
                            disc += 1
                # for k in keys:
                #     try:
                #         el0 = p0[k]
                #     except KeyError:
                #         print "Header fields are different between the two files, exitting"
                #         exit(-1)
                #     try:
                #         el1 = p1[k]
                #     except KeyError:
                #         print "Header fields are different between the two files, exitting"
                #         exit(-1)

                #     if el0 != el1:
                #         print "Scan %d charge %d: field %s: f0=%s, f1=%s" % (sid, charge,
                #                                                              k, el0, el1)

    print "%d target discrepancies found" % disc

def load_psms(filename):
    """Load an identification file, the output of spectrum identification.

    Standard DRIP output fields:
    (1)Kind (2)Sid (3)Frames (4)Score (5)Peptide (6)Obs_Inserts	(7)Theo_Deletes	(8)Obs_peaks_scored	
    (9)Theo_peaks_used	(10)Sum_obs_intensities	(11)Sum_scored_mz_dist	(12)Charge

    Assume the tab-delimited file has the following fields:
    Kind, Score, Sid, Charge

    Todo: change Sid to scan
    """
    targets = {}
    decoys = {}

    f = open(filename)
    reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]

    p = reader[0]

    if "Sid" in p:
        sidKey = "Sid"
    else:
        sidKey = "Scan"

    reader.sort(key = lambda r: int(r[sidKey]))

    if sidKey not in p:
        print "Expected field %s not present in file %s, exitting" % (sidKey, filename)
        exit(-1)
    if "Kind" not in p:
        print "Expected field Kind not present in file %s, exitting" % filename
        exit(-1)
    if "Score" not in p:
        print "Expected field Score not present in file %s, exitting" % filename
        exit(-1)
    if "Charge" not in p:
        print "Expected field Charge not present in file %s, exitting" % filename
        exit(-1)

    num_psms = 0

    keys = []
    psmKeys = set(["Kind", sidKey, "Charge", "Score", "Peptide"])
    for i in p:
        if i not in psmKeys:
            keys.append(i)

    for sid, rows in itertools.groupby(reader, lambda r: int(r[sidKey])):
        for p in rows:
            kind = p["Kind"].lower()
            try:
                c = int(p["Charge"])
            except TypeError:
                print "Could not convert charge %s for sid %d to int, exitting" % (p["Charge"], sid)
                exit(-1)

            el = {}
            if keys:
                for k in keys:
                    el[k] = p[k]

            el = PSM(p["Peptide"],
                     float(p["Score"]),
                     int(p[sidKey]),
                     p["Kind"],
                     int(p["Charge"]),
                     el)

            if kind=='t' or kind=='target':
                if (sid, c) in targets:
                    targets[sid, c].append(el)
                else:
                    targets[sid, c] = [el]
                num_psms += 1
            elif kind=='d' or kind=='decoy':
                if (sid, c) in decoys:
                    decoys[sid, c].append(el)
                else:
                    decoys[sid, c] = [el]
                num_psms += 1
    f.close()
    return targets, decoys, num_psms

def load_psm_library(filename):
    """Load high-confidence PSMs from a tab-delimited file with fields: Peptide, Scan, Charge
    """
    psms = {}
    f = open(filename)
    reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]
    reader.sort(key = lambda r: int(r["Scan"]))

    p = reader[0]
    if "Scan" not in p:
        print "Expected field Scan not present in file %s, exitting" % filename
        exit(-1)
    if "Peptide" not in p:
        print "Expected field Peptide not present in file %s, exitting" % filename
        exit(-1)
    if "Charge" not in p:
        print "Expected field Charge not present in file %s, exitting" % filename
        exit(-1)

    num_psms = 0

    keys = []
    psmKeys = set(["Scan", "Charge", "Peptide"])
    for i in p:
        if i not in psmKeys:
            keys.append(i)

    for sid, rows in itertools.groupby(reader, lambda r: int(r["Scan"])):
        for p in rows:
            try:
                c = int(p["Charge"])
            except TypeError:
                print "Could not convert charge %s for scan %d to int, exitting" % (p["Charge"], sid)
                exit(-1)

            el = {}
            if keys:
                for k in keys:
                    el[k] = p[k]

            el = PSM(p["Peptide"],
                     0.0,
                     int(p["Scan"]),
                     "t",
                     c,
                     el)
            if (sid, c) in psms:
                psms[sid, c].append(el)
            else:
                psms[sid, c] = [el]
            num_psms += 1

    f.close()
    return psms, num_psms

def load_drip_output(filename):
    """ Load all PSMs output by DRIP, or any tab-delimited output of a mass-spec experiment with field "Scan" to denote
        the spectrum identification number
        Todo: add a parser to load PSMs from a DRIP run, returning each PSM as an instance of the 
        dripPSM class

        Normal DRIP header fields:
        Kind, Sid or Scan, Frames, Score, Peptide, Obs_Inserts, Theo_Deletes, Obs_peaks_scored, Theo_peaks_used, Sum_obs_intensities, Sum_scored_mz_dist, Charge
    """

    all_psms = {} # dictionary with keys (x,y) where x = spectrum scan number, y = spectrum charge
    with open(filename, 'r') as f:
        reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]
    l = reader[0]
    if "Sid" not in l and "Scan" not in l:
        raise ValueError("No Scan/Sid field, exitting")
    if "Charge" not in l and "Ch" not in l:
        raise ValueError("No Charge field, exitting")
    if "Sid" in l:
        sidKey = "Sid"
    else:
        sidKey = "Scan"

    if "Charge" in l:
        chargeKey = "Charge"
    elif "Ch" in l:
        chargeKey = "Ch"
    else:
        print "Error: no Charge or Ch field in search identification file, exitting"
        exit(-1)

    psmKeys = set(["Kind", sidKey, chargeKey, "Score", "Peptide"])
    keys = []
    for k in l:
        if k not in psmKeys:
            keys.append(k)
            
    dripKeys = ["Kind", sidKey, "Frames", "Score", "Peptide", "Obs_Inserts", "Theo_Deletes", 
                "Obs_peaks_scored", "Theo_peaks_used", "Sum_obs_intensities", "Sum_scored_mz_dist", 
                chargeKey]

    for i, l in enumerate(reader):
        try:
            sid = int(l[sidKey])
        except ValueError:
            print "Could not convert scan number %s on line %d to int, exitting" % (l[sidKey], i+1)
        try:
            charge = int(l[chargeKey])
        except ValueError:
            print "Could not convert charge %s on line %d to int, exitting" % (l[chargeKey], i+1)

        if (sid,charge) not in all_psms:
            all_psms[sid,charge] = []

        el = {}
        for k in keys:
            el[k] = l[k]

        try:
            el = PSM(l["Peptide"],
                     float(l["Score"]),
                     int(l[sidKey]),
                     l["Kind"],
                     int(l[chargeKey]),
                     el)
        except KeyError:
            print "Standard PSM field not encountered, exitting"
            exit(-1)

        all_psms[sid,charge].append(el)

    return all_psms

def load_drip_psms(filename, training_decoys_by_charge, psms_by_charge):
    """ Load all PSMs output by DRIP, or any tab-delimited output of a mass-spec experiment with field "Scan" to denote
        the spectrum identification number
        Todo: add a parser to load PSMs from a DRIP run, returning each PSM as an instance of the 
        dripPSM class

        Normal DRIP header fields:
        Kind, Sid or Scan, Frames, Score, Peptide, Obs_Inserts, Theo_Deletes, Obs_peaks_scored, Theo_peaks_used, Sum_obs_intensities, Sum_scored_mz_dist, Charge
    """
    with open(filename, 'r') as f:
        reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]
    l = reader[0]
    if "Sid" not in l and "Scan" not in l:
        raise ValueError("No Scan/Sid field, exitting")
    if "Charge" not in l and "Ch" not in l:
        raise ValueError("No Charge field, exitting")
    if "Sid" in l:
        sidKey = "Sid"
    else:
        sidKey = "Scan"

    if "Charge" in l:
        chargeKey = "Charge"
    elif "Ch" in l:
        chargeKey = "Ch"
    else:
        print "Error: no Charge or Ch field in search identification file, exitting"
        exit(-1)

    psmKeys = set(["Kind", sidKey, chargeKey, "Score", "Peptide"])
    keys = []
    for k in l:
        if k not in psmKeys:
            keys.append(k)
            
    dripKeys = ["Kind", sidKey, "Frames", "Score", "Peptide", "Obs_Inserts", "Theo_Deletes", 
                "Obs_peaks_scored", "Theo_peaks_used", "Sum_obs_intensities", "Sum_scored_mz_dist", 
                chargeKey]

    for i, l in enumerate(reader):

        kind = l['Kind']

        try:
            sid = int(l[sidKey])
        except ValueError:
            print "Could not convert scan number %s on line %d to int, exitting" % (l[sidKey], i+1)
        try:
            charge = int(l[chargeKey])
        except ValueError:
            print "Could not convert charge %s on line %d to int, exitting" % (l[chargeKey], i+1)

        el = {}
        for k in keys:
            el[k] = l[k]
        try:
            el = PSM(l["Peptide"],
                     float(l["Score"]),
                     int(l[sidKey]),
                     l["Kind"],
                     int(l[chargeKey]),
                     el)
        except KeyError:
            print "Standard PSM field not encountered, exitting"
            exit(-1)

        if kind == 'r':
            # recalibration decoys, only need scores
            if charge in training_decoys_by_charge:
                training_decoys_by_charge[charge].append(float(l["Score"]))
            else:
                training_decoys_by_charge[charge] = [float(l["Score"])]
        else:
            if charge in psms_by_charge:
                psms_by_charge[charge].append(el)
            else:
                psms_by_charge[charge] = [el]

def load_pin_file(filename):
    """ Load all PSMs and features from a percolator PIN file, or any tab-delimited output of a mass-spec experiment with field "Scan" to denote
        the spectrum identification number
        Todo: add a parser to load PSMs from a DRIP run, returning each PSM as an instance of the 
        dripPSM class

        Normal tide with DRIP features
        SpecId	Label	ScanNr	lnrSp	deltLCn	deltCn	score	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	enzN	enzC	enzInt	lnNumSP	dm	absdM	insertions	deletions	peaksScoredA	theoPeaksUsedA	SumScoredIntensities	SumScoredMzDist	Peptide	Proteins
    """

    targets = {}
    decoys = {}
    with open(filename, 'r') as f:
        reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]
    l = reader[0]
    if "ScanNr" not in l:
        raise ValueError("No ScanNr field, exitting")
    if "Charge1" not in l:
        raise ValueError("No Charge1 field, exitting")

    # spectrum identification key for PIN files
    sidKey = "ScanNr" # note that this typically denotes retention time
    
    numPeps = 0

    maxCharge = 1
    chargeKeys = set([])
    # look at score key and charge keys
    scoreKey = ''
    for i in l:
        m = i.lower()
        if m == 'score':
            scoreKey = i
        if m[:-1]=='charge':
            chargeKeys.add(i)
            maxCharge = max(maxCharge, int(m[-1]))

    if not scoreKey:
        for i in l:
            if i.lower() == 'xcorr':
                scoreKey = i            

    # fields we have to keep track of
    psmKeys = set(["SpecId", "Label", sidKey, scoreKey, "Peptide", "Proteins"])
    keys = []
    for k in l:
        if k not in psmKeys and k not in chargeKeys:
            keys.append(k)
            
    for i, l in enumerate(reader):
        try:
            sid = int(l[sidKey])
        except ValueError:
            print "Could not convert scan number %s on line %d to int, exitting" % (l[sidKey], i+1)

        charge = 0
        # look for current PSM, encoded as a one-hot vector
        for c in chargeKeys:
            try:
                charge = int(l[c])
            except ValueError:
                print "Could not convert charge %s on line %d to int, exitting" % (l[c], i+1)

            if charge:
                charge = int(c[-1])
                break

        assert charge > 0, "No charge denoted with value 1 for PSM on line %d, exitting" % (i+1)

        el = {}
        for k in keys:
            el[k] = l[k]

        if l["Label"] == '1':
            kind = 't'
        elif l["Label"] == '-1':
            kind = 'd'
        else:
            print "Error: encountered label value %s, can only be -1 or 1, exitting" % l["Label"]
            exit(-1)

        try:
            el = PSM(l["Peptide"],
                     float(l[scoreKey]),
                     int(l[sidKey]),
                     kind,
                     charge,
                     el,
                     l["Proteins"],
                     l["SpecId"])
        except KeyError:
            print "Standard PSM field not encountered, exitting"
            exit(-1)

        if kind == 't':
            if (sid,charge) not in targets:
                targets[sid,charge] = []
            targets[sid,charge].append(el)
            numPeps += 1
        elif kind == 'd':
            if (sid,charge) not in decoys:
                decoys[sid,charge] = []
            decoys[sid,charge].append(el)
            numPeps += 1

    return targets,decoys,numPeps

def load_pin_return_dict(filename):
    """ Load all PSMs and features from a percolator PIN file, or any tab-delimited output of a mass-spec experiment with field "Scan" to denote
        the spectrum identification number
        Todo: add a parser to load PSMs from a DRIP run, returning each PSM as an instance of the 
        dripPSM class

        Normal tide with DRIP features
        SpecId	Label	ScanNr	lnrSp	deltLCn	deltCn	score	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	enzN	enzC	enzInt	lnNumSP	dm	absdM	insertions	deletions	peaksScoredA	theoPeaksUsedA	SumScoredIntensities	SumScoredMzDist	Peptide	Proteins
    """

    targets = {}
    decoys = {}
    with open(filename, 'r') as f:
        reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]
    l = reader[0]
    if "ScanNr" not in l:
        raise ValueError("No ScanNr field, exitting")
    if "Charge1" not in l:
        raise ValueError("No Charge1 field, exitting")

    # spectrum identification key for PIN files
    sidKey = "ScanNr" # note that this typically denotes retention time

    numPeps = 0

    maxCharge = 1
    chargeKeys = set([])
    # look at score key and charge keys
    scoreKey = ''
    for i in l:
        m = i.lower()
        if m == 'score':
            scoreKey = i
        if m[:-1]=='charge':
            chargeKeys.add(i)
            maxCharge = max(maxCharge, int(m[-1]))

    if not scoreKey:
        for i in l:
            if i.lower() == 'xcorr':
                scoreKey = i            

    # fields we have to keep track of
    psmKeys = set(["SpecId", "Label", sidKey, scoreKey, "Peptide", "Proteins"])
    keys = list(set(l.iterkeys()) - psmKeys)
    for i, l in enumerate(reader):
        try:
            sid = int(l[sidKey])
        except ValueError:
            print "Could not convert scan number %s on line %d to int, exitting" % (l[sidKey], i+1)

        charge = 0
        # look for current PSM, encoded as a one-hot vector
        for c in chargeKeys:
            try:
                charge = int(l[c])
            except ValueError:
                print "Could not convert charge %s on line %d to int, exitting" % (l[c], i+1)

            if charge:
                charge = int(c[-1])
                break

        assert charge > 0, "No charge denoted with value 1 for PSM on line %d, exitting" % (i+1)

        el = {}
        for k in keys:
            el[k] = l[k]

        if l["Label"] == '1':
            kind = 't'
        elif l["Label"] == '-1':
            kind = 'd'
        else:
            print "Error: encountered label value %s, can only be -1 or 1, exitting" % l["Label"]
            exit(-1)

        try:
            el = PSM(l["Peptide"],
                     float(l[scoreKey]),
                     int(l[sidKey]),
                     kind,
                     charge,
                     el,
                     l["Proteins"],
                     l["SpecId"])
        except KeyError:
            print "Standard PSM field not encountered, exitting"
            exit(-1)

        if kind == 't':
            targets[el.scan, el.peptide] = el # hash key for PSM class is (scan, peptide string), so
                                              # we shouldn't get collisions without adding charge as a key
            numPeps += 1
        elif kind == 'd':
            decoys[el.scan, el.peptide] = el # hash key for PSM class is (scan, peptide string), so
                                             # we shouldn't get collisions without adding charge as a key
            numPeps += 1

    return targets,decoys

def load_ident_all_charge_decoys(filename, decoy_db, sids = None):
    """Load an all decoy psms identification file, the output psms of a dataset.
       Note that this can be quite a large file, and care is taken to help ensure
       all decoys will fit in memory; this is done by mapping all peptide sequences
       to the integer index in their respective database

    Arguments:
        filename: Name of the tab-separated file, with Kind, Sid, Peptide,
            and Score fields.
        decoy_db: the decoy database, in crux digest output format of:
                  peptide_string\tpeptide_mass\n
        sids: Set of sids to load. If None, load all records.
        charge: whether the experiment involved multiple charges

    Returns:
        all_decoys, a list of tuples d such that:
        d[0] = decoy peptide string length
        d[1] = spectrum id
        d[2] = score
        d[3] = charge(if input charge=True)
    """
    all_decoys = []
    f = open(filename)
    reader = csv.DictReader(f, delimiter = '\t', skipinitialspace = True)

    # load target and decoy databases
    decoy_peptides = load_peptides(decoy_db)

    decoyf = lambda r: r["Kind"] == "d"
    for sid, rows in itertools.groupby(reader, lambda r: int(r["Sid"])):
        if not sids or sid in sids:
            records = list(rows)
            try:
                for bdecoy in itertools.ifilter(decoyf, records):
                    # save memory by storing peptides as ints; note, this temporarily bolsters current memory usage
                    try:
                        decoy_db_ind = decoy_peptides[bdecoy["Peptide"]]
                    except:
                        print "Decoy peptide %s not in database %s, skipping" % (bdecoy["Peptide"], decoy_db)
                        continue

                    all_decoys.append( float(bdecoy["Score"]) )
            except ValueError:
                print 'Record %d in %s is bad' % (reader.line_num, f.name)
                print 'Headers in bad ident file: %s' % reader.fieldnames
                raise

    f.close()
    return sorted(list(set(all_decoys)))

def load_ident_all_psms(filename, target_db, decoy_db,
                        sids = None, charge = False):
    """Load an identification file, the output of spectrum identification.

    Arguments:
        filename: Name of the tab-separated file, with Kind, Sid, Peptide,
            and Score fields.
        sids: Set of sids to load. If None, load all records.

    Returns:
        (targets, decoys): each is a list of records of the form (s, p, f, c)
        where s = spectrum id, an integer; p = best peptide, an instance of
        protein.peptide.Peptide; f = score, a floating point score for
        the peptide spectrum match; and c = charge, an integer returning the 
        identified charge.


        All the lines with Kind == 't' will be placed in the targets list;
        all the lines with Kind == 'd' will be placed in the decoys list.

    Post-conditions:
        1. len(targets) == len(decoys).
        2. all(t[0] == d[0] for t, d in zip(targets, decoy). The entries in
           targets and decoys are sorted by spectrum id.

    """
    all_psms = []
    f = open(filename)
    reader = csv.DictReader(f, delimiter = '\t', skipinitialspace = True)

    # load target and decoy databases
    target_peptides = load_peptides(target_db)
    decoy_peptides = load_peptides(decoy_db)

    targetf = lambda r: r["Kind"] == "t"
    decoyf = lambda r: r["Kind"] == "d"
    scoref = lambda r: float(r["Score"])
    for sid, rows in itertools.groupby(reader, lambda r: int(r["Sid"])):
        if not sids or sid in sids:
            records = list(rows)
            try:
                for btarget in itertools.ifilter(targetf, records):
                    # save memory by storing peptides as ints; note, this temporarily bolsters current memory usage
                    try:
                        target_db_ind = target_peptides[btarget["Peptide"]]
                    except:
                        print "Target peptide %s not in database %s, exitting" % (btarget["Peptide"], target_db)
                        exit(1)

                    # calculate target length
                    target_length = len(btarget["Peptide"])
                    if charge:
                        all_psms.append( [ True, target_length, sid, target_db_ind,
                                         float(btarget["Score"]), int(btarget["Charge"])])
                    else:
                        all_psms.append( [True, target_length, sid, target_db_ind,
                                         float(btarget["Score"])])

                for bdecoy in itertools.ifilter(decoyf, records):
                    try:
                        decoy_db_ind = decoy_peptides[bdecoy["Peptide"]]
                    except:
                        print "Decoy peptide %s not in database %s, skipping" % (bdecoy["Peptide"], decoy_db)
                        continue

                    # calculate decoy length
                    decoy_length = len(bdecoy["Peptide"])

                    if charge:
                        all_psms.append( [False, decoy_length, sid, decoy_db_ind,
                                        float(bdecoy["Score"]), int(bdecoy["Charge"])])
                    else:
                        all_psms.append( [False, decoy_length, sid, decoy_db_ind,
                                        float(bdecoy["Score"])])


            except ValueError:
                print 'Record %d in %s is bad' % (reader.line_num, f.name)
                print 'Headers in bad ident file: %s' % reader.fieldnames
                raise

    f.close()

    # sort all peptides by length for partition rank normalization
    pep_length = lambda r: r[1]
    all_psms.sort(key = pep_length)
    return all_psms
