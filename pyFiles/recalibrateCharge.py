#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

import itertools
import math
import numpy
# import scipy.interpolate
####### how to use linear interpolation:
####### y_interp = scipy.interpolate.interp1d(x,y)
####### new_range = range(10)
####### y_interp(new_range)

def write_drip_recal(output_file, targets, decoys):
    print "Writing new ident file %s" % (output_file)
    outfile = open(output_file, "w")
    outfile.write('Kind\tScan\tScore\tPeptide\tObs_Inserts\tTheo_Deletes\tObs_peaks_scored\tTheo_peaks_used\tSum_obs_intensities\tSum_scored_mz_dist\tCharge\tFlanking_nterm\tFlanking_cterm\n')
    keys = ['Kind', 'Scan', 'Score', 'Peptide', 'Obs_Inserts', 
            'Theo_Deletes', 'Obs_peaks_scored', 'Theo_peaks_used', 
            'Sum_obs_intensities', 'Sum_scored_mz_dist', 'Charge', 
            'Flanking_nterm', 'Flanking_cterm']
    for t,d in zip(targets,decoys):
        for k in keys:
            if k=='Kind':
                outfile.write('t\t')
            elif k=='Scan':
                outfile.write('%d\t' % t.scan)
            elif k=='Score':
                outfile.write('%f\t' % t.score)
            elif k=='Peptide':
                outfile.write('%s\t' % t.peptide)
            elif k=='Charge':
                outfile.write('%d\t' % t.charge)
            else:
                outfile.write('%s\t' % t.other[k])
        outfile.write('\n')

        for k in keys:
            if k=='Kind':
                outfile.write('d\t')
            elif k=='Scan':
                outfile.write('%d\t' % d.scan)
            elif k=='Score':
                outfile.write('%f\t' % d.score)
            elif k=='Peptide':
                outfile.write('%s\t' % d.peptide)
            elif k=='Charge':
                outfile.write('%d\t' % d.charge)
            else:
                outfile.write('%s\t' % d.other[k])
        outfile.write('\n')

    outfile.close()

def merge_idents(t,d,charges,topPsms = 1):
    # t and d are lists of list, where each list, say t[i] and d[i], correspond
    # to charge run charges[i].  for each x in t[i]: 
    # x[0] = sid
    # x[1] = peptide index
    # x[2] = peptide score
    if "Sid" in l:
        sidKey = "Sid"
    elif "Scan" in l:
        sidKey = "Scan"

    identSids = []
    # dictionaries of target and decoy PSMs, where the keys are spectrum ids
    targetPsms = {}
    decoyPsms = {}
    for c, tb,db in zip(charges, t,d):
        for tbb, dbb in zip(tb,db):
            # first check target PSMs
            if tbb[0] not in targetPsms:
                targetPsms[tbb[0]] = (tbb[1],tbb[2], c)
            else: # take higher scoring PSM
                if targetPsms[tbb[0]][1] < tbb[2]:
                    targetPsms[tbb[0]] = (tbb[1],tbb[2], c)
            # now decoys
            if dbb[0] not in decoyPsms:
                decoyPsms[dbb[0]] = (dbb[1],dbb[2], c)
            else: # take higher scoring PSM
                if decoyPsms[dbb[0]][1] < dbb[2]:
                    decoyPsms[tbb[0]] = (dbb[1],dbb[2], c)

    assert len(targetPsms)==len(decoyPsms)
    return targetPsms, decoyPsms

def returnTopPsms(psms, topPsms = 1):
    """
    psms is a tuple of PSM objects
    """
    targets = []
    decoys = []
    # resort psms by sid, take max psm per sid
    isTarget = lambda r: r.kind == 't'
    isDecoy = lambda r: r.kind == 'd'
    sidField = lambda r: r.scan
    scoref = lambda r: r.score
    # iterate through all PSMs per spectrum
    psms.sort(key = sidField)
    # Find max PSM per spectrum
    for sid, all_rows in itertools.groupby(psms, sidField):
        # print "Find max psm for sid=%d" % sid
        curr_psms = list(all_rows)
        # top_decoy = max(itertools.ifilter(isDecoy, curr_psms), key = scoref)
        # targets.append( ( sid, top_target[3], top_target[4] ) )
        # decoys.append( ( sid, top_decoy[3], top_decoy[4] ) )

        curr_psms.sort(key = scoref, reverse = True)
        for i,t in enumerate(itertools.ifilter(isTarget, curr_psms)):
            targets.append(t)
            if i >= ( topPsms - 1 ):
                break
        for i,d in enumerate(itertools.ifilter(isDecoy, curr_psms)):
            decoys.append(d)
            if i >= ( topPsms - 1 ):
                break

    return targets, decoys

def renorm_partition(partition, dist, percentile = 0.99):
    """ Fit scores of PSMs in partition to the distribution dist using linear interpolation
    Inputs:
        dist = set of two lists such that: 
               dist[0] = x-values of renormalized distribution, dist[1]=y-values
        partition = list of all target and decoy scores, a single
        equivalence class under the predefined equivalance relation(i.e. 
        peptide length, peptide charge, etc.).  Currently, the only
        supported/assumed equivalence relation is the length of the peptide.
        
    Post-conditions:
        Scores of PSMs in partition are fit to the distribution in dist

        fields: 0) isTarget, boolean
                1) peptide length, integer
                2) spectrum id, integer
                3) peptide index number, integer
                4) Score, float
    """
    ################# old
    # fields: 0) isTarget, boolean
    #         1) peptide length, integer
    #         2) spectrum id, integer
    #         3) peptide index number, integer
    #         4) Score, float
    ################# new code base: partition is a tuple of PSM objects

    partition.sort(key = lambda r: r.score)
    # loop through both the distribution and the partition
    part_ind = 0
    dist_ind = 0
    dist_x = dist[0]
    dist_y = dist[1]
    num_dist_els = len(dist_x)
    num_partition_els = len(partition)
    last_dist_ind = 0
    # check whether we must interpolate from the left tail first
    if partition[part_ind].score < dist_x[dist_ind]:
        # print "(%d,%d), (%d,%d), (%f,%f)" % (part_ind, num_partition_els, dist_ind, num_dist_els, partition[part_ind][scoreIdx], dist_x[dist_ind])
        dist_ind = int(math.floor((1.0-percentile)*num_dist_els))
        # make sure y intercepts are not equal
        while dist_x[dist_ind] == dist_x[last_dist_ind]:
            dist_ind += 1

        ub = dist_x[dist_ind]
        lb = dist_x[last_dist_ind]
        ub_y = dist_y[dist_ind]
        lb_y = dist_y[last_dist_ind]
        m = (ub_y-lb_y)/(ub-lb)
        while part_ind < num_partition_els and partition[part_ind].score <= ub:
            partition[part_ind].score = m*(partition[part_ind].score-lb)+lb_y
            part_ind += 1

        last_dist_ind = dist_ind

    # # reset last_dist_ind
    # last_dist_ind = 0
    while part_ind < num_partition_els:
        # print "(%d,%d), (%d,%d), (%f,%f)" % (part_ind, num_partition_els, dist_ind, num_dist_els, partition[part_ind][scoreIdx], dist_x[dist_ind])
        curr_score = partition[part_ind].score
        # bracket the current score
        while dist_ind < (num_dist_els-1) and dist_x[dist_ind] <= curr_score:
            dist_ind += 1
        # now we know current scoring interval it is [dist[last_dist_ind], dist[dist_ind]]
        # calculate slope
        ub = dist_x[dist_ind]
        lb = dist_x[last_dist_ind]
        ub_y = dist_y[dist_ind]
        lb_y = dist_y[last_dist_ind]
        # print "(%d,%d), (%d,%d), (%d,%d), (%f,%f)>%f" % (part_ind, num_partition_els, dist_ind, num_dist_els, last_dist_ind,dist_ind,ub, lb, curr_score)
        m = (ub_y-lb_y)/(ub-lb)
        # keep track of last dist_ind
        last_dist_ind = dist_ind
        # at the edge of the distribution, 
        # finish reranking all remaining partition elements
        if dist_ind == (num_dist_els-1):
            # linearly interpolate using the 99th percentile score
            last_dist_ind = int(math.floor(percentile*num_dist_els))
            ub = dist_x[dist_ind]
            lb = dist_x[last_dist_ind]
            ub_y = 1.0
            lb_y = dist_y[last_dist_ind]
            if ub_y==lb_y:
                m = 0.0
            else:
                m = (ub_y-lb_y)/(ub-lb)
            # m = (ub_y-lb_y)/(ub-lb)
            # print "ub_y=%f,lb_y=%f,ub=%f,lb=%f" % (ub_y,lb_y,ub,lb)
            # m = math.exp(math.log(ub_y-lb_y) - math.log(ub-lb))
            while part_ind < num_partition_els:
                partition[part_ind].score = m*(partition[part_ind].score-lb)+lb_y
                part_ind += 1
        else: # still moving down interpolated distribution         
            while part_ind < num_partition_els and partition[part_ind].score <= ub:
                partition[part_ind].score = m*(partition[part_ind].score-lb)+lb_y
                part_ind += 1

def standardNorm_partition(partition, dist):
    """ Fit scores of PSMs in partition to the distribution dist using linear interpolation
    Inputs:
        dist = set of two lists such that: 
               dist[0] = x-values of renormalized distribution, dist[1]=y-values
        partition = list of all target and decoy scores, a single
        equivalence class under the predefined equivalance relation(i.e. 
        peptide length, peptide charge, etc.).  Currently, the only
        supported/assumed equivalence relation is the length of the peptide.
        
    Post-conditions:
        Scores of PSMs in partition are fit to the distribution in dist

        fields: 0) isTarget, boolean
                1) peptide length, integer
                2) spectrum id, integer
                3) peptide index number, integer
                4) Score, float
    """
    ################# old
    # fields: 0) isTarget, boolean
    #         1) peptide length, integer
    #         2) spectrum id, integer
    #         3) peptide index number, integer
    #         4) Score, float
    ################# new code base: partition is a tuple of PSM objects
    m = numpy.mean(dist)
    v = numpy.std(dist)

    for ind in range(len(partition)):
        partition[ind].score = (partition[ind].score - m) / v

def decoy_charge_partition_dist(ch_scores):
    """ Calculate rank normalization distributions per charge run (partition) to fit
        scores to
    """
    ch_scores.sort()
    return [[x for x in ch_scores],[float(ind+1)/float(len(ch_scores)) for ind in range(len(ch_scores))]]

def recalibrate_charge_psms(output_file, psms_by_charge, 
                            training_decoys_by_charge,
                            percentile,
                            topPsms = 1):
    """Load an all psms file(for a dataset), calculate the rank normalization distribution per partition
       function(currently, the only supported partition function is peptide length)

    Arguments:
        filename: Name of the tab-separated file, with Kind, Sid, Peptide,
            and Score fields.
        all_psms_file: training data, where decoys are independent of testing target and decoy sets
        output_file: output ident file
    """
    
    # calculate distributions per partition(charge)
    psms = []
    for c in training_decoys_by_charge:
        standardNorm_partition(psms_by_charge[c], training_decoys_by_charge[c])
        # renorm_partition(psms_by_charge[c], 
        #                  decoy_charge_partition_dist(training_decoys_by_charge[c]), 
        #                  percentile)
        psms += psms_by_charge[c]

    tTop, dTop = returnTopPsms(psms, topPsms)
    write_drip_recal(output_file, tTop, dTop)
    # # also, write ch2 and ch3 psms to file
    # for c,tcurr,dcurr in zip(charges,t,d):
    #     ch_outfile = '%s-chargeRecalibrated-ch%d-ident.txt' % (basefile, c)
    #     psmUtils.write_merged_ident_list(ch_outfile, tcurr, dcurr, [c]*len(tcurr), [c]*len(tcurr), test_target_db, test_decoy_db)
