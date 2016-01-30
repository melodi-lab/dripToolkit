#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0
# Command line parsing utilities.

from __future__ import with_statement

__authors__ = ['John Halloran <halloj3@uw.edu>' ]

import re
import itertools

# Regular expressions to parse gmtkKernel output
################## regular expression, compile down for efficiency
lineVals=re.compile('\S+')
seg_pattern=re.compile('Segment (?P<segment>\d+) : (?P<frame>\d+) frames, (?P<uframe>\d+) usable frames, log\(PE\) = (?P<score>\S+) ', re.MULTILINE)
# num_means_pattern=re.compile('% (?P<numMeans>\d+) means \(there may be fewer due to -objsNotToUTilize\)')
# num_covars_pattern=re.compile('% (?P<numCovars>\d+) covars \(there may be fewer due to -objsNotToUTilize\)')
# num_components_pattern=re.compile('% (?P<numComponents>\d+) components \(there may be fewer due to -objsNotToUTilize\)')
# num_mdcpts_pattern=re.compile('% (?P<numMDCPTs>\d+) mdCpts \(there may be fewer due to -objsNotToUTilize\)')
component_pattern=re.compile('% (?P<component>\S+) nextMeans\.len (?P<meanDim>\d+)  nextDiagCovars\.len (?P<covarDim>\d+)')
####################### updated pattern for latest gmtk trunk, pulled 2014-08-26, foun din /homes/halloj3/gmtk_WORK_20140826
# dpmf_pattern=re.compile('% Dense1DPMF (?P<dpmfName>\S+):.*\.\.\. nextPmf\[i\] \.\.\..*(?P<dpmfVals>(\S+ )*)')
dpmf_pattern=re.compile('% Dense1DPMF (?P<dpmfName>\S+):  \.\.\. nextPmf\[i\] \.\.\. (?P<dpmfVals>(\S+ )*)')
mean_vals_pattern=re.compile('% MeanVector (?P<meanName>\S+):.*\.\.\. nextMeans\[i\] \.\.\..*\n(?P<meanVals>(\S+ )*)', re.MULTILINE)
covar_vals_pattern=re.compile('% DiagCovarVector (?P<covarName>\S+):.*\.\.\. nextDiagCovars\[i\] \.\.\..*\n(?P<covarVals>(\S+ )*)', re.MULTILINE)
cpt_vals_pattern=re.compile('% MDCPT (?P<cptName>\S+):.*\.\.\. nextMdcpt\[i\] \.\.\.\n(?P<cptVals>(\S+ )*)', re.MULTILINE)


############## gmtkMean class
class gmtkMean(object):
    """ GMTK mean class definition
    """

    def __init__(self, name = '', dimension = 1, 
                 vals = [], index = 0):

        assert dimension > 0, "Mean %s with dimension $d, must be > 0.  Exitting" % (name, dimension)
        assert vals and len(vals)==dimension, "Mean %s has dimension %d but %d values supplied, exitting." % (name, dimension, len(vals))

        self.name = name
        self.dimension = dimension
        self.vals = vals
        self.index = index

    def __cmp__(self,other):
        if self.name == other.name:
            if len(self.vals) == len(other.vals):
                if sum([abs(i-j) for i,j in zip(self.vals, other.vals)]) == 0.0:
                    return 1
        return 0

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        k = "%sDim%d:" % (self.name, self.dimension)
        for v in vals:
            k = k + ' ' + str(v)
        return k

def parse_kernel_output(kernel_output):
    # Parse the output of gmtkKernel
    try:
        f = open(kernel_output, "r")
    except IOError:
        print "Could not open file %s for reading, exitting" % kernel_output
        exit(-1)

    log = f.read()
    f.close()
             
    cpts = []
    means = []
    covars = []
    components = []
    dpmfs = []
    for match, subSection in itertools.izip(seg_pattern.finditer(log), 
                                            re.split('Segment', log)[1:]):
        segment = int(match.group('segment'))
        curr_frames = int(match.group('frame'))
        curr_used_frames = int(match.group('uframe'))
        probEvidence = float(match.group('score'))
        # for now, assume there is only one set of parameters, i.e. only one segment is run
        for cpt in cpt_vals_pattern.finditer(subSection):
            cptVals = []
            for cptVal in cpt.group('cptVals').split():
                cptVals.append(float(cptVal))
            cptName = cpt.group('cptName')
            # assert len(cptVals)==int(cpt.group('cptDim')), "Num parsed cpt %s values not equal to cpt dimensionality" % (cptName)
            cpts.append([cptName, cptVals])
        for component in component_pattern.finditer(subSection):
            meanCard = int(component.group('meanDim'))
            covarCard = int(component.group('covarDim'))
            components.append([component.group('component'), meanCard, covarCard])
        for num, mean in enumerate(mean_vals_pattern.finditer(subSection)):
            # assert len(mean.groups())==2, "Parsed means should have 2 groups"
            meanName = mean.group('meanName')
            meanVals = []
            for meanVal in mean.group('meanVals').split():
                meanVals.append(float(meanVal))
            assert len(meanVals)==components[num][1], "Num parsed mean %s values not equal to mean dimensionality" % (meanName)
            means.append([meanName, meanVals])
        for num, covar in enumerate(covar_vals_pattern.finditer(subSection)):
            # assert len(covar.groups())==2, "Parsed covariances should have 2 groups"
            covarName = covar.group('covarName')
            covarVals = []
            for covarVal in covar.group('covarVals').split():
                covarVals.append(float(covarVal))
            assert len(covarVals)==components[num][1], "Num parsed covar %s values not equal to covar dimensionality" % (covarName)
            covars.append([covarName, covarVals])
        for num, dpmf in enumerate(dpmf_pattern.finditer(subSection)):
            dpmfName = dpmf.group('dpmfName')
            dpmfVals = []
            for dpmfVal in dpmf.group('dpmfVals').split():
                dpmfVals.append(float(dpmfVal))
            # assert len(dpmfVals)==components[num][1], "Num parsed dpmf %s values not equal to dpmf dimensionality" % (dpmfName)
            dpmfs.append([dpmfName, dpmfVals])

    return means, covars, cpts, probEvidence, dpmfs

def chomp(allLines, totalLinesRead):
    # Eat all empty lines and comments, return parsed contents of first nonempty line
    # Also, discard any comments
    totalLines = len(allLines)
    l = allLines[totalLinesRead]
    split_line = lineVals.findall(l)
    lines_read = 1
    totalLinesRead += 1
    # Read all blank lines
    while(not split_line and totalLinesRead < totalLines):
        l = allLines[totalLinesRead]
        lines_read += 1
        totalLinesRead += 1
        split_line = lineVals.findall(l)
        if split_line: # skip any comments
            if split_line[0] == '%':
                split_line = []

    # search for comment in split_line
    if split_line:
        del_els = -1
        for i, s in enumerate(split_line):
            if s == '%':
                del_els = i
                break
        if del_els > -1:
            del split_line[del_els:]

    return split_line, lines_read

def parse_elements(allLines, totalLinesRead):
    # Function to read in dpmfs, means, and covariances
    totalLines = len(allLines)
    el_vector = []
    currLinesRead = 0
    num_els_line, lines_read = chomp(allLines, totalLinesRead)
    totalLinesRead += lines_read
    if len(num_els_line) > 1:
        print("Error: line %d should only contain number of proceeding elements" % totalLinesRead)
        exit(-1)
    else:
        num_els = int(num_els_line[0])
    for i in range(num_els): # read in all dpmfs
        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
        totalLinesRead += lines_read
        # el syntax: el_num[[:space:]]name[[:space:]]dimensionality[[:space:]]value(s)[[:space:]]
        el = []
        curr_ind = -1
        el_name = []
        el_dim = -1
        el_vals = []
        run = 1
        if curr_parsed: # make sure nonempty
            while run:
                for parsed in curr_parsed:
                    if curr_ind==-1:
                        try: # make it's a valid integer
                            curr_ind = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int element ind" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif len(el_name)==0:
                        el_name = parsed
                        # append to current element
                        el.append(el_name)
                    elif el_dim==-1:
                        try: # make sure it's a valid integer
                            el_dim = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int dimension" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                        # append to current element
                        el.append(el_dim)
                    else: # these must be the values
                        try: # make sure it's a valid float
                            curr_el = float(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid float" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                        el_vals.append(curr_el)
                # we must have at least read the current element index
                # decide whether we should continue reading lines
                if not el_name or el_dim < 0 or not el_vals:
                    curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                    totalLinesRead += lines_read
                else:
                    if len(el_vals) > el_dim:
                        print "Element %s with dimension %d, %d values supplied.  Exitting" % (el_name, el_dim, len(el_vals))
                        exit(-1)
                    elif len(el_vals) < el_dim:
                        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                        totalLinesRead += lines_read
                    else: # we're finished reading this element
                        run = 0                    
                        el.append(el_vals)
        el_vector.append(el)

    return el_vector, totalLinesRead

def parse_dense_cpts(allLines, totalLinesRead):
    # Function to read dense cpts, which have a slightly different syntax from dpmfs, 
    # 3 
    # 0 
    # stateTransitionProbs 
    # 1 % number parents
    # 150 3 % cardinalities
    # 8.65986587535374275e-01 1.28407679516182732e-01 5.60573294844298690e-03 
    # 8.74558545443081425e-01 6.83680641931374705e-02 5.70733903637811388e-02 
    # 8.50054112623554237e-01 9.99639249175895694e-02 4.99819624588561801e-02 
    # 6.35770836036747333e-01 3.62948543106983879e-01 1.28062085626879404e-03 
    # 8.42096460399923585e-01 1.54741669404714177e-01 3.16187019536216200e-03 
    # 6.50469963770039961e-01 2.33020024153178851e-01 1.16510012076781203e-01 
    # 4.08982325991748552e-01 5.86900408980366439e-01 4.11726502788500870e-03 
    # ...
    # output:
    # cpt[0] - name, cpt[1] - dimension(s), cpt[2] - flattened cpt
    el_vector = []
    num_els_line, lines_read = chomp(allLines, totalLinesRead)
    totalLinesRead += lines_read
    if len(num_els_line) > 1:
        print("Error: line %d should only contain number of proceeding elements" % totalLinesRead)
        exit(-1)
    else:
        num_els = int(num_els_line[0])
    for i in range(num_els): # read in all dpmfs
        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
        totalLinesRead += lines_read
        # el syntax: el_num[[:space:]]name[[:space:]]num_parents[[:space:]]dimensionalities[[:space:]]value(s)[[:space:]]
        el = []
        curr_ind = -1
        el_name = []
        el_num_parents = -1
        el_dim = []
        el_vals = []
        total_dim = 1
        run = 1
        if curr_parsed: # make sure nonempty
            while run:
                for parsed in curr_parsed:
                    if parsed == 'DirichletConst':
                        break
                    if curr_ind==-1:
                        try: # make sure it's a valid integer
                            curr_ind = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int element ind" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif len(el_name)==0:
                        el_name = parsed
                    elif el_num_parents==-1:
                        try: # make sure it's a valid integer
                            el_num_parents = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid number of parents" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif len(el_dim) < (el_num_parents+1):
                        try: # make sure it's a valid integer
                            curr_el_dim = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int dimension" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                        # append to current element
                        el_dim.append(curr_el_dim)
                        total_dim *= curr_el_dim
                    else: # these must be the values
                        try: # make sure it's a valid float
                            curr_el = float(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid float" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                        el_vals.append(curr_el)
                # we must have at least read the current element index
                # decide whether we should continue reading lines
                if not el_name or el_num_parents < 0 or len(el_dim) < (el_num_parents+1) or not el_vals:
                    curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                    totalLinesRead += lines_read
                else:
                    if len(el_vals) > total_dim:
                        print "Element %s with product of parents cardinalities %d, %d values supplied.  Exitting" % (el_name, total_dim, len(el_vals))
                        exit(-1)
                    elif len(el_vals) < total_dim: # keep reading
                        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                        totalLinesRead += lines_read
                    else: # we're finished reading this element
                        run = 0          
                        # append elements
                        el.append(el_name)
                        el.append(el_dim)
                        el.append(el_vals)
        el_vector.append(el)

    return el_vector, totalLinesRead

def parse_gaussian_components(allLines, totalLinesRead):
    # Function to read Gaussian components
    # example:
    # 150 
    # 0 
    # 39 
    # 0 gc0 
    # mean0 var0
    # output:
    # gaussian_component[0]: name, gaussian_component[1]: dimension, 
    # gaussian_component[2]: type
    # gaussian_component[3]: mean
    # gaussian_component[4]: covariance
    el_vector = []
    num_els_line, lines_read = chomp(allLines, totalLinesRead)
    totalLinesRead += lines_read
    if len(num_els_line) > 1:
        print("Error: line %d should only contain number of proceeding elements" % totalLinesRead)
        exit(-1)
    else:
        num_els = int(num_els_line[0])
    for i in range(num_els): # read in all Gaussian Components
        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
        totalLinesRead += lines_read
        # el syntax: el_num[[:space:]]name[[:space:]]dimensionality[[:space:]]value(s)[[:space:]]
        el = []
        curr_ind = -1
        el_name = []
        el_dim = -1
        el_type = -1
        el_mean = []
        el_covar = []
        run = 1
        if curr_parsed: # make sure nonempty
            while run:
                for parsed in curr_parsed:
                    if curr_ind==-1:
                        try: # make it's a valid integer
                            curr_ind = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int element ind" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif el_dim==-1:
                        try: # make sure it's a valid integer
                            el_dim = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int dimension" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif el_type == -1:
                        try: # make sure it's a valid integer
                            el_type = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int component type" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif len(el_name)==0:
                        el_name = parsed
                    elif len(el_mean)==0: # this must be mean 
                        el_mean = parsed
                    else: # this must be the covariance
                        el_covar = parsed
                # we must have at least read the current element index
                # decide whether we should continue reading lines
                if not el_name or not el_dim or not el_mean or not el_covar:
                    curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                    totalLinesRead += lines_read
                else:
                    run = 0  
                    el.append(el_name)
                    el.append(el_dim)
                    el.append(el_type)
                    el.append(el_mean)
                    el.append(el_covar)
        el_vector.append(el)
    return el_vector, totalLinesRead

def parse_gaussian_mixtures(allLines, totalLinesRead):
    # Function to read Gaussian components
    # example:
    # 150 
    # 0 
    # 39 
    # gm0 
    # 64 dpmf0 
    # gc0 gc0_cl5 gc0_cl4 gc0_cl4_cl0 gc0_cl3 gc0_cl3_cl1 gc0_cl3_cl0 gc0_cl3_cl0_cl0 gc0_cl2 gc0_cl2_cl2 gc0_cl2_cl1 gc0_cl2_cl1_cl0 gc0_cl2_cl0 gc0_cl2_cl0_cl1 gc0_cl2_cl0_cl0 gc0_cl2_cl0_cl0_cl0 gc0_cl1 gc0_cl1_cl3 gc0_cl1_cl2 gc0_cl1_cl2_cl0 gc0_cl1_cl1 gc0_cl1_cl1_cl1 gc0_cl1_cl1_cl0 gc0_cl1_cl1_cl0_cl0 gc0_cl1_cl0 gc0_cl1_cl0_cl2 gc0_cl1_cl0_cl1 gc0_cl1_cl0_cl1_cl0 gc0_cl1_cl0_cl0 gc0_cl1_cl0_cl0_cl1 gc0_cl1_cl0_cl0_cl0 gc0_cl1_cl0_cl0_cl0_cl0 gc0_cl0 gc0_cl0_cl4 gc0_cl0_cl3 gc0_cl0_cl3_cl0 gc0_cl0_cl2 gc0_cl0_cl2_cl1 gc0_cl0_cl2_cl0 gc0_cl0_cl2_cl0_cl0 gc0_cl0_cl1 gc0_cl0_cl1_cl2 gc0_cl0_cl1_cl1 gc0_cl0_cl1_cl1_cl0 gc0_cl0_cl1_cl0 gc0_cl0_cl1_cl0_cl1 gc0_cl0_cl1_cl0_cl0 gc0_cl0_cl1_cl0_cl0_cl0 gc0_cl0_cl0 gc0_cl0_cl0_cl3 gc0_cl0_cl0_cl2 gc0_cl0_cl0_cl2_cl0 gc0_cl0_cl0_cl1 gc0_cl0_cl0_cl1_cl1 gc0_cl0_cl0_cl1_cl0 gc0_cl0_cl0_cl1_cl0_cl0 gc0_cl0_cl0_cl0 gc0_cl0_cl0_cl0_cl2 gc0_cl0_cl0_cl0_cl1 gc0_cl0_cl0_cl0_cl1_cl0 gc0_cl0_cl0_cl0_cl0 gc0_cl0_cl0_cl0_cl0_cl1 gc0_cl0_cl0_cl0_cl0_cl0 gc0_cl0_cl0_cl0_cl0_cl0_cl0 
    # output:
    # gaussian_component parameters: 
    # gaussian_mixture[0]: name, gaussian_mixture[1]: mixture dim, gaussian_mixture[2]: num components, gaussian_mixture[3]: dpmf name, gaussian_mixture[4]: components
    el_vector = []
    num_els_line, lines_read = chomp(allLines, totalLinesRead)
    totalLinesRead += lines_read
    if len(num_els_line) > 1:
        print("Error: line %d should only contain number of proceeding elements" % totalLinesRead)
        exit(-1)
    else:
        num_els = int(num_els_line[0])
    for i in range(num_els): # read in all dpmfs
        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
        totalLinesRead += lines_read
        # el syntax: el_num[[:space:]]name[[:space:]]dimensionality[[:space:]]dpmf[[:space:]]value(s)[[:space:]]
        el = []
        curr_ind = -1
        el_name = []
        el_dim = -1
        el_dpmf = []
        el_num_comps = -1
        el_vals = []
        run = 1
        if curr_parsed: # make sure nonempty
            while run:
                for parsed in curr_parsed:
                    if curr_ind==-1:
                        try: # make it's a valid integer
                            curr_ind = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int element ind" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif el_dim==-1:
                        try: # make sure it's a valid integer
                            el_dim = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid int dimension" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif len(el_name)==0:
                        el_name = parsed
                        # append to current element
                    elif el_num_comps==-1:
                        try: # make sure it's a valid integer
                            el_num_comps = int(parsed)
                        except ValueError:
                            print("Error near line %d: %s not valid number of components" % 
                                  (totalLinesRead,parsed))
                            exit(-1)
                    elif len(el_dpmf)==0:
                        el_dpmf = parsed
                    else: # these must be the mixture components
                        el_vals.append(parsed)
                # we must have at least read the current element index
                # decide whether we should continue reading lines
                if el_dim < 0 or not el_name or el_num_comps < 0 or not el_dpmf or not el_vals:
                    curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                    totalLinesRead += lines_read
                else:
                    if len(el_vals) > el_num_comps:
                        print "Element %s with dimension %d, %d values supplied.  Exitting" % (el_name, el_dim, len(el_vals))
                        exit(-1)
                    elif len(el_vals) < el_num_comps:
                        curr_parsed, lines_read = chomp(allLines, totalLinesRead)
                        totalLinesRead += lines_read
                    else: # we're finished reading this element
                        run = 0                    
                        el.append(el_name)
                        el.append(el_dim)
                        el.append(el_num_comps)
                        el.append(el_dpmf)
                        el.append(el_vals)
        el_vector.append(el)

    return el_vector, totalLinesRead

def write_params(filename,mixtures,
                 components,dpmfs,dense_cpts,means,covars):
    # write params file
    output = open(filename, "w")
    # write dpmfs
    output.write("\n% dense PMFs\n")
    output.write("%d\n" % len(dpmfs))
    for i, dpmf in enumerate(dpmfs):
        output.write("%d\n" % (i))
        # dpmf parameters: 
        # dpmf[0]: name, dpmf[1]: dimension, dpmf[2]: value
        output.write("%s %d" % (dpmf[0],dpmf[1]))
        if len(dpmf[2])!=dpmf[1]:
            print("Error: dpmf %s has less elements than dimenions" % dpmf[0])
            output.close()
            exit(-1)
        for el in dpmf[2]:
            output.write(" %.16e" % (el))
        output.write("\n")

    print >> output, """
% sparse PMFs
0
"""
    # write means
    output.write("% means\n")
    output.write("%d\n" % len(means))
    for i, mean in enumerate(means):
        output.write("%d\n" % (i))
        # mean parameters: 
        # mean[0]: name, mean[1]: dimension, mean[2]: value
        output.write("%s %d" % (mean[0],mean[1]))
        if len(mean[2])!=mean[1]:
            print("Error: mean %s has less elements than dimenions" % mean[0])
            output.close()
            exit(-1)
        for el in mean[2]:
            output.write(" %.16e" % (el))
        output.write("\n")
    # write covars
    output.write("\n% diagonal covariance matrices\n")
    output.write("%d\n" % len(covars))
    for i, covar in enumerate(covars):
        output.write("%d\n" % (i))
        # covar parameters: 
        # covar[0]: name, covar[1]: dimension, covar[2]: value
        output.write("%s %d" % (covar[0],covar[1]))
        if len(covar[2])!=covar[1]:
            print("Error: covariance %s has less elements than dimenions" % covar[0])
            output.close()
            exit(-1)
        for el in covar[2]:
            output.write(" %.16e" % (el))
        output.write("\n")
        
    print >> output, """
% dlink matrices
0

% real valued N_i x M_i matrices
0

% Dense CPTs
"""
    output.write("%d\n" % len(dense_cpts))
    for i,cpt in enumerate(dense_cpts):
        # cpt syntax:
        # cpt[0] - name, cpt[1] - dimension(s), cpt[2] - flattened cpt
        output.write("%d\n%s\n%d\n" % (i,cpt[0], len(cpt[1])-1))
        for dim in cpt[1]:
            output.write("%d " % dim)
        output.write("\n")
        for val in cpt[2]:
            output.write("%.16e " % val)
        output.write("\n")

    # write components
    output.write("\n% Components\n")
    output.write("%d\n" % len(components))
    for i, gaussian_component in enumerate(components):
        output.write("%d\n" % (i))
        # gaussian_component parameters: 
        # gaussian_component[0]: name, gaussian_component[1]: dimension, 
        # gaussian_component[2]: type
        # gaussian_component[3]: mean
        # gaussian_component[4]: covariance
        # assume scalar gaussian for now
        output.write("%d\n%d " % (gaussian_component[1],gaussian_component[2])) # dimension and gaussian type
        output.write("%s\n%s %s\n" % (gaussian_component[0],gaussian_component[3],gaussian_component[4]))

    # write mixtures
    output.write("\n% Mixtures of components\n")
    output.write("%d\n" % len(mixtures))
    for i, gaussian_mixture in enumerate(mixtures):
        output.write("%d\n" % (i))
        # gaussian_mixture parameters: 
        # gaussian_mixture[0]: name, gaussian_mixture[1]: mixture dim, gaussian_mixture[2]: num components, gaussian_mixture[3]: dpmf name, gaussian_mixture[4]: components
        # assume scalar gaussian for now
        output.write("%d\n" % gaussian_mixture[1]) # dimension
        output.write("%s\n%d %s\n" % (gaussian_mixture[0],gaussian_mixture[2],gaussian_mixture[3]))
        for comp in gaussian_mixture[4]:
            output.write("%s " % comp)
        output.write("\n")

    print >> output, """
% Gaussian Switching Mixtures of Gaussians
0

% Logistic-Regression-based Switching Mixutres of Gaussians
0

% MLP-based Switching Mixtures of Gaussians
0

"""

def parse_params(params_file):
    totalLinesRead = 0
    # open params file
    params = open(params_file, "r")
    log = params.readlines()
    params.close()
    ##### comment lines to anchor where datatypes are
    ### dense PMFs
    dpmf_l0 = '% dense PMFs\n'
    dpmf_l1 = '% dense probability mass functions\n'
    ### sparse PMFs
    spmf_l0 = '% sparse PMFs\n'
    spmf_l1 = '% sparse probability mass functions\n'
    ### means
    means_l0 = '% means\n'
    ### covariances
    covars_l0 = '% diagonal covariance matrices\n'
    covars_l1 = '% vars\n'
    ### dlink matrices
    dlink_l0 = '% dlink matrices\n'
    ### weight matrices
    wmat_l0 = '% weight matrices\n'
    wmat_l1 = '% real valued N_i x M_i matrices\n'
    ### cpts
    dcpt_l0 = '% Dense CPTs\n'
    dcpt_l1 = '% dense conditional probability tables\n'
    ### gaussian components
    gc_l0 = '% Components\n'
    gc_l1 = '% Gaussian components\n'
    ### gaussian mixtures
    gm_l0 = '% Mixtures of components\n'
    gm_l1 = '% Gaussian components\n'
    while totalLinesRead < len(log): # read line by line
        l = log[totalLinesRead]
        totalLinesRead+=1
        # read means
        if l == means_l0:
            means, totalLinesRead = parse_elements(log, totalLinesRead)
        # read dpmfs
        if l == dpmf_l0 or l == dpmf_l1:
            dpmfs,totalLinesRead = parse_elements(log, totalLinesRead)
        # read Gaussian Covariances
        if l == covars_l0 or l == covars_l1:
            covars,totalLinesRead = parse_elements(log, totalLinesRead)
        # read gaussian components
        if l == gc_l0 or l == gc_l1:
            gaussian_components,totalLinesRead = parse_gaussian_components(log, totalLinesRead)
        # read gaussian mixtures
        if l == gm_l0 or l == gm_l1:
            gaussian_mixtures,totalLinesRead = parse_gaussian_mixtures(log, totalLinesRead)
        # read Dense CPTs
        if l == dcpt_l0 or l == dcpt_l1:
            dense_cpts,totalLinesRead = parse_dense_cpts(log, totalLinesRead)
    return gaussian_mixtures, gaussian_components, dpmfs, means, covars, dense_cpts
