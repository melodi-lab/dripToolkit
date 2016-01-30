#!/usr/bin/env python
#
# Written by John Halloran <halloj3@uw.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0

"""Parsers needed to generate inputs for visualizations:

- Merging shards into a single .out file: merge_shards
- Loading identification files: load_ident

"""

import csv
import glob
import itertools
import os
# import progress

from peptide import Peptide

def load_ident(filename, sids = None, charge = False, by_matched = False):
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
        2. all(t[0] == d[0] for t, d in zip(targets, decoy). Each spectrum has
           exactly one record in targets, and one in decoys. The entries in
           targets and decoys are sorted by spectrum id.

    """
    targets = [ ]
    decoys = [ ]

    f = open(filename)
    # reader = csv.DictReader(f, delimiter = '\t', skipinitialspace = True)
    reader = [l for l in csv.DictReader(f, delimiter = '\t', skipinitialspace = True)]

    reader.sort(key = lambda r: int(r["Sid"]))

    targetf = lambda r: r["Kind"] == "t"
    decoyf = lambda r: r["Kind"] == "d"
    scoref = lambda r: float(r["Score"])
    for sid, rows in itertools.groupby(reader, lambda r: int(r["Sid"])):
        if not sids or sid in sids:
            records = list(rows)
            try:
                btarget = max(itertools.ifilter(targetf, records), key = scoref)
                bdecoy = max(itertools.ifilter(decoyf, records), key = scoref)
                if by_matched:
                    if charge:
                        targets.append( (sid, Peptide(btarget["Peptide"]),
                                         int(btarget["Score"]), int(btarget["Charge"]), int(btarget["NumBY"])))
                        decoys.append( (sid, Peptide(bdecoy["Peptide"]),
                                        int(bdecoy["Score"]), int(bdecoy["Charge"]), int(bdecoy["NumBY"])))
                    else:
                        targets.append( (sid, Peptide(btarget["Peptide"]),
                                         int(btarget["Score"]), int(btarget["NumBY"])))
                        decoys.append( (sid, Peptide(bdecoy["Peptide"]),
                                        int(bdecoy["Score"]), int(bdecoy["NumBY"])) )
                else:
                    if charge:
                        targets.append( (sid, Peptide(btarget["Peptide"]),
                                         float(btarget["Score"]), int(btarget["Charge"])))
                        decoys.append( (sid, Peptide(bdecoy["Peptide"]),
                                        float(bdecoy["Score"]), int(bdecoy["Charge"])))
                    else:
                        targets.append( (sid, Peptide(btarget["Peptide"]),
                                         float(btarget["Score"])))
                        decoys.append( (sid, Peptide(bdecoy["Peptide"]),
                                        float(bdecoy["Score"])))


            except ValueError:
                print records
                print 'Record %d in %s is bad' % (sid, f.name)
                raise

    f.close()

    spectrum = lambda r: r[0]
    targets.sort(key = spectrum)
    decoys.sort(key = spectrum)
    if len(targets) != len(decoys):
        raise Exception( ('Input file %s is bad, len(targets) = %d, '
                          'len(decoys) = %d' % (filename, len(targets), len(decoys))) )

    return targets, decoys

def load_shards(directory, progressbar = False):
    """Collect the best target and decoy identifications from the shards.

    The results of parallel testing are a collection of shard-#.out files,
    where # = shard number. Each file contains all the scored peptide-spectrum
    matches for a few spectra, in the following format:

    <Kind> <Sid> <Peptide> <Score>

    where

    <Kind> = 'target' | 'decoy'.
    <Sid> is the scan id / spectrum id.
    <Peptide> is a Peptide object.
    <Score> is the score assigned by the DBN model.

    Arguments:
       directory: Where the shard-#.out files are located. The files must
          match the regex 'shard*.out' to be loaded.

    Returns:
       Two lists, targets and decoys. Each list consists of triplets of the
       form (sid, peptide, score), where peptide is a protein.peptide.Peptide
       object, sid is an int, and score is a float. The peptide is the best
       target or decoy peptide for the given spectrum.

    """
    targets = [ ]
    decoys = [ ]
    shard_files = glob.glob(os.path.join(directory, 'shard-*.out'))

    pb = None
    if progressbar:
        widgets = [ progress.Percentage(), progress.Bar(), progress.ETA() ]
        pb = progress.ProgressBar(widgets = widgets, maxval = len(shard_files))
        pb.start()

    nShardsLoaded = 0
    for fn in shard_files:
        nShardsLoaded = nShardsLoaded + 1
        if pb: pb.update(nShardsLoaded)
        if os.stat(fn).st_size == 0:
            print '%s is empty, skipping' % fn
            continue

        tlist = [ ]
        dlist = [ ]
        try:
            for t, d in _load_shard(fn):
                tlist.append(t)
                dlist.append(d)
        except (ValueError, IOError):
            print 'Some exception in %s...skipping that shard' % fn
            continue
        else:
            targets.extend(tlist)
            decoys.extend(dlist)

    spectrum = lambda r: r[0]
    targets.sort(key = spectrum)
    decoys.sort(key = spectrum)
    return (targets, decoys)

def _load_shard(filename):
    """Loads an individual shard-#.out file. See load_shards."""

    fieldnames = ('Kind', 'Sid', 'Charge', 'Peptide', 'Score')
    reader = csv.DictReader(open(filename), delimiter = '\t',
                            fieldnames = fieldnames, skipinitialspace = True)

    for sid, rows in itertools.groupby(reader, lambda r: int(r["Sid"])):
        records = list(rows)
        trecords = itertools.ifilter(lambda r: r["Kind"] == "target", records)
        drecords = itertools.ifilter(lambda r: r["Kind"] == "decoy", records)
        t = max(trecords, key = lambda r: float(r["Score"]))
        d = max(drecords, key = lambda r: float(r["Score"]))

        assert(t["Kind"] == "target")
        assert(d["Kind"] == "decoy")
        tr = (sid, str(Peptide(t["Peptide"])), float(t["Score"]), int(t["Charge"]))
        dr = (sid, str(Peptide(d["Peptide"])), float(d["Score"]), int(d["Charge"]) )
        yield (tr, dr)
