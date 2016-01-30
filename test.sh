#!/bin/bash
#
# Written by John Halloran <halloj3@ee.washington.edu>
#
# Copyright (C) 2016 John Halloran
# Licensed under the Open Software License version 3.0
# See COPYING or http://opensource.org/licenses/OSL-3.0


function trainTest {
    echo "Training DRIP"
    python -OO dripTrain.py \
    	--psm-library data/riptideTrainingData/strict-orbitrap.psm \
    	--spectra data/riptideTrainingData/strict-orbitrap.ms2 \
    	--output-mean-file dripLearned.means \
    	--output-covar-file dripLearned.covars \
    	--mods-spec 'C+57.0214'

    echo "Digesting protein database"
    python -OO dripDigest.py \
    	--min-length 6 \
    	--fasta data/yeast.fasta \
    	--enzyme 'trypsin/p' \
    	--monoisotopic-precursor true \
    	--missed-cleavages 0 \
    	--digestion 'full-digest'

    echo "Searching spectra"
    python -OO dripSearch.py \
    	--digest-dir 'dripDigest-output' \
    	--precursor-window 3.0 \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
    	--num-threads 8 \
    	--top-match 1 \
    	--spectra data/test.ms2 \
    	--output dripSearch-test-output
}

function trainTestRecalibrate {
    echo "Training DRIP"
    python -OO dripTrain.py \
    	--psm-library data/riptideTrainingData/strict-orbitrap.psm \
    	--spectra data/riptideTrainingData/strict-orbitrap.ms2 \
    	--output-mean-file dripLearned.means \
    	--output-covar-file dripLearned.covars \
    	--mods-spec 'C+57.0214'

    echo "Digesting protein database"
    python -OO dripDigest.py \
	--recalibrate True \
    	--min-length 6 \
    	--fasta data/yeast.fasta \
    	--enzyme 'trypsin/p' \
    	--monoisotopic-precursor true \
    	--missed-cleavages 0 \
    	--digestion 'full-digest'

    echo "Searching spectra"
    python -OO dripSearch.py \
    	--digest-dir 'dripDigest-output' \
    	--precursor-window 3.0 \
    	--learned-means dripLearned.means \
    	--learned-covars dripLearned.covars \
    	--num-threads 8 \
    	--top-match 1 \
    	--spectra data/test.ms2 \
    	--output dripSearch-test-output
}

function clusterTest {
    if [ ! $DRIPTOOLKIT ]
    then
	echo "Please set DRIPTOOLKIT environment variable a directory containing the built toolkit."
    else

	echo "Digesting protein database"
	python -OO dripDigest.py \
    	    --min-length 6 \
    	    --fasta data/yeast.fasta \
    	    --enzyme 'trypsin/p' \
    	    --monoisotopic-precursor true \
    	    --missed-cleavages 0 \
    	    --digestion 'full-digest'

	echo "Creating cluster jobs"
	python -OO dripSearch.py \
    	    --digest-dir 'dripDigest-output' \
    	    --precursor-window 3.0 \
    	    --learned-means dripLearned.means \
    	    --learned-covars dripLearned.covars \
    	    --num-threads 8 \
	    --num-cluster-jobs 4 \
    	    --top-match 1 \
    	    --spectra data/test.ms2 \
    	    --output dripSearch-clusterTest-output

	echo "Running created jobs (locally)"
	for j in encode/*.sh
	do
    	    echo $j
    	    ./$j
	done

	echo "Merging cluster results"
	python -OO dripSearch.py --merge-cluster-results True \
    	    --logDir log \
    	    --output dripSearch-clusterTest-output
    fi
}

trainTest
# trainTestRecalibrate
# clusterTest
