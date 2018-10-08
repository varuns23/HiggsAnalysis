#!/bin/bash
cd /afs/cern.ch/user/r/rcoelhol/work/public/CMSSW_8_0_12/src/HiggsAnalysis/GBRLikelihood/macros/eregdata
eval `scramv1 runtime -sh`
root -q -l -b -x eregtraining_data.C++g\(\)