#!/bin/bash

dest=/Users/frankmeier/Documents/psi/lambda/anLambdaBlifetime/notes/AN-12-044/trunk/img_detperf_reso

#cp -v {doAnalysis01plots*_,$dest/}MCtruthCt3dB0.pdf
cp -v {doAnalysis01plots*_,$dest/}MCtruthCt3dTruthB0.pdf
cp -v {doAnalysis01plots*_,$dest/}MCtruthCt3dCutTruthB0.pdf
#cp -v {doAnalysis01plots*_,$dest/}MCtruthd3dB0.pdf
cp -v {doAnalysis01plots*_,$dest/}MCtruthd3dTruthB0.pdf
cp -v {doAnalysis01plots*_,$dest/}MCtruthd3dCutTruthB0.pdf
#cp -v {doAnalysis01plots*_,$dest/}MCtruthpB0.pdf
cp -v {doAnalysis01plots*_,$dest/}MCtruthpTruthB0.pdf
cp -v {doAnalysis01plots*_,$dest/}MCtruthpCutTruthB0.pdf

cp -v {doAnalysis01plots*_,$dest/}TtruthVsT3dprofB0.pdf
cp -v {doAnalysis01plots*_,$dest/}TtruthVsT3dprof_cutsB0.pdf
cp -v {doAnalysis01plots*_,$dest/}DtruthVsD3dprofB0.pdf
cp -v {doAnalysis01plots*_,$dest/}DtruthVsD3dprof_cutsB0.pdf
cp -v {doAnalysis01plots*_,$dest/}PtruthVsP3dprofB0.pdf
cp -v {doAnalysis01plots*_,$dest/}PtruthVsP3dprof_cutsB0.pdf

