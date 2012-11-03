#!/bin/bash

dest=/Users/frankmeier/Documents/psi/lambda/anLambdaBlifetime/notes/AN-12-044/trunk/img_detperf_reso

#transferlist="MCtruthCt3dTruthB0 MCtruthCt3dCutTruthB0 MCtruthd3dTruthB0 MCtruthd3dCutTruthB0 MCtruthpTruthB0 MCtruthpCutTruthB0 TtruthVsT3dprofB0 TtruthVsT3dprof_cutsB0 DtruthVsD3dprofB0 DtruthVsD3dprof_cutsB0 PtruthVsP3dprofB0 PtruthVsP3dprof_cutsB0"
#transferlist="MCtruthCt3dTruthLb MCtruthCt3dCutTruthLb MCtruthd3dTruthLb MCtruthd3dCutTruthLb MCtruthpTruthLb MCtruthpCutTruthLb TtruthVsT3dprofLb TtruthVsT3dprof_cutsLb DtruthVsD3dprofLb DtruthVsD3dprof_cutsLb PtruthVsP3dprofLb PtruthVsP3dprof_cutsLb"

#transferlist="MCtruthCt3dTruth MCtruthCt3dCutTruth MCtruthd3dTruth MCtruthd3dCutTruth MCtruthpTruth MCtruthpCutTruth TtruthVsT3dprof TtruthVsT3dprof_cuts DtruthVsD3dprof DtruthVsD3dprof_cuts PtruthVsP3dprof PtruthVsP3dprof_cuts"
transferlist="Ct3dECut TtruthVsT3dprofVar TtruthVsPtprofVar TtruthVsEtaprofVar TtruthVsT3dprofVar_cuts TtruthVsPtprofVar_cuts TtruthVsEtaprofVar_cuts tEVsT3dprofVar tEVsPtprofVar tEVsetaprofVar tEVsT3dprofVar_cuts tEVsPtprofVar_cuts tEVsetaprofVar_cuts"
transferlist="MCtruthCt3dPull MCtruthCt3dPullCutTruth"


for i in $transferlist
do
    echo $i
    cp -v {doAnalysis01plots*_,$dest/}${i}B0.pdf
    cp -v {doAnalysis01plots*_,$dest/}${i}Lb.pdf
done

#cp -v {doAnalysis01plots*_,$dest/}MCtruthCt3dB0.pdf

#cp -v {doAnalysis01plots*_,$dest/}MCtruthCt3dLb.pdf

