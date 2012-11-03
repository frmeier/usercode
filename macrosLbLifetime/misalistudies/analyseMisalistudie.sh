#!/bin/bash

path="../../data"

list="vrt_r575_nomisali vrt_r578_radialEpsilon_minus vrt_r579_radialEpsilon_plus vrt_r584_telescopeEpsilon_minus vrt_r585_telescopeEpsilon_plus vrt_r573_layerRotEpsilon_minus vrt_r574_layerRotEpsilon_plus vrt_r592_bowingEpsilon_minus vrt_r570_bowingEpsilon_plus vrt_r588_zExpEpsilon_minus vrt_r589_zExpEpsilon_plus vrt_r586_twistEpsilon_minus vrt_r587_twistEpsilon_plus vrt_r571_ellipticalEpsilon_minus vrt_r572_ellipticalEpsilon_plus vrt_r582_skewEpsilon_minus2 vrt_r583_skewEpsilon_plus2 vrt_r590_saggitaEpsilon_minus vrt_r591_saggitaEpsilon_plus vrt_r576_radialEpsilon_layer1minus2 vrt_r577_radialEpsilon_layer1plus2"

outfile="fit.out"
echo "Start" > $outfile

for i in $list;
do
    filename=$path/${i}_lbbar.root
    root -b -q "Fit_2D.C(\"$filename\")" >> $outfile
done


