#!/bin/bash
set -e

vcf=$1
output=`basename $vcf .vcf`
stargazer_folder='/Users/lamnguyen/Works/pharmcat/bin/stargazer-grc38-v.2.0.2/'
stargazer="${stargazer_folder}/stargazer"
genes=('abcb1' 'cacna1s' 'cftr' 'cyp1a1' 'cyp1a2' 'cyp1b1' 'cyp2a6' 'cyp2a13' 
        'cyp2b6' 'cyp2c8' 'cyp2c9' 'cyp2c19' 'cyp2d6' 'cyp2e1' 'cyp2f1' 'cyp2j2' 'cyp2r1'
        'cyp2s1' 'cyp2w1' 'cyp3a4' 'cyp3a5' 'cyp3a7' 'cyp3a43' 'cyp4a11' 'cyp4a22' 'cyp4b1' 'cyp4f2'
        'cyp17a1' 'cyp19a1' 'cyp26a1' 'dpyd' 'g6pd' 'gstm1' 'gstp1' 'ifnl3' 'nat1' 'nat2' 'nudt15'
        'por' 'ptgis' 'ryr1' 'slc15a2' 'slc22a2' 'slco1b1' 'slco1b3' 'slco2b1' 'sult1a1' 'tbxas1'
        'tpmt' 'ugt1a1' 'ugt1a4' 'ugt2b7' 'ugt2b15' 'ugt2b17' 'vkorc1' 'xpc' '2c_cluster' 'abcg2')

for gene in ${genes[@]};
do
cd $stargazer_folder
python $stargazer -d chip -a grc38 -t $gene -o $output -i $vcf
done