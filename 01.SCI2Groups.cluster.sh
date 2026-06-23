#!/bin/bash
gene="../../../2.DEG/all_exp_no_duplicate_genesymbol.xls"
WD=$(pwd)
g1="Healthy"
g2="ADAL"
###########################################################

rm -rf  01_cluster_sci2group && mkdir 01_cluster_sci2group && cd 01_cluster_sci2group
/media/nfs/nfs02/liudy/biosoft/mambaforge/envs/R4_1_3/bin/Rscript $WD/cluster_analysis_sci2group.r \
$gene \
$g1 $g2
cd ..

###########################################################