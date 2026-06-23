#!/bin/bash
gene="../../../2.DEG/all_exp_no_duplicate_genesymbol.xls"
WD=$(pwd)
g1="Healthy"
g2="ADAL"

###########################################################

rm -rf  01_cluster && mkdir 01_cluster && cd 01_cluster

/media/nfs/nfs02/liudy/biosoft/mambaforge/envs/R4_1_3/bin/Rscript $WD/cluster_analysis.r \
$gene \
$g1 $g2
cd ..

###########################################################


