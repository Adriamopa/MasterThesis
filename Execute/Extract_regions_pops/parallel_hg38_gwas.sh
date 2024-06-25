#!/bin/bash

wd="$(pwd)"

mkdir -p $wd/../VCFs/GWAS_VCFs/Input_VCFs/

parallel -j16 -u "bash $wd/Extract_regions_pops/run_hg38_gwas.sh" < $wd/../VCFs/30X/Invs2Imp
