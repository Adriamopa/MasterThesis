#!/bin/bash

wd="$(pwd)"

mkdir -p $wd/Results_ScoreInvHap
mkdir -p $wd/Results_ScoreInvHap/Inputs

parallel -j16 -u "bash $wd/scoreInvHap/Inputmaker_ScoreInvHap.sh" < $wd/../VCFs/30X/Invs2Imp
