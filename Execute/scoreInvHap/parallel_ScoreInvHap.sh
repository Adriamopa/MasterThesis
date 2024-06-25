#!/bin/bash
## DON'T RUN THIS SCRIPT, run run_ScoreInvHap.sh instead
wd="$(pwd)"

mkdir -p $wd/Results_ScoreInvHap
mkdir -p $wd/Results_ScoreInvHap/Outputs

parallel -j4 -u "bash $wd/scoreInvHap/run_ScoreInvHap.sh" < $wd/../VCFs/30X/Invs2Imp
