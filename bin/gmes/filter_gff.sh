#!/bin/bash
# ==============================================================
# Tomas Bruna
#
# This script filters out gff entries with specified score
# lower than the desired threshold.
#
# If a score name is supplied, the script looks for it in gff's 9th column.
# Without a score name, script uses gff's 6th column.
#
# Usage ./filter_gff input_file.gff threshold [score_name]
# ==============================================================

if [ $# -lt 2 ]; then
    echo "Usage: $0 input_file.gff threshold [score_name]"
    exit 1
fi

if [ -z $3 ]; then
    awk -v threshold=$2 '
    {
        if ($6 + 0 >= threshold + 0) {
            print $0;
        }
    }
    ' $1
else
    awk -v threshold=$2 \
    -v score_name=$3 '
    {
        regex = score_name"=([^;]+);";
        match($0, regex);
        score = substr($0, RSTART + length(score_name) + 1, RLENGTH - length(score_name) - 1);
        if (score + 0 >= threshold + 0) {
            print $0;
        }
    }
    ' $1
fi
