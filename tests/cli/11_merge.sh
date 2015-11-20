#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses merge ../_output/query_split_found1.fa ../_output/query_split_notfound1.fa > ../_output/merge1.fa

ulysses merge -i ../_output/query_split_found2.fa ../_output/query_split_notfound2.fa > ../_output/merge2.fa

ulysses merge ../_output/query_split_found3.fa ../_output/query_split_notfound3.fa > ../_output/merge3.fa

ulysses merge -i ../_data/35268_notfound_inter.fa.tmp ../_data/35268_notfound_inter.fa.tmp > ../_output/merge4.fa

dff=`diff ../_data/35268_notfound_inter.fa.tmp ../_output/merge4.fa | wc -l`
if [ $dff -gt 0 ]; then
    exit 1
fi


