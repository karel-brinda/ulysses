#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses query_and_split -t 4 -r ../_output/1.bf ../_data/query.fa ../_output/query_split_found1.fa ../_output/query_split_notfound1.fa

ulysses query_and_split -t 4 -r ../_output/1_fr.bf ../_data/query.fa ../_output/query_split_found2.fa ../_output/query_split_notfound2.fa

ulysses query_and_split -t 4 -r ../_output/1_fr.bf ../_output/query_split_found2.fa ../_output/query_split_found3.fa ../_output/query_split_notfound3.fa

