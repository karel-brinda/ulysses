#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses query_and_split -r ../_output/1.bf ../_data/query.fa ../_output/query_split_found1.txt ../_output/query_split_notfound1.txt

ulysses query_and_split -r ../_output/1_fr.bf ../_data/query.fa ../_output/query_split_found2.txt ../_output/query_split_notfound2.txt
