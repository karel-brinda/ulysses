#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

cat ../_data/1.fa | ulysses create_many -t -a 200000 -h 10  -m ../_data/map_file.tab -F /dev/fd/0  ../_output/ > ../_output/1.bf.stats.txt

cat ../_data/1.fa | ulysses create_many -t -e ../_output/1.bf -a 200000 -h 10  -m ../_data/map_file2.tab -F /dev/fd/0  ../_output/ > ../_output/exclude_stats.txt

cat ../_data/1.fa | ulysses create_many -t -n ../_output/1.bf -a 200000 -h 10  -m ../_data/map_file3.tab -F /dev/fd/0  ../_output/ > ../_output/include_stats.txt
