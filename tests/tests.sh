#! /usr/bin/env bash

set -ex

../ulysses create 1.fa _1.bf
../ulysses create -s '###-##-###-###-#' -r 1.fa _1_fr.bf
../ulysses create 2.fa _2.bf
../ulysses bitwise -a _1.bf _2.bf _and.bf
../ulysses bitwise -t -o _1.bf _2.bf _or.bf > _or_stats.txt
../ulysses bitwise -x _1.bf _2.bf _xor.bf
../ulysses hamming _1.bf _2.bf > _hamming.txt
../ulysses symmdiffmat _1.bf _2.bf > _symmdiffmat.txt
../ulysses dump _1.bf > _dump1.txt
../ulysses dump -x _1.bf > _dump2.txt
../ulysses dump -b _1.bf > _dump3.txt
../ulysses stats _and.bf > _stats.txt
../ulysses query -r _1.bf query.fa > _query.txt
../ulysses query -r _1_fr.bf query.fa > _query_fr.txt
../ulysses shrink -f 4 _1.bf _1_shrink4.bf
../ulysses shrink -f 2 _1_shrink4.bf _1_shrink4.bf
../ulysses bitwise -s -t -a _1.bf _2.bf _and_shrink.bf > _and_stats.txt
../ulysses stats _1.bf > _stats_shrink.txt
../ulysses stats _1_shrink4.bf >> _stats_shrink.txt
cat 1.fa | ../ulysses create_many -t -a 200000 -h 10  -m map_file.tab -F /dev/fd/0  ./
cat 1.fa | ../ulysses create_many -t -e _1.bf -a 200000 -h 10  -m map_file2.tab -F /dev/fd/0  ./
