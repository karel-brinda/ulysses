#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses bitwise -a ../_output/1.bf ../_output/2.bf ../_output/and.bf

ulysses bitwise -t -o ../_output/1.bf ../_output/2.bf ../_output/or.bf > ../_output/or_stats.txt

ulysses bitwise -x ../_output/1.bf ../_output/2.bf ../_output/xor.bf

ulysses bitwise -s -t -a ../_output/1.bf ../_output/2.bf ../_output/and_shrink.bf > ../_output/and_stats.txt
