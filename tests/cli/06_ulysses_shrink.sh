#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses shrink -f 4 ../_output/1.bf ../_output/1_shrink4.bf

ulysses shrink -f 2 ../_output/1_shrink4.bf ../_output/1_shrink4.bf
