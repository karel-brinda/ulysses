#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses dump ../_output/1.bf > ../_output/dump1.txt

ulysses dump -x ../_output/1.bf > ../_output/dump2.txt

ulysses dump -b ../_output/1.bf > ../_output/dump3.txt
