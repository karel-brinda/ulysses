#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

../../ulysses hamming ../_output/1.bf ../_output/2.bf > ../_output/hamming.txt
