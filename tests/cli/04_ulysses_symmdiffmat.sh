#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

ulysses symmdiffmat ../_output/1.bf ../_output/2.bf > ../_output/symmdiffmat.txt
