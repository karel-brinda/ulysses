#! /usr/bin/env bash

set -eu
set -o pipefail

cd "$(dirname "$0")"

BIN=$(cd ../bin; pwd)
export PATH=$PATH:$BIN

ulysses_tests
