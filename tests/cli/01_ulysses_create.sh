#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

../../ulysses create ../_data/1.fa ../_output/1.bf

../../ulysses create -s '###-##-###-###-#' -r ../_data/1.fa ../_output/1_fr.bf

../../ulysses create ../_data/2.fa ../_output/2.bf

