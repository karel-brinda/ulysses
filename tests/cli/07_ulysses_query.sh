#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

../../ulysses query -r ../_output/1.bf ../_data/query.fa > ../_output/query.txt

../../ulysses query -r ../_output/1_fr.bf ../_data/query.fa > ../_output/query_fr.txt
