#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

../../ulysses stats ../_output/and.bf > ../_output/stats.txt

../../ulysses stats ../_output/1.bf > ../_output/stats_shrink.txt

../../ulysses stats ../_output/1_shrink4.bf >> ../_output/stats_shrink.txt
