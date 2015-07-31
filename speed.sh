#! /usr/bin/env bash

if [ ! -f chr21.fa ];
then
	curl -o chr21.fa.gz http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
	gzip -d chr21.fa.gz
fi

cmake . && make
touch gmon.out
time (./ulysses create -a 100000000 ./chr21.fa ./bf.bf) > $1 
cp gmon.out gmon.out.create
time (./ulysses bitwise -x bf.bf bf.bf xor.bf) >> $1 
cp gmon.out gmon.out.bitwise
time (./ulysses symmdiffmat bf.bf bf.bf xor.bf) >> $1 
cp gmon.out gmon.out.symmdiffmat
cat $1

