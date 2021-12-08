#!/bin/bash

if [[ $# -lt 1 ]]; then
    exit -1
fi
in_fa=$(readlink -f $1)
name=$(basename $(basename $in_fa .fa) .fasta)

CURL="curl -s --proxy http://markov.bch.msu.edu:9999"

pwd=$(pwd)
tmpdir=$TMPDIR/mmseqs2.$$
mkdir $tmpdir
cd $tmpdir
ln -sf $in_fa .

tg=$($CURL -F q=@$in_fa -F mode=env https://a3m.mmseqs.com/ticket/msa | jq -r '.id')
echo "Running ..." $tg

status=$($CURL https://a3m.mmseqs.com/ticket/$tg | jq -r '.status')
while [[ "$status" == "RUNNING" || "$status" == "PENDING" ]];
do
    status=$($CURL https://a3m.mmseqs.com/ticket/$tg | jq -r '.status')
    sleep 10
done

if [ "$status" == "COMPLETE" ]; then
    $CURL https://a3m.mmseqs.com/result/download/$tg > $name.mmseqs2.tgz
    tar xzf $name.mmseqs2.tgz
    cat $in_fa > tmp.a3m
    tail -n+3 uniref.a3m >> tmp.a3m
    tail -n+3 bfd.mgnify30.metaeuk30.smag30.a3m >> tmp.a3m
    tr -d '\000' < tmp.a3m > $name.mmseqs2.a3m
else
    echo "FAILED"
fi

cd $pwd
mv $tmpdir/$name.mmseqs2.a3m .
rm -rf $tmpdir
