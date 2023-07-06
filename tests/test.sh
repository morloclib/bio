#!/usr/bin/env bash

for testdir in test-*
do
    echo $testdir
    cd $testdir
    make
done
