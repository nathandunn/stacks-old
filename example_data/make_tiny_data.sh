#!/bin/sh
TARGET=$HOME/Desktop/stacks_female_only
rm -rf stacks_tut
tar xfz female_only_stacks_tut.tgz
rm -rf $TARGET
mv stacks_tut $TARGET
cp popmap_stacks_tut.txt $TARGET/popmap

