#!/bin/sh
TARGET=$HOME/Desktop/stacks_tut
rm -rf stacks
tar xfz stacks_tut.tar.gz
rm -rf $TARGET
mv stacks $TARGET
#cp popmap_stacks_tut.txt /tmp/stacks_tut/popmap

