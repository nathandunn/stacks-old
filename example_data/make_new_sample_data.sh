#!/bin/sh
TARGET=$HOME/Desktop/stacks_new_tut
rm -rf stacks
tar xfz stacks_new_tut.tgz
rm -rf $TARGET
mv stacks $TARGET
#cp popmap_stacks_tut.txt /tmp/stacks_tut/popmap

