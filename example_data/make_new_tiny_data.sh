#!/bin/sh
TARGET=$HOME/Desktop/new_tut_tiny
rm -rf stacks_tiny
tar xfz stacks_new_tut_tiny.tgz
rm -rf $TARGET
mv stacks_tiny $TARGET
#cp popmap_stacks_tut.txt $TARGET/popmap
#rm -f $TARGET/popmap

