#!/bin/sh
TARGET=$HOME/Desktop/pop_stacks1
rm -rf pop_stacks_tiny
tar xfz pop_stacks_tiny.tgz
rm -rf $TARGET
mv pop_stacks_tiny $TARGET
#cp popmap_stacks_tut.txt $TARGET/popmap
#rm -f $TARGET/popmap

