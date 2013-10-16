#!/bin/sh
TARGET=$HOME/Desktop/pop_stacks_tiniest
rm -rf pop_stacks_tiniest
tar xfz pop_stacks_tiniest.tgz
rm -rf $TARGET
mv pop_stacks_tiniest $TARGET
#cp popmap_stacks_tut.txt $TARGET/popmap
#rm -f $TARGET/popmap

