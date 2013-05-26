#!/bin/sh
gunzip pop_stacks_example.tar.gz

TARGET=$HOME/Desktop/stacks_large
tar xf pop_stacks_example.tar
mv -f $TARGET $TARGET.old
#rm -rf $HOME/Desktop/stacks_large
mv nathan_stakcs $TARGET
cp popmap $TARGET/popmap

