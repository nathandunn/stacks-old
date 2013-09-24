#!/bin/sh
TARGET=$HOME/Desktop/stacks_tut
rm -rf stacks
tar xfz stacks_tut.tar.gz
rm -rf $TARGET
mv stacks $TARGET
rm -f $TARGET/progeny_[0-9][0-9]*

