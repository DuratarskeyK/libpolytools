#!/bin/bash

#build script for libpolytools

#change this
basedir=/home/lichak/libpolytools
#this
clp_include=/home/lichak/usr/include/coin
#and this
clp_lib=/home/lichak/usr/lib

gcc -I$basedir/include -Wall -lgmp -DGMP -DSILENTLRS -fPIC -c *.c
g++ -I$basedir/include -I$clp_include -Wall -lgmp -DGMP -fPIC -c *.cpp

g++ -L$clp_lib -lClp -lCoinUtils -lbz2 -lz -lm -lgmp -z defs -shared -Wl,-soname,libpolytools.so -Wl,--rpath,$clp_lib -o libpolytools.so lrslib.o lrsgmp.o facet_cnt.o adjacency.o
