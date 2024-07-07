#!/bin/bash
mkdir build
gfortran -c -Jbuild -o build/integral.o modules/integral.f90
gfortran -c -Jbuild -o build/simulation.o modules/simulation.f90
gfortran -c -Jbuild -o build/qua.o modules/qua.f90   
gfortran -c -Jbuild -o build/mle.o modules/mle.f90     
gfortran -c -Jbuild -o build/aic.o modules/aic.f90 
mkdir bin
gfortran -o bin/main.exe -Ibuild src/main.f90 build/*.o
./bin/main
rm build/* bin/*
rmdir build bin