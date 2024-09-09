# Parameter estimation and model selection for stochastic differential equations for biological growth

Hay 2 formas para ejecutar la aplicación de consola
1. Ejecutando el script **execute.sh** ejecutando el siguiente comando en consola:

	    bash execute.sh

2. Ejecutar ejecutar los siguientes comandos en consola:

	    mkdir  build
	    gfortran  -c  -Jbuild  -o  build/rand_seed.o  modules/rand_seed.f90
	    gfortran  -c  -Jbuild  -o  build/integral.o  modules/integral.f90
	    gfortran  -c  -Jbuild  -o  build/simulation.o  modules/simulation.f90
	    gfortran  -c  -Jbuild  -o  build/qua.o  modules/qua.f90
	    gfortran  -c  -Jbuild  -o  build/mle.o  modules/mle.f90
	    gfortran  -c  -Jbuild  -o  build/aic.o  modules/aic.f90
	    gfortran  -c  -Jbuild  -o  build/menu.o  modules/menu.f90
	    mkdir  bin
	    gfortran  -o  bin/main.exe  -Ibuild  src/main.f90  build/*.o
	    ./bin/main
	    rm  build/*  bin/*
	    rmdir  build  bin
