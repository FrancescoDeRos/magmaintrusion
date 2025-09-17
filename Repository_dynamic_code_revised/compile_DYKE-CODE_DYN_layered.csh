#!/bin/csh 

#~ set FFLAGS = '-g -fbounds-check -ffpe-trap=[invalid,zero,overflow,underflow,precision] -finit-real=snan -Wall'
set FFLAGS = '-g -fbacktrace -Wall -fcheck=all'
set MODPATH=./MODULES_F90/
set MAINPATH=./MAIN/
gfortran -O3 -c ${MODPATH}DISL2D.f90 ${MODPATH}CRACK2D.f90 ${MODPATH}EXTERNAL_FIELD.f90 ${MODPATH}BELEMENT_4FE-STRESS_DYN-CRACK_layered.f90 ${MODPATH}COMP_FIELD_4FE-STRESS.f90 ${MODPATH}PROPAGATION_4FE-STRESS_DYN-CRACK_layered.f90 ${MAINPATH}MAIN_PROPAGATION_DYN-CRACK_layered.f90
rm *.mod           
gfortran -O3 -o ${MAINPATH}DYKE-SIMULATION_DYN-CRACK_layered *.o -llapack
rm *.o




