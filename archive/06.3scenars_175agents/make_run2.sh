#!/bin/sh
mpicc noeud.c -o mpi -lm
mpirun -np 15  ./mpi
