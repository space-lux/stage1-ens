#!/bin/sh
mpicc noeud.c -o mpi -lm
mpirun -np 175 ./mpi
