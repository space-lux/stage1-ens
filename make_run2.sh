#!/bin/sh
gcc noeud.c -lm -lpthread -o noeud -Ofast -march=native
./noeud
