#!/bin/sh
gcc noeud.c -lm -lpthread -o noeud -Ofast -march=native
./noeud
#wget "https://smsapi.free-mobile.fr/sendmsg?user= &pass= &msg=Calcul%20terminé%20!" -q --spider
