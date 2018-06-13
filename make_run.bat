gcc noeud.c -o noeud.exe -lm -lpthread -std=gnu99 -Ofast -march=native
noeud.exe
rem notif par SMS (manque ID&pass)
rem certutil -urlcache -split -f "https://smsapi.free-mobile.fr/sendmsg?user= &pass= &msg=Calcul%20terminé%20!" sms.html
