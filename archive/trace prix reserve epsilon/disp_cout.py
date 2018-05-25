#!/usr/bin/python3
import matplotlib.pyplot as plt

f=open("cout.csv","r")
fp=open("price.csv","r")

cout=[]
prix=[]
reserve=[]
dif=[]
sqdif=[]

for ligne in f:
    p=ligne.split(";")[1:]
    cout=cout+[[float(i) for i in p]]
    reserve=reserve+[float(ligne.split(";")[0])]
    
for ligne in fp:
	p=ligne.split(";")[1:]
	prix.append(float(p[1])-float(p[0]))
	

nnoeuds=15 #int(round(len(p)/2))

for r in range(len(reserve)):
	d=0
	sq=[]
	for noeud in range(nnoeuds):
		d+=(cout[r][noeud+nnoeuds]-cout[r][noeud])
		sq.append((cout[r][noeud+nnoeuds]-cout[r][noeud])**2)
	dif.append(d)
	sqdif.append(sq)

#print(reserve)

##t=list(range(len(prix)))
plt.plot(reserve,dif)
##plt.xlabel("Itérations")
##plt.ylabel("u")
plt.title("dif=f(reserve)")
plt.figure()

plt.plot(reserve,sqdif)
plt.title("sqdif=f(reserve)")
plt.figure()

plt.plot(reserve,prix)
plt.title("diff_prix=f(reserve)")

plt.show()

f.close()
fp.close()
