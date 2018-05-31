#!/usr/bin/python3
import matplotlib.pyplot as plt

f=open("cout.csv","r")
fp=open("price.csv","r")

cout=[]
prix=[]
ecart_type=[]
e_ts=[]
dif=[]
us=[]
u2s=[]
#sqdif=[]

for ligne in f:
    p=ligne.split(";")[1:]
    cout=cout+[[float(i) for i in p]]
    ecart_type.append(float(ligne.split(";")[0]))
    
for ligne in fp:
	p=ligne.split(";")[1:]
	prix.append(float(p[1]))
	

nnoeuds=15 #int(round(len(p)/2))
e_t=0
uavg=0
uavg2=0
i=0

for r in range(len(ecart_type)):
	d=0
	#sq=[]
	for noeud in range(nnoeuds):
		d+=(cout[r][noeud+nnoeuds]-cout[r][noeud])
		#sq.append((cout[r][noeud+nnoeuds]-cout[r][noeud])**2)
	dif.append(d)
	#sqdif.append(sq)
	e_t=ecart_type[r]
	uavg+=prix[r]
	uavg2+=prix[r]**2
	i+=1
	if r<len(ecart_type)-1:
		if e_t!=ecart_type[r+1]:
			e_ts.append(e_t)
			us.append(uavg/i)
			u2s.append((uavg2-uavg*uavg/i)/i)
			i=0
			uavg=0
	

#print(reserve)

##t=list(range(len(prix)))
"""plt.plot(ecart_type,dif)
##plt.xlabel("Itérations")
##plt.ylabel("u")
plt.title("dif=f(reserve)")
plt.figure()"""

"""plt.plot(reserve,sqdif)
plt.title("sqdif=f(reserve)")
plt.figure()"""

plt.plot(e_ts,us)
plt.title("prix=f(ecart_type)")
plt.figure()

plt.plot(e_ts,u2s)
plt.title("variance_prix=f(ecart_type)")

plt.show()

f.close()
fp.close()
