#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np

f2=open("2e_marche.csv","r")
fp=open("prevision.csv","r")
fs=open("simple.csv","r")
fr=open("reserve.csv","r")

cout_social=[[],[],[],[]]
cout_social_t=[]
prix=np.array([[],[],[],[]])
prix_t=[]
ecart_type=[]
ecart_type_t=[]
reserve=[]
prix_ecart_type=np.array()
cout_ecart_type=np.array()



for ligne in fs:
	p=np.asfarray(ligne.split(";"))
	cout_social[0].append(np.sum(p[1:]))
	prix[0].append(p[0])

for ligne in fr:
	p=np.asfarray(ligne.split(";"))
	reserve.append(p[0])
	cout_social[1].append(np.sum(p[2:]))
	prix[1].append(p[1])
	
for ligne in fp:
	p=np.asfarray(ligne.split(";"))
	cout_social[2].append(np.sum(p[1:]))
	prix[2].append(p[0])
	
for ligne in f2:
	p=np.asfarray(ligne.split(";"))
	ecart_type_t.append(p[0])
	cout_social_t.append(np.sum(p[2:]))
	prix_t.append(p[1])

n_et=0
while ecart_type_t[n_et]==ecart_type_t[0]:
	n_et+=1

ecart_type=ecart_type_t[::n_et]
cout_social[3] = np.mean(np.split(np.array(cout_social_t),
				np.arange(0, len(cout_social_t), n_et)), axis=1)
cout_ecart_type = np.std(np.split(np.array(cout_social_t),
				np.arange(0, len(cout_social_t), n_et)), axis=1)
prix[3] = np.mean(np.split(np.array(prix_t),
				np.arange(0, len(prix_t), n_et)), axis=1)
prix_ecart_type = np.std(np.split(np.array(prix_t),
				np.arange(0, len(prix_t), n_et)), axis=1)

plt.plot(reserve,prix[1],legend="1er prix avec réserve")
plt.plot(ecart_type,[prix[0][0] for i in range(len(prix[1]))],label="1er prix sans réserve")
plt.plot(ecart_type,[prix[2][0] for i in range(len(prix[1]))],label="1er prix avec estimation des incertitudes")
plt.legend()
plt.title("prix=f(reserve)")
plt.figure()

plt.plot(reserve,cout_social[1],legend="1er coût social avec réserve")
plt.plot(ecart_type,[cout_social[0][0] for i in range(len(cout_social[1]))],label="1er coût sans réserve")
plt.plot(ecart_type,[cout_social[2][0] for i in range(len(cout_social[1]))],label="1er coût avec estimation des incertitudes")
plt.legend()
plt.title("cout social=f(reserve)")
plt.figure()

plt.plot(ecart_type,prix[3],label="prix moyen")
plt.plot(ecart_type,prix[3]+prix_ecart_type,label="prix moyen+écart type prix")
plt.plot(ecart_type,prix[3]-prix_ecart_type,label="prix moyen-écart type prix")
plt.title("prix=f(ecart_type)")
plt.legend()
plt.figure()

plt.plot(ecart_type,cout_social[3],label="cout social moyen")
plt.plot(ecart_type,cout_social[3]+cout_ecart_type,label="coût moyen+écart type coût")
plt.plot(ecart_type,cout_social[3]-cout_ecart_type,label="coût moyen-écart type coût")
plt.title("cout social=f(ecart_type)")
plt.legend()
#plt.figure()

plt.show()

f2.close()
fp.close()
fr.close()
fs.close()
