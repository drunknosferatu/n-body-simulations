import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
start=time.time()
fig=plt.figure()
fig2=plt.figure()
ax=fig.add_subplot(projection='3d')
ax2=fig2.add_subplot()
def gera_part(sis,v):
    x=np.random.rand()/((7/4)**(1/3))
    y=np.random.rand()/((7/4)**(1/3))
    z=np.random.rand()/((7/4)**(1/3))
    vx=np.random.rand()*0
    vy=np.random.rand()*0
    vz=np.random.rand()*0
    part=[x,y,z]
    vpart=[vx,vy,vz]
    sis.append(part)
    v.append(vpart)
sis=[]
d=[]
aux=[]
q=[]
n=4000
nprint=n
mass=1.6605*1e-19
massb=1.6605*1e-27
charge=1.6*1e-19
m=[]
cont_el=0
t=0.0005
nb=200
j=0
v=[]#questionar se seria necessario gerar velocidades randômicas para as partículas
for i in range(n):
    gera_part(sis,v)
    if(cont_el!=nb):
        q.append(charge)
        cont_el+=1
        m.append(massb)
    else:
        m.append(mass)
        q.append(0)
v=np.array(v)
m=np.reshape(np.array(m),(n,1))
q=np.reshape(np.array(q),(n,1))
v-=np.mean(v*m,0)/np.mean(m)
auxplot=np.array(sis)
ax.scatter(np.reshape(auxplot[:,:-2],n),np.reshape(auxplot[:,1:-1],n),np.reshape(auxplot[:,2:],n),c='b',marker='D',s=m)
ael=[]
G=6*1e-11
c=9*1e9
controle=0
tsum=0
teste=[]
tf=[]
while(controle==0 and j<1e5):
    if(j):
        print(t)
        anorm=np.sqrt(np.sum(a**2,axis=1))
        controle=np.argwhere(anorm<=5e-30).shape[0]
        t=0.005/anorm[np.argmax(anorm)]
        tsum+=t
        tf.append(tsum)
        teste.append(anorm[400])
        sis=sis.tolist()
    aux=[]
    aux2=[]
    aux3=[]
    d=[]
    for i in range(n):
        aux=np.reshape(n*sis[i],(n,3))
        aux=np.array(sis)-aux
        d.append(aux)
    d=np.array(d)
    dnorm=d**2
    for i in range(n):
        aux2.append(np.sqrt(np.sum(dnorm[i],axis=1)))
    dnorm=np.array(aux2)
    for i in range(n):
        dnorm[i][i]=1
    args_coup=np.argwhere(dnorm<=1e-3)
    filt=args_coup.shape[0]
    for i in range(math.ceil(filt/2)):
        args_coup=np.delete(args_coup,np.argwhere(args_coup==[args_coup[i][1],args_coup[i][0]])[0][0],0)
    if(filt!=0):
        print(args_coup)
        for i in range(math.ceil(filt/2)):#neste trecho as partículas são acopladas, porém elas podem interagir de outra forma, perguntar ao vitor.
            auxd=[]
            indic=args_coup[i][0]
            flag_m=2
            flag_d=2
            if(i<args_coup.shape[0]-1):
                while(args_coup[i+1][0]==indic):
                    if(m[args_coup[i+1][1]][0]>m[args_coup[i][1]][0]):
                        flag_m=1
                    elif(m[args_coup[i+1][1]][0]<m[args_coup[i][1]][0]):
                        flag_m=0
                    if(dnorm[indic][args_coup[i][1]]>dnorm[indic][args_coup[i+1][1]]):
                        flag_d=1
                    elif(dnorm[indic][args_coup[i+1][1]]>dnorm[indic][args_coup[i][1]]):
                        flag_d=0
                    if (flag_d==1 and flag_m==1):
                        args_coup=np.delete(args_coup,i,0)
                    if (flag_d==0 and flag_m==0):
                        args_coup=np.delete(args_coup,i+1,0)
                    if (flag_d==0 and flag_m==1):
                        if(m[args_coup[i+1][1]][0]>dnorm[indic][args_coup[i][1]]**(-2)):
                            args_coup=np.delete(args_coup,i,0)
                        else:
                            args_coup=np.delete(args_coup,i+1,0)
                    if (flag_d==1 and flag_m==0):
                        if(m[args_coup[i][1]][0]>dnorm[indic][args_coup[i+1][1]]**(-2)):
                            args_coup=np.delete(args_coup,i+1,0)
                        else:
                            args_coup=np.delete(args_coup,i,0)
            print(args_coup)
            for l in range(3):#limpeza de excessos
                v[args_coup[i][1]][l]=(m[args_coup[i][1]]*(v[args_coup[i][1]][l]+a[args_coup[i][1]][l]*t)+m[indic]*(v[indic][l]+a[indic][l]*t))/(m[args_coup[i][1]]+m[indic])
            v=np.delete(v,indic,0)
            sis=np.delete(sis,indic,0).tolist()
            m[args_coup[i][1]]+=m[indic]
            q[args_coup[i][1]]+=q[indic]
            m=np.delete(m,indic,0)
            q=np.delete(q,indic,0)
            d=np.delete(d,indic,0)
            dnorm=np.delete(dnorm,indic,0)
            dnorm=np.delete(dnorm,indic,1)
            for k in range(n-1):
                auxd.append(np.delete(d[k],indic,0))
            d=np.array(auxd)
            n-=1
            args_coup-=1
    d=np.reshape(d,(n**2,3))
    dnorm=np.reshape(dnorm,(n**2))
    for i in range(3):
        aux3.append(np.reshape(d[ : ,i:1+i],(n**2))/dnorm**3)
    aux3=np.array(aux3)
    d=np.reshape(aux3.transpose(),(n,n,3))
    a=[]
    for i in range(n):
        a.append(np.reshape(d[i].transpose()@(G*m),3)-np.reshape(d[i].transpose()@(c*q*q[i][0]/m[i][0]),3))
    a=np.array(a)
    v+=a*t/2
    sis=np.array(sis)+v*t
    v+=a*t/2
    j+=1
print(m)
ax.scatter(np.reshape(sis[:,:-2],n),np.reshape(sis[:,1:-1],n),np.reshape(sis[:,2:],n),c='r',marker='o')
ax2.plot(tf,teste)
print(time.time()-start)
plt.show()

