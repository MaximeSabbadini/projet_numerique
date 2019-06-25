""" UE Projet intégrateur Numérique """



###########################################################################
# Modélisation de la trajectoire d'un objet naturel ou artificiel au      #
#                        voisinage d'une planète                          #
###########################################################################
                        
#________________________Maxime Sabbadini, 2019___________________________#                                                   




""" Importation des librairies """

import numpy as np 
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D


""" Déclaration des constantes """

G=6.67408e-11 #SI
M=5.972e24#Masse du corps central
rayon= 6371e3  #Rayon du corps central (m)

  #Distance terre-soleil (m) A enlever

""" Déclaration des conditions initiales """

ti=0      #Temps initial (s)
tf=4000  #Temps final (s)
dt=10    #Temps entre chaque points (s)
n=int((tf-ti)/dt)  #Nombre de points pris dans l'intervalle de temps

v0x=0          #m/s
v0y=7667+3170  
v0z=0

x0=6779e3     #m
y0=0
z0=0

""" Définition des fonctions """

def trajectoire_euler(r, v, ti, tf, dt, n):    #Fonction utilisant la méthode d'Euler
    v_x=v[0,0]                                            #Inutilisée 
    v_y=v[0,1]
    v_z=v[0,2]
    x=r[0,0]
    y=r[0,1]
    z=r[0,2]
    

    for i in range (1,n):
        re=np.sqrt(x**2+y**2+z**2)
        
        x+=dt*v_x
        y+=dt*v_y
        z+=dt*v_z
        
        r[i,0]=x
        r[i,1]=y
        r[i,2]=z
        
        v_x+=dt*(-G*M/re**3*x)
        v_y+=dt*(-G*M/re**3*y)
        v_z+=dt*(-G*M/re**3*z)
        
        v[i,0]=v_x
        v[i,1]=v_y        
        v[i,2]=v_z
        
    return v, r


r2 = np.zeros((3,n),dtype = float)   #Tableaux pour la méthode de RK4
v2 = np.zeros((3,n),dtype = float)

##Conditions initiales

r2[0,0] = x0
r2[1,0] = y0
r2[2,0] = z0
v2[0,0] = v0x
v2[1,0] = v0y
v2[2,0] = v0z

""" Fonction RK4 """

def f(r, v):
    x = r[0]
    y = r[1]
    z = r[2]
    re=np.sqrt(x**2+y**2+z**2)         #Définition du vecteur à 6 composantes
    vect = np.zeros(6)
    vect[0] = v[0]    
    vect[1] = v[1]
    vect[2] = v[2]
    vect[3] = -G*M/re**3*x
    vect[4] = -G*M/re**3*y
    vect[5] = -G*M/re**3*z
    return vect

def rk4(r, v, dt, n):
    for i in range(0,n-1):
    
        k1 = f(r[:,i],v[:,i])                         
        k2 = f(r[:,i]+ k1[0:3]*dt*0.5,v[:,i]+ k1[3:6]*dt*0.5)  #Intégration par RK4
        k3 = f(r[:,i]+ k2[0:3]*dt*0.5,v[:,i]+ k2[3:6]*dt*0.5)
        k4 = f(r[:,i]+ k3[0:3]*dt,v[:,i]+ k3[3:6]*dt)
        k = (k1 + 2*k2 + 2*k3 + k4)/6
        r[:,i+1] = r[:,i] + dt*k[0:3]
        v[:,i+1] = v[:,i] + dt*k[3:6]
        
    return r
 
L=rk4(r2,v2,dt,n)  #Appel de la fonction 

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]  #Tracé du corps central
x = rayon*np.cos(u)*np.sin(v)
y = rayon*np.sin(u)*np.sin(v)
z = rayon*np.cos(v)


##Tracé de/des trajectoire(s) 

fig=pl.figure(figsize=(8,5))
ax=fig.add_subplot(111, projection='3d')
ax.plot(L[0,:], L[1,:], L[2,:], marker='+', ls='', color='r',label='Trajectoire') #Tracé de la planète
ax.plot_wireframe(x, y, z, color="b", label='Planète')    #Tracé du corps central
pl.legend()
pl.show()







