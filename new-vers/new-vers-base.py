import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import time
from random import random
#check for negative numbers



class Node:
    def __init__(self):
        self.node1 = None
        self.node2 = None
        self.node3 = None
        self.node4 = None
        self.node5 = None
        self.node6 = None
        self.node7 = None
        self.node8 = None
        self.mass = 0
        self.charge = 0
        self.leng = 0
        self.center = 0
    

    def insert(self, sis, masses, bounds, charges):
        call_strings = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8"]
        middle = (bounds[1] - bounds[0]) / 2
        aux = np.concatenate(([bounds[0]], [bounds[0] + middle]), axis=0) 
        aux2 = np.concatenate(([bounds[0] + middle], [bounds[1]]), axis=0)
        new_bounds = [aux2,
                np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1),
                np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1),
                np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1),
                np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1),
                np.concatenate((aux[:,:-2],aux2[:,1:-1], aux[:,-1:]), axis=1),
                np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1),
                aux]
        for count, cond in enumerate(new_bounds):
            index = np.argwhere(np.all(sis>new_bounds[count][0], axis=1) & np.all(sis<new_bounds[count][1], axis=1)).flatten()
            if index.size > 0 and sis.shape[0] > 1:
                setattr(self, call_strings[count], Node())
                getattr(self, call_strings[count]).insert(sis[index],  masses[index], new_bounds[count], charges[index])
                setattr(self, 'mass', self.mass + getattr(self, call_strings[count]).mass)
                setattr(self, 'charge', self.charge + getattr(self, call_strings[count]).charge)
                setattr(self, 'center', self.center + getattr(self, call_strings[count]).center)
                setattr(self, 'leng', abs(bounds[0,0] - bounds[1,0]))
        if sis.shape[0] == 1:
            setattr(self, 'mass', np.sum(masses))
            setattr(self, 'charge', np.sum(charges))
            setattr(self, 'leng', abs(bounds[0,0] - bounds[1,0]))
            setattr(self, 'center', np.sum(masses[:,None] * sis, axis=0))



    def calc_forces(self, sis, charges, masses, accel):
        call_strings = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8"]
        aux_dist =  self.center / self. mass - sis
        aux_dist = np.around(aux_dist, decimals=15)
        ind_exclud = np.argwhere(np.any((aux_dist) != 0, axis=1)).flatten()
        dist = np.linalg.norm(aux_dist, axis=1)
        dist = np.around(dist, decimals=15)
        if ind_exclud.shape[0] != sis.shape[0]:
            oriented_dist = aux_dist[ind_exclud] / (dist[ind_exclud]**3)[:,None]
            accel[ind_exclud] += 6.7e-11 * self.mass * oriented_dist - 9e9 * charges[ind_exclud,None] * self.charge * oriented_dist / masses[ind_exclud,None]
        else:
            oriented_dist = aux_dist / (dist**3)[:,None]
            ind_theta = np.argwhere(self.leng / dist < 1).flatten()
            if ind_theta.shape[0] > 0:
                accel[ind_theta] += 6.7e-11 * self.mass * oriented_dist[ind_theta] - 9e9 * charges[ind_theta,None] * self.charge * oriented_dist[ind_theta] / masses[ind_theta,None]
            ind_theta = np.argwhere(self.leng / dist >= 1).flatten()
            if ind_theta.shape[0] > 0:
                for call in call_strings:
                     if getattr(self,call):
                        accel[ind_theta] = getattr(self, call).calc_forces(sis[ind_theta], charges[ind_theta], masses[ind_theta], accel[ind_theta])
        return accel



def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis) - 0.1]), np.array(3*[np.amax(sis) + 0.1])])

def calc_energ(sis, v, masses, charge):
    x = sis[:,0:1]
    y = sis[:,1:2]
    z = sis[:,2:3]
    distx = x.T - x
    disty = y.T - y
    distz = z.T - z
    inv_r = np.sqrt(distx**2 + disty**2 + distz**2)
    args = np.argwhere(inv_r >0)
    inv_r[args] = 1.0/inv_r[args]
    remake_sis =  np.array(np.where(dist!=0))
    args = np.argwhere(dist!=0)
    return np.dot(0.5 * masses, np.linalg.norm(v, axis = 1)**2) + np.sum(9e9 * charge**2 / remake_sis) + 6.7e-11 * np.sum(np.sum(np.triu(-(masses*masses.T)*inv_r,1)))


fig3d = plt.figure()
ax3d = fig3d.add_subplot(projection = '3d')

t = 0.005
control = None
n = 1024000
nBar = 51200 
BarMass = 1.6605*1e-27
nonBarMass = 1.6605*1e-15
charge = 1.6*1e-19

sis = np.random.randn(n, 3) * ((4.5) ** (1/3))
v = np.zeros((n,3))
charges = np.concatenate((nBar * [charge], np.zeros((n-nBar))), axis=None)
masses = np.concatenate((nBar * [BarMass], (n-nBar) * [nonBarMass]), axis=None)

auxt = start = time()
while not control:
    accel = np.zeros((n, 3))
    bounds = calc_bounds(sis)
    root = Node()
    root.insert(sis, masses, bounds, charges)
    root.calc_forces(sis, charges, masses, accel)
    v += accel *  t / 2
    sis += v * t
    v += accel * t / 2
    anorm = np.linalg.norm(accel, axis=1)
    amax = np.amax(anorm)
    print(amax)
    control = np.all(anorm <= 5e-5)
    t =  0.005 / amax if 0.005 / amax <= 0.5 else 0.5
    print(auxt - time())
    auxt = time()

print(time()-start)
ax3d.scatter(np.reshape(sis[:,:-2],n),np.reshape(sis[:,1:-1],n),np.reshape(sis[:,2:],n),c='r',marker='o')

plt.show()
