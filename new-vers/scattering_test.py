import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import time
#check for negative numbers
class Node:

    def __init__(self, masses, sis, charges, leng):
        self.node1 = None
        self.node2 = None
        self.node3 = None
        self.node4 = None
        self.node5 = None
        self.node6 = None
        self.node7 = None
        self.node8 = None
        self.mass = np.sum(masses)
        self.charge = np.sum(charges)
        self.leng = leng
        self.center = np.sum(masses[:,None] * sis, axis=0)
    

    


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
                setattr(self, call_strings[count], Node(masses[index], sis[index], charges[index], abs(new_bounds[count][0,0]-new_bounds[count][1,0])))
                getattr(self, call_strings[count]).insert(sis[index], masses[index], new_bounds[count], charges[index])

    



    def calc_forces(self, sis, charges, masses, accel):
        call_strings = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8"]
        aux_dist =  self.center / self. mass - sis
        aux_dist = np.around(aux_dist, decimals=15)
        ind_exclud = np.argwhere(np.any((aux_dist) != 0, axis=1)).flatten()
        dist = np.sqrt(np.sum((aux_dist)**2, axis=1))
        dist = np.around(dist, decimals=15)
        print(dist)
        if ind_exclud.shape[0] != sis.shape[0]:
            oriented_dist = aux_dist[ind_exclud] / (dist[ind_exclud]**3)[:,None]
            accel[ind_exclud] += 0*self.mass * oriented_dist - 9*1e9 * charges[ind_exclud,None] * self.charge * oriented_dist / masses[ind_exclud,None]
        else:
            oriented_dist = aux_dist / (dist**3)[:,None]
            ind_theta = np.argwhere(self.leng / dist < 0.5).flatten()
            if ind_theta.shape[0] > 0:
                accel[ind_theta] += 0*self.mass * oriented_dist[ind_theta] - 9*1e9 * charges[ind_theta,None] * self.charge * oriented_dist[ind_theta] / masses[ind_theta,None]
            ind_theta = np.argwhere(self.leng / dist > 0.5).flatten()
            if ind_theta.shape[0] > 0:
                for call in call_strings:
                     if getattr(self,call):
                        accel[ind_theta] = getattr(self, call).calc_forces(sis[ind_theta], charges[ind_theta], masses[ind_theta], accel[ind_theta])
        return accel




def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis) - 0.1]), np.array(3*[np.amax(sis) + 0.1])])


fig = plt.figure()
fig2 = plt.figure()
ax = fig.add_subplot()
ax2 = fig2.add_subplot()
xplot = []
yplot = []
xplot2 = []
yplot2 = []
tplot = []
dplot = []

t = 0.00005
control = None
n = 2
nBar = 2
BarMass = 1.6605*1e-27
nonBarMass = 1.6605*1e-15
charge = 1.6*1e-19


sis = np.array([[1.0,0.2,0],[0,0.3,0]])
v = np.array([[-1.6,0,0],[0,0,0]])
charges = np.concatenate((nBar * [charge], np.zeros((n-nBar))), axis=None)
masses = np.concatenate((nBar * [BarMass], (n-nBar) * [nonBarMass]), axis=None)

j=0
start = time()
while not control:
    accel = np.zeros((n, 3))
    bounds = calc_bounds(sis)
    root = Node(masses, sis, charges, abs(bounds[0,0]-bounds[1,0]))
    root.insert(sis, masses, bounds, charges) 
    xplot.append(sis[0,0])
    yplot.append(sis[0,1])
    xplot2.append(sis[1,0])
    yplot2.append(sis[1,1])
    root.calc_forces(sis, charges, masses, accel)
    v += accel *  t / 2
    sis += v * t
    v += accel * t / 2
    anorm = np.sqrt(np.sum(accel**2, axis=1))
    control = np.all(anorm <= 67*1e-3)
    j+=1
    t = 5*1e-5 / np.amax(anorm)
print(time()-start)
ax.plot(xplot,yplot)
ax.plot(xplot2,yplot2)
plt.show()
