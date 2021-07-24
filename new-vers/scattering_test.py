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
            accel[ind_exclud] += 0*self.mass * oriented_dist - 9*1e9 * charges[ind_exclud,None] * self.charge * oriented_dist / masses[ind_exclud,None]
        else:
            oriented_dist = aux_dist / (dist**3)[:,None]
            ind_theta = np.argwhere(self.leng / dist < 1).flatten()
            if ind_theta.shape[0] > 0:
                accel[ind_theta] += 0*self.mass * oriented_dist[ind_theta] - 9*1e9 * charges[ind_theta,None] * self.charge * oriented_dist[ind_theta] / masses[ind_theta,None]
            ind_theta = np.argwhere(self.leng / dist >= 1).flatten()
            if ind_theta.shape[0] > 0:
                for call in call_strings:
                     if getattr(self,call):
                         accel[ind_theta] = getattr(self, call).calc_forces(sis[ind_theta], charges[ind_theta], masses[ind_theta], accel[ind_theta])
        return accel




def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis) - 0.1]), np.array(3*[np.amax(sis) + 0.1])])


err = []
tplot = []
part_1_plot = []
part_2_plot = []
fig1 = plt.figure()
ax_err = fig1.add_subplot()
fig2 = plt.figure()
ax_mov = fig2.add_subplot()


t = 0.00005
control = None
n = 2
nBar = 2
BarMass = 1.6605*1e-27
nonBarMass = 1.6605*1e-15
charge = 1.6*1e-19

sis = np.array([[1.0,-0.2,0.0],[-5.0,0.0,0.0]])
v = np.array([[0.0,0.0,0.0],[1.4,0.0,0.0]])
charges = np.concatenate((nBar * [charge], np.zeros((n-nBar))), axis=None)
masses = np.concatenate(([BarMass],[BarMass] ), axis=None)
E0 = np.sum(0.5 * BarMass * np.linalg.norm(v, axis = 1)**2) + 9e9 * charge**2 / np.linalg.norm(sis[0] - sis[1])



j=0
start = time()
while not control:
    root = Node()
    accel = np.zeros((n, 3))
    bounds = calc_bounds(sis)
    root.insert(sis, masses, bounds, charges)
    root.calc_forces(sis, charges, masses, accel)
    v += accel *  t / 2
    sis += v * t
    v += accel * t / 2
    tplot.append(time()-start)
    err.append(abs(100*((E0 - np.sum(0.5 * BarMass * np.linalg.norm(v, axis = 1)**2) - 9e9 * charge**2 / np.linalg.norm(sis[0] - sis[1])) / E0)))
    part_1_plot.append(sis[0,:-1].flatten())
    part_2_plot.append(sis[1,:-1].flatten())
    anorm = np.linalg.norm(accel, axis=1)
    control = np.all(anorm <= 1e-7)
    j+=1
    t = 5e-4 / np.amax(anorm) if 5e-4 / np.amax(anorm) <= 0.5 else 0.5

v -= accel*t/2
r = v[1] - (v[0] * masses[0]  + v[1] * masses[1]) / root.mass
r = np.delete(r, 2)
cos = np.dot(r, [1,0])
sin = np.dot(r, [0,1])
tg_real = sis[1,1] / sis[1,0]
tg_calc = sin / (cos + 0.7)
print(100*(tg_calc-tg_real)/tg_real)

part_1_plot = np.array(part_1_plot).T
part_2_plot = np.array(part_2_plot).T
ax_err.set_xlabel('Tempo transcorrido (s)')
ax_err.set_ylabel('Erro percentual de energia (J)')
ax_err.plot(tplot, err)
plt.xlim(-2, 5)
plt.ylim(-2, 2)
ax_mov.set_xlabel('x')
ax_mov.set_ylabel('y')
ax_mov.plot(part_1_plot[0], part_1_plot[1], label = 'Partícula incidente')
ax_mov.plot(part_2_plot[0], part_2_plot[1], label = 'Partícula alvo')
ax_mov.legend()
plt.show()
