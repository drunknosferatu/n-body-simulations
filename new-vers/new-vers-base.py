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
        aux_dist = sis - self.center / self. mass
        ind_exclud = np.argwhere(np.any((aux_dist) != 0, axis=1)).flatten()       
        dist = np.sqrt(np.sum((aux_dist)**2, axis=1) + 1e-4)

        ind_theta = np.argwhere(self.leng / dist[ind_exclud] < 1).flatten()
        oriented_dist = aux_dist[ind_exclud] / (dist[ind_exclud]**3)[:,None]
        
        if ind_theta.shape[0] > 0:
            accel[ind_theta] +=  self.mass * oriented_dist[ind_theta] - 1e6 * charges[ind_theta,None] * self.charge * oriented_dist[ind_theta] / masses[ind_theta,None]
        ind_theta = np.argwhere(self.leng / dist[ind_exclud] > 1).flatten()
        if ind_theta.shape[0] > 0:
            for call in call_strings:
                if getattr(self,call):
                    getattr(self, call).calc_forces(sis[ind_theta], charges[ind_theta], masses[ind_theta], accel)


    def print_tree(self):
        print(self.mass)
        if self.node1:
            print("1")
            self.node1.print_tree()
        if self.node2:
            print("2")
            self.node2.print_tree()
        if self.node3:
            print("3")
            self.node3.print_tree()
        if self.node4:
            print("4")
            self.node4.print_tree()
        if self.node5: 
            print("5")
            self.node5.print_tree()
        if self.node6: 
            print("6")
            self.node6.print_tree()
        if self.node7:
            print("7")
            self.node7.print_tree()
        if self.node8:
            print("8")
            self.node8.print_tree()
        print("fim") 




def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis) - 0.1]), np.array(3*[np.amax(sis) + 0.1])])


fig = plt.figure()
ax = fig.add_subplot(projection = '3d')
t = 0.00005
control = None
n = 100000
nBar = 0
BarMass = 1
nonBarMass = 1
charge = 1.6 

sis = np.random.randn(n, 3) * 10
v = np.random.rand(n,3)
charges = np.concatenate((nBar * [charge], np.zeros((n-nBar))), axis=None)
masses = np.concatenate((nBar * [BarMass], (n-nBar) * [nonBarMass]), axis=None)

start = time()
while not control:
    accel = np.zeros((n, 3))
    bounds = calc_bounds(sis)
    root = Node(masses, sis, charges, abs(bounds[0,0]-bounds[1,0]))
    root.insert(sis, masses, bounds, charges) 
    print(t)
    root.calc_forces(sis, charges, masses, accel)
    v += accel *  t / 2
    sis += v * t
    v += accel * t / 2
    control = True #np.any(abs(accel) <= 1e-6)
    t = 0.00005 / np.amax(abs(accel))
#root.print_tree()i
np.lsq_min(
print(time()-start)
ax.scatter(np.reshape(sis[:,:-2],n),np.reshape(sis[:,1:-1],n),np.reshape(sis[:,2:],n),c='r',marker='o')
#plt.show()
