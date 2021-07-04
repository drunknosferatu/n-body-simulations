import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#check for negative numbers
class Node:

    def __init__(self, mass, pos, charge):
        self.node1 = None
        self.node2 = None
        self.node3 = None
        self.node4 = None
        self.node5 = None
        self.node6 = None
        self.node7 = None
        self.node8 = None
        self.mass = mass
        self.charge = charge
        self.center = mass*pos
    
    def cond_checking_support(self, pos, aux):
        cond = pos > aux[1]
        index = 0
        for count, state in enumerate(reversed(cond)):
            if state == False:
                index += 2**count
        return index


    def insert(self, pos, mass, bounds, charge):
        
        self.mass += mass
        self.center += mass * pos
        self.charge += charge
        call_strings = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8"]
        nodes = [self.node1, self.node2, self.node3, self.node4, self.node5, self.node6, self.node7, self.node8]
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
        if not np.any(nodes):
            index = self.cond_checking_support((self.center - mass * pos) / (self.mass - mass),  aux)
            setattr(self, call_strings[index],  Node(self.mass - mass, (self.center - mass * pos) / (self.mass - mass), self.charge-charge))
        index = self.cond_checking_support(pos, aux)
        if not getattr(self, call_strings[index]):
            setattr(self, call_strings[index], Node(mass, pos, charge))
            return
        else:
            getattr(self, call_strings[index]).insert(pos, mass, new_bounds[index], charge)
            return

    def calc_forces(self, pos, bounds, charge, mass, accel):
        middle = (bounds[1] - bounds[0]) / 2
        aux = np.concatenate(([bounds[0]], [bounds[0] + middle]), axis=0) 
        aux2 = np.concatenate(([bounds[0] +middle], [bounds[1]]), axis=0) 
        dist = np.sqrt(np.sum((pos - self.center / self.mass )**2) + 1e-2)
        if dist == 0:
            return
        leng = abs(bounds[0][0] - bounds[1][0])
        oriented_dist = (pos - self.center/self.mass) / (dist**3)
        new_bounds = [aux2,
                np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1),
                np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1),
                np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1),
                np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1),
                np.concatenate((aux[:,:-2],aux2[:,1:-1], aux[:,-1:]), axis=1),
                np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1),
                aux]
        if leng / dist < 1:
            accel +=  self.mass * oriented_dist + charge * self.charge * oriented_dist * 1e5 / mass
            return
        else:
            if self.node1:
                self.node1.calc_forces(pos, new_bounds[0], charge, mass, accel)
            if self.node2: 
                self.node2.calc_forces(pos, new_bounds[1], charge, mass, accel)
            if self.node3:
                self.node3.calc_forces(pos, new_bounds[2], charge, mass, accel)
            if self.node4: 
                self.node4.calc_forces(pos, new_bounds[3], charge, mass, accel)
            if self.node5: 
                self.node5.calc_forces(pos, new_bounds[4], charge, mass, accel)
            if self.node6: 
                self.node6.calc_forces(pos, new_bounds[5], charge, mass, accel)
            if self.node7: 
                self.node7.calc_forces(pos, new_bounds[6], charge, mass, accel)
            if self.node8: 
                self.node8.calc_forces(pos, new_bounds[7], charge, mass, accel)

    def print_tree(self):
        print(self.center)
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
        
def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis)]), np.array(3*[np.amax(sis)])])



fig = plt.figure()
ax = fig.add_subplot(projection = '3d')

t = 0.005
control = None
n = 100000
nBar = 0
BarMass = 1
nonBarMass = 1
charge = 1.6 


sis = np.random.rand(n, 3)
v = np.random.rand(n,3)
accel = np.zeros((n, 3))
charges = np.concatenate((nBar * [charge], np.zeros((n-nBar))), axis=None)
masses = np.concatenate((nBar * [BarMass], (n-nBar) * [nonBarMass]), axis=None)

while not control:
    root = Node(masses[0], sis[0], charges[0])
    bounds = calc_bounds(sis)
    for i in range(1, n):
        root.insert(sis[i], masses[i], bounds, charges[i])
    
    for i in range(1):
        root.calc_forces(sis[i], bounds, charges[i], masses[i], accel[i])
    v += accel *  t / 2
    sis += v * t
    v += accel * t / 2
    print(t)
    control = np.any(abs(accel) <= 1e-10)
    t = 0.005/ np.amax(abs(accel))
#root.print_tree()
#print(sis)
ax.scatter(np.reshape(sis[:,:-2],n),np.reshape(sis[:,1:-1],n),np.reshape(sis[:,2:],n),c='r',marker='o')
plt.show()
