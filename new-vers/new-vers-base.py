import numpy as np


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
        self.center = 0
    
    def cond_checking_support(self, pos, mass, new_bounds, sis, aux, nodes):
        cond = pos > aux[1]
        index = 0
        for count, state in enumerate(reversed(cond)):
            if state == False:
                index += 2**count
        nodes[index] = Node()
        nodes[index].insert(pos, mass, new_bounds[index], sis)


    def insert(self, pos, mass, bounds, sis):
        
        self.mass += mass
        self.center += mass*pos

        if self.mass == mass:
            return
        nodes = [self.node1, self.node2, self.node3, self.node4, self.node5, self.node6, self.node7, self.node8]

        middle = (bounds[1] - bounds[0]) / 2
        aux = np.concatenate(([bounds[0]], [bounds[0] + middle]), axis=0) 
        aux2 = np.concatenate(([bounds[0] + middle], [bounds[1]]), axis=0)
        index = 0
        new_bounds = [aux2,
                np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1),
                np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1),
                np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1),
                np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1),
                np.concatenate((aux[:,:-2],aux2[:,1:-1], aux[:,-1:]), axis=1),
                np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1),
                aux]
        
        if not np.any(nodes):
            self.cond_checking_support((self.center - mass * pos) / (self.mass - mass), self.mass - mass, new_bounds, sis, aux, nodes)
        self.cond_checking_support(pos, mass, new_bounds, sis, aux, nodes)

    def build_tree(self, sis, masses, bounds):
        i=0
        for pos, mass in zip(sis, masses):
            self.insert(pos, mass, bounds, sis)

    def calc_forces(self, pos, bounds, accel):

        aux = np.concatenate(([bounds[0]], [bounds[1]/2]), axis=0) 
        aux2 = np.concatenate(([bounds[1] / 2], [bounds[1]]), axis=0) 
        dist = np.sqrt(np.sum((pos - self.center / self.mass)**2))
        leng = abs(bounds[0][0] - bounds[1][0])
        
        if leng / dist < 1:
            accel += (self.mass * (pos - self.center/self.mass)) / (dist**3 )
        
        else:
            if self.node1:
                if self.node1.mass != 0:
                    self.node1.calc_forces(pos, aux2, accel)
            if self.node2:
                if self.node2.mass != 0:
                    self.node2.calc_forces(pos, np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1), accel)             
            if self.node3:
                if self.node3.mass != 0:
                    self.node3.calc_forces(pos, np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1), accel)
            if self.node4:
                if self.node4.mass != 0:
                    self.node4.calc_forces(pos, np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1), accel)
            if self.node5:
                if self.node5.mass != 0:
                    self.node5.calc_forces(pos, np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1), accel)
            if self.node6:
                if self.node6.mass != 0:
                    self.node6.calc_forces(pos, np.concatenate((aux[:,:-2], aux2[:,1:-1], aux[:,-1:]), axis=1), accel)
            if self.node7:
                if self.node7.mass != 0:
                    self.node7.calc_forces(pos, np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1), accel)
            if self.node8:
                if self.node8.mass != 0:
                    self.node8.calc_forces(pos, aux, accel)

def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis)-2]), np.array(3*[np.amax(sis)+2])])

root = Node()
n = 100000
nNonBar = 9
sis = np.random.rand(n, 3)*10
bounds = calc_bounds(sis)
v = np.random.rand(3)
barMass = 1
nonBarMass = 1
masses = np.concatenate((nNonBar*[nonBarMass], (n-nNonBar)*[barMass]), axis=None)
root.build_tree(sis, masses, bounds)
print("finished")
accel = np.random.rand(n, 3)
for i in range(n):
    root.calc_forces(sis[i], bounds, accel[i])

