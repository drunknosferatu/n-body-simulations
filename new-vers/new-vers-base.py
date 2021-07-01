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
    
    def insert_and_build(self, pos, mass, bounds, sis):
        flag = np.argwhere(np.all(sis > bounds[0], axis=1) & np.all(sis < bounds[1], axis=1)).shape[0]
        self.mass += mass
        self.center += mass*pos
        if flag == 1:
            return
        aux = np.concatenate(([bounds[0]], [bounds[1]/2]), axis=0) 
        aux2 = np.concatenate(([bounds[1]/2], [bounds[1]]), axis=0) 
        if pos[0] > bounds[1,0]/2:
            if pos[1] > bounds[1,1]/2:
                if pos[2] > bounds[1,2]/2:
                        self.node1 = Node()
                        self.node1.insert(pos, mass, aux2, sis)
                else:
                        self.node2 = Node()
                        self.node2.insert(pos, mass, np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1), sis)
                   
            elif pos[2] > bounds[1,2]/2:
                        self.node3 = Node()
                        self.node3.insert(pos, mass, np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1), sis)
            else:
                        self.node4 = Node()
                        self.node4.insert(pos, mass, np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1), sis)
        elif pos[1] > bounds[1,1]/2:
            if pos[2] > bounds[1,2]:
                        self.node5 = Node()
                        self.node5.insert(pos, mass, np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1), sis)
            else:               
                        self.node6 = Node()
                        self.node6.insert(pos, mass, np.concatenate((aux[:,:-2], aux2[:,1:-1], aux[:,-1:]), axis=1), sis)
        elif pos[2] > bounds[1,2]/2:
                        self.node7 = Node()
                        self.node7.insert(pos, mass, np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1), sis)
        else:
                        self.node8 = Node()
                        self.node8.insert(pos, mass, aux, sis)


n = 1000
nNonBar = 900
bounds = np.array([1])
sis = ((7/4)**1/3)*np.random.rand(n, 3)
v = np.random.rand(n, 3)
barMass = 1
nonBarMass = 1
k = 1
G = 1
masses = np.concatenate((nNonBar*[nonBarMass], (n-nNonBar)*[barMass]), axis=None)

