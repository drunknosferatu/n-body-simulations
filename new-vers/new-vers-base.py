import numpy as np

class Node:

    def __init__(self, center, mass):
        self.node1 = None
        self.node2 = None
        self.node3 = None
        self.node4 = None
        self.node5 = None
        self.node6 = None
        self.node7 = None
        self.node8 = None
        self.mass = mass
        self.center = center
    
    def insert(self, center, mass, bounds):
        aux = np.concatenate(([bounds[0]], [bounds[1]/2]), axis=0) 
        aux2 = np.concatenate(([bounds[1]/2], [bounds[1]]), axis=0) 
        if center[0] > bounds[1][0]/2:
            if center[1] > bounds[1][1]/2:
                if center[2] > bounds[1][2]/2:
                    if self.node1 == None:
                        self.node1 = Node(center, mass)
                    else:
                        self.node1.insert(center, mass, aux2)
                else:
                     if self.node2 == None:
                        self.node2 = Node(center, mass)
                    else:
                        self.node2.insert(center, mass, np.concatenate((aux2.T[:][ :-1], aux.T[:][-1:])).T
                   
            elif center[2] > bounds[1][2]/2:
                    if self.node3 == None:
                        self.node3 = Node(center, mass)
                    else:
                        self.node3.insert(center, mass, np.concatenate())
            else:
                    if self.node4 == None:
                        self.node4 = Node(center)
                    else:
                        self.node4.insert(center, mass, np.concatenate())
        elif center[1] > bounds[1][1]/2:
            if center[2] > bounds[1][2]:
                    if self.node5 == None:
                        self.node5 = Node(center, mass)
                    else:
                        self.node5.insert(center, mass, np.concatenate())
            else:
                    if self.node6 == None:
                        self.node6 = Node(center, mass)
                    else:
                        self.node6.insert(center, mass, np.concatenate())
        elif center[2] > bounds[1][2]/2:
                    if self.node7 == None:
                        self.node7 = Node(center, mass)
                    else:
                        self.node7.insert(center, mass, np.concatenate())
        else:
                    if self.node8 == None:
                        self.node8 = Node(center, mass)
                    else:
                        self.node8.insert(center, mass, np.concatenate())
n = 1000
nNonBar = 900
bounds = np.array()
sis = ((7/4)**1/3)*np.random.rand(n, 3)
v = np.random.rand(n, 3)
barMass = 1
nonBarMass = 1
k = 1
G = 1
masses = np.concatenate((nNonBar*[nonBarMass], (n-nNonBar)*[barMass]), axis=None)

