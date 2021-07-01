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
    
    def insert(self, pos, mass, bounds, sis):
        
        flag = np.argwhere(np.all(sis > bounds[0], axis=1) & np.all(sis < bounds[1], axis=1)).shape[0]
        if flag == 0:
            return
        
        self.mass += mass
        self.center += mass*pos

        if flag == 1:
            return

        middle = (bounds[1] - bounds[0]) / 2
        aux = np.concatenate(([bounds[0]], [bounds[0] + middle]), axis=0) 
        aux2 = np.concatenate(([bounds[0] + middle], [bounds[1]]), axis=0)

        if pos[0] > aux[1,0]:
            if pos[1] > aux[1,1]:
                if pos[2] > aux[1,2]:

                    self.node1 = Node()
                    self.node1.insert(pos, mass, aux2, sis)
                    return
                
                else:

                    self.node2 = Node()
                    self.node2.insert(pos, mass, np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1), sis)
                    return
                   
            elif pos[2] > aux[1,2]:

                self.node3 = Node()
                self.node3.insert(pos, mass, np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1), sis)
                return
            
            else:

                self.node4 = Node()
                self.node4.insert(pos, mass, np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1), sis)
                return
        
        elif pos[1] > aux[1,1]:
            if pos[2] > aux[1,2]:
           
                self.node5 = Node()
                self.node5.insert(pos, mass, np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1), sis)
                return
            
            else:
             
                self.node6 = Node()
                self.node6.insert(pos, mass, np.concatenate((aux[:,:-2], aux2[:,1:-1], aux[:,-1:]), axis=1), sis)
                return
        
        elif pos[2] > aux[1,2]:
        
            self.node7 = Node()
            self.node7.insert(pos, mass, np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1), sis)
            return
                    
        else:
      
            self.node8 = Node()
            self.node8.insert(pos, mass, aux, sis)
            return

    def build_tree(self, sis, masses, bounds):
        for pos, mass in zip(sis, masses):
            self.insert(pos, mass, bounds, sis)

    def calc_forces(self):
        if self == None: 
            return

        aux = np.concatenate(([bounds[0]], [bounds[1]/2]), axis=0) 
        aux2 = np.concatenate(([bounds[1] / 2], [bounds[1]]), axis=0) 
        dist = np.sqrt(np.sum((pos - self.center / self.mass)**2))
        leng = abs(bounds[0][0] - bounds[1][0])
        
        if leng / dist < 1:
            accel += (self.mass * (pos - center)) / (dist**3 )
        
        else:
            self.node1.calc_forces(pos, mass, aux2, sis)
            self.node2.calc_forces(pos, mass, np.concatenate((aux2[:,:-1], aux[:,-1:]), axis=1), sis)             
            self.node3.calc_forces(pos, mass, np.concatenate((aux2[:,:-2], aux[:,1:-1], aux2[:,-1:]), axis=1), sis)
            self.node4.calc_forces(pos, mass, np.concatenate((aux2[:,:-2], aux[:,-2:]), axis=1), sis)
            self.node5.calc_forces(pos, mass, np.concatenate((aux[:,:-2], aux2[:,-2:]), axis=1), sis)
            self.node6.calc_forces(pos, mass, np.concatenate((aux[:,:-2], aux2[:,1:-1], aux[:,-1:]), axis=1), sis)
            self.node7.calc_forces(pos, mass, np.concatenate((aux[:,:-1], aux2[:,-1:]), axis=1), sis)
            self.node8.calc_forces(pos, mass, aux, sis)

def calc_bounds(sis):
    return np.array([np.array(3*[np.amin(sis)-10]), np.array(3*[np.amax(sis)+10])])

root = Node()
n = 10
nNonBar = 9
sis = np.random.rand(n, 3)
bounds = calc_bounds(sis)
v = np.random.rand(n, 3)
barMass = 1
nonBarMass = 1
masses = np.concatenate((nNonBar*[nonBarMass], (n-nNonBar)*[barMass]), axis=None)
root.build_tree(sis, masses, bounds)
print(root.mass)
print(root.center/root.mass)

