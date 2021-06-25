import numpy as np
arr=np.random.rand(1,3)
arr2=np.random.rand(1,3)
arr4=np.concatenate((arr,arr2),axis=0)
arr5=np.concatenate((arr2,arr),axis=0)
print(arr4)
print(arr5)
arr3=np.concatenate((arr4.T[:][:-1],arr5.T[:][-1:]), axis=0).T
print(arr3)
