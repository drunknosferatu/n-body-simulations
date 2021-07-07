import numpy as np
arr=np.array([[2, 2]])
arr2=np.array([[1, 2]])
arr4=[None, None, None]
arr5=np.concatenate((arr2,arr2,arr),axis=0)
arr2=[1,1]
arr3=[3,3]
#arr3=np.concatenate((arr4[:,:-1],arr5[:,-1:]), axis=1)
print(arr5)
print(np.any(arr5>0, axis=1))
