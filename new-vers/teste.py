import numpy as np
arr=np.array([1, 2, 3])
arr2=np.array([1, 3, 3])
arr4=[None, None, None]
arr5=np.concatenate((arr2,arr),axis=0)
#arr3=np.concatenate((arr4[:,:-1],arr5[:,-1:]), axis=1)
if not np.any(arr4): 
    print(np.any(arr4))
    print(type(arr4))

