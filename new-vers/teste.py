import numpy as np
arr=np.random.rand(1,3)
arr2=np.random.rand(1,3)
arr4=np.concatenate((arr,arr2),axis=0)
arr5=np.concatenate((arr2,arr),axis=0)
#arr3=np.concatenate((arr4[:,:-1],arr5[:,-1:]), axis=1)
print(arr4)
print(np.argwhere(np.all(arr4>[0,0,0],axis=1)&np.all(arr4<[0.7,0.7,0.7],axis=1)).shape[0])
