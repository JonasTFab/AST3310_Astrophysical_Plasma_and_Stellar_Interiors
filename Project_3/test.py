import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


x = 300
y = 100
x_len = np.linspace(1,x,x)
y_len = np.linspace(1,y,y)

mat = np.zeros((y,x))
zx = signal.gaussian(x,std=int(x/10))
zy = signal.gaussian(y,std=int(y/10))

for i in range(x):
    for j in range(y):
        mat[j,i] = zx[i]*zy[j]


X,Y = np.meshgrid(x_len,y_len)


#plt.contourf(X,Y,mat)
#plt.colorbar()
#plt.show()


mat_test = np.zeros((12,4))
for i in range(len(mat_test[:,0])):
    for j in range(len(mat_test [0,:])):
        mat_test[i,j] = np.random.randint(1,6)
print(mat_test, "\n")

mat_test=np.roll(mat_test, 1, axis=0)
print(mat_test)
