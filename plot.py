import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize = (10,8))
ax = plt.axes(projection='3d')

t = np.arange(0, 101, 1)
x = np.arange(0, 121, 1)

T, X = np.meshgrid(t, x)

Z = np.genfromtxt("build/output.csv", delimiter=",", dtype=np.double)
print(Z.shape)

surf = ax.plot_surface(T, X, Z, cmap=plt.cm.cividis)
#surf = ax.plot_surface(T, X, Z, cmap=plt.cm.coolwarm)

ax.set_xlabel('t', labelpad = 20)
ax.set_ylabel('x', labelpad = 20)
ax.set_zlabel('z', labelpad = 20)

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()