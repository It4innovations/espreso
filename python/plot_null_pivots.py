import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


#def randrange(n, vmin, vmax):
#    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure(0)
plt.clf()
#ax = fig.add_subplot(111, projection='3d')
ax = fig.gca(projection='3d')

ne = np.array([1,1,1])*3

Lx = 30
Ly = 30
Lz = 30
x = np.linspace(0,1,ne[0]+1)*Lx;y = np.linspace(0,1,ne[1]+1)*Ly;z = np.linspace(0,1,ne[2]+1)*Lz

#Z,Y,X = np.meshgrid(z,y,x)
#xx = X.ravel();yy = Y.ravel();zz = Z.ravel()
nPoints = (ne[0]+1)*(ne[1]+1)*(ne[2]+1)
xx=np.zeros(nPoints)
yy=np.zeros(nPoints)
zz=np.zeros(nPoints)
cnt=0
for k in range(ne[2]+1):
    for j in range(ne[1]+1):
        for i in range(ne[0]+1):
            xx[cnt]=x[i]
            yy[cnt]=y[j]
            zz[cnt]=z[k]
            label = str(cnt)
            ax.text(xx[cnt], yy[cnt], zz[cnt], label,fontsize=8)            
            
            cnt+=1


 


#null_pivots = np.array([12, 39, 115, 131, 134, 192])-1
null_pivots = np.array([39, 46, 47, 60, 147, 148])-1   
null_nodes = np.int64(np.floor(null_pivots/3))

x_np = xx[null_nodes];y_np = yy[null_nodes];z_np = zz[null_nodes]
u_np = x_np.copy();v_np = y_np.copy();w_np = z_np.copy()

del1 = max([xx.max()-xx.min(),yy.max()-yy.min(),zz.max()-zz.min()])
len1 = 0.1*del1
arr = null_pivots-3*null_nodes
for i in range(null_pivots.shape[0]):
    if arr[i]==0:
        u_np[i]+=len1
    if arr[i]==1:
        v_np[i]+=len1
    if arr[i]==2:
        w_np[i]+=len1

print(null_pivots-3*null_nodes)

for i in range(null_pivots.shape[0]):
    ax.plot([x_np[i], u_np[i]], [y_np[i],v_np[i]],zs=[z_np[i],w_np[i]],c='r')
 
ax.scatter(xx,yy,zz, c='k', marker='.')
ax.scatter(x_np,y_np,z_np,c='r',marker='o',s=35)
#ax.scatter(u_np,v_np,w_np,c='r',marker='.',s=35)



ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

ax.axis('equal')






plt.show()


#from mpl_toolkits.mplot3d import axes3d
#import matplotlib.pyplot as plt
#import numpy as np
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
#x, y, z = np.meshgrid(np.arange(-0.8, 1, 0.2),
#                      np.arange(-0.8, 1, 0.2),
#                      np.arange(-0.8, 1, 0.8))
#
#u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
#v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
#w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
#     np.sin(np.pi * z))
#
#ax.quiver(x, y, z, u, v, w, length=0.1)
#
#plt.show()
