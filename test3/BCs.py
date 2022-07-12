import numpy as np

def xderiv2(arr,dx=None):
	if dx==None:
		dx=np.divide(1.0,arr.shape[0]-1)
	nx=arr.shape[0]
	ny=arr.shape[1]
	hidx=np.divide(0.5,dx)
	darr=hidx*(np.roll(arr,-1,axis=0)-np.roll(arr,+1,axis=0))
	darr[0,:]=hidx*(-3.*arr[0,:]+4.*arr[1,:]-arr[2,:])
	darr[nx-1,:]=hidx*(3.*arr[nx-1,:]-4.*arr[nx-2,:]+arr[nx-3,:])
	return darr

def yderiv2(arr,dy=None):
	if dy==None:
		dy=np.divide(1.0,arr.shape[1]-1)	
	nx=arr.shape[0]
	ny=arr.shape[1]
	hidy=np.divide(0.5,dy)
	darr=hidy*(np.roll(arr,-1,axis=1)-np.roll(arr,1,axis=1))
	darr[:,0]=hidy*(-3.*arr[:,0]+4.*arr[:,1]-arr[:,2])
	darr[:,ny-1]=hidy*(3.*arr[:,ny-1]-4.*arr[:,ny-2]+arr[:,ny-3])
	return darr

dim=64
x = np.divide(np.linspace(0, dim-1, dim),dim)
y = np.divide(np.linspace(0, dim-1, dim),dim)
yv, xv = np.meshgrid(x, y)

kk = 2*np.pi
lambda0 = 1.0
ll = np.sqrt(np.power(kk,2) - np.power(lambda0,2))

bx = ll*np.sin(kk*xv)*np.cosh(ll*(1.))
by = lambda0*np.sin(kk*xv)*np.sinh(ll*(1.))
bz = kk*np.cos(kk*xv)*np.sinh(ll*(1.))

dxby = xderiv2(by)/dim
dybx = yderiv2(bx)/dim

jz = dxby - dybx
alpha = (dim-1)*np.divide(jz,bz)
alpha = alpha*0+np.median(alpha)
sig=alpha*0.01
bx = np.divide(bx,np.abs(bz).max())
by = np.divide(by,np.abs(bz).max())
bz = np.divide(bz,np.abs(bz).max())

np.transpose(bz).astype('float32').tofile('./data/bz0.dat')
np.transpose(alpha).astype('float32').tofile('./data/alpha0.dat')
np.transpose(sig).astype('float32').tofile('./data/sig0.dat')



