import numpy as np

dim=64
dimx=64
dimy=64
dimz=64
tagname='bxyz'
f    = open('./correct_results/B_0002_0032.dat','rb')
data = np.fromfile(f, dtype=np.double, count=dim*dim*dim*3)
f.close()

bxyz = np.transpose(np.reshape(data,(3,dimz,dimy,dimx)),(3,2,1,0))
bxyz = bxyz.astype(np.float)

f = open("./correct_results/bxyz.vtk", "w+")
f.write('# vtk DataFile Version 2.0\n')
f.write('Volume example\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_POINTS\n')
f.write("DIMENSIONS "+"{0:5d}".format(int(dimx))+' '+"{0:5d}".format(int(dimy))+' '+"{0:5d}".format(int(dimz))+"\n")
f.write("ASPECT_RATIO "+"{:.1f}".format(1.0)+' '+"{:.1f}".format(1.0)+' '+"{:.1f}".format(1.0)+"\n")
f.write("ORIGIN "+"{:1d}".format(0)+' '+"{:1d}".format(0)+' '+"{:1d}".format(0)+"\n")
f.write("POINT_DATA "+"{0:10d}".format(int(dimx*dimy*dimz))+"\n")
f.write("VECTORS "+tagname+" FLOAT\n")

for k in range(dimz):
	for j in range(dimy):
		for i in range(dimx):
			string="{:.3f}".format(bxyz[i,j,k,0])+' '+"{:.3f}".format(bxyz[i,j,k,1])+' '+"{:.3f}".format(bxyz[i,j,k,2])+"\n"
			f.write(string)
f.close()

# f.write(bytearray(bxyz.reshape([int(dimx*dimy*dimz),3])))
