import pandas as pd
import numpy as np
try:
	import f90nml
except:
	print('Python package f90nml is needed, please install it.')

import glob,os
try:
	from pyevtk.hl import gridToVTK
except:
	print("For 3D visualization, the python package pyevtk.hl is needed to write a VTK format file, please check and install it.")


class cfitio():
	def __init__(self,par=None,quiet=False):
		nml    = f90nml.read(par)
		self.resultDir= nml['filename_par']['OutFileName']
		self.dimx = nml['cal_par']['dimx']
		self.dimy = nml['cal_par']['dimy']
		self.dimz = nml['cal_par']['dimz']
		self.ncycle= nml['cal_par']['ncycle']
		self.nloop= nml['cal_par']['nloop']
		self.savNum= nml['cal_par']['savNum']
		self.totalDim=self.dimx*self.dimy*self.dimz
		self.quiet=quiet
		self.x = np.divide(np.linspace(0, self.dimx-1, self.dimz),self.dimx)
		self.y = np.divide(np.linspace(0, self.dimy-1, self.dimy),self.dimx)
		self.z = np.divide(np.linspace(0, self.dimz-1, self.dimz),self.dimx)
		self.vectorMessage='In the 3D vector data array, axis=0 is the x dir, axis=1 is the y dir, axis=2 is the z dir, axis=3 is the 3 components of the vector field.'
		self.scalarMessage='In the 3D scalar data array, axis=0 is the x dir, axis=1 is the y dir, axis=2 is the z dir the scalar field.'

		# self.xv, self.yv, self.zv = np.meshgrid(x, y, z, indexing='ij')
	def log(self):
		filename=os.path.join(self.resultDir,'cfit_log.csv')
		try:
			log=pd.read_csv(filename)
		except:
			print('File '+filename+' not found, please check path and file name.')
			log=0
		return log
	def B3D(self,cycle=0,loop=0):
		filename=os.path.join(self.resultDir,'B_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
		try:
			if ~quiet:
				print(self.vectorMessage)
			with open(filename,'rb') as f:
				data = np.fromfile(f, dtype=np.double, count=self.totalDim*3)
			bxyz = np.transpose(np.reshape(data,(3,self.dimz,self.dimy,self.dimx)),(3,2,1,0))
		except:
			print('File '+filename+' not found, please check the cycle and loop number.')
			bxyz=0
		return bxyz
	def Bc3D(self,cycle=0,loop=0):
		filename=os.path.join(self.resultDir,'Bc_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
		try:
			if ~quiet:
				print(self.vectorMessage)
			with open(filename,'rb') as f:
				data = np.fromfile(f, dtype=np.double, count=self.totalDim*3)
			bc = np.transpose(np.reshape(data,(3,self.dimz,self.dimy,self.dimx)),(3,2,1,0))
		except:
			print('File '+filename+' not found, please check the cycle and loop number.')
			bc=0
		return bc
	def J3D(self,cycle=0,loop=0):
		filename=os.path.join(self.resultDir,'Jc_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
		try:
			if ~quiet:
				print(self.vectorMessage)
			with open(filename,'rb') as f:
				data = np.fromfile(f, dtype=np.double, count=self.totalDim*3)
			jxyz = np.transpose(np.reshape(data,(3,self.dimz,self.dimy,self.dimx)),(3,2,1,0))
		except:
			print('File '+filename+' not found, please check the cycle and loop number.')
			jxyz=0
		return jxyz
	def A3D(self,cycle=0,loop=0):
		if (~self.quiet):
			print('A reminder that the file A_0000_0000.dat contains the vector potential for the initial potential field.')
			print('Other A files contain only the vector potential for the current carrying field Bc.')
		filename=os.path.join(self.resultDir,'A_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
		try:
			if ~quiet:
				print(self.vectorMessage)
			with open(filename,'rb') as f:
				data = np.fromfile(f, dtype=np.double, count=self.totalDim*3)
			Axyz = np.transpose(np.reshape(data,(3,self.dimz,self.dimy,self.dimx)),(3,2,1,0))
		except:
			print('File '+filename+' not found, please check the cycle and loop number.')
			Axyz=0
		return Axyz
	def Alpha3D(self,cycle=0,loop=0):
		filename=os.path.join(self.resultDir,'Alpha3d_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
		try:
			if ~quiet:
				print(self.scalarMessage)
			with open(filename,'rb') as f:
				data = np.fromfile(f, dtype=np.double, count=self.totalDim)
			alpha = np.transpose(np.reshape(data,(self.dimz,self.dimy,self.dimx)),(2,1,0))
		except:
			print('File '+filename+' not found, please check the cycle and loop number.')
			alpha=0
		return alpha
	def B3D_toVtk(self,cycle=0,loop=0):
		try:
			if(~self.quiet):
				print("Vtr File will be stored in "+self.resultDir)
			filename=os.path.join(self.resultDir,'B_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
			bxyz=self.B3D(cycle=cycle,loop=loop)
			vtkname=filename.replace('.dat','')
			gridToVTK(vtkname, self.x, self.y, self.z, pointData = {"B" : (bxyz[...,0],bxyz[...,1],bxyz[...,2])})
		except:
			print('Cannot read the .dat file or write the .vtr file.')
	def J3D_toVtk(self,cycle=0,loop=0):
		try:
			if(~self.quiet):
				print("Vtr File will be stored in "+self.resultDir)
			filename=os.path.join(self.resultDir,'Jc_'+str(cycle).zfill(4)+'_'+str(loop).zfill(4)+'.dat')
			jxyz=self.J3D(cycle=cycle,loop=loop)
			vtkname=filename.replace('.dat','')
			gridToVTK(vtkname, self.x, self.y, self.z, pointData = {"J" : (bxyz[...,0],bxyz[...,1],bxyz[...,2])})
		except:
			print('Cannot read the .dat file or write the .vtr file.')