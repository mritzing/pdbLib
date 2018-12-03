import glob
import pdb
import numpy as np



class atomClass: 
	""" Class used to store information about individual atoms 
		Constructor values
		param idNum: ID number
		param x: X coordinate
		param y: Y coordinate
		param z: Z coordinate
		param temp: temperature
		param molType: Molecule type (C , N , O, etc..)
	"""
	def __init__(self, idNum, x, y, z, temp ,molType):
		self.idNum = idNum
		self.x=  x
		self.y = y
		self.z = z
		self.temp = temp
		self.molType = molType

	def getID(self):
		return self.idNum


	def getAtomLoc(self):
		""" Returns [x,y,z] list of atoms coordinates """
		return [double(self.x), double(self.y), double(self.z)]



class connectClass: 
	""" Class used to store information about conect lines in PDB file 
		atomDict = atomClass objects from the ID numbers in connect line
		ex  : CONECT    4    3    5    7 
		atomDict would contain key and associated atom items of values 4, 3, 5, 7 
	"""

	def __init__(self, atomDict):
		"""Dictionary of atomClass objects"""
		self.atomDict = atomDict

	def getDict(self):
		return self.atomDict


	def getLocationArray(self):
		""" Returns numpy array of positions of atoms within object 
			return ex: 
			np.array([[46.900, -34.882, 16.524],
					  [47.305, -39.387,  13.282],
					  [49.426, -32.298, 10.930]])  
		"""
		arr = np.array()
		for key, atom in enumerate(self.atomDict):
			arr.append(atom.getLoc)
		return arr


