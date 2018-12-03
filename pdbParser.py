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

	# https://stackoverflow.com/questions/20184992/finding-3d-distances-using-an-inbuilt-function-in-python
	def getDist(self, atom_):
		""" Returns distance between two points
			param atom_ : atomClass object being compared
		"""
		squared_dist = np.sum(self.getAtomLoc()**2 + atom_.getAtomLoc()**2)
		return np.sqrt(squared_dist)

	#https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates
	def getAngle(self, centerPoint, endPoint):
		"""
			Returns angle between 3 points with one end being location of this object
			param centerPoint : atomClass object of the center point
			param endPoint : atomClass object of the second end point
		"""
		ba = self.getAtomLoc() - centerPoint.getAtomLoc()
		bc = endPoint.getAtomLoc() - centerPoint.getAtomLoc()
		
		cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		angle = np.arccos(cosine_angle)
		return(np.degree(angle))


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

