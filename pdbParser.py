import numpy as np
from atomicMass import atomic_mass
from collections import Counter
import rmsd
""" Classes useful in parsing pdb files down to individual components """

class atomClass: 
	""" Class used to store information about individual atoms 
		Constructor values
		param idNum: ID number
		param x: X coordinate
		param y: Y coordinate
		param z: Z coordinate
		param temp: temperature
		param ele: Element type (C , N , O, etc..)
	"""
	def __init__(self, idNum, x, y, z, temp ,ele):
		self.idNum = idNum
		self.x=  x
		self.y = y
		self.z = z
		self.temp = temp
		self.ele = ele

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
		""" Returns angle between 3 points with one end being location of this object
			param centerPoint : atomClass object of the center point
			param endPoint : atomClass object of the second end point
		"""
		ba = self.getAtomLoc() - centerPoint.getAtomLoc()
		bc = endPoint.getAtomLoc() - centerPoint.getAtomLoc()

		cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		angle = np.arccos(cosine_angle)
		return(np.degree(angle))

	def get_mol_mass(self):
		return atomicMass.atomic_mass[self.ele]

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

class chainClass: 
	""" Chain stores atomDict and connectList of elements between END markers in file
	"""
	def __init__(self):
		self.atomDict = {}
		self.connectList = []
		self.atomCount = Counter()

	#TODO move down to actual parser class once I start tackling files w/ more than one chain
	def process(self, fileName):
		""" Main function for processing files
			Stores associated items atomDict and connectList
			Counts number of elements to atomCount object
		"""
		for line in open(files[0], 'r'):
			lineArr = line.split()
			if lineArr[0] == "HETATM":
				self.atomDict[int(lineArr[1])] = atomClass(lineArr[1], lineArr[6], lineArr[7], lineArr[8], lineArr[10], lineArr[-1])
				self.atomCount[lineArr[-1]] += 1
			elif lineArr[0] == "CONECT":
				self.connectDict = {};
				for atomNum in lineArr[1:]:
					connectDict[int(atomNum)] = self.atomDict[int(atomNum)]
				self.connectList.append(connectClass(self.connectDict))

	def getMakeup(self):
		return(self.atomCount)


	def getMolMass(self):
		""" Returns mass of elements in object
		"""
		totalMass = 0
		for element in self.totals.keys():
			elementMass = atomicMass.atomic_mass[element]*self.totals[element]
			totalMass = totalMass + elementMass
		return totalMass

	#TODO decide what to store as an object
	def getPosArray(self):
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

	def getAtoms(self):
		return atomDict

	def getConnects(self):
		return connectList

	def chainRMSD(self, chain_):
		""" Calculates RMSD using https://github.com/charnley/rmsd
			param chain_ : chain being compared
		"""
		return(rmsd.rmsd(self.getPosArray(), chain_.getPosArray()))

