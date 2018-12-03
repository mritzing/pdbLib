import numpy as np
from atomicMass import atomic_mass
from collections import Counter
import rmsd
from itertools import groupby, chain
import pdb
""" Classes useful in parsing pdb files down to individual components """
#big TODO make tests
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
		return [float(self.x), float(self.y), float(self.z)]

	# https://stackoverflow.com/questions/20184992/finding-3d-distances-using-an-inbuilt-function-in-python
	def getDist(self, atom_):
		""" Returns distance between two points
			param atom_ : atomClass object being compared
		"""
		p1 = np.array(self.getAtomLoc())
		p2 = np.array(atom_.getAtomLoc())
		dist = np.linalg.norm(p1- p2)
		return (dist)

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
		self.arr = None

	def getDict(self):
		return self.atomDict


	def getConnectPos(self):
		""" Returns numpy array of positions of atoms within object 
			return ex: 
			np.array([[46.900, -34.882, 16.524],
					  [47.305, -39.387,  13.282],
					  [49.426, -32.298, 10.930]])  
		"""
		if self.arr is None:
			posList = [];
			for key in self.atomDict:
				posList.append((self.atomDict[key].getAtomLoc()))
			self.arr = np.vstack(posList)
		return self.arr


	def connectRMSD(self, connect_):
		return(rmsd.rmsd(self.getConnectPos(), connect_.getConnectPos()))

class chainClass: 
	""" Chain stores atomDict and connectList of elements between END markers in file
	"""
	def __init__(self, text):
		self.atomDict = {}
		self.connectList = []
		self.atomCount = Counter()
		self.text = text
		self.processLines()
		self.arr = None


	def processLines(self):
		""" Main function for processing files
			Stores associated items atomDict and connectList
			Counts number of elements to atomCount object
		"""
		for line in self.text.split('\n'):
			lineArr = line.split()
			if lineArr != []:
				if lineArr[0] == "HETATM":
					self.atomDict[int(lineArr[1])] = atomClass(lineArr[1], lineArr[6], lineArr[7], lineArr[8], lineArr[10], lineArr[-1])
					self.atomCount[lineArr[-1]] += 1
				elif lineArr[0] == "CONECT":
					connectDict = {};
					for atomNum in lineArr[1:]:
						connectDict[int(atomNum)] = self.atomDict[int(atomNum)]
					self.connectList.append(connectClass(connectDict))

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
		if self.arr is None:
			posList = [];
			for key in self.atomDict:
				posList.append((self.atomDict[key].getAtomLoc()))
			self.arr = np.vstack(posList)
		return self.arr


	def getAtoms(self):
		return self.atomDict

	def getConnects(self):
		return self.connectList

	def chainRMSD(self, chain_):
		""" Calculates RMSD using https://github.com/charnley/rmsd
			param chain_ : chain being compared
		"""
		return(rmsd.rmsd(self.getPosArray(), chain_.getPosArray()))





class pdbParser():
	def __init__(self, fileName):
		self.chainList=[]
		self.fileName = fileName
		self.process()

	#TODO check if multiple chains can get processed
	def process(self):
		""" Main function for processing files
			Stores associated items atomDict and connectList
			Counts number of elements to atomCount object
		"""
		file = open(self.fileName,'r').read().split('END')
		for textBlock in file:
			self.chainList.append(chainClass(textBlock))

	def getChains(self):
		return(self.chainList)

	def getFileName(self):
		return(self.fileName)

	# Mainly a testing function # 
	def getTotalAtoms(self):
		total = 0
		for chain in self.chainList:
			total = total + len(chain.getAtoms())
		return total

# testing function
if __name__ == "__main__":
	print("Testing")