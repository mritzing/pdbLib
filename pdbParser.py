import glob
import pdb
import numpy as np

class atomClass: 
	def __init__(self, idNum, x, y, z, temp ,molType):
		self.idNum = idNum
		self.x= x
		self.y = y
		self.z = z
		self.temp = temp
		self.molType = molType

	def getID(self):
		return self.idNum

	#use numpy array to get rmsd
	#https://github.com/charnley/rmsd/blob/master/example.py
	def getLoc(self):
		return np.array[double(x), double(y), double(z)]



class connectClass: 
	def __init__(self, atomList):
		self.atomList = atomList

	def getList():
		return atomList

atomDict = {}
connectList = []
files = glob.glob('*nspose1.pdb')
for line in open(files[0], 'r'):
	lineArr = line.split()
	if lineArr[0] == "HETATM":
		atomDict[int(lineArr[1])] = atomClass(lineArr[1], lineArr[6], lineArr[7], lineArr[8], lineArr[10], lineArr[-1])
	elif lineArr[0] == "CONECT":
		connectDict = {};
		for atomNum in lineArr[1:]:
			connectDict[int(atomNum)] = atomDict[int(atomNum)]
		connectList.append(connectClass(connectDict))
else :
	print(len(atomDict))
	print(len(connectList))
	print("END")