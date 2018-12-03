from pdbParser import atomClass, connectClass
import glob
import pdb


if __name__ == "__main__":
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
		print(connectList)
		print("END")