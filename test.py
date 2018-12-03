import pdbParser
import glob
import pdb


if __name__ == "__main__":
	files = glob.glob('*nspose1.pdb')
	print (files)
	pdbs = []
	for file in files:
		pdbs.append(pdbParser.pdbParser(file))
	#Find max distance
	startChain = pdbs[0].getChains()[0]
	startAtoms = startChain.getAtoms()
	endChain = pdbs[-1].getChains()[0]
	endAtoms = endChain.getAtoms()
	maxDist = 0

	posArray = endChain.getPosArray()
	xPos = posArray[:,0]
	yPos = posArray[:,1]
	zPos = posArray[:,2]


	for key in startAtoms:
		distance = startAtoms[key].getDist(endAtoms[key])
		if distance > maxDist:
			maxDist = distance
			maxKey = key
	print(maxDist)