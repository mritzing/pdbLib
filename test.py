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

	distances = []
	for key in startAtoms:
		distances.append(startAtoms[key].getDist(endAtoms[key]))

#	rmsdArr = []

#	startConnects = startChain.getConnects()
#	endConnects = endChain.getConnects()

#	for inx, item in enumerate(startConnects):
#		rmsdArr.append(item.connectRMSD(endConnects[inx]))
	posStart = startChain.getPosArray()
	xPosStart = posStart[:,0]
	yPosStart = posStart[:,1]
	zPosStart = posStart[:,2]

	pdb.set_trace()