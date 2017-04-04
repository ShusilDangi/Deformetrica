import csv
import xml.etree.ElementTree as ET
import numpy as np
import sklearn.metrics.pairwise as kernel
import os


def pathToFile(ROOT_DIR,FILENAME):
	filePath = []
	for dirpath, dirnames, filenames in os.walk(ROOT_DIR):
		for filename in filenames:
			if filename.lower() == FILENAME.lower():
				filePath.append(os.path.join(dirpath,filename))
	return filePath


def readCP(fileName):
	CP = []
	with open(fileName, 'rb') as csvfile:
		filereader = csv.reader(csvfile, delimiter=' ', quotechar='|')
		for row in filereader:
			arr = []
			for char in row:
				arr.append(float(char))
			CP.append(arr)
	return np.asarray(CP)


def readMOM(fileName):
	MOM = []
	with open(fileName, 'rb') as csvfile:
		filereader = csv.reader(csvfile, delimiter=' ', quotechar='|')
		for row in filereader:
			arr = []
			for char in row:
				arr.append(float(char))
			MOM.append(arr)
	MOM = np.asarray(MOM)
	return [MOM]


def readMOMAtlas(fileName):
	MOMS = []
	with open(fileName, 'rb') as csvfile:
		filereader = csv.reader(csvfile, delimiter=' ', quotechar='|')
		firstLine = filereader.next()
		nSubjects = int(firstLine[0])
		nCtrlPoints = int(firstLine[1])
		for i in range(nSubjects):
			line = filereader.next()
			MOMSubject = []
			for j in range(nCtrlPoints):
				line = filereader.next()
				arr = []
				for char in line:
					arr.append(float(char))
				MOMSubject.append(arr)
			MOMS.append(np.asarray(MOMSubject))
	return MOMS


def readKernelWidth(fileName):
	tree = ET.parse(fileName)
	root = tree.getroot()
	kernelWidth = float(root.find('kernel-width').text)
	return kernelWidth


def addWeightsToXML(fileName,weights):
	tree = ET.parse(fileName)
	root = tree.getroot()
	subElement = ET.Element('regularity-weights')
	subElement.text = weights
	root.append(subElement)
	tree.write(fileName)


def gaussianKernel(CP,kernelWidth):
	K = kernel.rbf_kernel(CP,CP,1/(kernelWidth**2))
	return K


def moment2Velocity(K,MOM):
	v = np.dot(K,MOM)
	return v


def velocity2Moment(K,v):
	MOM = np.dot(np.linalg.pinv(K),v)
	return MOM


def geodesicDist(K,MOM):
	v = moment2Velocity(K,MOM)
	dist = np.sum(np.multiply(v,MOM))
	return dist


def dotProduct(m1,m2,K):
	prod = np.sum(np.multiply(m1,np.dot(K,m2)))
	return prod


def gaussianWeight(dist,var):
	return(np.exp(-dist/(2*var)))


def eucledianAngle(v1,v2):
	angle = np.sum(np.multiply(v1,v2))/(np.sum(np.multiply(v1,v1))*np.sum(np.multiply(v2,v2)))
	return angle


