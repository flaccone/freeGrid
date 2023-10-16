# Project: FREEGRID 
# updated 20230530
# This file computes buildability metrics of a mesh
# required Python >= 3.7 (64 bit), numpy, pymeshlab

import os
import sys
import csv
from email import header
from encodings import utf_8
import pymeshlab
import numpy as np
import math
import statistics
from itertools import permutations
import argparse

class mesh:
    def __init__(self, meshPath):
        self.path = meshPath
        self.computeBasicAttr()
        self.computeEdges()
        self.computeStarVE()
    
    def computeBasicAttr(self):
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(self.path)
        # vertices
        self.vn = ms.current_mesh().vertex_number()
        self.vm = ms.current_mesh().vertex_matrix()
        self.vNorm = ms.current_mesh().vertex_normal_matrix()
        # faces
        self.fm = ms.current_mesh().polygonal_face_list()
        self.fn = len(self.fm)
        # dimensions
        self.bb = ms.current_mesh().bounding_box()
        self.out_dict = ms.get_geometric_measures()
        # edges
        # avg_el = self.out_dict['avg_edge_length']
        # tot_el = self.out_dict['total_edge_length']
        # self.en = math.ceil(tot_el/avg_el)

    def printData(self):
        print('Verts number: ', self.vn)
        print('pFace number: ', self.fn)
        print('Edges number: ', self.en)

    def computeEdges(self):
        faceMatrix = self.fm
        edgeMa = []
        originFace = []
        for i in range (len(faceMatrix)):
            face = faceMatrix[i]
            for ii in range (np.shape(face)[0]):
                if face[ii] < face[ii-1]:
                    edge = [ face[ii], face[ii-1] ]
                else:
                    edge = [ face[ii-1], face[ii] ]
                edgeMa.append(edge)  # all edges (with duplicates)
                originFace.append(i) # corresponding face number
        edgeMatrix = np.array(edgeMa)
        em = np.unique(edgeMatrix, axis=0)
        
        # assert(self.en == np.shape(em)[0])
        self.en = np.shape(em)[0]

        # assebly face adjacency matrix per edge (one or two face number per edge)
        en = np.shape(em)[0] 
        adj = []
        for e in range(int(en)):
            currEdge = em[e]
            faces = []
            currEdgelist = currEdge.tolist()
            # a = edgeMa.index(currEdgelist)
            dupId = [i for i,edgeMa in enumerate(edgeMa) if edgeMa==currEdgelist]
            faces = [ originFace[i] for i in dupId]
            adj.append(faces)
        self.em = em
        self.adj = adj 

    def computeStarVE(self):
        #for each vertex find the incident edges (this star is unordered)
        em = self.em
        starVE = []
        for v in range(self.vn):
            #find num vert in edges and write a list
            rowId, colId = np.where(v==em)
            # print(rowId)
            starVE.append(rowId)
        self.starVE = starVE
   

def computeBins(nparray, binSize, plot=False, roundPrecision = 0):
    # computes the number of non-null bins of a distribution (array) centered on the mode value for a given binSize  
    mode = statistics.mode(np.round(nparray,roundPrecision))
    mean = statistics.mean(nparray)
    devs = statistics.stdev(nparray)
    nLeftBins = abs((mode - 0.5*binSize - nparray.min() ) / binSize) + 1
    nRightBins =abs((nparray.max() - mode + 0.5*binSize ) / binSize) + 1
    range0 = mode - 0.5*binSize - int(nLeftBins)* binSize
    range1 = mode + 0.5*binSize + int(nRightBins)*binSize
    bin = np.histogram(nparray, bins = np.arange(range0, range1, binSize), density=False)[0]
    fullBin = bin[bin != 0]
    numFullBins = len(fullBin)
    
    return mode, numFullBins, mean, devs


def projectEdgeOnPlane(pv0, pv1, v0norm):
    eArray = np.append(pv0, pv1, axis=0).reshape(2,3)
    e = np.diff(eArray, axis=0).squeeze()
    ap = np.dot(e, v0norm)
    pv1proj = pv1 - np.dot(ap, v0norm)

    eproj = np.append(pv0, pv1proj, axis=0).reshape(2,3)
    edproj = np.diff(eproj, axis=0).squeeze()
    unitEdproj = edproj / np.linalg.norm(pv0-pv1proj)
    dist = ap
    return unitEdproj, dist


def bestFittingPlane(verts):
    #SVD best fitting plane
    x, y, z = (( []) for i in range(3) )
    for i in range(len(verts)):
        x.append(verts[i][0])
        y.append(verts[i][1])
        z.append(verts[i][2])
    v = np.array([np.array(x), np.array(y), np.array(z)])

    # subtract out the centroid and take the SVD
    centroid = np.mean(v, axis=1, keepdims=True)
    svd = np.linalg.svd(v - centroid)

    # Extract the left singular vectors
    left = svd[0]
    # the corresponding left singular vector is the normal vector of the best-fitting plane
    normFit = left[:, -1]
    return normFit, centroid.squeeze()


def computeBeamSimilarity(mesh, plot_bool=False):
    en = mesh.en    
    vm = mesh.vm
    em = mesh.em

    #compute beam lenghts
    beamLenghts = []
    for i in range (int(en)):
        v0id = em[i,0]
        v1id = em[i,1]
        v0 = np.array(vm[v0id])
        # print('v0: ', v0)
        v1 = np.array(vm[v1id])
        # print('v1: ', v1)
        l= np.linalg.norm(v1 - v0)
        beamLenghts.append(l)
    bl = np.array(beamLenghts)
    bl2 = np.sort(bl[bl>0.1])
    totLen =  np.sum(bl)
    totLenSB = np.sum(bl2)
    print('Total length of the beams ' + str( (np.around(totLen,3)) ) )
    print('Total length of the beams (excluding small edges) ' + str( (np.around(totLenSB,3)) ) )
   
    mode, numFullBins, mean, stdev = computeBins(bl2,0.005,plot_bool,roundPrecision=1)
    beamStatsRow = [mode, numFullBins, mean, stdev]

    return beamStatsRow


def pointsBestFit(ptSet1, ptSet2, perm): 
    # based on the paper "Least-Squares Rigid Motion Using SVD"
    # the point sets have same num of nodes = same valence
    # all points in a set lie on a unit sphere
    valence = len(ptSet1)
    distanceMetrics = []

    for i in range(len(perm)):
        currentPerm = perm[i]
        XM = ptSet1[currentPerm,:].transpose() # vertices to be aligned
        YM = ptSet2.transpose() # centroid
        WM = np.identity(valence)

        SM = np.matmul( np.matmul(XM, WM), YM.transpose() )

        U, S, VT = np.linalg.svd(SM)

        MM = np.eye(U.shape[1])
        MM[MM.shape[0]-1, MM.shape[1]-1] = np.linalg.det( np.matmul(VT.transpose(), U.transpose()) )

        RM = np.matmul( np.matmul(VT.transpose(), MM) , U.transpose() )

        newXM = np.matmul(RM, XM)
        rotatedptSet1 = newXM.transpose()
        
        norma = np.linalg.norm(rotatedptSet1 - ptSet2,axis=1)
        distanceMetrics.append(math.sqrt( (np.sum(np.square(norma))) /valence))
 
    minDistanceMetrics = min(distanceMetrics)
    return minDistanceMetrics


def pointSetBestFit(pointSetList, centroid, perm):
    # computes the pointBestFit for a list of point sets and a centroid
    # the point sets have same num of nodes = same valence
    dList = []
    for i in range(len(pointSetList)):
        dList.append(pointsBestFit(pointSetList[i], centroid, perm) )
    return dList


def nodeClustering(pointSetList):
    # clustering using "The farthest point strategy for progressive image sampling"

    delta = 0.01 # convergence tolerance (in the paper Liu et al. (Eng. Str. 2023) is set to 0.001)
    # the point sets have same num of nodes = same valence
    centroids = []
    counter = 0
    centroids.append( pointSetList[0] )

    # computing permutations
    valence = centroids[0].shape[0]
    num = list(range(0,valence))
    perm = list(permutations(num))
    
    firstDist = pointSetBestFit(pointSetList, centroids[counter], perm )
    minValue = max(firstDist)

    if minValue < delta:
        numberOfClusters = 1
    else:
        dist = np.array(firstDist).reshape(1,len(firstDist))
        minDistFromClusters = dist.copy()   # initialization
        minValue = max(firstDist)           # initialization
        while minValue > delta:
            counter = counter + 1
            newId = np.unravel_index(np.argmax(minDistFromClusters, axis=None), dist.shape)
            centroids.append( pointSetList[newId[1]] )
            addedDist = pointSetBestFit(pointSetList, centroids[counter], perm )
            dist = np.row_stack ((dist,np.array(addedDist)))
            minDistFromClusters = np.min(dist,axis=0)
            minValue = np.max(minDistFromClusters)
        numberOfClusters = counter + 1

    return numberOfClusters


def computeNodeSimilarity(mesh, plot_bool=False):
    # metrics adapted from Liu et al. (Eng. Str. 2023)
    vm = mesh.vm
    vn = mesh.vNorm
    em = mesh.em
    st = mesh.starVE
    valenceVector = []
    pointSetList = []
    for i in range(len(vm)):
        valenceVector.append(len (st[i]) )
        starEdges = em[ st[i] ]
        # flip indices to have vm id as v0s
        idV1 = np.extract(i!=starEdges, starEdges)
        v1s = vm[idV1]
        v0s = np.full_like(v1s,vm[i])

        # normalize edges
        vEdges = (v1s - v0s)
        norms =  np.linalg.norm(vEdges, axis=1).reshape(len(st[i]),1)
        pointSetList.append(vEdges/norms)

    # define at least one cluster per node valence, the number of clusters is increased if needed
    valenceArray = np.array(valenceVector)
    clustersPerValence = np.unique(valenceArray)
    totNumOfClusters = 0
    for j in range(len(clustersPerValence)):
        rowsId = np.where(valenceArray==clustersPerValence[j])
        ids = rowsId[0].tolist()
        pointSetForClustering = [pointSetList[i] for i in ids]
        totNumOfClusters = totNumOfClusters + nodeClustering(pointSetForClustering)
        print("Number of joint clusters for valence " + str(clustersPerValence[j]) +  ": " + str(nodeClustering(pointSetForClustering)))
    
    return totNumOfClusters


def computePlanarity(mesh, binDimension, plot_bool=False):
    fm = mesh.fm
    vm = mesh.vm
    plan = []
    for i in range(int(mesh.fn)):
        pVerts = [ vm[f] for f in fm[i]]
        if len(pVerts) < 4:
            plan.append(0)
        else:
            #svd best fitting plane
            bfpNorm, centroid = bestFittingPlane(pVerts)
            avgDist = []
            halfPerimeter = []
            for e in range(len(pVerts)):
                edge, dist = projectEdgeOnPlane(centroid,pVerts[e], bfpNorm)
                #avg of dist
                avgDist.append(abs(dist)/len(pVerts))
                #half perimeter of the original face
                halfPerimeter.append(0.5 * np.linalg.norm(pVerts[e]-pVerts[e-1]) )
            #planarity
            plan.append( sum(avgDist)/sum(halfPerimeter) )

    planArray = np.array(plan)
    mode, numBins, mean, stdev = computeBins(planArray, binDimension, plot_bool, roundPrecision=2)
    planarity = [mode, numBins, mean, stdev]

    return planarity




################ MAIN ################

# parser for bool argument
parser = argparse.ArgumentParser(
    prog= "freeGridBuildability",
    description= "Computes the buildability metrics given an input mesh (.obj or .ply) and a number of cross section (int)"
    )
parser.add_argument('dir', help='The mesh path (.obj or .ply)')
parser.add_argument('--s', type=int, required=False, help='The number of cross sections (int)')
args = parser.parse_args()
print('Mesh path: ' + str(args.dir))
if args.s != None:
    print('Number of cross sections: ' + str(args.s))
else:
    print('Number of cross sections: not provided')

# check correct mesh format
meshPath = args.dir
if  not (str(meshPath).endswith(".obj") or str(meshPath).endswith(".ply") ):
    print("Please check the mesh format (.obj or .ply)")
    sys.exit(1)
filename = os.path.split(meshPath)[1]
name = os.path.splitext(filename)[0]

# output path == location of the .py file
outPath = os.path.dirname(os.path.realpath(__file__))
print(outPath)

#csv creation
c = open(os.path.join(str(outPath), str(name) + '_bMetrics.csv'), 'w', newline='', encoding='utf_8')
writer = csv.writer(c)
header = ['mesh name', '1 + \ bar(\Delta_0) [-]' , '#(N_0) [-]', '#(J_0) [-]', '1 + \ tilde(l_0) [-]', '#(C_0) [-]' ]
writer.writerow(header)

# load the mesh
m = mesh(str(meshPath))
m.printData()

# doublecheck input
try:
    m
except NameError:
    sys.exit('Please provide a mesh (.obj or .ply)')


# BUILDABILITY PERFORMANCE METRICS
# planarity
planarityRow = computePlanarity(m,0.01)
planarityMetric = np.around(1 + planarityRow[2],3)
print('Face planarity metric: ', planarityMetric)

# number of joints
vn = m.vn
print('Joint number: ', vn)

# node similarity according to Liu et al. (Eng. Str. 2023)
nSimilarityMetric = computeNodeSimilarity(m)
print('Joint similarity metric: ', nSimilarityMetric)

# beam lengths std
beamRow = computeBeamSimilarity(m)
lengthSimilarityMetric = np.around( 1+ ( beamRow[3] / beamRow[2] ), 3)
print('Beams similarity metric: ', lengthSimilarityMetric)

# number of cross sections used
crossSectionsMetric = args.s
print('Number of cross-sections: ', crossSectionsMetric)

# all results
caseRow = [str(name)] + [planarityMetric] + [vn] + [nSimilarityMetric] + [lengthSimilarityMetric] + [crossSectionsMetric]
print("***************")
print("RESULTS: ")
print(caseRow)

writer.writerow(caseRow)
c.close()
