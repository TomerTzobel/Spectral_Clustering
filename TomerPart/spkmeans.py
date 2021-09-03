import numpy as np
import pandas as pd
import math
import sys
import spkmeansmodule as spkm

def get_init_centroids(points, k):
    np.random.seed(0)
    u1 = np.random.choice(len(points))
    z = 1
    centroids = np.zeros((k,), dtype=np.int32)
    centroids[0] = u1
    MinDisForPoint = [math.inf for i in range(len(points))]
    while(z < k):
        sum_D = float(0)
        distances = np.zeros((len(points)))
        for i, point in enumerate(points):
            newMinPoint = eval_distance(point, points[(centroids[z-1])])
            MinDisForPoint[i] = min(MinDisForPoint[i],newMinPoint)
            distances[i] = MinDisForPoint[i]
            sum_D += MinDisForPoint[i]
        for i in range(len(distances)):
            distances[i] = distances[i]/sum_D
        centroids[z] = np.random.choice(len(points), p=distances)
        z += 1
    return centroids

def eval_distance(a, b):
    norm = np.linalg.norm(a-b)
    return norm**2

ERR_MSG = "An Error Has Occured"
#main
args_arr = sys.argv[1:]
k = int(args_arr[0])
goal = args_arr[1]
filename = args_arr[2]

if(goal != "spk"):
    spkm.fit(goal, filename)

if(goal == "spk"):
    pointsMatrix = np.array(spkm.get_normalized_matrix_wrapper(k, filename))
    pointsNumber, dimension = pointsMatrix.shape
    k = dimension # in case k was 0
    init_centroids = get_init_centroids(pointsMatrix, k)
    listPoints = pointsMatrix.flatten().tolist()
    spkm.kmeans_pp(
        listPoints, init_centroids.tolist(), dimension, pointsNumber, k)
