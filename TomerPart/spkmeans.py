import numpy as np
import pandas as pd
import math
import sys
import nsc_wrapper

def kmeans_pp(points, k):
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
max_iter = 300
k = int(args_arr[0])
goal = args_arr[1]
filename = args_arr[2]

if(goal != "spk"):
    nsc_wrapper.fit(k, goal, filename)

try:
    data = open(file, "r")
except:
    raise Exception(ERR_MSG)
table = pd.read_csv(data, header=None)
data.close()
pointsMatrix = table.to_numpy()
#more validations:
pointsNumber, dimension = pointsMatrix.shape
if not (0 < k < pointsNumber):
    raise Exception(ERR_MSG)
if (pointsNumber == 0):
    raise Exception(ERR_MSG)
# kmeans plus plus:
initCent = kmeans_pp(pointsMatrix, k)
joined_string = ",".join([str(element) for element in initCent.tolist()])
print(joined_string)

listPoints = pointsMatrix.flatten().tolist()
listCent = []
for cent in initCent:
    listCent.extend(pointsMatrix[cent].tolist())
centroidsFinal = mykmeanssp.fit(
    listPoints, listCent, initCent.tolist(), dimension, pointsNumber, k, max_iter)
# print result:
joined_string = ",".join([str(element) for element in initCent.tolist()])
print(joined_string)
for i in range(len(centroidsFinal)):
    for j in range(dimension):
        print((np.round_(centroidsFinal[i][j], decimals=4)), end="")
        if j < dimension - 1:
            print(',', end="")
    print("")
