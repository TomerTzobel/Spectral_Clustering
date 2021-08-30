import numpy as np
import pandas as pd
import math
import sys
import mykmeanssp


def validateArgs(args):
    if len(args) == 3:
        if not args[0].isnumeric():
            raise Exception("1st arg is not a positive number")
    elif len(args) == 4:
        if not args[0].isnumeric():
            raise Exception("1st arg is not a positive number")
        elif not args[1].isnumeric():
            raise Exception("2nd arg is not a positive number")
    else:
        raise Exception("Invalid amount of args")


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

# main script
# args processing:
args_arr = sys.argv[1:]
validateArgs(args_arr)
k = int(args_arr[0])
if len(args_arr) == 4:
    max_iter = int(args_arr[1])
    filePath1 = args_arr[2]
    filePath2 = args_arr[3]
else:
    max_iter = 300
    filePath1 = args_arr[1]
    filePath2 = args_arr[2]
# reading and merging data files:
try:
    data1 = open(filePath1, "r")
    data2 = open(filePath2, "r")
except:
    raise Exception("Bad path provided")
table1 = pd.read_csv(data1, header=None)
table2 = pd.read_csv(data2, header=None)
data1.close()
data2.close()
mergeTabels = pd.merge(table1, table2, on=0).sort_values(0)
mergeTabels = mergeTabels.drop(columns=[0])  # delete the key value
pointsMatrix = mergeTabels.to_numpy()
#more validations:
pointsNumber, dimension = pointsMatrix.shape
if not (0 < k < pointsNumber):
    raise Exception("Invalid k")
if (pointsNumber == 0):
    raise Exception("Empty data")
# kmeans plus plus:
initCent = kmeans_pp(pointsMatrix, k)
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
