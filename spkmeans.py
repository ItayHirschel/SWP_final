import numpy as np
import sys
import mykmeanssp
import myspkmodule

EPS_KMEANS = 0
MAX_ITER_KMEANS = 300
MODES = ["wam", "ddg", "spk", "lnorm", "jacobi"] # legal work modes

# prints a given list as a vector
def printVector(list):
    printline = ""
    for index in range(len(list)):
        if (index > 0):
            printline = printline + ","
        printline = printline + ("%.4f" % list[index])
    print(printline)
# prints a list as a Xdim x Ydim Matrix 
def printMatrix(data, Xdim, Ydim):
    for index in range(Ydim):
        templst = data[index * Xdim : (index + 1) * Xdim]
        printVector(templst)

# reads the arguments from console
def read_ARGS(arg_lst):                                       
    try:
        if len(arg_lst) == 4:
            file = arg_lst[3]
            k = int(arg_lst[1])
            goal = arg_lst[2]
            if k < 0 or (goal not in MODES):
                raise Exception
            return file, k, goal
        else:
            raise Exception
        
        
    except:
        print("Invalid Input!")
        sys.exit()

def read_IN_FILE(IN_FILE):  #read the input file and save data to matrix of lists
    f = open(IN_FILE, 'r') # opens file
    rows = f.readlines() # creates stream of lines from the file
    f.close() # closes file
    points = []
    for row in rows: # translates the stream into points
        st_vec = row.split(",")                 #split string by ','            
        if len(row) > 0:                        #if row isn't empty
            vec =  [float(x) for x in st_vec]    #take vector made from list of strings
            points.append(vec)
    return np.array(points)

# chooses randomly the first centroid
def init(points):
    arr = [i for i in range(len(points))]
    l = np.random.choice(arr,1).tolist()
    ind = l[0]
    return points[ind]

# chooses random centroid based of probability vector prob
def randomCentroid(points, prob):
    arr = [i for i in range(len(points))]
    l = np.random.choice(arr, size = 1 ,replace= False, p = prob).tolist()
    ind = l[0]
    return points[ind]

# calculates the squared distance between the vectors x and y
def dist_sq(x,y):
    sum = 0
    i = 1  # the first coordinate is the index of the observation
    while(i < len(x)):
        sum += (x[i] - y[i]) ** 2
        i += 1
    return sum

# returns the index of the nearest centroid (between the i+1 existing centroids) to a given vector x
def find_min(x, centroids, i):
    min = -1
    for t in range(i + 1):
        d = dist_sq(centroids[t],x)
        if min == -1:
            min =d
        else:
            if d < min:
                min = d
    return min

# calculates probability of the point in place l to be chosen as the next centroid
def calcP(l,D,sum):
    num = D[l]/sum
    return num

# calculates the sum of the array
def calcSum(array):
    sum = 0
    for a in array:
        sum += a
    return sum  

# Kmeans++ algorithm from HW2:
def Kmeans(points, k, N, VecDim): 
    length = N
    if length < k:
        raise Exception
    

    D = [0 for i in range(length)]
    P = [0 for i in range(length)]
    
    centroids = [ [] for i in range(k)]
    centroids[0] = init(points)
    
    I = 0
    while I < k :                                      # 0 <= i <= k-1
        for l_1 in range(length):                            # 0 <= l <= N - 1
            D[l_1] = find_min(points[l_1],centroids, I)
        
        s = calcSum(D)
        for l_2 in range(length):                            # 0 <= l <= N - 1
            P[l_2] = calcP(l_2,D,s)
        I += 1
        if I == k:
            break
        centroids[I] = randomCentroid(points,P) #chooses the next random centroid

    printline = ""
    for t in range(k): 
        if (t > 0):
            printline = printline + ","

        printline = printline + str(int(centroids[t][0])) 
    print(printline) # prints the chosen indices of the observations
    
    # adjusting the arrays to be used in C
    CENTROID_arr = np.array(centroids)[:,1:]
    POINTS_arr = np.array(points)[:,1:]
    points_lst = POINTS_arr.flatten().tolist()
    centroids_lst = CENTROID_arr.flatten().tolist()
    
    # applying the KMeans++ algorithm
    RESULT_lst = mykmeanssp.fit(points_lst, centroids_lst, N, VecDim, k, MAX_ITER_KMEANS, EPS_KMEANS)

    # printing the final centroids
    for index in range(k):
        templst = RESULT_lst[index * VecDim : (index + 1) * VecDim]
        printVector(templst)

# "unflatten" lists
def UnFlat(flat_lst, VecDim, N):
    lst = []
    for index in range(N):
        templst = flat_lst[index * VecDim : (index + 1) * VecDim]
        lst.append(templst)
    return lst

# the main algorithm
def MAIN(FILE, GOAL, K):
    np.random.seed(0)
    POINTS = read_IN_FILE(FILE)
    N = POINTS.shape[0]
    D = POINTS.shape[1]
    
    # checks if the input is a square matrix if the goal is jacobi
    if GOAL == "jacobi" and N != D: 
        raise Exception


    POINTS_lst = POINTS.flatten().tolist() # adjusts points to be used in C
    Vecs, Vals, largestKnorma, K = myspkmodule.spkFit(GOAL, K, D, N, POINTS_lst) # applying SPKMeans algorithm
    # printing results if goal is not spk
    if GOAL not in ["jacobi", "spk"]:
        printMatrix(Vecs, N, N)
    elif GOAL == "jacobi":
        printVector(Vals)
        printMatrix(Vecs, N, N)
    # if goal is spk, apply KMeans++ algorithm:
    elif GOAL == "spk":
        # treating each row as a point in R^k, and passing the points, in the right format, to KMeans++ :
        obs = UnFlat(largestKnorma, K, N) 
        indexed_obs = []
        
        for i in range(N):   # adding index to each observation
            tmp = [i] + obs[i]
            indexed_obs.append(tmp)
            
        
            
        Kmeans(indexed_obs, K, N, K) # Kmeans++ 


# the main function:
if __name__ == "__main__":
    
    arg_lst = sys.argv
    FILE, K, GOAL = read_ARGS(arg_lst) # reads arguments from console

    try: # apply the main algorithm 
       MAIN(FILE, GOAL, K)
    except: # case of error
        print("An Error Has Occured")
        sys.exit()