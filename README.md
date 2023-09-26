# SWP_final
The project implements a version the normalized spectral clustering algorithm.
Normalized Spectral Clustering We present the Normalized Spectral Clustering algorithm:

# The Normalized Spectral Clustering Algorithm:
## 1 : Form the weighted adjacency matrix $W$ from $X$
Assume we are given the input matrix $X$ then we define the Weighted Adjacency Matrix (WAM) $W$ to be:
$$ W_{i,j} = \exp{\left(-\frac{\|X_{i} - X_{j}\|}{2}\right)} $$
Where $ W_{i,j} $ is the element in $W$ on row i and column j, and $X_{i}$ is row i in matrix $X$ treated as a vector.
$\|X_{i} - X_{j}\|$ is the euclidean distance between the two vectors.

## 2 : Compute the normalized graph Laplacian $L_{norm}$
The diagonal degree matrix (DDM)  $D \in \mathbb{R}^{n \times n}$  is defined as $D = (d_{i,j})_{i,j=1,...,n} $ such that:

if i = j :

$ d_{i,j} = \sum_{z = 1}^{n} w_{iz} $

else :

$ d_{i,j} = 0 $

The normalized graph Laplacian $L_{norm} \in \mathbb{R}^{n \times n}$ is defined as:

$$ L_{norm} = I - D^{-\frac{1}{2}} W D^{-\frac{1}{2}} $$

The reason we are interested in Lnorm is that it has all eigenvalues $\lambda_1 \ge \dots \ge \lambda_n \ge 0 $ are real and
non-negative.

## 3 : Determine k and obtain the largest k eigenvectors $u_1, \dots , u_k$ of $L_{norm}$

First we obtain the eigenvalues of $L_{norm}$ by diagonising the matrix according to the Jacobi Algorithm.
Then, if we were not given $k$ as an argument, we obtain $k$ by using the heuristic:
$$ \delta_i = |\lambda_i - \lambda_{i+1}|$$ 
$$ k = \arg\max_{1 \le i \le \frac{n}{2}} \delta_i $$

## 4 : Let $U \in \mathbb{R}^{n \times k}$ be the matrix containing the vectors  $u_1, \dots , u_k$ as columns
We cut the matrix $L_{norm}$ such that only the k rightmost columns remain.
## 5 : Form the matrix $T \in \mathbb{R}^{n \times k}$ from $U$ by renormalizing each of U's rows to have unit length

## 6 : Treating each row of $T$ as a point in $\mathbb{R}^k$, cluster them into k clusters via the K-means algorithm

# Files :


## spkmeans.py
A Python interface for the code.


## spkmeans.c
A C interface for the code, that implements the Normalized Spectral Clustering algorithm.


## spkmeans.h
Header file.


## spkmeansmodule.c 
The Python C API wrapper.


## setup.py
The C API set up file.


## comp.sh
The compilation script.

## kmeans.c
Implements the K-Means algorithm.

# Setup Instructions 
- Run the setup.py file by the command :\
"python3 setup.py build_ext --inplace"
- Compile the C file by the command:\
"bash comp.sh"

# Running instructions
requirements:
- numpy

After setting up, you can run the project by two means:

## Python Interface :
Run the command : \
"python3 spkmeans.py \<arg1\> \<arg2\> \<arg3\>"

- \<arg1\> : The k argument (must be given)
- \<arg2\> : Mode 
    - spk : Perform full spectral kmeans
    - wam : Calculate the Weighted Adjacency Matrix
    - ddg : Calculate Diagonal Degree Matrix
    - lnorm : Calculate and output the Normalized Graph Laplacian
    - jacobi : Calculate and output the eigenvalues and eigenvectors
- \<arg3\> : Path to input file, which contains:
    - For **jacobi** a symmetric matrix.
    - **Else** the data points in a csv format.

**Outputs:**
* **spk** : The first line will be the indices of the observations chosen by the K-means++
algorithm as the initial centroids. The second line onward will be the calculated final centroids from the K-means algorithm, separated by a comma, such that each centroid is in a line of its own.
* **jacobi** : Case of 'Jacobi': The first line will be the eigenvalues, second line onward will be the corresponding eigenvectors (printed as columns).
* **else** : output the required matrix separated by a comma, such that each row is in a line of its own.

## C Interface
Run the command : \
"./spkmeans \<arg1\> \<arg2\>"
- \<arg1\> : Mode 
    - wam : Calculate the Weighted Adjacency Matrix
    - ddg : Calculate Diagonal Degree Matrix
    - lnorm : Calculate and output the Normalized Graph Laplacian
    - jacobi : Calculate and output the eigenvalues and eigenvectors
- \<arg2\> : Path to input file, which contains:
    - For **jacobi** a symmetric matrix.
    - **Else** the data points in a csv format.

**Outputs:**
* **jacobi** : Case of 'Jacobi': The first line will be the eigenvalues, second line onward will be the corresponding eigenvectors (printed as columns).
* **else** : output the required matrix separated by a comma, such that each row is in a line of its own.