import numpy as np
import numpy.random as rng

rng.seed(0)

# Number of points
N = 2001

# Timestamps and images
t = np.linspace(0.0, 100.0, N)
image = np.zeros(N, dtype="int64")
for i in range(0, N):
    image[i] = i % 4

# The time delays
tau = np.array([0.0, 10.0, 20.0, 30.0])
t2 = t.copy()
for i in range(0, N):
    t2[i] -= tau[image[i]]

# Covariance matrix
C = np.zeros((N, N))
C = np.matrix(C)
for i in range(0, N):
    for j in range(0, N):
        delta_t = t2[i] - t2[j]
        C[i, j] = np.exp(-np.abs(delta_t) / 50.0)
for i in range(0, N):
    C[i, j] += 0.1**2

# Cholesky decomposition
L = np.linalg.cholesky(C)

# Normals
n = np.matrix(rng.randn(N)).T

# Data
y = (L*n).T

data = np.empty((N, 4))
data[:,0] = t
data[:,1] = y
data[:,2] = 0.1
data[:,3] = image

np.savetxt("data.txt", data)

