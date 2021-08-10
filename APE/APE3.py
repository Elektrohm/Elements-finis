from numpy import *
from numpy.linalg import solve

A = array([[1, -1, 0],
           [-1, 2, -1],
           [0, -1, 1]])
b = array([[1/12], [1/2], [5/12]])

U = solve(A, (0.5**2)*b)
print(U)
