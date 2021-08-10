from numpy import *
import matplotlib.pyplot as plt
import matplotlib


A = array([[2, -1], [-1, 2]])
b = array([1, 1])
xinit = array([1.5, 1.0])
Xmin = -4.5
Xmax = 3.5
Ymin = -2.5
Ymax = 3.5
eps = 0.0001

# GRADIENT CONJUGUE :
x = xinit.copy()
conjugate = [(x[0], x[1])]
i = 0
delta = 1
imax = 2

r = b - A @ x
d = r.copy()

while i < imax and delta > eps:
     delta = r @ r
     s = A @ d
     alpha = delta / (r @ s)
     x = x + alpha * d
     r = r - alpha * s
     beta = (r @ r) / delta
     d = r + beta * d
     conjugate.append((x[0], x[1]))
     i += 1
     print(f"x = {x}")
     print(" = Iteration %4d : %14.7e" % (i, delta))

conjugate = array(conjugate)

matplotlib.rcParams['toolbar'] = 'None'
myColorMap = matplotlib.cm.jet

n = 30

X, Y = meshgrid(linspace(Xmin, Xmax, n), linspace(Ymin, Ymax, n))
J = zeros((n, n))
for i in range(n):
     for j in range(n):
          x = array([X[i, j], Y[i, j]])
          J[i, j] = (x @ A @ x / 2 - b @ x)

plt.figure("Conjugate gradients - steepest descent")
plt.contourf(X, Y, J, 12, cmap=myColorMap)
plt.axis('equal')
plt.axis('off')
plt.plot(conjugate[:, 0], conjugate[:, 1], '-or')
plt.grid()
plt.show()