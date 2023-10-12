import numpy as np
import math
import matplotlib.pyplot as plt

def f1(X):
    return math.sqrt(1 - X**2)
def f2(X):
    return math.tan(X)

def NewX(X):
    return math.atan(math.sqrt((1 - X**2)))
 
def Solution(X0, Epsilon):
    Xk = X0
    Xk1 = NewX(Xk)
    while (math.fabs(Xk1 - Xk) > Epsilon):
        Xk = Xk1
        Xk1 = NewX(Xk)
    return Xk1

plt.figure(figsize = (10, 10))

X = np.arange(-1, 1, 0.001)
Y = [f1(i) for i in X]
plt.plot(X, Y, 'y')

Y = [-i for i in Y]
plt.plot(X, Y, 'y')

Y = [f2(i) for i in X]
plt.plot(X, Y, 'g')

plt.grid()
plt.show()

Epsilon = 1e-6
X0 = 0.6
X = Solution(X0, Epsilon)
print("Ответ x =", format(X, '.2f'), ", y =", format(math.tan(X),  '.2f'))