import numpy as np
import math
import matplotlib.pyplot as plt

def f1(X):
    return (1.03 + X**2) / X
def f2(X):
    return math.sqrt(1.98 + 2 * X**3)

def NextIter(X, Y):
    X1 = ((Y**2 - 1.98) / 2)**(1.0 / 3)
    Y1 = X + 1.03 / X
    return X1, Y1

def Solution(X0, Y0, eps):
    Xk, Yk   = X0, Y0
    Xk1, Yk1 = NextIter(Xk, Yk)
    cnt = 0
    while (math.sqrt((Xk1 - Xk)**2 + (Yk1 - Yk)**2) > eps):
        cnt += 1
        Xk, Yk   = Xk1, Yk1
        Xk1, Yk1 = NextIter(Xk, Yk)
    print("Число итераций:", cnt)
    return Xk1, Yk1


plt.figure(figsize = (10, 10))

X = np.arange(-0.99, 2, 0.001)
Y = [f2(i) for i in X]
plt.plot(X, Y, 'y')

Y = [-i for i in Y]
plt.plot(X, Y, 'y')

X1 = np.arange(-1.5, -0.1, 0.001)
X2 = np.arange(0.1, 2, 0.001)
Y = [f1(i) for i in X1]
plt.plot(X1, Y, 'g')
Y = [f1(i) for i in X2]
plt.plot(X2, Y, 'g')

plt.grid()
plt.show()

eps = 1e-4
x0, y0 = 1.0, 2.0
x, y = Solution(x0, y0, eps)
print("Ответ: x =", format(x, '.3f'), ", y =", format(y,  '.3f'))

