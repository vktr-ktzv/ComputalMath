import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sp

def Xn(x1, x2, N):
    X = []
    h = (x2 - x1) / N
    for i in range(1, N, 1):
        X.append(x1 + i * h) 
    return X

def eigenvalues(x1, x2, N):
    L = []
    X = x2 - x1
    h = X * 1.0 / N
    for i in range(1, N, 1):
        L.append(4.0 / (h * h) * math.sin(math.pi * i * h / 2.0 / X)**2)
    return L              

def eigenfunction(x1, x2, N, Xn):
    F = []
    X = x2 - x1
    for i in range(1, N, 1):
        F.append(math.sqrt(2.0 / X) * math.sin(math.pi * i * Xn / X))
    return F                    
                
def G(X):
    G = []
    for x in X:
        G.append(x**3)
    return G
            
def Solution(x1, x2, X, Ck_, Lambdas):
    N = len(X) + 1
    U = []
    for i in range(1, N, 1):
        EigenF = eigenfunction(x1, x2, N, X[i - 1])
        U.append(0)
        for j in range(1, N, 1):
            U[i - 1] += Ck_[j - 1] / Lambdas[j - 1] * EigenF[j - 1]
    return U

    
def main():
    x0 = 0
    x1 = 1
    N = 100
    X = Xn(x0, x1, N)
    Lambdas = eigenvalues(x0, x1, N)
    
    Matrix = []
    for i in range(1, N, 1):
        Matrix.append(eigenfunction(x0, x1, N, X[i - 1]))
    
    g = G(X)
    Ck_ = np.linalg.solve(Matrix, g)
    u = Solution(x0, x1, X, Ck_, Lambdas)
   
    plt.figure(figsize = (15, 8))
    plt.title("Решение задачи методом Фурье:")

    plt.plot(X, u, 'r', label = "u(x)")

    plt.legend()
    plt.grid()
    plt.show()
    
if __name__ == '__main__':
    main()