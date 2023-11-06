import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.linalg import norm

def GetMatrix(N): # A
    Arr = np.zeros(shape = (N, N))
    
    for i in range(N):
        if (i > 0):
            Arr[i][i - 1] = -1
        if (i < (N - 1)):
            Arr[i][i + 1] = -1
            
    for i in range(N):
        Arr[i][i] = 2 + i**2 / N**2
        
    Arr[0][N - 1] = -1
    Arr[N - 1][0] = -1
    
    return Arr

def GetRightPart(N): # f 
    f = np.zeros(shape = (N, 1))
    
    for i in range(N):
        f[i, 0] = (1 + N**2 * math.sin(math.pi / N)**2) * math.sin(2 * math.pi * i / N)
        
    return f

def GershCircles(A): # return (cntr, R)
    Circles = []
    
    for i in range(len(A)):
        cntr = A[i, i]
        r = 0.0
        for j in range(len(A)):
            if  j != i:
                r += math.fabs(A[i, j])
                
        Circles.append((cntr, r))
        
    return Circles

def KrilovsMatrix(A, Y): #  Yk = A**k * Y 
    N = len(A)
    MtxKr = []
    Yk = Y
    
    for k in range(N, 0, -1):
        Yk = np.dot(A, Yk)
        MtxKr.append(np.reshape(Yk, (1, N)))
        
    Mtx = np.zeros(shape = (N, N))
    
    for i in range(N):
        for j in range(N):
            Mtx[i, j] = MtxKr[i][0][j]
            
    return np.transpose(Mtx)

def NewIterPoly(Xk, P): # X_(k+1) = (X_(k)**N * P[0] + ... + X_(k)**2 * P[N - 2] + P[N]) / (- P[N - 1]) 
    N = len(P) - 1
    Xk1 = P[N] / (- P[N - 1])
    for j in range(2, N + 1):
        Xk1 += P[N - j] * Xk**j / (- P[N - 1])
    return Xk1
    
def GetEigsByMPI(A, P, Eps):
    Eigs = []
    Circles = GershCircles(A)
    print (Circles)
    N = len(A)
    
    for i in range(len(A)):
        Xk = Circles[i][0] - Circles[i][1] / 2
        Xk1 = NewIterPoly(Xk, P)
        
        while (math.fabs(Xk - Xk1) > Eps):
            Xk = Xk1
            Xk1 = NewIterPoly(Xk, P)
        
        P1 = [1, -Xk1]
        P = np.polydiv(P, P1)[0]
        Eigs.append(Xk1)
    
    return Eigs

def GetOptTauByKrilov(A, Eps):
    N = len(A)
    Y = np.ones(shape = (N, 1))
    Mat = KrilovsMatrix(A, Y)
    P = np.linalg.solve(Mat, Y)
    P1 = []
    
    for i in range(len(P) - 1, -1, -1):
        P1.append(P[i, 0])
        
    P1.append(-1)
    Eigs = GetEigsByMPI(A, P1, Eps)
    return 2.0 / (max(Eigs) + min(Eigs))

def GetOptTauByLinalg(A):
    Eigs, Vectors = np.linalg.eig(A)
    return 2.0 / (max(Eigs) + min(Eigs))
    
def GetNextXinMPI(Xk, A, f, Tau):
    Xk1 = (np.eye(len(A)) - Tau * np.matrix(A)) * Xk + Tau * f
    return Xk1

def SolveByMPI(A, f, Tau, Eps):
    Xk = np.ones(shape = (len(A), 1))
    Xk1 = GetNextXinMPI(Xk, A, f, Tau)
    Norm = norm((Xk1 - Xk), ord=2)
    Count = 0
    Norms = []
    
    while (Norm > Eps):
        Norms.append(Norm)
        Count += 1
        Xk = Xk1
        Xk1 = GetNextXinMPI(Xk, A, f, Tau)
        Norm = norm((Xk1 - Xk), ord=2)    
    return Xk1, Norms

def PrintGraph(Nrm):
    plt.figure(figsize = (10, 10))
    plt.xlabel('Шаг итерации')
    plt.ylabel('Норма невязки')
    plt.plot(np.arange(len(Nrm)), Nrm, 'r')
    plt.grid()
    plt.show()

####################################
Eps = 1e-4
N = 6
A = GetMatrix(N)
print(" A:\n", A)
f = GetRightPart(N)
print("f:\n", f)
X0 = np.linalg.solve(A, f)
print("\n решение np.linalg \n X0 =\n", X0)

X, Norms = SolveByMPI(A, f, 0.4, Eps)
dX = norm(X - X0, ord=2)
print("\n\n", "Произвольное (Tau =", 0.4, "\nX:\n", X, \
        "\n\nНорма разницы ||X - X0||:", format(dX, '.10f'))
PrintGraph(Norms)

KrilTAu = GetOptTauByKrilov(A, Eps)
X, Norms = SolveByMPI(A, f, KrilTAu , Eps)
dX = norm(X - X0, ord=2)
print("\n\n", "Метод Крылова Tau =", KrilTAu, "\nX:\n", X, \
        "\n\nНорма разницы ||X - X0||:", format(dX, '.10f'))
PrintGraph(Norms)

LinalgTAu = GetOptTauByLinalg(A)
X, Norms = SolveByMPI(A, f,LinalgTAu , Eps)
dX = norm(X - X0, ord=2)
print("\n\n", "Метод Крылова Tau =", LinalgTAu, "\nX:\n", X, \
        "\n\nНорма разницы ||X - X0||:", format(dX, '.10f'))
PrintGraph(Norms)