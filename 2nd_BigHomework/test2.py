import scipy
from scipy import integrate
import numpy as np
import math
import matplotlib.pyplot as plt

def GetFunction(X : float) -> float:
    return math.log(100.0 - X) / (10.0 - math.sqrt(X))

def Getlk(t : float, k : int, ArgValues : list) -> float: #возвращает lk(t)
    # ArgValues - массив t
    n = len(ArgValues)
    lk = 1
    for j in range(n):
        if k != j:
            lk *= ((t - ArgValues[j]) / (ArgValues[k] - ArgValues[j]))  
        else: 1    
    return lk

def GetLagrangePoly(t : float, ArgValues : list, FunctionValues : list) -> float: #возвращает значение L(t)
    n = len(ArgValues)
    Value = 0
    for k in range(n):
        Value +=  (Getlk(t, k, ArgValues) * FunctionValues[k])
    return Value

def GetLagrangeValues(Args : list, ArgValues : list, FunctionValues : list) -> list:
    # Возвращает массив значений полинома Лагранжа для каждого значения аргумента из Args
    # ArgValues - массив аргумента t 
    # FunctionValues - массив значений функции в точках ArgValues
    Values = []
    for Arg in Args:
        Values.append(GetLagrangePoly(Arg, ArgValues, FunctionValues))
    return Values

def SimpsonIntegral(FunctionsVal : list, a : float, b : float) -> float:
    # Считает интеграл от поточечно заданной функции на отрезке [a, b] методом Симпсона
    Int = 0
    N = len(FunctionsVal)
    k = int(N / 2)
    h = (b - a) / N
    for i in range(1, k):
        F1 = FunctionsVal[2*i]
        F2 = FunctionsVal[2*i - 1]
        F3 = FunctionsVal[2*i - 2]
        Int += h / 3.0 * (F1 + 4 * F2 + F3)
    return Int

def NewtonIteration(Xk : float, P : float, P1 : float) -> float:
    # X_{k+1} = X_{k} - P(Xk) / Р1(Xk)
    # P - значение полинома Лежандра для Xk
    # P1 - значение производной полинома Лежандра для Xk
    return Xk - P / P1

def GetLejandrPoly(X : float, N : int) -> float:
    # Возвращает N-ый полином Лежандра для X
    if (N == 0):
        return 1
    if (N == 1):
        return X
    return (2.0 * N + 1) * X * GetLejandrPoly(X, N - 1) / (N + 1) - N * GetLejandrPoly(X, N - 2) / (N + 1)

def GetLejandrDerriative(X : float, N : int) -> float:
    # Возвращает первую производную N-ного полинома Лежандра для X
    return N * (GetLejandrPoly(X, N - 1) - X * GetLejandrPoly(X, N)) / (1 - X * X)

def GetLejandrZeroes(N : int) -> list:
    # Возвращает N нулей полинома Лежандра
    #, вычисленные итеративно по методу Ньютона с начальным приближением: X0 = cos(pi(4i - 1)/(4N + 2))
    Zeros = []
    Epsilon = 1e-3
    for i in range(1, N + 1):
        Xk = math.cos(math.pi * (4 * i - 1) / (4 * N + 2))
        Xk1 = NewtonIteration(Xk, GetLejandrPoly(Xk, N), GetLejandrDerriative(Xk, N))
        while (abs(Xk - Xk1) > Epsilon):
            Xk = Xk1
            Xk1 = NewtonIteration(Xk, GetLejandrPoly(Xk, N), GetLejandrDerriative(Xk, N))
        Zeros.append(Xk1)
    return Zeros

def FromMinusOneToOneInterval(Start : float, Stop : float, Vars : list) -> list:
    # Делает замену переменных в узлах квадратуры 
    HalfSum = (Start + Stop) / 2.0
    HalfDiff = (Stop - Start) / 2.0
    NewInterval = [HalfSum + HalfDiff * T for T in Vars]
    return NewInterval

def CalculateIntegralUsingGaussQuadrature(F : "function", a : float, b : float, N : int) -> float:
    # Считает интеграл от F на отрезке [a, b] методом квадратур Гаусса по N узлам
    Int = 0
    NodesT = GetLejandrZeroes(N)
    NodesX = FromMinusOneToOneInterval(a + 1e-1, b, NodesT)
    for k in range(1, N + 1):
        # Интегрирование методом Симпсона по 1000 точкам
        Args = np.arange(a, b, 0.001)
        BaseLagranValues = [Getlk(i, k - 1, NodesX) for i in Args]
        Ck = SimpsonIntegral(BaseLagranValues, a, b)
        Fk = F(NodesX[k - 1])
        Int += Ck * Fk
    return Int
    


Start_ = 0
End_ = 10

ArrayN = np.arange(2, 10)
Errors =[]
Errors2 = []

Exact = scipy.integrate.quad(GetFunction, Start_, End_)
print("Точное решение и его ошибка (I, delta I) =", Exact)

for N in ArrayN:
    Gauss = CalculateIntegralUsingGaussQuadrature(GetFunction, Start_, End_, N)
    FuncValues = [GetFunction(i) for i in range (Start_, End_)]
    Simpson = SimpsonIntegral(FuncValues, 0, 10)
    Errors2.append(abs(Gauss - Simpson)* 100/Simpson)
    print("\nN =", N, ": I =", format(Gauss, '.10f'))
    Errors.append(abs(Gauss - Exact[0]) * 100 / Exact[0])

plt.figure(figsize = (10, 10))
plt.title("Кривая зависимости относительной ошибки интегрирования от количества узлов")
plt.scatter(ArrayN, Errors2, marker = '^', color = 'g')
Args = np.arange(min(ArrayN), max(ArrayN), 0.01)
NewtonVals = GetLagrangeValues(Args, ArrayN, Errors2)
plt.plot(Args, NewtonVals, color = 'r')
plt.xlabel('Количество узлов n')
plt.ylabel('Ошибка, %')
plt.grid()

plt.show()

print(max(Errors2))
