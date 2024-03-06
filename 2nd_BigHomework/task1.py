import numpy as np
import math
import matplotlib.pyplot as plt

def MakeStepPoints(Start : float, Stop : float, N : int) -> float:
    # Возвращает массив из N точек c одинаковым интервалом
    Points = []
    Step = (Stop - Start) / (N - 1)
    for PointNum in range(N):
        Points.append(Start + PointNum * Step)
    return Points

def GetFunction(X : float) -> float: # возвращает заданную функцию
    return 1.0 / (1 + 25 * X * X)

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

##############################################

def GetDividedDifference(k : int, n : int, ArgValues : list, FunctionValues : list) -> float:
    # Возвращает разделенную разность 
    # ArgValues - массив значений t
    # FunctionValues - массив значений функции 
    if (k == n):
        return FunctionValues[0]
    F2 = GetDividedDifference(k + 1, n, ArgValues, FunctionValues[1:])
    F1 = GetDividedDifference(k, n - 1, ArgValues, FunctionValues[:-1])
    t2 = ArgValues[n]
    t1 = ArgValues[k]
    DivDiff = (F2 - F1) / (t2 - t1)
    return DivDiff

def GetNewtonPoly(ArgValues : list, FunctionValues : list) -> list:
    # Возвращает коэффициенты полинома Ньютона, являющиеся разделенными разностями
    n = len(ArgValues)
    DivDiffs = []
    for k in range(n):
        F = GetDividedDifference(0, k, ArgValues, FunctionValues)
        DivDiffs.append(F)
    return DivDiffs
    
def getNewtonPolyVal(t : float, ArgValues : list, NewtonPol : list) -> float:
    # Возвращает значение полинома Ньютона в точке t
    Result = 0
    n = len(ArgValues)
    for k in range(n):
        Mul = 1
        for i in range(k):
            Mul *= (t - ArgValues[i])
        Result += Mul * NewtonPol[k]
    return Result

def GetNewtonValues(Args : list, ArgValues : list, FunctionValues : list) -> list:
    # Возвращает масcив значений полинома Ньютона 
    NewtonValues = []
    NewtonPol = GetNewtonPoly(ArgValues, FunctionValues)
    for Arg in Args:
        NewtonValues.append(getNewtonPolyVal(Arg, ArgValues, NewtonPol))
    return NewtonValues

def GetChebZeros(Start : float, Stop : float, N : int) -> float:
    # Возвращает массив нулей полинома Чебышева
    Zeros = []
    HalfSum = (Start + Stop) / 2.0
    HalfDiff = (Stop - Start) / 2.0
    for ZeroNum in range(1, N + 1):
        Zero = HalfSum + HalfDiff * math.cos((2 * ZeroNum - 1) * math.pi / (2 * N))
        Zeros.append(Zero)
    return Zeros

#########################################

Start_ = -1
End_ = 1

plt.figure(figsize = (10, 10))
plt.title("Интерполяционный многочлен Лагранжа при разных n")

Args = np.arange(-1, 1.01, 0.01)
Vals = [GetFunction(Arg) for Arg in Args]
plt.plot(Args, Vals, 'r', label = "График исходной функции")

ArrayN = [4, 6, 10]

for N in ArrayN:
    ArgValues = MakeStepPoints(Start_, End_, N)
    FunctionValues = []
    for Arg in ArgValues:
        FunctionValues.append(GetFunction(Arg))
    plt.scatter(ArgValues, FunctionValues, marker = "^")
    Args = np.arange(-1, 1.01, 0.01)
    LagrangeVals = GetLagrangeValues(Args, ArgValues, FunctionValues)
    NameGraph = "Лагранж n =" + str(N)
    plt.plot(Args, LagrangeVals, label = NameGraph)
    
    Errors_ = np.zeros(201)
    for i in range(201):
        Errors_[i] = (Vals[i] - LagrangeVals[i])
    print(max(Errors_))

plt.grid()
plt.legend()


plt.figure(figsize = (10, 10))
plt.title("Интерполяционный многочлен Ньютона с узлами в нулях полинома Чебышева при разных n")

Args = np.arange(-1, 1.01, 0.01)
Vals = [GetFunction(Arg) for Arg in Args]
plt.plot(Args, Vals, 'r', label = "График исходной функции")

ArrayN = [4, 6, 10]

for N in ArrayN:
    NewtonArgValues = GetChebZeros(Start_, End_, N)
    FunctionValues = []
    for Arg in NewtonArgValues:
        FunctionValues.append(GetFunction(Arg))
    plt.scatter(NewtonArgValues, FunctionValues, marker = "^")
    Args = np.arange(-1, 1.01, 0.01)
    NewtonVals = GetNewtonValues(Args, NewtonArgValues, FunctionValues)
    NameGraph = "Ньютон n =" + str(N)
    plt.plot(Args, NewtonVals, label = NameGraph)
    
    Errors_ = np.zeros(201)
    for i in range(201):
        Errors_[i] = (Vals[i] - NewtonVals[i])
    print(max(Errors_))

plt.grid()
plt.legend()
plt.show()




    

