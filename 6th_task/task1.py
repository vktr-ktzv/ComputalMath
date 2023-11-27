import numpy as np 

a = 0.0
b = 3.0
n = 10000 

def Sipson(f, a: float, b: float, n: int) -> float:
    assert a < b
    
    h = (b - a) / float(n)
    ans = 0.0

    for i in range(0, n, 2):
        x_i = a + i * h
        x_j = x_i + h
        x_k = x_i + 2 * h

        f_i = f(x_i)
        f_j = f(x_j)
        f_k = f(x_k)

        ans += (f_i + 4.0 * f_j + f_k)

    return (ans / 3.0) * h

def Trapezoid(f, a: float, b: float, n: int) -> float:
    assert a < b
    
    h = (b - a) / float(n)
    ans = 0.0
    
    for i in range (0, n):
        x_i = a + i * h
        f_i = f(x_i)
        f_i1 = f(x_i + h)
        
        ans += (f_i + f_i1)
        
    
    return (ans / 2.0) * h 

def f_1(x: float) -> float: 
    return np.sin(100.0 * x) * np.exp(-x ** 2.0) * np.cos(2.0 * x)



TrapAns = Trapezoid(f_1, a, b, n)
SipsonAns = Sipson(f_1, a, b, n)

print("Trapezoid:", TrapAns)
print("Sipson:", SipsonAns)
print("|Simpson - Trapezoid|:", abs(SipsonAns - TrapAns))