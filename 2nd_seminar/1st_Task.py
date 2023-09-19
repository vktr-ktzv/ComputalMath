import matplotlib as plt
import numpy as np
from sympy import diff, symbols, cos, sin

def comb(n, k) :
    if k ==0:
        return 1
    
    kk = 1
    nn = n 
    for i in range(1, k) :
        nn *= n - i 
        kk *= i + 1    
    return nn // kk
    

def dev(n):
    ans = 0
    h = 1e-5
    for k in range (0, n+1):
        print('comb',comb(n, k))
        print('sin', sin(k*h))
        print('k',(-1)**k, k)

        ans += comb(n, k)*sin(k*h)*((-1)**k)
        
    ans = ans/(h**n)
    return ans

def fact(n):
    tmp = 1
    while n > 1:
        tmp *= n
        n -= 1
    return tmp
        


def func(n):
    ans = 0
    for i in (0, n):
        ans += dev(n)*(0.5**n)/fact(n) 
    return ans 

    

du = 1e-3
u = np.sin(0.5)
eps = u
n = 0

print(dev(5))

