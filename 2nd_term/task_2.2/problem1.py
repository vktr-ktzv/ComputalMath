import math
import numpy as np 
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

A = [[1/2 - np.sqrt(3) / 6, 0], [np.sqrt(3) / 3, - 1 / 2 - np.sqrt(3) / 2]]
c = [1 / 2 - np.sqrt(3) / 6, - 1 / 2 - np.sqrt(3) / 6]
b = [1 + np.sqrt(3) / 6, - np.sqrt(3) / 6]
e = [1, 1]

def create_plot(x, y, title = ''):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(title)
    plt.grid()
    plt.show()

# dv/dt = -sin(u) = f  v(0) = 1
# du/dt = v = g u(0) = 0
# functions = [-sin(u), v]
def f(t, v, u):
    return - np.sin(u)

def g(t, v, u):
    return v



def R(z):
    nominator = np.eye(2) - np.dot(z, A) + np.dot(z, np.dot(e, b))
    denominator = np.eye(2) - np.dot(z, A)

    return np.linalg.det(nominator) / np.linalg.det(denominator)

def make_data():
    xgrid = np.arange(-6, 1, 0.01)
    ygrid = np.arange(-2, 2, 0.01)
    coordinates = []
    for x_ in xgrid:
        for y_ in ygrid:
            coordinates.append([x_, y_, abs(R(x_ + 1j * y_))])

    return coordinates

def make_graph():
    coord = make_data()
    x_temp, y_temp = [] ,[]
    
    for t in coord:
        #print(t)
        if(t[2] <= 1):
            x_temp.append(t[0])
            y_temp.append(t[1])
    

    plt.scatter(x_temp, y_temp)
    plt.title('Область устойчивости метода Рунге-Кутта 2-го порядка (R(z) <= 1)')
    plt.show()

def count_kq_for_RK(f, g, t0, v0, u0, h = 1e-3):
    k = [f(t0 + c[0] * h, v0, u0),
         f(t0 + c[1] * h, v0 + h * f(t0 + c[0] * h, v0, u0), u0 + h * f(t0 + c[0] * h, v0, u0))]
    q = [g(t0 + c[0] * h, v0, u0),
         g(t0 + c[1] * h, v0 + h * g(t0 + c[0] * h, v0, u0), u0 + h * g(t0 + c[0] * h, v0, u0))]
    
    DeltaH = 2 * h
    while(DeltaH > h):
        temp_k = [f(t0 + c[0] * h, v0 + h * (A[0][0] * k[0] + A[0][1] * k[1]), u0 + h * (A[0][0] * k[0] + A[0][1] * k[1])),
                  f(t0 + c[1] * h, v0 + h * (A[1][0] * k[0] + A[1][1] * k[1]), u0 + h * (A[1][0] * k[0] + A[1][1] * k[1]))]
        temp_q = [g(t0 + c[0] * h, v0 + h * (A[0][0] * k[0] + A[0][1] * k[1]), u0 + h * (A[0][0] * k[0] + A[0][1] * k[1])),
                  g(t0 + c[1] * h, v0 + h * (A[1][0] * k[0] + A[1][1] * k[1]), u0 + h * (A[1][0] * k[0] + A[1][1] * k[1]))]
        
        DeltaH = max(abs(k[0] - temp_k[0]), abs(k[1] - temp_k[1]), abs(q[0] - temp_q[0]), abs(q[1] - temp_q[1]))

        k = temp_k
        q = temp_q
        

    return k, q

def Runge_Khutta(f, g, v0, u0, t0, t1, h = 1e-3):
    t_list = np.arange(t0, t1, h)
    v_list = [v0]
    u_list = [u0]

    for i in range(len(t_list) - 1):
        k, q = count_kq_for_RK(f, g, t_list[i], v_list[i], u_list[i], h)
        v_list.append(v_list[i] + h * (b[0] * k[0] + b[1] * k[1]))
        u_list.append(u_list[i] + h * (b[0] * q[0] + b[1] * q[1]))

    
    create_plot(t_list, u_list)


def main() :
    make_graph()
    Runge_Khutta(f, g, 1, 0, 0, 4 * np.pi, h = 1e-3)

    return 0

if (__name__ == "__main__"):
    main()