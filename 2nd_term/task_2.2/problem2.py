import math
import numpy as np 
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# dy1/dt = y1 - y1 * y2
# dy2/dt = -y2 + y1 * y2
# y1(0) = 2, y2(0) = 2

def f(t, y1, y2):
    return y1 - y1 * y2
def g(t, y1, y2):
    return -y2 + y1 * y2

def create_plot(x, y, title = ''):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("t")
    plt.ylabel("y")
    plt.title(title)
    plt.grid()
    plt.show()

def main():
    h = 1e-3
    t_list = np.arange(0, 10, h)
    y1_list = [2]
    y2_list = [2]
    # u1 = h * 23/12 * f0
    y1_list.append(h * (23 / 12) * f(t_list[0], y1_list[0], y2_list[0])) 
    y2_list.append(h * (23 / 12) * g(t_list[0], y1_list[0], y2_list[0]))

    # u2 = u0 + h * ((23/12) * f1 - 16/12 * f0)
    y1_list.append(y1_list[0] + h * ((23 / 12) * f(t_list[1], y1_list[1], y2_list[1]) - \
                                     (16/12) * f(t_list[0], y1_list[0], y2_list[0])))
    y2_list.append(y2_list[0] + h * ((23 / 12) * g(t_list[1], y1_list[1], y2_list[1]) - \
                                     (16/12) * g(t_list[0], y1_list[0], y2_list[0])))

    for i in range(3, len(t_list)):
        # u_i = u_(i-2) + h * (23/12 * f_(i-1) - 16/12 * f_(i-2) + 5/12 * f_(i - 3))
        y1_list.append(y1_list[i - 2] + h * ((23/12) * f(t_list[i-1], y1_list[i-1], y2_list[i-1]) - \
                                             (16/12) * f(t_list[i-2], y1_list[i-2], y2_list[i-2]) + \
                                             (5/12) * f(t_list[i-3], y1_list[i-3], y2_list[i-3])))
        y2_list.append(y2_list[i - 2] + h * ((23/12) * g(t_list[i-1], y1_list[i-1], y2_list[i-1]) - \
                                             (16/12) * g(t_list[i-2], y1_list[i-2], y2_list[i-2]) + \
                                             (5/12) * g(t_list[i-3], y1_list[i-3], y2_list[i-3])))

    create_plot(t_list, y1_list, 'y1(t)')
    create_plot(t_list, y2_list, 'y2(t)')

    return 0

if (__name__ == "__main__"):
    main()
