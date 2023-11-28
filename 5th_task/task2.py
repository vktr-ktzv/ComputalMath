import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def getl_k(xs: np.array, k: int) -> np.array:
    poly = np.array([1.0])
    div = 1.0
    
    for n, x in enumerate(xs):
        if n == k:
            continue
            
        poly = np.polymul(poly, np.array([1.0, -x]))
        div *= (xs[k] - xs[n])

    poly /= div
    return poly

def LangrangePolynom(xs: np.array, ys: np.array) -> np.array:
    finalPoly = np.array([0.0])
    for n, y in enumerate(ys):
        finalPoly = np.polyadd(finalPoly, getl_k(xs, n) * y) 
    return finalPoly

def MakePlot(axes, xval: np.array, yval: np.array, poly: np.array):
    points = 100000
    xval_linspace = np.linspace(min(xval), max(xval), points)

    axes.scatter(xval, yval, color="red")
    axes.plot(xval_linspace, np.polyval(poly, xval_linspace))

    axes.get_xaxis().set_visible(True)
    axes.get_yaxis().set_visible(True)

    axes.grid(True, which="both")
    axes.set_xlabel("$x$")
    axes.set_ylabel("$y$")

fig, ax = plt.subplots(1, 1, figsize=(10, 6))

xs = np.array([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ys = np.array([0.02, 0.079, 0.175, 0.303, 0.459, 0.638, 0.831, 1.03, 1.23, 1.42])
poly_y_t = LangrangePolynom(xs, ys)

MakePlot(ax, xs, ys, poly_y_t)

plt.show()