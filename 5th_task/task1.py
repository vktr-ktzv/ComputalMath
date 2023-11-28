import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def DiffN(xval: np.array, yval: np.array) -> float: #разделенные разности
    n = len(xval)

    if n == 1:
        return yval[0]
    elif n == 2:
        x1, x2 = xval
        y1, y2 = yval
        return (y2 - y1) / (x2 - x1)

    x0 = xval[0]
    xn = xval[-1]
    
    return (DiffN(xval[1:], yval[1:]) - DiffN(xval[0:-1], yval[0:-1])) / (xn - x0)

def NewtonsInterpolation(xval: np.array, yval: np.array) -> np.array:
    curPolynom = np.array([1.0])
    finalPolynom = np.array([0.0])

    for n, x in enumerate(xval):
        finalPolynom = np.polyadd(finalPolynom, DiffN(xval[0:n + 1], yval[0:n + 1]) * curPolynom)
        multiplePol = np.array([1.0, -x])
        curPolynom = np.polymul(curPolynom, multiplePol)
    
    return finalPolynom

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
    
    
fig, axes = plt.subplots(1, 1, figsize=(10, 6))

xval = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
yval = np.array([1.0, 0.8, 0.5, 0.307, 0.2, 0.137, 0.1, 0.075, 0.06, 0.047, 0.039])
polynom = NewtonsInterpolation(xval, yval)

MakePlot(axes, xval, yval, polynom)
plt.show()
print(polynom, " - коэффициенты полинома")