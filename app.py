import numpy as np
import math

import matplotlib as mpl
import matplotlib.pyplot as plt


# ----------- CONSTANTS -----------
# Initial x, y values
Y = 2.0
X = 0.0
# Grid parameters
ε = 0.00000001
h = 0.0001
a = 1000
# Oscillating parameters
A = 0.9
ω = 10
# Time parameters
START = 0
END = 200
# ----------------------------------


class vector:
    y = 0
    x = 0
    t = 0


class support_vector:
    y = 0
    x = 0


class forced_oscillator:
    v1 = support_vector()
    v2 = support_vector()
    v3 = vector()

    A = [[0] * 2, [0] * 2]
    F = [0] * 2
    
    def __init__(self):
        self.v1.x = X
        self.v1.y = Y

        self.v2.x = X
        self.v2.y = Y 

        self.v3.x = X
        self.v3.y = Y 
        self.v3.t = START


    def get_v(self):
        return self.v3


    def set_v(self, y, x):
        self.v1.x = x
        self.v1.y = y


    def set_time(self, t):
        self.v3.t = t


    def is_equel_epsilon(self):
        return ((self.v3.y - self.v2.y) > ε) or ((self.v3.x - self.v2.x) > ε)


    def count_newton(self):
        self.count_inv_matrix()
        self.count_func()

        self.v3.y = self.v2.y - self.A[0][0]*self.F[0] - self.A[0][1]*self.F[1]
        self.v3.x = self.v2.x - self.A[1][0]*self.F[0] - self.A[1][1]*self.F[1]
    

    def count_inv_matrix(self):
        self.v2.y = self.v3.y
        self.v2.x = self.v3.x

        det = h*a*(self.v2.y**2 - 1) + 1 + a*h*h
        self.A[0][0] = (1/det)*(-1 - a*h*(self.v2.y**2 - 1))
        self.A[0][1] = (1/det)*(a*h)
        self.A[1][0] = (1/det)*(-h)
        self.A[1][1] = (1/det)*(-1)


    def count_func(self):
        self.F[0] = float(self.v1.y + a*h*(-((self.v2.y**3)/3 - self.v2.y) + self.v2.x) - self.v2.y)
        self.F[1] = float(self.v1.x - h*self.v2.y + h*A*math.cos(ω*self.v3.t) - self.v2.x)  
        


if __name__ == "__main__":
    f_o = forced_oscillator()
    x, y, t = [], [], []

    while (f_o.get_v().t < END):
        while True:
            f_o.count_newton()
            if f_o.is_equel_epsilon:
                break

        f_o.set_v(f_o.get_v().y, f_o.get_v().x)
        f_o.set_time(f_o.get_v().t + h)

        x.append(f_o.get_v().x)
        y.append(f_o.get_v().y)
        t.append(f_o.get_v().t)

    plt.plot(y, x)
    plt.show() # Plot (x,y)
    plt.plot(t, y)
    plt.show() # Plot (y,t)
    plt.plot(t, x)
    plt.show() # Plot (x, t)