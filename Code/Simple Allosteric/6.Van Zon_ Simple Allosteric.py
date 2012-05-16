from matplotlib.pylab import plot, show
from scipy.integrate import odeint
from numpy import arange, zeros, size
from time import sleep

kpf = 13.6
kdps = 0.908
f6 = kdps
b0 = kdps
kAf = 3.45e7

t = arange(0,96,0.1)
#structure of y: concentrations C, AC, C', A, molar
y0 = [0.58,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.012]

def rate(y,t0):
  C = y[0:7]
  AC = y[7:13]
  Cp = y[13:20]
  A = y[20]
  
  yp = zeros(size(y))
  for i in range(0,6):
    kAb = 10**(i+1)
    if i > 0:
      yp[i] = kpf*AC[i-1] + kAb*AC[i] - kAf*A*C[i]
      yp[13+i] = kdps*(Cp[i+1]-Cp[i])
    yp[i+7] = -kpf*AC[i] - kAb*AC[i] + kAf*A*C[i]
    yp[20] += -kAf*A*C[i] + (kAb+kpf)*AC[i]
  yp[0] = 10*AC[0] - kAf*A*C[0] + b0*Cp[0]
  yp[6] = kpf*AC[5]-f6*C[6]
  yp[13] = kdps*Cp[1] - b0*Cp[0]
  yp[19] = -kdps*Cp[6] + f6*C[6]
  
  #input pathway
  if abs(t0 - 30) < 0.1 :
    yp[20] += 0.1
  
  return yp
  
y = odeint(rate, y0, t)
y = y.transpose()
C = y[0:7]
AC = y[7:13]
Cp = y[13:20]
A = y[20]
  
plot(t,C[6])
show()
