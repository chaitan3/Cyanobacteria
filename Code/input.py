from matplotlib.pyplot import plot, show, gca, fill_between, xlabel, ylabel
from scipy.integrate import odeint
from numpy import arange, zeros, size
from time import sleep

#Todo: KaiBC feedback
#Problems
#Too high generation
#Include Quinone reaction for input?
#Tweaking parameters?

#Verification
#kaiA matching - cikA paper
#kaiC matching - rhythms paper
#response to a dark pulse
#back to dark from light

#Van Zon Allosteric model constants
kpf = 13.6
kdps = 0.908
f6 = kdps
b0 = kdps
kAf = 3.45e7

#kaiA degradation rate
kaiAd =4e-3
KkaiAd = 0.001
#kaiA generation
kaiAg = 9e-3
KkaiAg = 0.14
#Dark deactivation of kaiA
lAd = 6e-3
KlAd = 0.012

#input pathway, every 12 hrs
#l = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1]
l = [1,0,1,1,1,1,1,1,1,1,1,1,0,1,0,1]

#Simulation time
hours = 192

t = arange(0,hours,0.1)
#structure of y: concentrations C, AC, C', A, molar
#y0 = [0.58,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.012]
y0 = [  2.08958651e-03,   1.22275346e-02,   1.37539199e-01,   2.16654131e-01,
   6.30963569e-02,   2.94735721e-02,   4.74654668e-05,   3.02402796e-03,
   3.64195579e-03,   4.58027901e-03,   7.30463031e-04,   2.13005592e-05,
   9.95118144e-07,   4.29134314e-02,   3.00661718e-02,   1.84397044e-02,
   9.63880083e-03,   4.13111399e-03,   1.36711226e-03,   3.16798254e-04,
   9.78653245e-07]

def rate(y,t0):
  
 #Extract the variables
  C = y[0:7]
  AC = y[7:13]
  Cp = y[13:20]
  A = y[20]
  
  #initialise the derivative and assign values acc the Van Zon model
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
  
  
  #Transcription
  #kaiA generation
  yp[20] += kaiAg*C[6]/(C[6]+ KkaiAg)
  
  #kaiA degradation
  yp[20] -= kaiAd*A/(A+KkaiAd)
  
  #Light input
  i = int(t0/12);
  if i == hours/12:
    i -= 1
  light = l[i]
  
  #Dark quinone degradation of kaiA
  if light == 0:
    yp[20] -= lAd*A/(A+KlAd)
  
  return yp

#Solve the ode  
y = odeint(rate, y0, t)
y = y.transpose()
#Extract the variables
C = y[0:7]
AC = y[7:13]
Cp = y[13:20]
A = y[20]

#Plot the 12 hour grid
gca().set_xticks(range(12,hours+1,12))
gca().grid(True,linestyle='-',color='0.75')
gca().yaxis.grid(False)

#Plot the values to Check
#plot(t,C[6])
plot(t,A+AC[0]+AC[1]+AC[2]+AC[3]+AC[4]+AC[5]) 
#savefig('1.png')

#Get the axis limits
li = gca().axis()

#Plot the LD cycles
for i in range(0, len(l)):
	if l[i] == 0:
		fill_between([12*i, 12*(i+1)], [1, 1], alpha=0.1)

#restore the axis limits
gca().axis(li)
xlabel('hours')
ylabel('concentration')

show()
