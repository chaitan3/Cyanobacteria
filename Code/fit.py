from matplotlib.pyplot import plot, show, gca, fill_between, xlabel, ylabel
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from numpy import arange, zeros, size, array
from time import sleep

#structure of y: concentrations C, AC, C', A, molar
#y0 = [0.58,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.012]
y0 = [  2.08958651e-03,   1.22275346e-02,   1.37539199e-01,   2.16654131e-01,
   6.30963569e-02,   2.94735721e-02,   4.74654668e-05,   3.02402796e-03,
   3.64195579e-03,   4.58027901e-03,   7.30463031e-04,   2.13005592e-05,
   9.95118144e-07,   4.29134314e-02,   3.00661718e-02,   1.84397044e-02,
   9.63880083e-03,   4.13111399e-03,   1.36711226e-03,   3.16798254e-04,
   9.78653245e-07]

def model(ti, kpf, kdps, kAf):
  
  #Van Zon Allosteric model constants
  #kpf = 13.6
  #kdps = 0.908
  f6 = kdps
  b0 = kdps
  #kAf = 3.45e7
  print kpf, kdps, kAf

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
    
    return yp

  #Solve the ode  
  yc = odeint(rate, y0, ti)
  yc = yc.transpose()

  return yc[6]

ydata = array([  4.74654668e-05,   4.46648372e-05,   4.21604873e-05,
         3.99305985e-05,   3.79555773e-05,   3.62178640e-05,
         3.47017723e-05,   3.33933466e-05,   3.22803153e-05,
         3.13520020e-05,   3.05992590e-05,   3.00144371e-05,
         2.95913552e-05,   2.93252923e-05,   2.92129947e-05,
         2.92526973e-05,   2.94441593e-05,   2.97887146e-05,
         3.02893355e-05,   3.09506998e-05,   3.17792642e-05,
         3.27833339e-05,   3.39731335e-05,   3.53608223e-05,
         3.69604878e-05,   3.87880855e-05,   4.08612629e-05,
         4.31991708e-05,   4.58221214e-05,   4.87512318e-05,
         5.20080563e-05,   5.56142884e-05,   5.95916369e-05,
         6.39619677e-05,   6.87477909e-05,   7.39731290e-05,
         7.96647276e-05,   8.58536038e-05,   9.25768147e-05,
         9.98792628e-05,   1.07815549e-04,   1.16452382e-04,
         1.25870547e-04,   1.36167527e-04,   1.47459879e-04,
         1.59886529e-04,   1.73612627e-04,   1.88834004e-04,
         2.05782834e-04,   2.24734206e-04,   2.46015120e-04,
         2.70014181e-04,   2.97195018e-04,   3.28109438e-04,
         3.63419741e-04,   4.03918256e-04,   4.50554933e-04,
         5.04465937e-04,   5.67008077e-04,   6.39797774e-04,
         7.24740758e-04,   8.24065037e-04,   9.40337711e-04,
         1.07646169e-03,   1.23564420e-03,   1.42132182e-03,
         1.63705112e-03,   1.88638703e-03,   2.17272401e-03,
         2.49917196e-03,   2.86847313e-03,   3.28297871e-03,
         3.74469005e-03,   4.25535933e-03,   4.81661494e-03,
         5.43011461e-03,   6.09768831e-03,   6.82147512e-03,
         7.60403686e-03,   8.44846105e-03,   9.35844679e-03,
         1.03383836e-02,   1.13934222e-02,   1.25295504e-02,
         1.37536638e-02,   1.50736466e-02,   1.64984514e-02,
         1.80381843e-02,   1.97041783e-02,   2.15090546e-02,
         2.34667511e-02,   2.55924944e-02,   2.79026806e-02,
         3.04146237e-02,   3.31461049e-02,   3.61146410e-02,
         3.93364315e-02,   4.28248998e-02,   4.65889309e-02,
         5.06305507e-02,   5.49428067e-02,   5.95077186e-02,
         6.42950850e-02,   6.92623258e-02,   7.43559461e-02,
         7.95143712e-02,   8.46718030e-02,   8.97624810e-02,
         9.47246737e-02,   9.95036464e-02,   1.04053651e-01,
         1.08338689e-01,   1.12332327e-01,   1.16016812e-01,
         1.19381784e-01,   1.22422937e-01,   1.25140449e-01,
         1.27537733e-01,   1.29620254e-01,   1.31394553e-01,
         1.32867454e-01,   1.34045419e-01,   1.34934087e-01,
         1.35537859e-01,   1.35859644e-01,   1.35900649e-01,
         1.35660239e-01,   1.35135874e-01,   1.34323104e-01,
         1.33215661e-01,   1.31805671e-01,   1.30084056e-01,
         1.28041194e-01,   1.25667893e-01,   1.22956847e-01,
         1.19904563e-01,   1.16513780e-01,   1.12796038e-01,
         1.08774166e-01,   1.04483769e-01,   9.99733046e-02,
         9.53021797e-02,   9.05367942e-02,   8.57457709e-02,
         8.09947806e-02,   7.63422863e-02,   7.18365372e-02,
         6.75142946e-02,   6.34006278e-02,   5.95096830e-02,
         5.58459880e-02,   5.24059402e-02,   4.91792680e-02,
         4.61505737e-02,   4.33008784e-02,   4.06092014e-02,
         3.80543254e-02,   3.56167037e-02,   3.32803583e-02,
         3.10344171e-02,   2.88738360e-02,   2.67989866e-02,
         2.48141775e-02,   2.29255707e-02,   2.11391674e-02,
         1.94593323e-02,   1.78881170e-02,   1.64251943e-02,
         1.50682226e-02,   1.38133368e-02,   1.26556414e-02,
         1.15896263e-02,   1.06094960e-02,   9.70938102e-03,
         8.88351160e-03,   8.12631862e-03,   7.43249698e-03,
         6.79704168e-03,   6.21526426e-03,   5.68279564e-03,
         5.19557856e-03,   4.74986708e-03,   4.34219381e-03,
         3.96936462e-03,   3.62844465e-03,   3.31673795e-03,
         3.03176642e-03,   2.77125910e-03,   2.53312892e-03,
         2.31546164e-03,   2.11651439e-03,   1.93468538e-03,
         1.76850767e-03,   1.61663581e-03,   1.47784304e-03,
         1.35100934e-03,   1.23510618e-03,   1.12919330e-03,
         1.03241180e-03,   9.43976196e-04,   8.63168443e-04,
         7.89331855e-04,   7.21866253e-04,   6.60222941e-04,
         6.03900576e-04,   5.52441056e-04,   5.05425581e-04,
         4.62471382e-04,   4.23228763e-04,   3.87378276e-04,
         3.54627860e-04,   3.24710586e-04,   2.97382934e-04,
         2.72422285e-04,   2.49625329e-04,   2.28805883e-04,
         2.09794126e-04,   1.92435523e-04,   1.76588397e-04,
         1.62123287e-04,   1.48922169e-04,   1.36877358e-04,
         1.25890526e-04,   1.15871922e-04,   1.06739598e-04,
         9.84186944e-05,   9.08409755e-05,   8.39442465e-05,
         7.76718254e-05,   7.19720045e-05,   6.67976875e-05,
         6.21059881e-05,   5.78579032e-05,   5.40179455e-05,
         5.05538330e-05,   4.74362703e-05,   4.46387171e-05,
         4.21371470e-05,   3.99098440e-05,   3.79372465e-05,
         3.62017994e-05,   3.46878164e-05,   3.33813719e-05,
         3.22702089e-05,   3.13436605e-05,   3.05925980e-05,
         3.00093888e-05,   2.95878672e-05,   2.93233263e-05,
         2.92125259e-05,   2.92537153e-05,   2.94466686e-05,
         2.97927364e-05,   3.02949057e-05,   3.09578658e-05,
         3.17880926e-05,   3.27939333e-05,   3.39856202e-05,
         3.53752930e-05,   3.69770974e-05,   3.88069937e-05,
         4.08826552e-05,   4.32232407e-05,   4.58490735e-05,
         4.87812683e-05,   5.20413878e-05,   5.56511243e-05,
         5.96321883e-05,   6.40064548e-05,   6.87964309e-05,
         7.40261726e-05,   7.97224462e-05,   8.59163272e-05,
         9.26449365e-05,   9.99531724e-05,   1.07895827e-04,
         1.16539841e-04,   1.25966152e-04,   1.36272239e-04,
         1.47574898e-04,   1.60013350e-04,   1.73752992e-04,
         1.88989985e-04,   2.05956814e-04,   2.24929279e-04,
         2.46234949e-04,   2.70262851e-04,   2.97476247e-04,
         3.28430027e-04,   3.63786883e-04,   4.04340423e-04,
         4.51042373e-04,   5.05029347e-04,   5.67663209e-04,
         6.40561190e-04,   7.25632361e-04,   8.25108385e-04,
         9.41560533e-04,   1.07789250e-03,   1.23731346e-03,
         1.42326171e-03,   1.63930035e-03,   1.88897935e-03,
         2.17569096e-03,   2.50254265e-03,   2.87227220e-03,
         3.28722714e-03,   3.74940684e-03,   4.26055955e-03,
         4.82231491e-03,   5.43633096e-03,   6.10444136e-03,
         6.82878696e-03,   7.61193469e-03,   8.45697773e-03,
         9.36762121e-03,   1.03482612e-02,   1.14040581e-02,
         1.25410050e-02,   1.37660139e-02,   1.50869796e-02,
         1.65128448e-02,   1.80537436e-02,   1.97210241e-02,
         2.15273139e-02,   2.34865644e-02,   2.56140123e-02,
         2.79260664e-02,   3.04400456e-02,   3.31737404e-02,
         3.61446506e-02,   3.93689693e-02,   4.28600975e-02,
         4.66268159e-02,   5.06711180e-02,   5.49859382e-02,
         5.95531855e-02,   6.43425342e-02,   6.93112765e-02,
         7.44058405e-02,   7.95645642e-02,   8.47216569e-02,
         8.98113950e-02,   9.47720626e-02,   9.95490107e-02,
         1.04096597e-01,   1.08378916e-01,   1.12369624e-01,
         1.16051028e-01,   1.19412868e-01,   1.22450871e-01,
         1.25165254e-01,   1.27559459e-01,   1.29638961e-01,
         1.31410304e-01,   1.32880313e-01,   1.34055447e-01,
         1.34941319e-01,   1.35542334e-01,   1.35861394e-01,
         1.35899672e-01,   1.35656520e-01,   1.35129378e-01,
         1.34313776e-01,   1.33203434e-01,   1.31790461e-01,
         1.30065778e-01,   1.28019752e-01,   1.25643201e-01,
         1.22928839e-01,   1.19873234e-01,   1.16479179e-01,
         1.12758347e-01,   1.08733663e-01,   1.04440873e-01,
         9.99285547e-02,   9.52561611e-02,   9.04902039e-02,
         8.56992799e-02,   8.09490191e-02,   7.62977502e-02,
         7.17936294e-02,   6.74732868e-02,   6.33617411e-02,
         5.94730149e-02,   5.58115518e-02,   5.23736447e-02,
         4.91489721e-02,   4.61221113e-02,   4.32740506e-02,
         4.05838032e-02,   3.80301572e-02,   3.55935891e-02,
         3.32581645e-02,   3.10130660e-02,   2.88533069e-02,
         2.67792988e-02,   2.47953892e-02,   2.29077463e-02,
         2.11223571e-02,   1.94435652e-02,   1.78733989e-02,
         1.64115184e-02,   1.50555580e-02,   1.38016408e-02,
         1.26448639e-02,   1.15797125e-02,   1.06003867e-02,
         9.70101992e-03,   8.87584520e-03,   8.11929397e-03,
         7.42606416e-03,   6.79115308e-03,   6.20987499e-03,
         5.67786480e-03,   5.19107040e-03,   4.74574266e-03,
         4.33842127e-03,   3.96591793e-03,   3.62529756e-03,
         3.31386390e-03,   3.02914101e-03,   2.76885821e-03,
         2.53093338e-03,   2.31345826e-03,   2.11468452e-03,
         1.93300350e-03,   1.76697164e-03,   1.61523387e-03,
         1.47656535e-03,   1.34984434e-03,   1.23404457e-03,
         1.12822395e-03,   1.03152605e-03,   9.43168278e-04,
         8.62431564e-04,   7.88660017e-04,   7.21250945e-04,
         6.59658898e-04,   6.03384434e-04,   5.51969253e-04,
         5.04994273e-04,   4.62076974e-04,   4.22868057e-04,
         3.87048288e-04,   3.54325989e-04,   3.24434594e-04,
         2.97130577e-04,   2.72191462e-04,   2.49414018e-04,
         2.28612625e-04,   2.09617699e-04,   1.92274347e-04,
         1.76441128e-04,   1.61988877e-04,   1.48799709e-04,
         1.36765914e-04,   1.25789155e-04,   1.15779817e-04,
         1.06655925e-04,   9.83425330e-05,   9.07715452e-05,
         8.38810592e-05,   7.76144068e-05,   7.19198618e-05,
         6.67503726e-05,   6.20631343e-05,   5.78191631e-05,
         5.39829649e-05,   5.05223017e-05,   4.74079310e-05,
         4.46133379e-05,   4.21145138e-05,   3.98897418e-05,
         3.79195055e-05,   3.61862592e-05,   3.46743339e-05,
         3.33698188e-05,   3.22604719e-05,   3.13356476e-05])

#Simulation time
hours = 48

t = arange(0,hours,0.1)    

popt, pcov = curve_fit(model, t, ydata, p0=[13.6, 0.908, 3.45e7])
print popt
print pcov
