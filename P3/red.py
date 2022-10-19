import matplotlib.pyplot as plt
import math
import numpy as np

dt = 0.001
u = 0.0021
s = 4
Esyn = -1.92

#PD1
xPD1 = 1.605030702160478
yPD1 = -4.737031630427527
zPD1 = 3.308475838267355
xresPD1 = []
yresPD1 = []
zresPD1 = []

#PD2
xPD2 = -0.8237772862346795
yPD2 = -2.6984138562485125
zPD2 = 3.134926915912511
xresPD2 = []
yresPD2 = []
zresPD2 = []

#AB
xAB = 0.5804651879462951
yAB = -5.5582642862033795
zAB = 3.1741742459130817
xresAB = []
yresAB = []
zresAB = []

#LP
xLP = -0.945
yLP = -3.4594
zLP = 3.350
mslowLP = 0.5
xresLP = []
yresLP = []
zresLP = []

#PY
xPY = -0.945
yPY = -3.4594
zPY = 3.350
mslowPY = 0.5
xresPY = []
yresPY = []
zresPY = []

#regular
e = 3
#irregular
#e = 3.281

time = []
contador = 0
for i in range(int(4000/dt)):
    contador += 1
    milisegundo = int(i*dt)
    time.append(i*dt)

    #PD1
    xresPD1.append(xPD1)
    yresPD1.append(yPD1)
    zresPD1.append(zPD1)
    IPD1 = (xPD1 - xPD2)*0.332 +  (xPD1 - xAB)*0.325
    xNextPD1 = xPD1 + dt*(yPD1 + 3*xPD1**2 - xPD1**3 - zPD1 + e - IPD1)
    yNextPD1 = yPD1 + dt*(1 - 5*xPD1**2 - yPD1)
    zNextPD1 = zPD1 + dt*(-zPD1 + s*(xPD1 + 1.6))*u
    xPD1 = xNextPD1
    yPD1 = yNextPD1
    zPD1 = zNextPD1

    #PD2
    xresPD2.append(xPD2)
    yresPD2.append(yPD2)
    zresPD2.append(zPD2)
    IPD2 = (xPD2 - xPD1)*0.332 + (xPD2 - xAB)*0.548
    xNextPD2 = xPD2 + dt*(yPD2 + 3*xPD2**2 - xPD2**3 - zPD2 + e - IPD2)
    yNextPD2 = yPD2 + dt*(1 - 5*xPD2**2 - yPD2)
    zNextPD2 = zPD2 + dt*(-zPD2 + s*(xPD2 + 1.6))*u
    xPD2 = xNextPD2
    yPD2 = yNextPD2
    zPD2 = zNextPD2

    #AB
    xresAB.append(xAB)
    yresAB.append(yAB)
    zresAB.append(zAB)
    IAB = (xAB - xPD1)*0.325 + (xAB - xPD2)*0.548
    IAB += (0.585*(xAB - Esyn) / (1 + np.exp(0.44*(-1.66 - xLP))))
    xNextAB = xAB + dt*(yAB + 3*xAB**2 - xAB**3 - zAB + e - IAB)
    yNextAB = yAB + dt*(1 - 5*xAB**2 - yAB)
    zNextAB = zAB + dt*(-zAB + s*(xAB + 1.6))*u
    xAB = xNextAB
    yAB = yNextAB
    zAB = zNextAB

    #LP
    xresLP.append(xLP)
    yresLP.append(yLP)
    zresLP.append(zLP)
    mslowLP = mslowLP + dt*(((0.74*(1-mslowLP)) / (1 + np.exp(-1.74-xAB))) - (0.007*mslowLP))
    ILP = (0.029*mslowLP*(xLP - Esyn))
    ILP += (0.112*(xLP - Esyn) / (1 + np.exp(0.44*(-1.66 - xAB))))
    ILP += (0.186*(xLP - Esyn) / (1 + np.exp(0.44*(-1.66 - xPY))))
    xNextLP = xLP + dt*(yLP + 3*xLP**2 - xLP**3 - zLP + e - ILP)
    yNextLP = yLP + dt*(1 - 5*xLP**2 - yLP)
    zNextLP = zLP + dt*(-zLP + s*(xLP + 1.6))*u
    xLP = xNextLP
    yLP = yNextLP
    zLP = zNextLP

    #PY
    xresPY.append(xPY)
    yresPY.append(yPY)
    zresPY.append(zPY)
    mslowPY = mslowPY + dt*(((0.74*(1-mslowPY)) / (1 + np.exp(-1.74-xAB))) - (0.015*mslowPY))
    IPY = (0.032*mslowPY*(xPY - Esyn))
    IPY += (0.120*(xPY - Esyn) / (1 + np.exp(0.44*(-1.66 - xAB))))
    IPY += (0.241*(xPY - Esyn) / (1 + np.exp(0.44*(-1.66 - xLP))))
    xNextPY = xPY + dt*(yPY + 3*xPY**2 - xPY**3 - zPY + e - IPY)
    yNextPY = yPY + dt*(1 - 5*xPY**2 - yPY)
    zNextPY = zPY + dt*(-zPY + s*(xPY + 1.6))*u
    xPY = xNextPY
    yPY = yNextPY
    zPY = zNextPY

contador = int(contador*0.75)
# plt.plot(time, xresPD1, label='PD1')
# plt.plot(time, xresPD2, label='PD2')
plt.plot(time[contador:], xresAB[contador:], label ='AB')
plt.plot(time[contador:], xresLP[contador:], label='LP')
plt.plot(time[contador:], xresPY[contador:], label ='PY')
plt.legend(loc='upper right')
plt.title('Potencial de membrana / t', fontsize=24)
plt.ylabel('Potencial (mV)', fontsize=18)
plt.xlabel('Tiempo (ms)', fontsize=18)

plt.show()
