import matplotlib.pyplot as plt
import math

dt = 0.001
Cm = 1
VNa = 50
VK = -77
VL = -54.387
gNa = 120
gK = 36
gL = 0.3
v = -65
m = 0.053
h = 0.6
n = 0.318

alfa = 2
beta = 1
Esyn = 0
T = 0
r = 0
gsyn = 1
Isyn = 0

vres = []
mres = []
hres = []
nres = []
time = []
rres = []
Isynres = []

f = open("h-h-model.txt", "w")
f.write(f'{str(v)},{str(m)},{str(h)},{str(n)},0,{str(r)}')

for i in range(int(100/dt)):
    milisegundo = int(i*dt)
    if milisegundo == 30 or milisegundo == 40 or milisegundo == 50:
        T = 1
    else:
        T = 0

    vres.append(v)
    mres.append(m)
    hres.append(h)
    nres.append(n)
    time.append(i*dt)
    rres.append(r)
    Isynres.append(Isyn)

    aM = (0.1 * (-v-40)) / (math.exp((-v-40)/10) - 1)
    bM = 4 * math.exp((-v-65)/18)
    aH = 0.07 * math.exp((-v-65)/20)
    bH = 1 / (math.exp((-v-35)/10) + 1)
    aN = (0.01 * (-v-55)) / (math.exp((-v-55)/10) - 1)
    bN = 0.125 * math.exp((-v-65)/80)

    Isyn = gsyn * r * (v - Esyn)
    vNext = v + dt*(-gL*(v - VL) - gNa*(m**3)*h*(v - VNa) - gK*(n**4)*(v - VK) - Isyn)
    mNext = m + dt*(aM*(1-m) - bM*m)
    hNext = h + dt*(aH*(1-h) - bH*h)
    nNext = n + dt*(aN*(1-n) - bN*n)
    r = r + dt*((alfa*T) * (1-r) - (beta*r))

    v = vNext
    m = mNext
    h = hNext
    n = nNext

    f.write(f'{str(v)},{str(m)},{str(h)},{str(n)},{str(i*dt)},{str(r)}\n')
f.close()

fig, axes = plt.subplots(2,1)
axes[0].plot(time, vres, label='Potencial de membrana')
axes[0].plot(time, Isynres, label = 'Corriente sináptica')
axes[0].set_title('Potencial de membrana / t', fontsize=24, loc='left')
axes[0].legend(loc="upper right", fontsize=18)
axes[0].set_ylabel('Potencial (mV)', fontsize=18)
axes[0].set_xlabel('Tiempo (ms)', fontsize=18)

axes[1].plot(time, rres)
axes[1].set_title('r / t', fontsize=24, loc='left')
axes[1].set_ylabel('r', fontsize=24)
axes[1].set_xlabel('Tiempo (ms)', fontsize=18)
axes[1].set_ylim(ymin=0, ymax=1)

plt.show()
