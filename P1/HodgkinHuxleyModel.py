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

Iext = 0

vres = []
mres = []
hres = []
nres = []
time = []

f = open("h-h-model.txt", "w")
f.write(f'{str(v)},{str(m)},{str(h)},{str(n)},0')

for i in range(int(100/dt)):

    vres.append(v)
    mres.append(m)
    hres.append(h)
    nres.append(n)
    time.append(i*dt)

    aM = (0.1 * (-v-40)) / (math.exp((-v-40)/10) - 1)
    bM = 4 * math.exp((-v-65)/18)
    aH = 0.07 * math.exp((-v-65)/20)
    bH = 1 / (math.exp((-v-35)/10) + 1)
    aN = (0.01 * (-v-55)) / (math.exp((-v-55)/10) - 1)
    bN = 0.125 * math.exp((-v-65)/80)

    vNext = v + dt*(Iext - gL*(v - VL) - gNa*(m**3)*h*(v - VNa) - gK*(n**4)*(v - VK))
    mNext = m + dt*(aM*(1-m) - bM*m)
    hNext = h + dt*(aH*(1-h) - bH*h)
    nNext = n + dt*(aN*(1-n) - bN*n)

    v = vNext
    m = mNext
    h = hNext
    n = nNext

    f.write(f'{str(v)},{str(m)},{str(h)},{str(n)},{str(i*dt)}\n')

    if i*dt > 50:
        Iext = 0

f.close()


fig, axes = plt.subplots(2,1)
axes[0].plot(time, vres)
axes[0].set_title('Potencial de membrana / t')
axes[0].set_ylabel('Potencial (mV)')
axes[0].set_xlabel('Tiempo (ms)')
#axes[0].set_ylim(ymin=-85, ymax=50)

axes[1].plot(time, mres, label = "m")
axes[1].plot(time, hres, label = "h")
axes[1].plot(time, nres, label = "n")
axes[1].legend(loc="upper right")
axes[1].set_title('Variables de conductancia / t')
axes[1].set_ylabel('Conductancia')
axes[1].set_xlabel('Tiempo (ms)')

plt.show()
