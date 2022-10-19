import matplotlib.pyplot as plt
import math

dt = 0.001
u = 0.0021
s = 4
x = 1.559
y = -2.294
z = 3.084
contador = 0

#regular
#e = 3
#irregular
e = 3.281

xres = []
yres = []
zres = []
time = []

f = open("hr-model.txt", "w")
f.write(f'{str(x)},{str(y)},{str(z)},0')

for i in range(int(3000/dt)):
    milisegundo = int(i*dt)

    xres.append(x)
    yres.append(y)
    zres.append(z)
    time.append(i*dt)

    xNext = x + dt*(y + 3*x**2 - x**3 - z + e)
    yNext = y + dt*(1 - 5*x**2 - y)
    zNext = z + dt*(-z + s*(x + 1.6))*u

    x = xNext
    y = yNext
    z = zNext
    if contador == 100:
        contador = 0
        f.write(f'{str(x)},{str(y)},{str(z)},{str(i*dt)}\n')
    contador += 1
f.close()


plt.plot(time, xres)
plt.title('Potencial de membrana / t irregular', fontsize=24)
plt.ylabel('Potencial (mV)', fontsize=18)
plt.xlabel('Tiempo (ms)', fontsize=18)

plt.show()
