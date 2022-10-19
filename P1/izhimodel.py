import matplotlib.pyplot as plt

v = -65
I = 10
dt = 0.1
a = 0.02
b = 0.2
c = -50
d = 2
u = b * v
result = []
time = []

f = open("izhimodel"+str(dt)+".txt", "w")
f.write(f'{str(v)}\t0\n')
for i in range(6000):
    vNext = v + dt*(0.04*v*v + 5*v +140 -u + I)
    uNext = u + dt*a*(b*v - u)
    if vNext >= 30:
        vNext = c
        uNext = uNext + d
    v = vNext
    u = uNext

    f.write(f'{str(v)}\t{str(i*dt)}\n')
    result.append(v)
    time.append(i*dt)

f.close()

plt.plot(time, result)
plt.show()
