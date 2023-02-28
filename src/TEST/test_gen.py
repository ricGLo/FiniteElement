import numpy as np

N = 20
a = 0
b = 1

x = np.linspace(a, b, num = N)

with open("prueba.txt", 'w') as f:
    f.write(f"{N} {N-1} \n")

    for coord in x:
        f.write(f"{coord}\n")

    for n in range(N-1):
        f.write(f"{n} {n+1}\n")
        

