import numpy as np

for num in range(15):
    x = np.random.random(2**(num+1))
    with open("data" + str(num) + ".txt", "w") as file:
        for i in x:
            file.write(str(i))
            file.write("\n")

    xr = np.fft.fft(x)
    with open("expected" + str(num) + ".txt", "w") as file:
        for i in xr:
            file.write(str(i.real) + " " + str(i.imag))
            file.write("\n")


