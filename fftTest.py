import numpy as np
import matplotlib.pyplot as plt
# a = [[1,2,3,4,5,6,7,8],
#     [1,2,3,4,5,6,7,8]]
a = np.loadtxt("rhino.txt", delimiter=",")

f = a[:,0] + 1j * a[:,1]

plt.plot(f.real, f.imag)
# plt.show()

n = 512

dft = (1/n) * np.fft.fft(f)
# b = np.abs(dft)
# c = np.fft.fftshift(b)
# plt.figure()
# plt.plot(c.real)
# plt.show()

tdft = dft
# tdft[0] = 0

idft = n * np.fft.ifft(tdft)
plt.figure()
plt.plot(idft.real, idft.imag)
plt.show()

# b = np.fft.fftshift(b)
# c = n * np.fft.ifft2(b)
# c = np.fft.ifftshift(c)

a = [[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8],
[1,2,3,4,5,6,7,8]]

N = 8*8
d = (1/N) * np.fft.fft2(a)
i = N * np.fft.ifft2(d)
print(i)
