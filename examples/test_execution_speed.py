import time
import numpy as np
from pfapack.ctypes import pfaffian as cpfaffian
from pfapack.ctypes import pfaffian_batched as cpfaffian_batched

np.random.seed(0)

num_walkers = 10
num_shadows = 1

M = num_walkers * num_shadows 
N = 8 

A = np.random.randn(M, N, N) + 1j*np.random.randn(M, N, N)

A = A - np.transpose(A, (0,2,1))


start_time = time.time()
pfaffians = [cpfaffian(a) for a in A]
end_time = time.time()

print("Complex-valued matrices")

print(f"Time taken to compute {M} {N}x{N} Pfaffians (looped):     {(end_time - start_time)*1e6:14.2f} microseconds")

start_time = time.time()
pfaffians = cpfaffian_batched(A)
end_time = time.time()

print(f"Time taken to compute {M} {N}x{N} Pfaffians (batched):    {(end_time - start_time)*1e6:14.2f} microseconds")


start_time = time.time()
determinants = [np.linalg.det(a) for a in A]
end_time = time.time()

print(f"Time taken to compute {M} {N}x{N} determinants (looped):  {(end_time - start_time)*1e6:14.2f} microseconds")


start_time = time.time()
determinants = np.linalg.det(A)
end_time = time.time()

print(f"Time taken to compute {M} {N}x{N} determinants (batched): {(end_time - start_time)*1e6:14.2f} microseconds")
