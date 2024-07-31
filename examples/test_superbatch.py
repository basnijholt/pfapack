import time
import numpy as np
from pfapack.ctypes import pfaffian as cpfaffian
from pfapack.ctypes import pfaffian_batched as cpfaffian_batched
from pfapack.ctypes import pfaffian_batched_4d as cpfaffian_batched_4d

np.random.seed(0)

num_replicas = 5
num_walkers = 10
num_shadows = 1

M = num_walkers * num_shadows 
N = 8 

A = np.random.randn(num_replicas, M, N, N) + 1j*np.random.randn(num_replicas, M, N, N)
A = A - np.transpose(A, (0, 1, 3, 2))

print("Complex-valued matrices")

start_time = time.time()
pfaffians = [[cpfaffian(a) for a in replica] for replica in A]
end_time = time.time()
print(f"Time taken to compute {num_replicas}x{M} {N}x{N} Pfaffians (looped):     {(end_time - start_time)*1e6:14.2f} microseconds")

start_time = time.time()
pfaffians = [cpfaffian_batched(replica) for replica in A]
end_time = time.time()
print(f"Time taken to compute {num_replicas}x{M} {N}x{N} Pfaffians (batched 3D): {(end_time - start_time)*1e6:14.2f} microseconds")

start_time = time.time()
pfaffians = cpfaffian_batched_4d(A)
end_time = time.time()
print(f"Time taken to compute {num_replicas}x{M} {N}x{N} Pfaffians (batched 4D): {(end_time - start_time)*1e6:14.2f} microseconds")

start_time = time.time()
determinants = [[np.linalg.det(a) for a in replica] for replica in A]
end_time = time.time()
print(f"Time taken to compute {num_replicas}x{M} {N}x{N} determinants (looped):  {(end_time - start_time)*1e6:14.2f} microseconds")

start_time = time.time()
determinants = [np.linalg.det(replica) for replica in A]
end_time = time.time()
print(f"Time taken to compute {num_replicas}x{M} {N}x{N} determinants (batched 3D): {(end_time - start_time)*1e6:14.2f} microseconds")

# NumPy doesn't have a built-in 4D determinant function, so we'll use the 3D version in a loop
start_time = time.time()
determinants = np.array([np.linalg.det(replica) for replica in A])
end_time = time.time()
print(f"Time taken to compute {num_replicas}x{M} {N}x{N} determinants (batched 4D): {(end_time - start_time)*1e6:14.2f} microseconds")
