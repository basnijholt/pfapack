import time
import numpy as np
import matplotlib.pyplot as plt
from pfapack import pfaffian as pfa
from pfapack.ctypes import pfaffian as cpfaffian
from pfapack.ctypes import pfaffian_batched as cpfaffian_batched
from pfapack.ctypes import pfaffian_batched_4d as cpfaffian_batched_4d

def measure_performance(num_replicas, N, num_iterations=50):
    M = N + 1  # Number of matrices per replica
    total_matrices = num_replicas * M
    
    times = [[] for _ in range(4)]  # One list for each method
    results = [None] * 4  # To store results for verification

    for _ in range(num_iterations):
        A = np.random.randn(num_replicas, M, N, N) + 1j*np.random.randn(num_replicas, M, N, N)
        A = A - np.transpose(A, (0, 1, 3, 2))

        # 1. Fully looped (Python implementation)
        start_time = time.time()
        results[0] = [[pfa.pfaffian(a) for a in replica] for replica in A]
        times[0].append(time.time() - start_time)

        # 2. Fully looped (C implementation)
        start_time = time.time()
        results[1] = [[cpfaffian(a) for a in replica] for replica in A]
        times[1].append(time.time() - start_time)

        # 3. Partially batched (3D, loop over replicas)
        start_time = time.time()
        results[2] = [cpfaffian_batched(replica) for replica in A]
        times[2].append(time.time() - start_time)

        # 4. Fully batched (4D)
        start_time = time.time()
        results[3] = cpfaffian_batched_4d(A)
        times[3].append(time.time() - start_time)

        # Verify that all methods produce the same results
        for i in range(1, 4):
            np.testing.assert_allclose(np.array(results[0]), np.array(results[i]), rtol=1e-13, atol=1e-13)

    # Calculate average times
    avg_times = [np.mean(method_times) for method_times in times]
    return avg_times, total_matrices

# Set up the test
np.random.seed(0)
N = 8  # Size of each matrix
replica_counts = [1, 5, 10, 50, 100, 500, 1000]
methods = ['Python (looped)', 'C (looped)', 'C (batched 3D)', 'C (batched 4D)']

# Measure performance
results = [measure_performance(count, N) for count in replica_counts]
times, total_matrices = zip(*results)

# Plot results
plt.figure(figsize=(4.75, 3.25))
for i, method in enumerate(methods):
    method_times = [result[i] for result in times]
    plt.loglog(total_matrices, method_times, marker='o', label=method)

plt.xlabel('Total number of matrices')
plt.ylabel('Average time (seconds)')
plt.title(f'Pfaffian Calculation Performance ({N}x{N})')
plt.legend(fontsize='x-small')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.tight_layout()
plt.savefig('pfaffian_performance.png', dpi=600, bbox_inches='tight')
plt.close()

# Print results for the largest number of matrices
print(f"Performance for {total_matrices[-1]} matrices:")
print(f"Structure: {replica_counts[-1]} replicas x {N+1} matrices per replica")
print(f"Each matrix: {N}x{N} (complex-valued)")
print("\n{:<20} {:>15} {:>15}".format("Method", "Time (µs)", "Speedup"))
print("-" * 53)

python_time = times[-1][0]
for method, time in zip(methods, times[-1]):
    speedup = python_time / time
    print("{:<20} {:>15.2f} {:>15.2f}x".format(method, time*1e6, speedup))
