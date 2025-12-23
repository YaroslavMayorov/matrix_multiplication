# Matrix multiplication

## Operation Counts and Performance Analysis
### Theoretical and Observed Operation Counts

For a 4×4 matrix, the theoretical number of scalar operations for each algorithm, confirmed by explicit operation counting in the code, is as follows:

**Naive matrix multiplication**
- Multiplications: 64
- Additions: 48

**Strassen’s algorithm**
- Multiplications: 49
- Additions: 198

**Iterative Winograd algorithm**
- Multiplications: 48
- Additions: 144

The observed operation counts approximately matched the expected theoretical values.
<img width="879" height="574" alt="Screenshot 2025-12-23 at 19 22 26" src="https://github.com/user-attachments/assets/a7c05ba5-9fd6-4828-a299-249f07d34fb0" />
<img width="321" height="241" alt="Screenshot 2025-12-23 at 19 22 32" src="https://github.com/user-attachments/assets/4bc43aa0-b44a-4f49-9c20-8fb4ad86cadd" />

<img width="868" height="572" alt="Screenshot 2025-12-23 at 19 23 49" src="https://github.com/user-attachments/assets/338c5499-4811-4167-b213-706b9f2866d1" />
<img width="323" height="232" alt="Screenshot 2025-12-23 at 19 23 55" src="https://github.com/user-attachments/assets/408442e3-9f4f-48d0-bf0d-64e97e3bae5b" />


### Effect of Reducing Multiplications on Runtime

In experiments on matrix sizes ranging from 4×4 to 512×512, the following consistent performance ordering was observed:

Winograd → fastest,
Naive → intermediate,
Strassen → slowest.

Reducing the number of multiplications generally led to faster runtime, but only when the reduction was not offset by a large increase in additions and algorithmic overhead. The Winograd algorithm provided the best balance between fewer multiplications and a moderate increase in additions, resulting in superior performance across all tested data types and matrix structures.

<img width="868" height="572" alt="Screenshot 2025-12-23 at 19 24 33" src="https://github.com/user-attachments/assets/5f0223b4-ae59-4da3-9db8-13104ac0ee1e" />
<img width="864" height="564" alt="Screenshot 2025-12-23 at 19 24 42" src="https://github.com/user-attachments/assets/e7e98fab-f815-4268-819d-a65f6972ca87" />
<img width="864" height="571" alt="Screenshot 2025-12-23 at 19 24 53" src="https://github.com/user-attachments/assets/f3e0fa70-856a-44d1-8339-daf10dabd957" />
<img width="858" height="562" alt="Screenshot 2025-12-23 at 19 25 03" src="https://github.com/user-attachments/assets/40926b60-462d-4a42-a65f-33a0cc9b5884" />
<img width="852" height="553" alt="Screenshot 2025-12-23 at 19 25 20" src="https://github.com/user-attachments/assets/9ba7a446-1637-4fbb-91ea-0bbaf34c5fd8" />
<img width="850" height="553" alt="Screenshot 2025-12-23 at 19 25 28" src="https://github.com/user-attachments/assets/cebb4643-759a-4977-a0ce-16df0be5fe3f" />
<img width="851" height="553" alt="Screenshot 2025-12-23 at 19 25 34" src="https://github.com/user-attachments/assets/22966b13-c7cc-4fbd-9c86-45ae876bb829" />
<img width="851" height="554" alt="Screenshot 2025-12-23 at 19 25 41" src="https://github.com/user-attachments/assets/b2fb5efa-6bd0-4e24-a893-ba2804898c18" />
<img width="851" height="543" alt="Screenshot 2025-12-23 at 19 25 49" src="https://github.com/user-attachments/assets/f17800c0-86d0-44d6-aa81-596b7a373963" />


### Why Strassen Was Slower Despite Fewer Multiplications

Although Strassen’s algorithm reduces the number of multiplications, it performed worse than the naive algorithm in all tested cases due to several factors:

- **Large increase in additions**
For a 4×4 matrix, saving 15 multiplications comes at the cost of approximately 150 additional additions. While additions are cheaper than multiplications on modern CPUs/GPUs, they are not free. At small and medium matrix sizes, the cost of these extra additions outweighs the benefit of fewer multiplications.

- **Recursion and memory overhead**
Strassen relies on recursive calls and the creation of multiple temporary submatrices. This introduces:

- function call overhead,

- additional memory allocations,

- increased memory traffic.

For matrix sizes up to 512×512, these overheads dominate the theoretical asymptotic advantage.

Poor cache locality
The intermediate matrix combinations in Strassen disrupt linear memory access patterns, reducing cache efficiency and limiting vectorization opportunities.

As a result, Strassen typically becomes advantageous only for much larger matrix sizes than those tested.

### Why the Performance Gap Was Smallest for Complex Numbers

The smallest performance difference between Winograd and the naive algorithm was observed for complex-valued matrices. This is expected because a single complex multiplication already consists of multiple real multiplications and additions. Consequently, the relative benefit of saving a small number of complex multiplications is reduced, while the overhead of additional additions becomes more significant.

In other words, when the base arithmetic operation is already expensive, algorithmic optimizations at the scalar level yield a smaller relative performance gain.

## Memory Overhead and Access Patterns

The memory requirements of the three implemented algorithms differ significantly and directly affect their performance characteristics.

Memory usage by algorithm

**Naive multiplication**:

Uses approximately 3n^2 memory: matrices A, B, and the result C.
No additional scratch matrices are required.

**Iterative Winograd**:

Uses approximately 3n^2 + 2n memory: matrices A, B, C, plus two auxiliary vectors (`row_factors` and `col_factors`) of size n.

**Strassen:**

Uses approximately 7.75n^2 memory due to multiple intermediate submatrices and temporaries created at each recursion level.

Performs many allocations:

8 submatrices for splitting A and B,

7 temporary matrices M1…M7,

4 result blocks C11…C22,

plus recursive allocations at deeper levels.

<img width="697" height="559" alt="Screenshot 2025-12-23 at 19 26 59" src="https://github.com/user-attachments/assets/2b8988b6-1ac4-4383-a50b-db89ab26260d" />


### Impact on performance

The extra memory usage in Strassen leads to:

- Increased data movement: many intermediate matrices must be written to and read from memory.

- Reduced cache locality: data is accessed in multiple passes and through temporary buffers rather than streamed linearly.

- Higher allocation and copying overhead: especially costly inside recursive calls.

In contrast, naive and Winograd algorithms have simpler memory access patterns and much smaller working sets, allowing better cache utilization. As a result, a larger fraction of their operations is served from cache rather than main memory, which reduces memory latency and leads to consistently faster runtime compared to Strassen in the tested matrix sizes. 

### Observations on large matrices

For larger matrix sizes, runtime increased as expected due to the O(n^3) nature of all three algorithms. However, no abrupt slowdowns or non-linear jumps were observed in the timing plots: the growth remained smooth and approximately linear with respect to the expected asymptotic trend.


## Numerical stability

I compared Winograd and Strassen against the naive result as a reference.

For int and double cases (random/symmetric/identity), the observed error was essentially ~0 (within numerical noise / exactness where applicable).

For complex matrices, I observed non-zero errors that grow with n, and Strassen’s error is larger than Winograd’s

<img width="849" height="546" alt="Screenshot 2025-12-23 at 19 27 34" src="https://github.com/user-attachments/assets/d4a91125-6329-4cd4-aa92-f9ddb94c0a36" />


### Why Strassen is less stable (theory matching the experiment)

Strassen introduces more subtraction-heavy intermediate expressions like (B12−B22), (B21−B11), and combinations such as
C11=M1+M4-M5+M7, C22=M1−M2+M3+M6.

These intermediates can involve subtraction of similar-sized numbers, which causes catastrophic cancellation: significant digits are lost, and rounding error becomes relatively larger.

Those errors then propagate through recursion: at each level, you create more intermediate sums/differences, so the rounding noise can accumulate and amplify.

### Why complex numbers showed the issue more clearly

Complex multiplication and addition internally perform more floating-point operations (multiple real multiplications/additions per complex multiply), which increases the opportunities for rounding error.

## Hardware and Parallelism Impact

All custom implementations (naive, Winograd, and Strassen) were executed on the CPU without explicit parallelization or low-level hardware-specific optimizations

### Comparison with BLAS

When compared against an optimized BLAS implementation (CPU/GPU library matrix multiplication), the performance gap was dramatic:

BLAS was approximately 100× faster than both naive and Winograd.

BLAS was approximately 1000× faster than Strassen.

This gap highlights the dominant role of hardware-aware optimizations over pure algorithmic operation counts.

### Parallelism considerations

Naive multiplication is inherently well-suited for parallelism: each output element (or block) can be computed independently. This is why hardware vendors invest heavily in optimizing this kernel.

Winograd, while still structurally similar to naive multiplication, introduces extra precomputation steps but remains amenable to parallel execution.

Strassen, however, relies on recursive decomposition and complex intermediate combinations, which complicates efficient parallel scheduling and load balancing, especially on GPUs.

## Overall Takeaways

1. **Fewer multiplications do not automatically mean faster execution.**
Although Winograd and Strassen reduce the number of multiplications compared to the naive algorithm, only Winograd consistently translated this reduction into better runtime. Strassen, despite its lower asymptotic multiplication count, suffered from large constant factors due to extra additions, recursion, temporary matrices, and less regular memory access.

2. **Algorithmic simplicity matters in practice.**
The naive algorithm, while asymptotically inferior, benefits from a simple structure that aligns well with modern hardware. Its regular loops, predictable memory access, and low overhead make it surprisingly competitive, especially when compared to more complex recursive algorithms like Strassen.

3. **Practical performance is dominated by overhead and hardware effects.**
Memory usage, cache behavior, numerical stability, and constant factors dominated performance for all tested sizes. These effects prevented Strassen from becoming competitive within the tested range, despite its better theoretical complexity.

4. **Asymptotic gains may matter only at very large scales.**
Theoretical improvements such as lowering the exponent in the recursive complexity may become relevant for extremely large matrices, but only if implementations can efficiently control memory usage, numerical error, and data movement. Without such careful engineering, practical factors will continue to limit the usefulness of more complex algorithms.


Overall, this exercise provided a concrete understanding of why matrix multiplication is so hard to improve in practice. It showed that progress in this area is not just about reducing arithmetic counts, but about balancing algebraic cleverness with hardware realities, memory behavior, and numerical robustness.
