# DAA
Aisultan Nuriman 
SE 2430
Introduction
The program Assignment1Demo.java demonstrates the implementation and analysis of several classic algorithms from the field of Design and Analysis of Algorithms. It includes sorting algorithms (MergeSort and QuickSort), a deterministic selection algorithm (Median-of-Medians), and a geometric algorithm (Closest Pair of Points in 2D). Additionally, it integrates a custom Metrics system to measure performance aspects such as comparisons, allocations, and recursion depth.
The goal of this assignment is not only to implement these algorithms correctly, but also to evaluate and compare their efficiency using experimental metrics.
________________________________________
Implemented Components
1. Metrics Class
•	Tracks comparisons, allocations, and recursion depth using AtomicLong.
•	Provides methods enter() and exit() to update recursion depth dynamically.
•	Designed to be reusable across all algorithms.
•	Outputs a summary string (e.g., comps=..., allocs=..., depth=...).
2. Utilities
•	Contains helper methods:
o	swap(int[], i, j) for exchanging elements.
o	shuffle(int[]) for randomizing arrays.
o	partition(...) implementing Lomuto’s partition scheme.
•	Used in QuickSort and Deterministic Select.
3. MergeSort
•	Implements top-down MergeSort with:
o	Insertion sort cutoff for small subarrays (≤16 elements).
o	Buffer array for efficient merging.
•	Measures comparisons and allocations during merging and insertion.
•	Recursive calls are tracked with the Metrics system.
4. QuickSort
•	Implements randomized QuickSort:
o	Chooses a pivot randomly.
o	Uses tail recursion elimination by recursing only on the smaller partition and iterating on the larger.
o	Ensures recursion depth is limited.
•	Shuffles the array initially to avoid worst-case scenarios.
5. Deterministic Select (Median-of-Medians)
•	Finds the k-th smallest element in an array.
•	Uses Median-of-Medians (MoM5) strategy for pivot selection.
•	Guarantees linear time complexity in the worst case.
•	Includes:
o	pivotMedianOfMedians() to find a robust pivot.
o	Recursive partitioning until the target index is found.
•	Allocations and comparisons are tracked through Metrics.
6. Closest Pair of Points (2D)
•	Implements a divide-and-conquer solution:
o	Points sorted by x-coordinates.
o	Recursively computes minimum distances in left and right halves.
o	Merges results by checking a vertical strip around the median.
o	Maintains efficiency with auxiliary arrays and y-sorting.
•	Runs in O(n log n) time compared to brute-force O(n²).
________________________________________
Testing and Validation
The main() method runs multiple experiments:
1.	MergeSort
o	Sorts an array of 20,000 random integers.
o	Validates correctness by comparing with Arrays.sort.
o	Reports execution time and metrics.
2.	QuickSort
o	Similar validation as MergeSort.
o	Ensures correct ordering with randomized pivot selection.
3.	Deterministic Select
o	Finds the median element (k = n/2).
o	Verifies correctness by comparing with a fully sorted copy.
o	Reports execution time and the selected value.
4.	Closest Pair
o	Generates 20,000 random 2D points.
o	Computes the closest pair distance using divide-and-conquer.
o	Compares with brute-force results for small input sizes (n=300) to ensure correctness.
Conclusion
This Java program successfully demonstrates and validates multiple classical algorithms:
•	MergeSort shows stable, predictable performance with low recursion depth.
•	QuickSort achieves practical efficiency with randomized pivots and tail recursion optimization.
•	Deterministic Select illustrates worst-case linear-time selection with the Median-of-Medians method.
•	Closest Pair of Points highlights an elegant divide-and-conquer approach in computational geometry.
The inclusion of a custom Metrics system enhances the understanding of algorithm performance beyond runtime, allowing insight into internal operations such as comparisons, allocations, and recursion depth. This makes the assignment a comprehensive study in both correctness and efficiency of fundamental algorithms.

