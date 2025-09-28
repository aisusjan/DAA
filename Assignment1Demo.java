/*
Single-file Java program you can copy-paste into IntelliJ IDEA (or any Java IDE).
It implements:
 - MergeSort (with reusable buffer + small-n cutoff)
 - QuickSort (randomized pivot, recurse on smaller partition)
 - Deterministic Select (median-of-medians / MoM5)
 - Closest Pair (2D divide-and-conquer)
 - Simple Metrics (comparisons, allocations, recursion depth)
 - Small correctness checks and timing in main()

Save as Assignment1Demo.java and run (Java 17+ recommended).
*/

import java.util.*;
import java.util.concurrent.atomic.AtomicLong;

public class Assignment1Demo {

    /* --------------------- Metrics --------------------- */
    public static class Metrics {
        public final AtomicLong comparisons = new AtomicLong();
        public final AtomicLong allocations = new AtomicLong();
        public final AtomicLong maxRecDepth = new AtomicLong();
        private final ThreadLocal<Long> depth = ThreadLocal.withInitial(() -> 0L);

        public void enter() {
            long d = depth.get() + 1;
            depth.set(d);
            maxRecDepth.getAndUpdate(x -> Math.max(x, d));
        }

        public void exit() {
            depth.set(depth.get() - 1);
        }

        public void addComparison(long c) { comparisons.addAndGet(c); }
        public void addAllocation(long a) { allocations.addAndGet(a); }

        public void reset() {
            comparisons.set(0);
            allocations.set(0);
            maxRecDepth.set(0);
            depth.set(0L);
        }

        @Override
        public String toString() {
            return String.format("comps=%d allocs=%d depth=%d",
                    comparisons.get(), allocations.get(), maxRecDepth.get());
        }
    }

    /* --------------------- Utilities --------------------- */
    public static final class Utils {
        private static final Random RNG = new Random();

        private Utils() {}

        public static void swap(int[] a, int i, int j) {
            int t = a[i]; a[i] = a[j]; a[j] = t;
        }

        public static void shuffle(int[] a) {
            for (int i = a.length - 1; i > 0; --i) {
                int j = RNG.nextInt(i + 1);
                swap(a, i, j);
            }
        }

        // Lomuto partition: pivotIndex is index of pivot; returns final pivot index
        public static int partition(int[] a, int lo, int hi, int pivotIndex) {
            int pivot = a[pivotIndex];
            swap(a, pivotIndex, hi);
            int store = lo;
            for (int i = lo; i < hi; ++i) {
                if (a[i] < pivot) {
                    swap(a, i, store);
                    store++;
                }
            }
            swap(a, store, hi);
            return store;
        }
    }

    /* --------------------- MergeSort --------------------- */
    public static final class MergeSort {
        private MergeSort() {}

        public static void sort(int[] a, Metrics m) {
            if (a == null || a.length <= 1) return;
            int[] buf = new int[a.length];
            sort(a, buf, 0, a.length - 1, m);
        }

        private static void sort(int[] a, int[] buf, int lo, int hi, Metrics m) {
            if (lo >= hi) return;
            int n = hi - lo + 1;
            if (n <= 16) { // cutoff to insertion
                insertion(a, lo, hi, m);
                return;
            }
            int mid = lo + ((hi - lo) >> 1);
            m.enter();
            try {
                sort(a, buf, lo, mid, m);
                sort(a, buf, mid + 1, hi, m);
            } finally { m.exit(); }
            merge(a, buf, lo, mid, hi, m);
        }

        private static void merge(int[] a, int[] buf, int lo, int mid, int hi, Metrics m) {
            System.arraycopy(a, lo, buf, lo, hi - lo + 1);
            int i = lo, j = mid + 1, k = lo;
            while (i <= mid && j <= hi) {
                m.addComparison(1);
                if (buf[i] <= buf[j]) a[k++] = buf[i++];
                else a[k++] = buf[j++];
                m.addAllocation(1);
            }
            while (i <= mid) { a[k++] = buf[i++]; m.addAllocation(1); }
            while (j <= hi) { a[k++] = buf[j++]; m.addAllocation(1); }
        }

        private static void insertion(int[] a, int lo, int hi, Metrics m) {
            for (int i = lo + 1; i <= hi; ++i) {
                int key = a[i];
                int j = i - 1;
                while (j >= lo) {
                    m.addComparison(1);
                    if (a[j] > key) { a[j + 1] = a[j]; m.addAllocation(1); j--; }
                    else break;
                }
                a[j + 1] = key; m.addAllocation(1);
            }
        }
    }

    /* --------------------- QuickSort --------------------- */
    public static final class QuickSort {
        private static final Random RNG = new Random();

        private QuickSort() {}

        public static void sort(int[] a, Metrics m) {
            if (a == null || a.length <= 1) return;
            Utils.shuffle(a);
            sort(a, 0, a.length - 1, m);
        }

        // tail recursion elimination: loop and recurse on smaller side
        private static void sort(int[] a, int lo, int hi, Metrics m) {
            while (lo < hi) {
                int pivot = lo + RNG.nextInt(hi - lo + 1);
                int p = Utils.partition(a, lo, hi, pivot);
                int leftSize = p - lo;
                int rightSize = hi - p;
                m.addComparison(1);
                m.enter();
                try {
                    if (leftSize < rightSize) {
                        sort(a, lo, p - 1, m); // recurse smaller
                        lo = p + 1; // iterate on larger
                    } else {
                        sort(a, p + 1, hi, m);
                        hi = p - 1;
                    }
                } finally { m.exit(); }
            }
        }
    }

    /* --------------------- Deterministic Select (Median-of-Medians) --------------------- */
    public static final class DeterministicSelect {
        private DeterministicSelect() {}

        // Public convenience: select k-th smallest (0-based)
        public static int select(int[] a, int k, Metrics m) {
            if (a == null || k < 0 || k >= a.length) throw new IllegalArgumentException();
            return select(a, 0, a.length - 1, k, m);
        }

        // returns value of k-th smallest in a[lo..hi]
        private static int select(int[] a, int lo, int hi, int k, Metrics m) {
            while (true) {
                if (lo == hi) return a[lo];
                m.enter();
                try {
                    int pivotIndex = pivotMedianOfMedians(a, lo, hi, m);
                    int p = Utils.partition(a, lo, hi, pivotIndex);
                    if (k == p) return a[p];
                    else if (k < p) hi = p - 1;
                    else lo = p + 1;
                } finally { m.exit(); }
            }
        }

        // choose pivot index using median-of-medians (groups of 5)
        private static int pivotMedianOfMedians(int[] a, int lo, int hi, Metrics m) {
            int n = hi - lo + 1;
            if (n <= 5) {
                insertionSortRange(a, lo, hi, m);
                return lo + n / 2;
            }
            int numGroups = (n + 4) / 5;
            // move each group's median to the front region [lo .. lo+numGroups-1]
            for (int i = 0; i < numGroups; ++i) {
                int gLo = lo + i * 5;
                int gHi = Math.min(gLo + 4, hi);
                insertionSortRange(a, gLo, gHi, m);
                int medianIndex = gLo + (gHi - gLo) / 2;
                Utils.swap(a, lo + i, medianIndex);
                m.addAllocation(1);
            }
            // find median of the medians by selecting on the front region
            int medianOfMediansValue = selectValue(a, lo, lo + numGroups - 1, lo + numGroups / 2, m);
            // find an index of this value between lo..hi (safe because value is from that range)
            for (int i = lo; i <= hi; ++i) if (a[i] == medianOfMediansValue) return i;
            // fallback (shouldn't happen): return lo
            return lo;
        }

        // helper that returns value (k-th smallest) within given bounds; uses select loop
        private static int selectValue(int[] a, int lo, int hi, int kIndex, Metrics m) {
            // kIndex is absolute index in array we want (the earlier code uses lo + something),
            // but we want a rank relative to lo..hi -> translate:
            int relativeK = kIndex - lo;
            if (relativeK < 0 || relativeK > (hi - lo)) relativeK = (hi - lo) / 2; // safe fallback
            // Copy the subarray (small) to a temp array to avoid disturbing outer array structure during recursion
            // (This keeps things simpler and safe.)
            int len = hi - lo + 1;
            int[] tmp = new int[len];
            System.arraycopy(a, lo, tmp, 0, len);
            // Now perform deterministic select on tmp for relativeK
            Metrics tempM = new Metrics(); // local metrics to avoid polluting outer counters
            int val = selectInTemp(tmp, 0, len - 1, relativeK, tempM);
            // we can choose to add some of tempM counts to outer m to reflect work
            m.addComparison(tempM.comparisons.get());
            m.addAllocation(tempM.allocations.get());
            return val;
        }

        // internal select used for the small temp array; it's a straightforward recursive implementation
        private static int selectInTemp(int[] a, int lo, int hi, int k, Metrics m) {
            while (true) {
                if (lo == hi) return a[lo];
                int n = hi - lo + 1;
                if (n <= 5) {
                    insertionSortRange(a, lo, hi, m);
                    return a[lo + k];
                }
                int numGroups = (n + 4) / 5;
                for (int i = 0; i < numGroups; ++i) {
                    int gLo = lo + i * 5;
                    int gHi = Math.min(gLo + 4, hi);
                    insertionSortRange(a, gLo, gHi, m);
                    int medianIndex = gLo + (gHi - gLo) / 2;
                    // move median to front
                    int dest = lo + i;
                    int t = a[medianIndex]; a[medianIndex] = a[dest]; a[dest] = t;
                    m.addAllocation(1);
                }
                int mid = lo + (numGroups - 1) / 2;
                int pivotVal = selectInTemp(a, lo, lo + numGroups - 1, mid - lo, m); // median of medians value
                // find pivot index
                int pivotIndex = lo;
                for (int i = lo; i <= hi; ++i) if (a[i] == pivotVal) { pivotIndex = i; break; }
                int p = partitionInTemp(a, lo, hi, pivotIndex);
                int rank = p - lo;
                if (k == rank) return a[p];
                else if (k < rank) hi = p - 1;
                else { k = k - rank - 1; lo = p + 1; }
            }
        }

        private static int partitionInTemp(int[] a, int lo, int hi, int pivotIndex) {
            int pivot = a[pivotIndex];
            int t = a[pivotIndex]; a[pivotIndex] = a[hi]; a[hi] = t;
            int store = lo;
            for (int i = lo; i < hi; ++i) {
                if (a[i] < pivot) {
                    int tt = a[i]; a[i] = a[store]; a[store] = tt;
                    store++;
                }
            }
            int tt = a[store]; a[store] = a[hi]; a[hi] = tt;
            return store;
        }

        private static void insertionSortRange(int[] a, int lo, int hi, Metrics m) {
            for (int i = lo + 1; i <= hi; ++i) {
                int key = a[i];
                int j = i - 1;
                while (j >= lo) {
                    m.addComparison(1);
                    if (a[j] > key) { a[j + 1] = a[j]; m.addAllocation(1); j--; }
                    else break;
                }
                a[j + 1] = key; m.addAllocation(1);
            }
        }
    }

    /* --------------------- Closest Pair (2D) --------------------- */
    public static final class ClosestPair {
        private ClosestPair() {}

        // pts: double[n][2] where each row is {x,y}
        public static double closestPair(double[][] pts, Metrics m) {
            if (pts == null || pts.length < 2) return Double.POSITIVE_INFINITY;
            double[][] byX = pts.clone();
            Arrays.sort(byX, Comparator.comparingDouble(p -> p[0]));
            double[][] aux = new double[pts.length][2];
            return rec(byX, aux, 0, pts.length - 1, m);
        }

        private static double dist(double[] a, double[] b) {
            double dx = a[0] - b[0];
            double dy = a[1] - b[1];
            return Math.hypot(dx, dy);
        }

        private static double brute(double[][] a, int lo, int hi) {
            double best = Double.POSITIVE_INFINITY;
            for (int i = lo; i <= hi; ++i)
                for (int j = i + 1; j <= hi; ++j)
                    best = Math.min(best, dist(a[i], a[j]));
            return best;
        }

        private static double rec(double[][] byX, double[][] aux, int lo, int hi, Metrics m) {
            int n = hi - lo + 1;
            if (n <= 3) return brute(byX, lo, hi);
            int mid = lo + (hi - lo) / 2;
            double midx = byX[mid][0];
            m.enter();
            double d1;
            try { d1 = rec(byX, aux, lo, mid, m); }
            finally { m.exit(); }
            m.enter();
            double d2;
            try { d2 = rec(byX, aux, mid + 1, hi, m); }
            finally { m.exit(); }
            double d = Math.min(d1, d2);
            // collect points within strip
            int sc = 0;
            for (int i = lo; i <= hi; ++i) {
                if (Math.abs(byX[i][0] - midx) < d) aux[sc++] = byX[i];
            }
            Arrays.sort(aux, 0, sc, Comparator.comparingDouble(p -> p[1]));
            for (int i = 0; i < sc; ++i) {
                for (int j = i + 1; j < sc && (aux[j][1] - aux[i][1]) < d; ++j) {
                    d = Math.min(d, dist(aux[i], aux[j]));
                }
            }
            return d;
        }
    }

    /* --------------------- Main / Demo --------------------- */
    public static void main(String[] args) {
        Metrics m = new Metrics();

        // Small correctness test for sorting algorithms
        int n = 20000; // modest size so this runs quickly
        Random rnd = new Random(12345);
        int[] base = new int[n];
        for (int i = 0; i < n; ++i) base[i] = rnd.nextInt();

        // MergeSort
        int[] a1 = base.clone();
        m.reset();
        long t0 = System.nanoTime();
        MergeSort.sort(a1, m);
        long t1 = System.nanoTime();
        boolean okMerge = Arrays.equals(a1, sortedCopy(base));
        System.out.printf("MergeSort: n=%d time_ms=%.3f %s\n", n, (t1 - t0) / 1e6, okMerge ? "OK" : "FAIL");
        System.out.println("  " + m);

        // QuickSort
        int[] a2 = base.clone();
        m.reset();
        t0 = System.nanoTime();
        QuickSort.sort(a2, m);
        t1 = System.nanoTime();
        boolean okQuick = Arrays.equals(a2, sortedCopy(base));
        System.out.printf("QuickSort: n=%d time_ms=%.3f %s\n", n, (t1 - t0) / 1e6, okQuick ? "OK" : "FAIL");
        System.out.println("  " + m);

        // Deterministic Select (median)
        int[] a3 = base.clone();
        int k = n / 2;
        m.reset();
        t0 = System.nanoTime();
        int med = DeterministicSelect.select(a3, k, m);
        t1 = System.nanoTime();
        int[] sorted = sortedCopy(base);
        boolean okSelect = med == sorted[k];
        System.out.printf("DeterministicSelect: k=%d time_ms=%.3f %s (value=%d)\n", k, (t1 - t0) / 1e6, okSelect ? "OK" : "FAIL", med);
        System.out.println("  " + m);

        // Closest pair (random points)
        int pN = 20000;
        double[][] pts = new double[pN][2];
        Random rr = new Random(42);
        for (int i = 0; i < pN; ++i) { pts[i][0] = rr.nextDouble(); pts[i][1] = rr.nextDouble(); }
        m.reset();
        t0 = System.nanoTime();
        double d = ClosestPair.closestPair(pts, m);
        t1 = System.nanoTime();
        System.out.printf("ClosestPair: n=%d time_ms=%.3f distance=%.6g\n", pN, (t1 - t0) / 1e6, d);
        System.out.println("  " + m);

        // Extra small validation: ClosestPair vs brute for small n
        int small = 300;
        double[][] smallPts = new double[small][2];
        Random r2 = new Random(7);
        for (int i = 0; i < small; ++i) { smallPts[i][0] = r2.nextDouble(); smallPts[i][1] = r2.nextDouble(); }
        double brute = bruteForceClosest(smallPts);
        double fast = ClosestPair.closestPair(smallPts, m);
        System.out.printf("Closest pair small n=%d: brute=%.6g fast=%.6g (diff=%.6g)\n", small, brute, fast, Math.abs(brute - fast));
    }

    private static int[] sortedCopy(int[] a) {
        int[] b = a.clone();
        Arrays.sort(b);
        return b;
    }

    private static double bruteForceClosest(double[][] pts) {
        double best = Double.POSITIVE_INFINITY;
        int n = pts.length;
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
            double dx = pts[i][0] - pts[j][0];
            double dy = pts[i][1] - pts[j][1];
            best = Math.min(best, Math.hypot(dx, dy));
        }
        return best;
    }
}
