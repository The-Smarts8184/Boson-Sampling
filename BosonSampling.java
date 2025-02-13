import java.util.*;

public class BosonSampling {

    // Generate a random unitary matrix (for simplicity, using identity matrix for now)
    public static double[][] generateUnitaryMatrix(int n) {
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return matrix;
    }

    // Compute the permanent of a matrix (Ryser's algorithm)
    public static double computePermanent(double[][] matrix) {
        int n = matrix.length;
        int[] indices = new int[n];
        for (int i = 0; i < n; i++) indices[i] = i;

        double permanent = 0.0;
        int sign = 1;

        do {
            double product = 1.0;
            for (int i = 0; i < n; i++) {
                product *= matrix[i][indices[i]];
            }
            permanent += sign * product;
            sign = -sign;
        } while (nextPermutation(indices));

        return Math.abs(permanent);
    }

    // Helper to generate next permutation (lexicographic order)
    private static boolean nextPermutation(int[] array) {
        int i = array.length - 2;
        while (i >= 0 && array[i] >= array[i + 1]) i--;
        if (i < 0) return false;
        int j = array.length - 1;
        while (array[j] <= array[i]) j--;
        swap(array, i, j);
        reverse(array, i + 1, array.length - 1);
        return true;
    }

    private static void swap(int[] array, int i, int j) {
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    private static void reverse(int[] array, int start, int end) {
        while (start < end) swap(array, start++, end--);
    }

    public static void main(String[] args) {
        int n = 75; // Number of modes/photons
        double[][] unitaryMatrix = generateUnitaryMatrix(n);

        System.out.println("Unitary Matrix:");
        for (double[] row : unitaryMatrix) System.out.println(Arrays.toString(row));

        double permanent = computePermanent(unitaryMatrix);
        System.out.println("Permanent of the matrix: " + permanent);

        // Further steps: Generate photon states, apply unitary transformation, simulate sampling
    }
}
