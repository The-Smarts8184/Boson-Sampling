// Boson Sampling Simulation in Java

import java.util.*;

public class BosonSampling {

    // Generate a random unitary matrix using Gram-Schmidt process
    public static double[][] generateUnitaryMatrix(int n) {
        Random rand = new Random();
        double[][] matrix = new double[n][n];

        // Generate random complex vectors
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = rand.nextGaussian();
            }
        }

        // Perform Gram-Schmidt orthogonalization
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                double dot = 0;
                for (int k = 0; k < n; k++) {
                    dot += matrix[k][i] * matrix[k][j];
                }
                for (int k = 0; k < n; k++) {
                    matrix[k][i] -= dot * matrix[k][j];
                }
            }
            double norm = 0;
            for (int k = 0; k < n; k++) {
                norm += matrix[k][i] * matrix[k][i];
            }
            norm = Math.sqrt(norm);
            for (int k = 0; k < n; k++) {
                matrix[k][i] /= norm;
            }
        }

        return matrix;
    }


    // Generate initial photon states (Fock states)
    public static int[] generatePhotonState(int n, int photons) {
        int[] state = new int[n];
        for (int i = 0; i < photons; i++) {
            state[i % n]++;
        }
        return state;
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
        int n = 25; // Number of modes
        int photons = 25; // Number of photons

        double[][] unitaryMatrix = generateUnitaryMatrix(n);
        int[] photonState = generatePhotonState(n, photons);

        System.out.println("Initial Photon State: " + Arrays.toString(photonState));
        System.out.println("Unitary Matrix:");
        for (double[] row : unitaryMatrix) System.out.println(Arrays.toString(row));

        double permanent = computePermanent(unitaryMatrix);
        System.out.println("Permanent of the matrix: " + permanent);
    }
}
