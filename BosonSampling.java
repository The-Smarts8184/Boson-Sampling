// Boson Sampling Simulation in Java
import java.util.*;
import java.util.concurrent.*;

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

    // Apply the unitary transformation to the photon state
    public static double[] applyUnitary(double[][] unitary, int[] photonState) {
        int n = unitary.length;
        double[] outputState = new double[n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                outputState[i] += unitary[i][j] * photonState[j];
            }
        }

        return outputState;
    }


    // Compute the permanent of a matrix (Ryser's algorithm) using multiple threads
    // Compute the permanent using multi-threading
    public static double computePermanent(double[][] matrix) throws InterruptedException, ExecutionException {
        int n = matrix.length;
        int[] indices = new int[n];
        for (int i = 0; i < n; i++) indices[i] = i;

        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future<Double>> futures = new ArrayList<>();

        do {
            final int[] perm = Arrays.copyOf(indices, n);
            futures.add(executor.submit(() -> {
                double product = 1.0;
                for (int i = 0; i < n; i++) {
                    product *= matrix[i][perm[i]];
                }
                return product;
            }));
        } while (nextPermutation(indices));

        double permanent = 0.0;
        for (Future<Double> f : futures) permanent += f.get();
        executor.shutdown();

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

    public static void main(String[] args) throws InterruptedException, ExecutionException {
        System.out.println("Note: Run the program with increased heap size, e.g., java -Xmx4G BosonSampling");
        int n = 3, photons = 2;
        double[][] unitaryMatrix = generateUnitaryMatrix(n);
        int[] photonState = generatePhotonState(n, photons);

        System.out.println("Initial Photon State: " + Arrays.toString(photonState));
        System.out.println("Unitary Matrix:");
        for (double[] row : unitaryMatrix) System.out.println(Arrays.toString(row));

        double[] outputState = applyUnitary(unitaryMatrix, photonState);
        System.out.println("Output Photon State After Unitary Transformation: " + Arrays.toString(outputState));

        double permanent = computePermanent(unitaryMatrix);
        System.out.println("Permanent of the matrix: " + permanent);
    }
}
