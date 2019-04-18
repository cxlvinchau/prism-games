package prism;

import java.io.*;
import java.util.*;
import java.math.*;
import parser.ast.*;

public class Bisimulation {

	private static final double zero = 0.000001;

	public static void main(String[] args) {
		new Bisimulation().go(args);
	}

	public void go(String args[]) {
		try {
			// Initialize
			PrismLog mainLog = new PrismFileLog("stdout");
			Prism prism = new Prism(mainLog);
			prism.initialise();
			prism.setEngine(Prism.EXPLICIT);

			// Parse model
			ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			prism.loadPRISMModel(modulesFile);

			// Get model
			prism.buildModel();
			explicit.Model model = prism.getBuiltModelExplicit();
			explicit.MDPSparse mdp = (explicit.MDPSparse) model;

			// Compute action matrices
			int numStates = mdp.getNumStates();
			HashMap<String, double[][]> actionMatrices = new HashMap<String, double[][]>();

			for (int s = 0; s < mdp.getNumStates(); s++) {
				String action = mdp.getAction(s, 0).toString();
				for (Iterator<Map.Entry<Integer, Double>> iter = mdp.getTransitionsIterator(s, 0); iter.hasNext();) {
					Map.Entry<Integer, Double> entry = iter.next();
					if (actionMatrices.containsKey(action)) {
						actionMatrices.get(action)[s][entry.getKey()] = entry.getValue();
					} else {
						actionMatrices.put(action, new double[numStates][numStates]);
						actionMatrices.get(action)[s][entry.getKey()] = entry.getValue();
					}
				}
			}

			List<double[][]> matrices = new LinkedList<double[][]>();
			for (Map.Entry<String, double[][]> entry : actionMatrices.entrySet()) {
				matrices.add(entry.getValue());
				System.out.println("Matrix for " + entry.getKey());
				printMatrix(entry.getValue());
				System.out.println("");
			}
			
			performIteration(matrices);
			
			
			
			List<double[]> vectorCollection = new LinkedList<double[]>();
			double[] v1 = {1,1,1,1,1,1,1};
			double[] v2 = {0,0,0,1,0,0,0};
			double[] v3 = {0,0,1,0,0,0,0};
			double[] v4 = {1,1,0,0,1,1,1};
			vectorCollection.add(v1);
			vectorCollection.add(v4);
			vectorCollection.add(v3);
			//vectorCollection.add(v2);
			double[] vector = {0,0.5,0,0,1,0,0};
			System.out.println(isLinearlyIndependent(vectorCollection, v2));
			
			

		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e);
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e);
			System.exit(1);
		}
	}

	
	public static void performIteration(List<double[][]> matrices) {
		if(matrices.isEmpty()) return;
		
		int numStates = matrices.get(0).length;
		List<double[]> vectorCollection = new LinkedList<double[]>();
		double[] initial = new double[numStates];
		for(int i = 0; i < numStates; i++) {
			initial[i] = 1.0;
		}
		vectorCollection.add(initial);
		for(int i = 0; i < numStates; i++) {
			for(int v = 0; v < vectorCollection.size(); v++) {
				for(double[][] matrix : matrices) {
					double[] newVector = multMatrixVector(matrix, vectorCollection.get(v));
					if(isLinearlyIndependent(vectorCollection, newVector)) {
						printVector(newVector);
						vectorCollection.add(newVector);
					}
				}
			}
		}
		
		
		System.out.println("vectors: ");
		for(double[] vector : vectorCollection) {
			printVector(vector);
		}
		
	}
	
	private static void printVector(double[] vector) {
		for(int i = 0; i < vector.length; i++) {
			System.out.print(vector[i] + " ");
		}
		System.out.println("\n");
	}
	
	private static double[] multMatrixVector(double[][] matrix, double[] vector) {
		double[] result = new double[vector.length];
		for(int r = 0; r < matrix.length; r++) {
			double tmp = 0.0;
			for(int c = 0; c < matrix[r].length; c++) {
				tmp += matrix[r][c] * vector[c];
			}
			result[r] = tmp;
		}
		return result;
	}
	
	/*
	 * Check if vector is linearly independent of vectorCollection
	 *
	 */
	private static boolean isLinearlyIndependent(List<double[]> vectorCollection, double[] vector) {
		double[][] matrix = new double[vectorCollection.size()+1][vector.length];
		for (int r = 0; r < vectorCollection.size() + 1; r++) {
			for (int c = 0; c < vector.length; c++) {
				if (r == vectorCollection.size()) {
					matrix[r][c] = vector[c];
				} else {
					matrix[r][c] = vectorCollection.get(r)[c];
				}
			}
		}
		
		gaussElimination(matrix);
		
		for (int r = 0; r < matrix.length; r++) {
			boolean zeroRow = true;
			for (int c = 0; c < matrix[0].length; c++) {
				if (!isZero(matrix[r][c])) {
					zeroRow = false;
					break;
				}
			}
			if (zeroRow) {
				return false;
			}
		}

		return true;
	}

	private static void gaussElimination(double[][] matrix) {
		// Perform Gaussian elimination with matrix
		int c = 0;
		for (int r = 0; r < matrix.length; r++) {
			boolean terminate = false;
			while (r != matrix.length - 1 && isZero(matrix[r][c])) {
				for (int i = r + 1; i < matrix.length; i++) {
					if (!isZero(matrix[i][c])) {
						swapRow(r, i, matrix);
						break;
					}
					if (i == matrix.length - 1) {
						c++;
						if (c >= matrix[r].length) {
							terminate = true;
							break;
						}
					}
				}
				if (terminate) {
					break;
				}
			}
			if (terminate) {
				break;
			}
			for (int r1 = r + 1; r1 < matrix.length; r1++) {
				addRow(multiplyRow(-matrix[r1][c] / matrix[r][c], matrix[r]), matrix, r1);
			}
			c++;
			if (c >= matrix[r].length) {
				break;
			}
		}
	}

	private static void printMatrix(double[][] matrix) {
		// Print result
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				System.out.print(matrix[r][c] + " ");
			}
			System.out.print("\n");
		}
	}

	private static boolean isZero(double val) {
		return Math.abs(val) < zero;
	}

	private static void addRow(double[] vector, double[][] matrix, int row) {
		for (int c = 0; c < vector.length; c++) {
			matrix[row][c] += vector[c];
		}
	}

	private static double[] multiplyRow(double scalar, double[] row) {
		double[] tmp = row.clone();
		for (int c = 0; c < tmp.length; c++) {
			tmp[c] *= scalar;
		}
		return tmp;
	}

	private static void swapRow(int r1, int r2, double[][] matrix) {
		for (int c = 0; c < matrix[r1].length; c++) {
			double tmp = matrix[r1][c];
			matrix[r1][c] = matrix[r2][c];
			matrix[r2][c] = tmp;
		}
	}

}
