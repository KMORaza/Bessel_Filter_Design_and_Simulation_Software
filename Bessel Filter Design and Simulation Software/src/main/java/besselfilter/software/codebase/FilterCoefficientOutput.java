package besselfilter.software.codebase;

import java.util.Arrays;

public class FilterCoefficientOutput {
    public String getCoefficientOutput(String filterType, int order, double[] denominator, BesselFilterModel model) {
        StringBuilder output = new StringBuilder();

        // Transfer Function
        double[] numerator = model.getNumeratorCoefficients(filterType, order);
        output.append("Transfer Function Coefficients:\n");
        output.append("Numerator: ").append(Arrays.toString(numerator)).append("\n");
        output.append("Denominator: ").append(Arrays.toString(denominator)).append("\n");

        // State-Space Representation
        double[][][] stateSpace = model.getStateSpace(filterType, order, denominator);
        output.append("\nState-Space Representation:\n");
        output.append("A Matrix:\n").append(matrixToString(stateSpace[0])).append("\n");
        output.append("B Matrix:\n").append(matrixToString(stateSpace[1])).append("\n");
        output.append("C Matrix:\n").append(matrixToString(stateSpace[2])).append("\n");
        output.append("D Matrix:\n").append(matrixToString(stateSpace[3])).append("\n");

        // Second-Order Sections
        double[][] sos = model.getSecondOrderSections(denominator, filterType, order);
        output.append("\nSecond-Order Sections (SOS):\n");
        for (int i = 0; i < sos.length; i++) {
            output.append("Section ").append(i + 1).append(": ").append(Arrays.toString(sos[i])).append("\n");
        }

        return output.toString();
    }

    private String matrixToString(double[][] matrix) {
        StringBuilder sb = new StringBuilder();
        for (double[] row : matrix) {
            sb.append(Arrays.toString(row)).append("\n");
        }
        return sb.toString();
    }
}