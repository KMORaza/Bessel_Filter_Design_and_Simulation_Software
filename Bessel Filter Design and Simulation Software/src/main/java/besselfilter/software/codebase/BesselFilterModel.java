package besselfilter.software.codebase;

import java.util.ArrayList;
import java.util.Arrays;

public class BesselFilterModel {
    private StringBuilder designDetails;
    private ArrayList<Complex> poles;
    private ArrayList<Complex> zeros;

    public BesselFilterModel() {
        designDetails = new StringBuilder();
        poles = new ArrayList<>();
        zeros = new ArrayList<>();
    }

    public double[] designBesselFilter(String filterType, int order, double cutoffFreq, double centerFreq, boolean isDigital, double sampleRate) {
        designDetails.setLength(0);
        designDetails.append(String.format("Designing %s Bessel Filter\nOrder: %d\nCutoff Frequency: %.2f Hz\n", filterType, order, cutoffFreq));
        if (filterType.contains("Band")) {
            designDetails.append(String.format("Center Frequency: %.2f Hz\n", centerFreq));
        }
        if (isDigital) {
            designDetails.append(String.format("Digital IIR, Sampling Rate: %.2f Hz\n", sampleRate));
        }

        poles.clear();
        zeros.clear();
        double[] analogDenominator = designAnalogLowPassPrototype(order);
        designDetails.append("Analog Low-pass Prototype Denominator: ").append(Arrays.toString(analogDenominator)).append("\n");

        normalizeToUnityCutoff(analogDenominator, 1.0);

        double[] transformedDenominator = switch (filterType) {
            case "High-pass" -> lowPassToHighPass(analogDenominator, cutoffFreq);
            case "Band-pass" -> lowPassToBandPass(analogDenominator, cutoffFreq, centerFreq);
            case "Band-stop" -> lowPassToBandStop(analogDenominator, cutoffFreq, centerFreq);
            default -> analogDenominator;
        };
        designDetails.append("Transformed Denominator: ").append(Arrays.toString(transformedDenominator)).append("\n");

        double[] finalDenominator = isDigital ? bilinearTransform(transformedDenominator, sampleRate) : transformedDenominator;
        if (isDigital) {
            designDetails.append("Digital IIR Denominator: ").append(Arrays.toString(finalDenominator)).append("\n");
        }

        computeZeros(filterType, order, cutoffFreq, centerFreq, isDigital);
        return finalDenominator;
    }

    private double[] designAnalogLowPassPrototype(int order) {
        if (order < 1) {
            throw new IllegalArgumentException("Order must be at least 1.");
        }

        double[] coeffs = computeBesselPolynomial(order);
        designDetails.append("Bessel Polynomial Coefficients: ").append(Arrays.toString(coeffs)).append("\n");

        poles = findPolynomialRoots(coeffs, order);
        designDetails.append("Poles: ").append(poles.toString()).append("\n");

        return polesToDenominator(poles);
    }

    private double[] computeBesselPolynomial(int n) {
        double[] coeffs = new double[n + 1];
        coeffs[n] = 1;
        if (n == 0) return coeffs;

        double[] prev = new double[] {1};
        double[] curr = new double[] {1, 1};

        for (int k = 2; k <= n; k++) {
            double[] next = new double[k + 1];
            next[k] = 1;
            for (int i = 0; i < k; i++) {
                next[i] = curr[i] * (2 * k - 1) + (i == 0 ? 0 : prev[i - 1]);
            }
            prev = curr;
            curr = next;
        }
        return curr;
    }

    private ArrayList<Complex> findPolynomialRoots(double[] coeffs, int order) {
        ArrayList<Complex> roots = new ArrayList<>();
        Complex[][] knownPoles = {
                {new Complex(-1, 0)},
                {new Complex(-1.622, -0.837), new Complex(-1.622, 0.837)},
                {new Complex(-2.322, 0), new Complex(-1.419, -1.591), new Complex(-1.419, 1.591)},
                {new Complex(-2.896, -0.867), new Complex(-2.896, 0.867), new Complex(-1.242, -2.242), new Complex(-1.242, 2.242)}
        };
        if (order <= 4) {
            roots.addAll(Arrays.asList(knownPoles[order - 1]));
        } else {
            throw new UnsupportedOperationException("Orders > 4 require numerical root finding.");
        }
        return roots;
    }

    private double[] polesToDenominator(ArrayList<Complex> poles) {
        double[] denom = new double[] {1};
        for (Complex pole : poles) {
            double[] poly = pole.isReal() ?
                    new double[] {-pole.getReal(), 1} :
                    new double[] {pole.getReal() * pole.getReal() + pole.getImaginary() * pole.getImaginary(), -2 * pole.getReal(), 1};
            denom = convolve(denom, poly);
        }
        return denom;
    }

    private void normalizeToUnityCutoff(double[] coeffs, double cutoffFreq) {
        double factor = Math.pow(cutoffFreq, coeffs.length - 1);
        for (int i = 0; i < coeffs.length; i++) {
            coeffs[i] /= factor;
        }
    }

    private double[] lowPassToHighPass(double[] coeffs, double cutoffFreq) {
        ArrayList<Complex> hpPoles = new ArrayList<>();
        Complex cutoffSquared = new Complex(cutoffFreq * cutoffFreq, 0);
        for (Complex p : poles) {
            if (p.getReal() != 0 || p.getImaginary() != 0) {
                hpPoles.add(cutoffSquared.divide(p));
            }
        }
        poles = hpPoles;
        return polesToDenominator(poles);
    }

    private double[] lowPassToBandPass(double[] coeffs, double cutoffFreq, double centerFreq) {
        double bandwidth = cutoffFreq;
        double omega0 = 2 * Math.PI * centerFreq;
        Complex omega0Squared = new Complex(omega0 * omega0, 0);
        ArrayList<Complex> bpPoles = new ArrayList<>();
        for (Complex p : poles) {
            Complex term1 = p.multiply(new Complex(bandwidth, 0));
            Complex term2 = term1.multiply(term1).add(omega0Squared);
            Complex sqrtTerm = term2.sqrt();
            Complex two = new Complex(2, 0);
            bpPoles.add(term1.add(sqrtTerm).divide(two));
            bpPoles.add(term1.subtract(sqrtTerm).divide(two));
        }
        poles = bpPoles;
        return polesToDenominator(poles);
    }

    private double[] lowPassToBandStop(double[] coeffs, double cutoffFreq, double centerFreq) {
        double bandwidth = cutoffFreq;
        double omega0 = 2 * Math.PI * centerFreq;
        Complex omega0Squared = new Complex(omega0 * omega0, 0);
        ArrayList<Complex> bsPoles = new ArrayList<>();
        for (Complex p : poles) {
            Complex term1 = omega0Squared.divide(p.multiply(new Complex(bandwidth, 0)));
            Complex term2 = term1.multiply(term1).add(omega0Squared);
            Complex sqrtTerm = term2.sqrt();
            Complex two = new Complex(2, 0);
            bsPoles.add(term1.add(sqrtTerm).divide(two));
            bsPoles.add(term1.subtract(sqrtTerm).divide(two));
        }
        poles = bsPoles;
        return polesToDenominator(poles);
    }

    private double[] bilinearTransform(double[] analogDenom, double sampleRate) {
        double T = 1.0 / sampleRate;
        Complex c = new Complex(2.0 / T, 0);
        ArrayList<Complex> digitalPoles = new ArrayList<>();
        ArrayList<Complex> digitalZeros = new ArrayList<>();
        for (Complex p : poles) {
            Complex term = p.multiply(new Complex(T / 2.0, 0));
            Complex one = new Complex(1, 0);
            digitalPoles.add(one.add(term).divide(one.subtract(term)));
        }
        for (Complex z : zeros) {
            Complex term = z.multiply(new Complex(T / 2.0, 0));
            Complex one = new Complex(1, 0);
            digitalZeros.add(one.add(term).divide(one.subtract(term)));
        }
        poles = digitalPoles;
        zeros = digitalZeros;
        return polesToDenominator(poles);
    }

    private void computeZeros(String filterType, int order, double cutoffFreq, double centerFreq, boolean isDigital) {
        zeros.clear();
        switch (filterType) {
            case "Low-pass":
                // No finite zeros for low-pass Bessel filter
                break;
            case "High-pass":
                for (int i = 0; i < order; i++) {
                    zeros.add(new Complex(0, 0)); // Zeros at origin
                }
                break;
            case "Band-pass":
                double omega0 = 2 * Math.PI * centerFreq;
                for (int i = 0; i < order; i++) {
                    zeros.add(new Complex(0, omega0)); // Zeros at ±jω0
                    zeros.add(new Complex(0, -omega0));
                }
                break;
            case "Band-stop":
                for (int i = 0; i < order; i++) {
                    zeros.add(new Complex(0, 0)); // Zeros at origin
                }
                break;
        }
    }

    public ArrayList<Complex> getPoles() {
        return new ArrayList<>(poles);
    }

    public ArrayList<Complex> getZeros() {
        return new ArrayList<>(zeros);
    }

    public double[] getNumeratorCoefficients(String filterType, int order) {
        // Assuming unity gain numerator for simplicity
        double[] num = new double[order + 1];
        switch (filterType) {
            case "Low-pass":
                num[0] = 1.0; // Constant numerator
                break;
            case "High-pass":
                num[order] = 1.0; // s^order
                break;
            case "Band-pass":
                num[order] = 1.0; // s^order
                break;
            case "Band-stop":
                num[0] = 1.0; // Constant numerator
                num[order] = 1.0; // s^order
                break;
        }
        return num;
    }

    public double[][][] getStateSpace(String filterType, int order, double[] denominator) {
        // Companion form state-space representation
        int n = order;
        double[][] A = new double[n][n];
        double[][] B = new double[n][1];
        double[][] C = new double[1][n];
        double[][] D = new double[1][1];

        // A matrix: Companion form
        for (int i = 0; i < n - 1; i++) {
            A[i][i + 1] = 1.0;
        }
        for (int i = 0; i < n; i++) {
            A[n - 1][i] = -denominator[i];
        }

        // B matrix
        B[n - 1][0] = 1.0;

        // C matrix
        C[0][0] = 1.0;

        // D matrix
        D[0][0] = 0.0;

        return new double[][][]{A, B, C, D};
    }

    public double[][] getSecondOrderSections(double[] denominator, String filterType, int order) {
        // Determine effective order based on filter type
        int effectiveOrder = filterType.equals("Band-pass") || filterType.equals("Band-stop") ? 2 * order : order;
        int numSections = (effectiveOrder + 1) / 2; // Number of second-order sections
        double[][] sos = new double[numSections][6]; // [b0, b1, b2, a0, a1, a2]

        for (int i = 0; i < numSections; i++) {
            // Default coefficients
            sos[i][0] = 1.0; // b0
            sos[i][1] = 0.0; // b1
            sos[i][2] = 0.0; // b2
            sos[i][3] = 1.0; // a0
            sos[i][4] = 0.0; // a1
            sos[i][5] = 0.0; // a2

            // Assign coefficients only if poles are available
            if (i * 2 < poles.size()) {
                Complex pole1 = poles.get(i * 2);
                sos[i][4] = -2 * pole1.getReal(); // a1 = -2 * Re(pole1)
                if (i * 2 + 1 < poles.size()) {
                    Complex pole2 = poles.get(i * 2 + 1);
                    // a2 = Re(pole1)^2 + Im(pole1)^2 (for conjugate pairs)
                    sos[i][5] = pole1.getReal() * pole1.getReal() + pole1.getImaginary() * pole1.getImaginary();
                } else if (pole1.isReal()) {
                    // Handle single real pole (first-order section)
                    sos[i][5] = 0.0; // No quadratic term
                }
            }
        }
        return sos;
    }

    public int autoDetermineOrder(double cutoffFreq) {
        return Math.max(1, Math.min(10, (int) (Math.log10(cutoffFreq / 100) * 2)));
    }

    public double computeMagnitudeResponse(double[] coeffs, double freq, boolean isDigital, double sampleRate) {
        double omega = isDigital ? Math.tan(Math.PI * freq / sampleRate) : 2 * Math.PI * freq;
        Complex s = new Complex(0, omega);
        Complex denominator = evaluatePolynomial(coeffs, s);
        double magnitude = 1.0 / denominator.abs();
        if (Double.isNaN(magnitude) || Double.isInfinite(magnitude)) {
            System.out.println("Invalid magnitude at freq " + freq + ": " + magnitude);
            return 0.0;
        }
        return magnitude;
    }

    public double computePhaseResponse(double[] coeffs, double freq, boolean isDigital, double sampleRate) {
        double omega = isDigital ? Math.tan(Math.PI * freq / sampleRate) : 2 * Math.PI * freq;
        Complex s = new Complex(0, omega);
        Complex denominator = evaluatePolynomial(coeffs, s);
        double phase = -Math.toDegrees(Math.atan2(denominator.getImaginary(), denominator.getReal()));
        if (Double.isNaN(phase) || Double.isInfinite(phase)) {
            System.out.println("Invalid phase at freq " + freq + ": " + phase);
            return 0.0;
        }
        return phase;
    }

    public double computeGroupDelay(double[] coeffs, double freq, boolean isDigital, double sampleRate) {
        double delta = 1e-6;
        double omega = isDigital ? Math.tan(Math.PI * freq / sampleRate) : 2 * Math.PI * freq;
        double phase1 = computePhaseResponse(coeffs, freq, isDigital, sampleRate);
        double phase2 = computePhaseResponse(coeffs, freq + delta, isDigital, sampleRate);
        double groupDelay = -(phase2 - phase1) / (2 * Math.PI * delta);
        if (Double.isNaN(groupDelay) || Double.isInfinite(groupDelay)) {
            System.out.println("Invalid group delay at freq " + freq + ": " + groupDelay);
            return 0.0;
        }
        return groupDelay;
    }

    public double[] computeStepResponse(double[] coeffs, double sampleRate, int samples) {
        if (coeffs == null || coeffs.length < 1 || Double.isNaN(sampleRate) || sampleRate <= 0) {
            System.out.println("Invalid input for step response: coeffs=" + Arrays.toString(coeffs) + ", sampleRate=" + sampleRate);
            return new double[samples];
        }

        double[] response = new double[samples];
        double[] state = new double[coeffs.length - 1];
        double T = 1.0 / sampleRate;
        for (int n = 0; n < samples; n++) {
            double input = 1.0; // Step input
            double output = input;
            for (int i = 0; i < state.length; i++) {
                if (!Double.isFinite(coeffs[i]) || !Double.isFinite(state[i])) {
                    System.out.println("Invalid state or coefficient at sample " + n + ": coeff[" + i + "]=" + coeffs[i] + ", state[" + i + "]=" + state[i]);
                    return new double[samples]; // Return zeros if unstable
                }
                output -= coeffs[i] * state[i];
            }
            if (Double.isNaN(output) || Double.isInfinite(output)) {
                System.out.println("Invalid step response at sample " + n + ": " + output);
                response[n] = 0.0;
            } else {
                response[n] = output;
            }
            for (int i = state.length - 1; i > 0; i--) {
                state[i] = state[i - 1];
            }
            state[0] = output;
        }
        return response;
    }

    public double[] computeImpulseResponse(double[] coeffs, double sampleRate, int samples) {
        if (coeffs == null || coeffs.length < 1 || Double.isNaN(sampleRate) || sampleRate <= 0) {
            System.out.println("Invalid input for impulse response: coeffs=" + Arrays.toString(coeffs) + ", sampleRate=" + sampleRate);
            return new double[samples];
        }

        double[] response = new double[samples];
        double[] state = new double[coeffs.length - 1];
        double T = 1.0 / sampleRate;
        for (int n = 0; n < samples; n++) {
            double input = (n == 0) ? 1.0 : 0.0; // Impulse input
            double output = input;
            for (int i = 0; i < state.length; i++) {
                if (!Double.isFinite(coeffs[i]) || !Double.isFinite(state[i])) {
                    System.out.println("Invalid state or coefficient at sample " + n + ": coeff[" + i + "]=" + coeffs[i] + ", state[" + i + "]=" + state[i]);
                    return new double[samples]; // Return zeros if unstable
                }
                output -= coeffs[i] * state[i];
            }
            if (Double.isNaN(output) || Double.isInfinite(output)) {
                System.out.println("Invalid impulse response at sample " + n + ": " + output);
                response[n] = 0.0;
            } else {
                response[n] = output;
            }
            for (int i = state.length - 1; i > 0; i--) {
                state[i] = state[i - 1];
            }
            state[0] = output;
        }
        return response;
    }

    private double[] convolve(double[] a, double[] b) {
        double[] result = new double[a.length + b.length - 1];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b.length; j++) {
                result[i + j] += a[i] * b[j];
            }
        }
        return result;
    }

    private Complex evaluatePolynomial(double[] coeffs, Complex x) {
        Complex result = new Complex(0, 0);
        for (int i = coeffs.length - 1; i >= 0; i--) {
            result = result.multiply(x).add(new Complex(coeffs[i], 0));
        }
        return result;
    }

    public String getFilterDetails() {
        return designDetails.toString();
    }

    public static class Complex {
        private final double re;
        private final double im;

        public Complex(double real, double imaginary) {
            this.re = real;
            this.im = imaginary;
        }

        public double getReal() {
            return re;
        }

        public double getImaginary() {
            return im;
        }

        public Complex add(Complex other) {
            return new Complex(re + other.re, im + other.im);
        }

        public Complex subtract(Complex other) {
            return new Complex(re - other.re, im - other.im);
        }

        public Complex multiply(Complex other) {
            return new Complex(re * other.re - im * other.im, re * other.im + im * other.re);
        }

        public Complex multiply(double scalar) {
            return new Complex(re * scalar, im * scalar);
        }

        public Complex divide(Complex other) {
            double denom = other.re * other.re + other.im * other.im;
            if (denom == 0) return new Complex(0, 0);
            return new Complex((re * other.re + im * other.im) / denom, (im * other.re - re * other.im) / denom);
        }

        public double abs() {
            return Math.sqrt(re * re + im * im);
        }

        public Complex sqrt() {
            double r = Math.sqrt(abs());
            double theta = Math.atan2(im, re) / 2;
            return new Complex(r * Math.cos(theta), r * Math.sin(theta));
        }

        public boolean isReal() {
            return Math.abs(im) < 1e-10;
        }

        @Override
        public String toString() {
            return im >= 0 ? re + " + " + im + "i" : re + " - " + Math.abs(im) + "i";
        }
    }
}