package besselfilter.software.codebase;

import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;

public class StepAndImpulseResponse {
    public void plotStepResponse(LineChart<Number, Number> chart, double[] coefficients, double sampleRate, double plotDuration, BesselFilterModel model) {
        chart.getData().clear();
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Step Response");

        int samples = (int) (plotDuration * sampleRate);
        samples = Math.max(100, Math.min(samples, 1000)); // Bound samples for performance
        double[] response = model.computeStepResponse(coefficients, sampleRate, samples);
        double T = 1.0 / sampleRate;

        double maxAmplitude = 0;
        int validPoints = 0;
        for (int i = 0; i < samples; i++) {
            double value = response[i];
            if (Double.isFinite(value)) {
                series.getData().add(new XYChart.Data<>(i * T, value));
                maxAmplitude = Math.max(maxAmplitude, Math.abs(value));
                validPoints++;
            }
        }

        if (validPoints == 0) {
            System.out.println("No valid data points for step response plot");
            return;
        }

        chart.getData().add(series);
        NumberAxis yAxis = (NumberAxis) chart.getYAxis();
        yAxis.setAutoRanging(false);
        yAxis.setLowerBound(-maxAmplitude * 1.1);
        yAxis.setUpperBound(maxAmplitude * 1.1);
        yAxis.setTickUnit(maxAmplitude / 5);
        chart.setCreateSymbols(false);
        chart.requestLayout();
        System.out.println("Plotted step response: samples=" + samples + ", maxAmplitude=" + maxAmplitude);
    }

    public void plotImpulseResponse(LineChart<Number, Number> chart, double[] coefficients, double sampleRate, double plotDuration, BesselFilterModel model) {
        chart.getData().clear();
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Impulse Response");

        int samples = (int) (plotDuration * sampleRate);
        samples = Math.max(100, Math.min(samples, 1000)); // Bound samples for performance
        double[] response = model.computeImpulseResponse(coefficients, sampleRate, samples);
        double T = 1.0 / sampleRate;

        double maxAmplitude = 0;
        int validPoints = 0;
        for (int i = 0; i < samples; i++) {
            double value = response[i];
            if (Double.isFinite(value)) {
                series.getData().add(new XYChart.Data<>(i * T, value));
                maxAmplitude = Math.max(maxAmplitude, Math.abs(value));
                validPoints++;
            }
        }

        if (validPoints == 0) {
            System.out.println("No valid data points for impulse response plot");
            return;
        }

        chart.getData().add(series);
        NumberAxis yAxis = (NumberAxis) chart.getYAxis();
        yAxis.setAutoRanging(false);
        yAxis.setLowerBound(-maxAmplitude * 1.1);
        yAxis.setUpperBound(maxAmplitude * 1.1);
        yAxis.setTickUnit(maxAmplitude / 5);
        chart.setCreateSymbols(false);
        chart.requestLayout();
        System.out.println("Plotted impulse response: samples=" + samples + ", maxAmplitude=" + maxAmplitude);
    }
}