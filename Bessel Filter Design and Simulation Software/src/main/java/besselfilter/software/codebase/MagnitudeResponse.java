package besselfilter.software.codebase;

import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;

public class MagnitudeResponse {
    public void plotMagnitudeResponse(LineChart<Number, Number> chart, NumberAxis yAxis, boolean linearScale, double[] coefficients, boolean isDigital, double sampleRate, double cutoffFreq, BesselFilterModel model) {
        chart.getData().clear();
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Magnitude Response");

        double maxFreq = isDigital ? sampleRate / 2 : cutoffFreq * 10;
        int points = 500;
        System.out.println("Plotting magnitude response, maxFreq: " + maxFreq + ", points: " + points);
        for (int i = 0; i <= points; i++) {
            double freq = (maxFreq / points) * i;
            double magnitude = model.computeMagnitudeResponse(coefficients, freq, isDigital, sampleRate);
            double value = linearScale ? magnitude : 20 * Math.log10(Math.max(magnitude, 1e-10)); // Avoid log(0)
            if (Double.isFinite(value)) {
                series.getData().add(new XYChart.Data<>(freq, value));
            }
            if (i % 100 == 0) {
                System.out.println("Magnitude at freq " + freq + ": " + value);
            }
        }

        XYChart.Series<Number, Number> cutoffSeries = new XYChart.Series<>();
        cutoffSeries.setName("Cutoff Frequency");
        double cutoffMagnitude = model.computeMagnitudeResponse(coefficients, cutoffFreq, isDigital, sampleRate);
        double cutoffValue = linearScale ? cutoffMagnitude : 20 * Math.log10(Math.max(cutoffMagnitude, 1e-10));
        if (Double.isFinite(cutoffValue)) {
            cutoffSeries.getData().add(new XYChart.Data<>(cutoffFreq, cutoffValue));
            cutoffSeries.getData().forEach(data -> {
                if (data.getNode() != null) {
                    data.getNode().setStyle("-fx-shape: \"M0,0 L0,10\"; -fx-stroke: red; -fx-stroke-width: 2;");
                }
            });
        }
        System.out.println("Cutoff at freq " + cutoffFreq + ": " + cutoffValue);

        chart.getData().addAll(series, cutoffSeries);
        yAxis.setLabel(linearScale ? "Magnitude (Linear)" : "Magnitude (dB)");
        chart.requestLayout();
    }
}