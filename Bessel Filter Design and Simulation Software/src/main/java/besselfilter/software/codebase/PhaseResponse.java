package besselfilter.software.codebase;

import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart;

public class PhaseResponse {
    public void plotPhaseResponse(LineChart<Number, Number> chart, double[] coefficients, boolean isDigital, double sampleRate, double cutoffFreq, BesselFilterModel model) {
        chart.getData().clear();
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Phase Response");

        double maxFreq = isDigital ? sampleRate / 2 : cutoffFreq * 10;
        int points = 500;
        System.out.println("Plotting phase response, maxFreq: " + maxFreq + ", points: " + points);
        for (int i = 0; i <= points; i++) {
            double freq = (maxFreq / points) * i;
            double phase = model.computePhaseResponse(coefficients, freq, isDigital, sampleRate);
            if (Double.isFinite(phase)) {
                series.getData().add(new XYChart.Data<>(freq, phase));
            }
            if (i % 100 == 0) {
                System.out.println("Phase at freq " + freq + ": " + phase);
            }
        }

        chart.getData().add(series);
        chart.requestLayout();
    }

    public void plotGroupDelay(LineChart<Number, Number> chart, double[] coefficients, boolean isDigital, double sampleRate, double cutoffFreq, BesselFilterModel model) {
        chart.getData().clear();
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName("Group Delay");

        double maxFreq = isDigital ? sampleRate / 2 : cutoffFreq * 10;
        int points = 500;
        System.out.println("Plotting group delay, maxFreq: " + maxFreq + ", points: " + points);
        for (int i = 0; i <= points; i++) {
            double freq = (maxFreq / points) * i;
            double groupDelay = model.computeGroupDelay(coefficients, freq, isDigital, sampleRate);
            if (Double.isFinite(groupDelay)) {
                series.getData().add(new XYChart.Data<>(freq, groupDelay));
            }
            if (i % 100 == 0) {
                System.out.println("Group delay at freq " + freq + ": " + groupDelay);
            }
        }

        chart.getData().add(series);
        chart.requestLayout();
    }
}