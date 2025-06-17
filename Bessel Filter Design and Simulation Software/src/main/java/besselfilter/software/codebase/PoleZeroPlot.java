package besselfilter.software.codebase;

import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;

public class PoleZeroPlot {
    public void plotPoleZero(ScatterChart<Number, Number> chart, BesselFilterModel model, boolean isDigital) {
        chart.getData().clear();

        // Poles series (X)
        XYChart.Series<Number, Number> poleSeries = new XYChart.Series<>();
        poleSeries.setName("Poles");
        for (BesselFilterModel.Complex pole : model.getPoles()) {
            poleSeries.getData().add(new XYChart.Data<>(pole.getReal(), pole.getImaginary()));
        }

        // Zeros series (O)
        XYChart.Series<Number, Number> zeroSeries = new XYChart.Series<>();
        zeroSeries.setName("Zeros");
        for (BesselFilterModel.Complex zero : model.getZeros()) {
            zeroSeries.getData().add(new XYChart.Data<>(zero.getReal(), zero.getImaginary()));
        }

        // Unit circle for digital filters
        XYChart.Series<Number, Number> unitCircleSeries = new XYChart.Series<>();
        unitCircleSeries.setName("Unit Circle");
        if (isDigital) {
            int points = 100;
            for (int i = 0; i <= points; i++) {
                double theta = 2 * Math.PI * i / points;
                unitCircleSeries.getData().add(new XYChart.Data<>(Math.cos(theta), Math.sin(theta)));
            }
        }

        chart.getData().addAll(poleSeries, zeroSeries, unitCircleSeries);
        chart.requestLayout();
    }
}