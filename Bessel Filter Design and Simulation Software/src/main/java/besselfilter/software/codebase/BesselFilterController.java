package besselfilter.software.codebase;

import javafx.application.Platform;
import javafx.fxml.FXML;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.control.*;
import javafx.collections.FXCollections;

public class BesselFilterController {
    @FXML private ComboBox<String> filterTypeCombo;
    @FXML private TextField orderField, cutoffField, centerField, sampleRateField, plotDurationField;
    @FXML private CheckBox autoOrderCheck, digitalCheck, linearScaleCheck;
    @FXML private Button designButton;
    @FXML private LineChart<Number, Number> magnitudeChart, phaseChart, groupDelayChart, stepResponseChart, impulseResponseChart;
    @FXML private ScatterChart<Number, Number> poleZeroChart;
    @FXML private NumberAxis magnitudeYAxis;
    @FXML private TabPane responseTabs;
    @FXML private TextArea outputArea;

    private BesselFilterModel model;
    private MagnitudeResponse magnitudeResponse;
    private PhaseResponse phaseResponse;
    private StepAndImpulseResponse stepImpulseResponse;
    private PoleZeroPlot poleZeroPlot;
    private FilterCoefficientOutput coefficientOutput;

    @FXML
    public void initialize() {
        model = new BesselFilterModel();
        magnitudeResponse = new MagnitudeResponse();
        phaseResponse = new PhaseResponse();
        stepImpulseResponse = new StepAndImpulseResponse();
        poleZeroPlot = new PoleZeroPlot();
        coefficientOutput = new FilterCoefficientOutput();
        filterTypeCombo.setItems(FXCollections.observableArrayList("Low-pass", "High-pass", "Band-pass", "Band-stop"));
        filterTypeCombo.setValue("Low-pass");
        configureChart(magnitudeChart);
        configureChart(phaseChart);
        configureChart(groupDelayChart);
        configureChart(stepResponseChart);
        configureChart(impulseResponseChart);
        configureChart(poleZeroChart);
        linearScaleCheck.setSelected(false);
        magnitudeYAxis.setLabel("Magnitude (dB)");
        plotDurationField.setText("0.01"); // Default plot duration: 10ms
    }

    private void configureChart(LineChart<Number, Number> chart) {
        chart.setCreateSymbols(false);
        chart.setAnimated(false);
        chart.getXAxis().setAutoRanging(true);
        chart.getYAxis().setAutoRanging(true);
    }

    private void configureChart(ScatterChart<Number, Number> chart) {
        chart.setAnimated(false);
        chart.getXAxis().setAutoRanging(true);
        chart.getYAxis().setAutoRanging(true);
    }

    @FXML
    private void designFilter() {
        try {
            String filterType = filterTypeCombo.getValue();
            double cutoffFreq = Double.parseDouble(cutoffField.getText());
            double centerFreq = filterType.contains("Band") ? Double.parseDouble(centerField.getText()) : 0;
            int order = autoOrderCheck.isSelected() ? model.autoDetermineOrder(cutoffFreq) : Integer.parseInt(orderField.getText());
            boolean isDigital = digitalCheck.isSelected();
            double sampleRate = isDigital ? Double.parseDouble(sampleRateField.getText()) : 44100;
            double plotDuration = Double.parseDouble(plotDurationField.getText());

            if (plotDuration <= 0) {
                throw new NumberFormatException("Plot duration must be positive.");
            }

            double[] coefficients = model.designBesselFilter(filterType, order, cutoffFreq, centerFreq, isDigital, sampleRate);
            String result = model.getFilterDetails();
            String coefficientDetails = coefficientOutput.getCoefficientOutput(filterType, order, coefficients, model);
            outputArea.setText(result + "\n" + coefficientDetails);

            Platform.runLater(() -> plotAllResponses(filterType, order, coefficients, isDigital, sampleRate, cutoffFreq, plotDuration));
        } catch (NumberFormatException e) {
            outputArea.setText("Error: Invalid input. Please enter valid numbers.\n" + e.getMessage());
        }
    }

    @FXML
    private void toggleMagnitudeScale() {
        designFilter();
    }

    private void plotAllResponses(String filterType, int order, double[] coefficients, boolean isDigital, double sampleRate, double cutoffFreq, double plotDuration) {
        magnitudeChart.getData().clear();
        phaseChart.getData().clear();
        groupDelayChart.getData().clear();
        stepResponseChart.getData().clear();
        impulseResponseChart.getData().clear();
        poleZeroChart.getData().clear();

        magnitudeResponse.plotMagnitudeResponse(magnitudeChart, magnitudeYAxis, linearScaleCheck.isSelected(), coefficients, isDigital, sampleRate, cutoffFreq, model);
        phaseResponse.plotPhaseResponse(phaseChart, coefficients, isDigital, sampleRate, cutoffFreq, model);
        phaseResponse.plotGroupDelay(groupDelayChart, coefficients, isDigital, sampleRate, cutoffFreq, model);
        stepImpulseResponse.plotStepResponse(stepResponseChart, coefficients, sampleRate, plotDuration, model);
        stepImpulseResponse.plotImpulseResponse(impulseResponseChart, coefficients, sampleRate, plotDuration, model);
        poleZeroPlot.plotPoleZero(poleZeroChart, model, isDigital);

        magnitudeChart.requestLayout();
        phaseChart.requestLayout();
        groupDelayChart.requestLayout();
        stepResponseChart.requestLayout();
        impulseResponseChart.requestLayout();
        poleZeroChart.requestLayout();
    }
}