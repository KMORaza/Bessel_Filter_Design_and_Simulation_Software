module besselfilter.software.codebase.besselfilter {
    requires javafx.controls;
    requires javafx.fxml;


    opens besselfilter.software.codebase to javafx.fxml;
    exports besselfilter.software.codebase;
}