import sys
import os
import subprocess
import xml.etree.ElementTree as ET
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QLineEdit, QVBoxLayout, QHBoxLayout, QPushButton, QRadioButton, QScrollArea, QTabWidget, QProgressBar, QGroupBox, QComboBox

from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtGui import QIcon

class MakeThread(QThread):
    make_finished = pyqtSignal(bool, str, str)
    make_progress = pyqtSignal(int)

    def __init__(self, window):
        super().__init__()
        self.window = window

    def run(self):
            # Apply the line edit values to the constants.h file
            self.window.apply_line_edit_values()


            # Execute the "make" command
            current_directory = os.getcwd()
            process = subprocess.Popen(["make"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=current_directory)
            stdout, stderr = process.communicate()

            if process.returncode == 0:
                success = True
                message = "Make successful"
            else:
                success = False
                message = "Make failed"

            self.make_finished.emit(success, stdout.decode("utf-8"), stderr.decode("utf-8"))

class AllCleanThread(QThread):
    allclean_finished = pyqtSignal(bool, str, str)

    def __init__(self):
        super().__init__()

    def run(self):
        # Execute the "make allclean" command
        current_directory = os.getcwd()
        process = subprocess.Popen(["make", "allclean"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=current_directory)
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            success = True
            message = "Make allclean successful"
        else:
            success = False
            message = "Make allclean failed"

        self.allclean_finished.emit(success, stdout.decode("utf-8"), stderr.decode("utf-8"))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Wave-Particle Interaction Simulation")
        self.setMinimumSize(800, 600)  # Set the minimum window size
        self.setMaximumSize(800, 600)  # Set the maximum window size

        current_directory = os.getcwd()
        app_directory = os.path.join(current_directory, "app")

        # Set app icon
        icon_path = os.path.join(app_directory, "resources", "icons", "icon.png")
        app_icon = QIcon(icon_path)
        self.setWindowIcon(app_icon)

        # Create the main widget and layout
        self.main_widget = QWidget(self)
        self.layout = QVBoxLayout(self.main_widget)

        # Add a tab widget
        self.tab_widget = QTabWidget(self.main_widget)


        # Add the first tab
        self.tab1_widget = QWidget(self)
        self.tab1_layout = QVBoxLayout(self.tab1_widget)

        # Create a scroll area for the first tab
        tab1_scroll_area = QScrollArea(self)
        tab1_scroll_area.setWidgetResizable(True)
        tab1_scroll_area.setWidget(self.tab1_widget)

        self.xml_filepath = os.path.join(current_directory, "app", "resources", "configuration", "constants.xml")
        # Read constant values from "constants.xml" file
        constants_by_namespace = self.read_config_file()
        
        # Dictionary with all the metric units and their multiplication factors
        self.metric_units = {"T":1000000000000, "G":1000000000, "M":1000000, "K":1000, "1":1, "m":0.001, "μ":0.000001, "n":0.000000001, "p":0.000000000001}

        # Add group boxes and layouts for each namespace
        for namespace, constants in constants_by_namespace.items():
            group_box = QGroupBox(namespace)
            group_layout = QVBoxLayout(group_box)

            for name, value, description, metric_unit, physical_unit in constants:
                label = QLabel(name, self.tab1_widget)
                label.setToolTip(description)  # Set the description as a tooltip
                lineedit = QLineEdit(value, self.tab1_widget)

                # Dropdown for Metric Unit
                combo_box_metric = QComboBox(group_box)
                for key in self.metric_units:
                    combo_box_metric.addItem(key)
                combo_box_metric.setCurrentText(metric_unit)  # Set default value from the XML document

                # Dropdown for Physical Unit
                lineedit_physical = QLineEdit()
                lineedit_physical.setReadOnly(True) # The physical value is fixed and comes from the XML document
                lineedit_physical.setText(physical_unit)  # Set default value
                lineedit_physical.setFixedWidth(100)  # Set the width of the QLineEdit
                lineedit_physical.setFixedHeight(30)  # Set the height of the QLineEdit

                layout = QHBoxLayout()
                layout.addWidget(label)
                layout.addWidget(lineedit)
                layout.addWidget(combo_box_metric)
                layout.addWidget(lineedit_physical)

                group_layout.addLayout(layout)

            self.tab1_layout.addWidget(group_box)

        # Add buttons for make, make allclean, and simulate
        make_button = QPushButton("Build", self.tab1_widget)
        make_button.clicked.connect(self.make_clicked)
        allclean_button = QPushButton("Clean", self.tab1_widget)
        allclean_button.clicked.connect(self.allclean_clicked)

        # Create a layout for the buttons and progress bar
        button_layout = QHBoxLayout()
        button_layout.addWidget(make_button)
        button_layout.addWidget(allclean_button)

        # Add the buttons layout to the main layout
        self.tab1_layout.addLayout(button_layout)

        # Add the progress bar
        self.progress_bar = QProgressBar(self.tab1_widget)
        self.tab1_layout.addWidget(self.progress_bar)

        # Add labels and line edits for simulation time
        no_wpi_label = QLabel("No WPI simulation time", self.tab1_widget)
        self.no_wpi_lineedit = QLineEdit(self.tab1_widget)
        wpi_label = QLabel("WPI simulation time", self.tab1_widget)
        self.wpi_lineedit = QLineEdit(self.tab1_widget)
        self.tab1_layout.addWidget(no_wpi_label)
        self.tab1_layout.addWidget(self.no_wpi_lineedit)
        self.tab1_layout.addWidget(wpi_label)
        self.tab1_layout.addWidget(self.wpi_lineedit)

        # Apply styles to line edits
        self.apply_line_edit_styles(self.no_wpi_lineedit)
        self.apply_line_edit_styles(self.wpi_lineedit)


        # Add radio button for simulation type
        self.simulation_radio = QRadioButton("Bell equations")
        self.tab1_layout.addWidget(self.simulation_radio)

        # Run simulation button
        simulate_button = QPushButton("Simulate", self.tab1_widget)
        simulate_button.clicked.connect(self.simulate_clicked)
        self.tab1_layout.addWidget(simulate_button)

        # Add the first tab to the tab widget
        self.tab_widget.addTab(tab1_scroll_area, "Simulation")
        self.layout.addWidget(self.tab_widget)

        # Set the central widget
        self.setCentralWidget(self.main_widget)

        # Create make and allclean threads
        self.make_thread = MakeThread(self)
        self.allclean_thread = AllCleanThread()

        # Connect signals from threads to update GUI
        self.make_thread.make_finished.connect(self.handle_make_finished)
        self.allclean_thread.allclean_finished.connect(self.handle_allclean_finished)
        self.make_thread.make_progress.connect(self.update_progress_bar)

    def apply_line_edit_styles(self, line_edit):
        line_edit.setStyleSheet("""
            QLineEdit {
                padding: 5px;
                border: 1px solid #ccc;
                border-radius: 3px;
                background-color: #fff;
                font-size: 14px;
                width: 400px; 
            }
            QLineEdit:focus {
                border: 1px solid #6e98e0;
                outline: none;
            }
        """)


    def read_config_file(self):
        constants_by_namespace = {}

        tree = ET.parse(self.xml_filepath)
        root = tree.getroot()

        for namespace_element in root:
            namespace = namespace_element.tag
            constants = []

            for element in namespace_element:
                name = element.find("name").text
                value = element.find("value").text
                description = element.find("description").text
                metric_unit = element.find("metric_unit").text
                physical_unit = element.find("physical_unit").text

                constants.append((name, value, description, metric_unit, physical_unit))

            constants_by_namespace[namespace] = constants

        return constants_by_namespace

    def handle_make_finished(self, success, stdout, stderr):
            if success:
                print("Make successful")
                self.progress_bar.setValue(100)
            else:
                print("Make failed")
            print("STDOUT:", stdout)
            print("STDERR:", stderr)

    def update_progress_bar(self, progress):
        self.progress_bar.setValue(progress)

    def make_clicked(self):
        # Reset the progress bar to 0%
        self.progress_bar.setValue(0)

        # Handle the make button click event
        self.make_thread.start()

    def apply_line_edit_values(self):
        edited_xml_filepath = os.path.join(self.xml_filepath.split("constants.xml")[0], "edited_constants.xml")
        tree = ET.parse(self.xml_filepath)
        root = tree.getroot()

        for namespace_element in root:
            for element in namespace_element:
                name = element.find("name").text
                line_edit_value = self.get_line_edit_value(name)
                # Edit new value in the new xml configuration file
                element.find("value").text = line_edit_value
                # The value is multiplied by the corresponding multiplication factor
                # That means that the new metric unit will be 1
                element.find("metric_unit").text = "1"

        tree.write(edited_xml_filepath)

    def get_line_edit_value(self, name):
        for group_box in self.tab1_widget.findChildren(QGroupBox):
            for layout in group_box.findChildren(QVBoxLayout):
                for index in range(layout.count()):
                    item = layout.itemAt(index)
                    if isinstance(item, QHBoxLayout):
                        label = item.itemAt(0).widget()
                        line_edit = item.itemAt(1).widget()
                        combobox = item.itemAt(2).widget()  
                        if isinstance(label, QLabel) and label.text() == name and isinstance(line_edit, QLineEdit) and isinstance(combobox, QComboBox):
                            # Typecast to float to multiply with the multiplicative factor
                            value = float(line_edit.text())
                            unit = combobox.currentText()
                            multiplication_factor = self.metric_units[unit]
                            new_value = value * multiplication_factor
                            # Typecast back to string
                            new_value = str(new_value)
                            return new_value

        return ""  # Return an empty string if the QLineEdit for the given name is not found


    def allclean_clicked(self):
        # Reset the progress bar to 0%
        self.progress_bar.setValue(0)
        # Handle the make allclean button click event
        self.allclean_thread.start()

    def handle_allclean_finished(self, success, stdout, stderr):
        if success:
            print("Make allclean successful")
        else:
            print("Make allclean failed")
        print("STDOUT:", stdout)
        print("STDERR:", stderr)

    def simulate_clicked(self):
        no_wpi_value = self.no_wpi_lineedit.text()
        wpi_value = self.wpi_lineedit.text()
        simulation_type = "-bell" if self.simulation_radio.isChecked() else ""

        # Execute the "./tracer" command with the provided arguments
        command = ["./tracer", no_wpi_value, wpi_value, simulation_type]
        current_directory = os.getcwd()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=current_directory)
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            success = True
            message = "Simulation successful"
        else:
            success = False
            message = "Simulation failed"

        print(message)
        print("STDOUT:", stdout.decode("utf-8"))
        print("STDERR:", stderr.decode("utf-8"))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())