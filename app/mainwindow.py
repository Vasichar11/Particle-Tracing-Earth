import sys
import os
import subprocess
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QLineEdit, QVBoxLayout, QHBoxLayout, QPushButton, QRadioButton, QScrollArea, QTabWidget, QProgressBar
from PyQt5.QtCore import QThread, pyqtSignal, QThread

class MakeThread(QThread):
    make_finished = pyqtSignal(bool, str, str)
    make_progress = pyqtSignal(int)

    def __init__(self, window):
        super().__init__()
        self.window = window

    def run(self):
        # Apply the line edit values to the constants.h file
        self.window.apply_line_edit_values()

        for progress in range(1, 101):
            self.make_progress.emit(progress)
            QThread.msleep(150)  # Sleep for 50 milliseconds (0.05 seconds)

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

        # Read constant values from "constants.h" file
        constant_values = self.read_constants_file("headers/constants.h")

        # Add labels and line edits for constant values
        self.constant_labels = []
        self.constant_lineedits = []
        # Add labels and line edits for constant values
        for label_text, default_value, comment in constant_values:
            label = QLabel(label_text, self.tab1_widget)
            label.setToolTip(comment)  # Set the comment as a tooltip
            lineedit = QLineEdit(default_value, self.tab1_widget)
            self.constant_labels.append(label)
            self.constant_lineedits.append(lineedit)

            # Create a horizontal layout for each constant value
            constant_layout = QHBoxLayout()
            constant_layout.addWidget(label)
            constant_layout.addWidget(lineedit)

            self.tab1_layout.addLayout(constant_layout)

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

        # Add radio button for simulation type
        self.simulation_radio = QRadioButton("Bell equations")
        self.tab1_layout.addWidget(self.simulation_radio)
        
        # Run simulation button
        simulate_button = QPushButton("Simulate", self.tab1_widget)
        simulate_button.clicked.connect(self.simulate_clicked)  
        self.tab1_layout.addWidget(simulate_button)

        # Add the first tab to the tab widget
        self.tab_widget.addTab(tab1_scroll_area, "Simulation")

        # Create the "Visualization" tab
        self.plots_tab = QWidget(self)
        self.plots_layout = QVBoxLayout(self.plots_tab)
        self.tab_widget.addTab(self.plots_tab, "Visualization")

        # Add the tab widget to the main layout
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

    def update_progress_bar(self, progress):
        self.progress_bar.setValue(progress)

    def read_constants_file(self, filename):
        constants_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), filename)
        constant_values = []
        with open(constants_path, "r") as file:
            for line in file:
                if line.strip().startswith("const"):
                    parts = line.strip().split("=")
                    if len(parts) == 2:
                        constant_name = parts[0].split()[-1].strip()
                        constant_value = parts[1].split(";")[0].strip()

                        comment_start_index = line.find("//")
                        if comment_start_index != -1:
                            comment = line[comment_start_index+2:].strip()
                        else:
                            comment = ""

                        constant_values.append((constant_name, constant_value, comment))
        return constant_values

    def make_clicked(self):
        # Handle the make button click event
        self.make_thread.start()

    def handle_make_finished(self, success, stdout, stderr):
        if success:
            print("Make successful")
            self.progress_bar.setValue(100)
        else:
            print("Make failed")
        print("STDOUT:", stdout)
        print("STDERR:", stderr)

    def apply_line_edit_values(self):
        constants_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "headers/constants.h")
        with open(constants_path, "r") as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if line.strip().startswith("const"):
                parts = line.strip().split("=")
                if len(parts) == 2:
                    constant_name = parts[0].split()[-1].strip()
                    constant_value = parts[1].split(";")[0].strip()

                    # Find the corresponding line edit for the constant
                    for j in range(len(self.constant_labels)):
                        if constant_name == self.constant_labels[j].text():
                            new_value = self.constant_lineedits[j].text()

                            # Replace the value in the line with the new value
                            lines[i] = line.replace(constant_value, new_value)

        with open(constants_path, "w") as file:
            file.writelines(lines)

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
