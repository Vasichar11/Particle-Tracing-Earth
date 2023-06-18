import sys
from PyQt5.QtWidgets import QApplication
from app.mainwindow import MainWindow


def tracer():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    tracer()