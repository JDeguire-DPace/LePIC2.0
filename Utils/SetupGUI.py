#!/usr/bin/env python3
import sys
import os
from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QFileDialog,
    QVBoxLayout, QMessageBox
)
from PySide6.QtGui import QCursor
from PySide6.QtCore import Qt

# --- Directories ---
START_DIR = os.getcwd()                                  # where the user runs the script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # where the script resides
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "PathSimulation.txt")  # always written next to SetupGUI.py
INPUT_DIR = os.path.join(SCRIPT_DIR, START_DIR +"/Inputs")
os.makedirs(INPUT_DIR, exist_ok=True)  

class LePICSetup(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LePIC – Simulation Setup")
        self.resize(420, 180)
        self.setWindowFlag(Qt.WindowStaysOnTopHint)
        self.center_on_screen()

        # --- UI ---
        self.label = QLabel("Is this a NEW simulation?", alignment=Qt.AlignCenter)
        self.label.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")

        self.btn_yes = QPushButton("Yes")
        self.btn_no = QPushButton("No")
        self.btn_yes.clicked.connect(self.new_simulation)
        self.btn_no.clicked.connect(self.close)

        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.btn_yes)
        layout.addWidget(self.btn_no)
        self.setLayout(layout)

    def center_on_screen(self):
        """Center the window on the primary display."""
        screen = QApplication.primaryScreen()
        screen_geometry = screen.availableGeometry()
        window_geometry = self.frameGeometry()
        window_geometry.moveCenter(screen_geometry.center())
        self.move(window_geometry.topLeft())

    def new_simulation(self):
        """Ask user for 4 files and save paths."""
        files = {}
        prompts = [
            ("Select conditions file", "Conditions"),
            ("Select particle data file", "Particles"),
            ("Select geometry file", "Geometry"),
            ("Select boundary file", "Boundary"),
            ("Select field file", "Field")
        ]

        for title, key in prompts:
            file_path, _ = QFileDialog.getOpenFileName(
                self, title, INPUT_DIR, "All Files (*.*)"
            )
            if not file_path:
                QMessageBox.warning(self, "LePIC", f"{key} file not selected. Aborting setup.")
                return
            files[key] = file_path

        # Write selected paths
        with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
            for k, v in files.items():
                f.write(f"{k} = {v}\n")

        QMessageBox.information(self, "LePIC", f"Paths saved to:\n{OUTPUT_FILE}")
        self.close()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = LePICSetup()
    window.show()
    sys.exit(app.exec())
