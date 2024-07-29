import sys
import numpy as np
import logging
import tifffile
import os
from typing import Tuple, List, Optional, Any, Dict

from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget,
                             QPushButton, QSlider, QLabel, QSpinBox, QDoubleSpinBox,
                             QComboBox, QFileDialog, QMessageBox, QCheckBox, QDockWidget, QSizePolicy, QTableWidget,
                             QTableWidgetItem, QDialog, QGridLayout, QTabWidget, QTextEdit, QAction, QFormLayout, QGroupBox, QScrollArea)
from PyQt5.QtCore import Qt, QSize
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from PyQt5.QtCore import QTimer, QEvent
from PyQt5.QtGui import QColor, QVector3D, QImage, QMouseEvent, QWheelEvent
import traceback

from matplotlib import pyplot as plt
from skimage.feature import blob_log
from skimage import measure

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from scipy.ndimage import distance_transform_edt, center_of_mass, binary_dilation, gaussian_filter
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtCore import pyqtSlot

from data_generation import DataGenerator
from blob_detection import BlobAnalyzer
from biological_simulation import BiologicalSimulator
from file_operations import ImporterFactory
#from volume_processor import VolumeProcessor

import OpenGL
print(f"PyQtGraph version: {pg.__version__}")
print(f"OpenGL version: {OpenGL.__version__}")


class BiologicalSimulationWidget(QWidget):
    simulationRequested = pyqtSignal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()

    def initUI(self):
        main_layout = QVBoxLayout(self)

        # Create a scroll area
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        main_layout.addWidget(scroll_area)

        # Create a widget to hold all the content
        content_widget = QWidget()
        scroll_area.setWidget(content_widget)
        content_layout = QVBoxLayout(content_widget)

        # Create tab widget
        self.tab_widget = QTabWidget()
        content_layout.addWidget(self.tab_widget)

        # Create tabs
        self.create_protein_tab()
        self.create_structure_tab()
        self.create_organelles_tab()
        self.create_calcium_tab()

        # Simulate Button
        self.simulate_button = QPushButton("Simulate")
        self.simulate_button.clicked.connect(self.requestSimulation)
        main_layout.addWidget(self.simulate_button)

        # Set the size policy to allow the widget to be resized
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Set a minimum size for the widget
        self.setMinimumSize(400, 600)

    def create_protein_tab(self):
        protein_tab = QWidget()
        layout = QVBoxLayout(protein_tab)

        # Protein Dynamics group
        protein_group = QGroupBox("Protein Dynamics")
        protein_layout = QFormLayout(protein_group)

        self.diffusion_checkbox = QCheckBox()
        self.diffusion_coefficient = QDoubleSpinBox()
        self.diffusion_coefficient.setRange(0, 100)
        self.diffusion_coefficient.setValue(1)

        self.transport_checkbox = QCheckBox()
        self.transport_velocity = QDoubleSpinBox()
        self.transport_velocity.setRange(-10, 10)
        self.transport_velocity.setValue(1)

        protein_layout.addRow("Protein Diffusion:", self.diffusion_checkbox)
        protein_layout.addRow("Diffusion Coefficient:", self.diffusion_coefficient)
        protein_layout.addRow("Active Transport:", self.transport_checkbox)
        protein_layout.addRow("Transport Velocity:", self.transport_velocity)

        layout.addWidget(protein_group)
        self.tab_widget.addTab(protein_tab, "Protein Dynamics")

    def create_structure_tab(self):
        structure_tab = QWidget()
        layout = QVBoxLayout(structure_tab)

        # Cell Structure group
        structure_group = QGroupBox("Cell Structure")
        structure_layout = QFormLayout(structure_group)

        self.cell_type_combo = QComboBox()
        self.cell_type_combo.addItems(['spherical', 'neuron', 'epithelial', 'muscle'])
        self.cell_type_combo.currentTextChanged.connect(self.toggle_neuron_options)

        self.cell_membrane_checkbox = QCheckBox("Cell Membrane")
        self.nucleus_checkbox = QCheckBox("Nucleus")
        self.er_checkbox = QCheckBox("Endoplasmic Reticulum")
        self.mitochondria_checkbox = QCheckBox("Mitochondria")
        self.cytoskeleton_checkbox = QCheckBox("Cytoskeleton")

        structure_layout.addRow("Cell Type:", self.cell_type_combo)
        structure_layout.addRow(self.cell_membrane_checkbox)
        structure_layout.addRow(self.nucleus_checkbox)
        structure_layout.addRow(self.er_checkbox)
        structure_layout.addRow(self.mitochondria_checkbox)
        structure_layout.addRow(self.cytoskeleton_checkbox)

        layout.addWidget(structure_group)

        # Cell Parameters group
        params_group = QGroupBox("Cell Parameters")
        params_layout = QFormLayout(params_group)

        self.cell_radius = QSpinBox()
        self.cell_radius.setRange(5, 50)
        self.cell_radius.setValue(20)

        self.membrane_thickness = QSpinBox()
        self.membrane_thickness.setRange(1, 5)
        self.membrane_thickness.setValue(1)

        self.nucleus_radius = QSpinBox()
        self.nucleus_radius.setRange(1, 20)
        self.nucleus_radius.setValue(5)

        self.pixel_size_x = QDoubleSpinBox()
        self.pixel_size_y = QDoubleSpinBox()
        self.pixel_size_z = QDoubleSpinBox()
        for spinbox in [self.pixel_size_x, self.pixel_size_y, self.pixel_size_z]:
            spinbox.setRange(0.1, 10)
            spinbox.setSingleStep(0.1)
            spinbox.setValue(1)

        self.er_density = QDoubleSpinBox()
        self.er_density.setRange(0.05, 0.2)
        self.er_density.setSingleStep(0.01)
        self.er_density.setValue(0.1)


        params_layout.addRow("Cell Radius:", self.cell_radius)
        params_layout.addRow("Membrane Thickness:", self.membrane_thickness)
        params_layout.addRow("Nucleus Radius:", self.nucleus_radius)
        params_layout.addRow("ER Density:", self.er_density)
        params_layout.addRow("Pixel Size X:", self.pixel_size_x)
        params_layout.addRow("Pixel Size Y:", self.pixel_size_y)
        params_layout.addRow("Pixel Size Z:", self.pixel_size_z)

        layout.addWidget(params_group)

        # Neuron-specific options
        self.neuron_options = QGroupBox("Neuron Options")
        neuron_layout = QFormLayout(self.neuron_options)

        self.soma_radius = QSpinBox()
        self.soma_radius.setRange(1, 20)
        self.soma_radius.setValue(5)

        self.axon_length = QSpinBox()
        self.axon_length.setRange(10, 100)
        self.axon_length.setValue(50)

        self.axon_width = QSpinBox()
        self.axon_width.setRange(1, 10)
        self.axon_width.setValue(2)

        self.num_dendrites = QSpinBox()
        self.num_dendrites.setRange(1, 10)
        self.num_dendrites.setValue(5)

        self.dendrite_length = QSpinBox()
        self.dendrite_length.setRange(5, 50)
        self.dendrite_length.setValue(25)

        neuron_layout.addRow("Soma Radius:", self.soma_radius)
        neuron_layout.addRow("Axon Length:", self.axon_length)
        neuron_layout.addRow("Axon Width:", self.axon_width)
        neuron_layout.addRow("Number of Dendrites:", self.num_dendrites)
        neuron_layout.addRow("Dendrite Length:", self.dendrite_length)

        layout.addWidget(self.neuron_options)
        self.neuron_options.setVisible(False)

        self.tab_widget.addTab(structure_tab, "Cell Structure")

    def create_organelles_tab(self):
        organelles_tab = QWidget()
        layout = QVBoxLayout(organelles_tab)

        # Mitochondria group
        mito_group = QGroupBox("Mitochondria")
        mito_layout = QFormLayout(mito_group)

        self.mito_count = QSpinBox()
        self.mito_count.setRange(10, 200)
        self.mito_count.setValue(50)

        self.mito_size_min = QSpinBox()
        self.mito_size_max = QSpinBox()
        self.mito_size_min.setRange(1, 10)
        self.mito_size_max.setRange(1, 10)
        self.mito_size_min.setValue(3)
        self.mito_size_max.setValue(8)

        mito_layout.addRow("Count:", self.mito_count)
        mito_layout.addRow("Min Size:", self.mito_size_min)
        mito_layout.addRow("Max Size:", self.mito_size_max)

        layout.addWidget(mito_group)

        # Cytoskeleton group
        cyto_group = QGroupBox("Cytoskeleton")
        cyto_layout = QFormLayout(cyto_group)

        self.actin_density = QDoubleSpinBox()
        self.actin_density.setRange(0.01, 0.1)
        self.actin_density.setSingleStep(0.01)
        self.actin_density.setValue(0.05)

        self.microtubule_density = QDoubleSpinBox()
        self.microtubule_density.setRange(0.01, 0.1)
        self.microtubule_density.setSingleStep(0.01)
        self.microtubule_density.setValue(0.02)

        cyto_layout.addRow("Actin Density:", self.actin_density)
        cyto_layout.addRow("Microtubule Density:", self.microtubule_density)

        layout.addWidget(cyto_group)

        self.tab_widget.addTab(organelles_tab, "Organelles")

    def create_calcium_tab(self):
        calcium_tab = QWidget()
        layout = QVBoxLayout(calcium_tab)

        # Calcium Signaling group
        calcium_group = QGroupBox("Calcium Signaling")
        calcium_layout = QFormLayout(calcium_group)

        self.calcium_combo = QComboBox()
        self.calcium_combo.addItems(['None', 'Blip', 'Puff', 'Wave'])

        self.calcium_intensity = QDoubleSpinBox()
        self.calcium_intensity.setRange(0, 1)
        self.calcium_intensity.setSingleStep(0.1)
        self.calcium_intensity.setValue(0.5)

        self.calcium_duration = QSpinBox()
        self.calcium_duration.setRange(1, 100)
        self.calcium_duration.setValue(10)

        calcium_layout.addRow("Signal Type:", self.calcium_combo)
        calcium_layout.addRow("Signal Intensity:", self.calcium_intensity)
        calcium_layout.addRow("Signal Duration:", self.calcium_duration)

        layout.addWidget(calcium_group)

        self.tab_widget.addTab(calcium_tab, "Calcium Signaling")

    def toggle_neuron_options(self, cell_type):
        self.neuron_options.setVisible(cell_type == 'neuron')

    def requestSimulation(self):
        params = {
            'protein_diffusion': {
                'enabled': self.diffusion_checkbox.isChecked(),
                'coefficient': self.diffusion_coefficient.value()
            },
            'active_transport': {
                'enabled': self.transport_checkbox.isChecked(),
                'velocity': self.transport_velocity.value()
            },
            'cellular_structures': {
                'cell_membrane': self.cell_membrane_checkbox.isChecked(),
                'nucleus': self.nucleus_checkbox.isChecked(),
                'er': self.er_checkbox.isChecked(),
                'mitochondria': self.mitochondria_checkbox.isChecked(),
                'cytoskeleton': self.cytoskeleton_checkbox.isChecked()
            },
            'cell_type': self.cell_type_combo.currentText(),
            'cell_radius': self.cell_radius.value(),
            'membrane_thickness': self.membrane_thickness.value(),
            'nucleus_radius': self.nucleus_radius.value(),
            'er_density': self.er_density.value(),  # Add this line
            'pixel_size': (self.pixel_size_x.value(), self.pixel_size_y.value(), self.pixel_size_z.value()),
            'mitochondria': {
                'count': self.mito_count.value(),
                'size_range': (self.mito_size_min.value(), self.mito_size_max.value())
            },
            'cytoskeleton': {
                'actin_density': self.actin_density.value(),
                'microtubule_density': self.microtubule_density.value()
            },
            'calcium_signal': {
                'type': self.calcium_combo.currentText(),
                'intensity': self.calcium_intensity.value(),
                'duration': self.calcium_duration.value()
            },
            'neuron': {
                'soma_radius': self.soma_radius.value() if self.cell_type_combo.currentText() == 'neuron' else None,
                'axon_length': self.axon_length.value() if self.cell_type_combo.currentText() == 'neuron' else None,
                'axon_width': self.axon_width.value() if self.cell_type_combo.currentText() == 'neuron' else None,
                'num_dendrites': self.num_dendrites.value() if self.cell_type_combo.currentText() == 'neuron' else None,
                'dendrite_length': self.dendrite_length.value() if self.cell_type_combo.currentText() == 'neuron' else None
            }
        }
        self.simulationRequested.emit(params)


class BiologicalSimulationWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Biological Simulation")
        self.simulation_widget = BiologicalSimulationWidget()
        self.setCentralWidget(self.simulation_widget)
        self.resize(400, 300)  # Set an initial size for the window


class VolumeProcessor:
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def apply_threshold(self, data, threshold):
        try:
            return np.where(data > threshold, data, 0)
        except Exception as e:
            self.logger.error(f"Error in applying threshold: {str(e)}")
            return data

    def apply_gaussian_filter(self, data, sigma):
        try:
            from scipy.ndimage import gaussian_filter
            return gaussian_filter(data, sigma)
        except ImportError:
            self.logger.error("SciPy not installed. Cannot apply Gaussian filter.")
            return data
        except Exception as e:
            self.logger.error(f"Error in applying Gaussian filter: {str(e)}")
            return data

    def calculate_statistics(self, data):
        try:
            stats = {
                'mean': np.mean(data),
                'std': np.std(data),
                'min': np.min(data),
                'max': np.max(data)
            }
            self.logger.info(f"Statistics calculated: {stats}")
            return stats
        except Exception as e:
            self.logger.error(f"Error in calculating statistics: {str(e)}")
            return {}



class BlobAnalysisDialog(QDialog):
    def __init__(self, blob_analyzer, parent=None):
        super().__init__(parent)
        self.blob_analyzer = blob_analyzer
        self.setWindowTitle("Blob Analysis Results")
        self.setGeometry(100, 100, 800, 600)
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        # Add time point selection
        self.timeComboBox = QComboBox()
        self.timeComboBox.addItems([str(t) for t in self.blob_analyzer.time_points])
        self.timeComboBox.currentIndexChanged.connect(self.updateAnalysis)
        layout.addWidget(self.timeComboBox)

        self.tabWidget = QTabWidget()
        layout.addWidget(self.tabWidget)

        # Add a close button
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.close)
        layout.addWidget(closeButton)

        self.setLayout(layout)

        self.updateAnalysis()

    def updateAnalysis(self):
        self.tabWidget.clear()
        time_point = int(self.timeComboBox.currentText())

        self.addDistanceAnalysisTab(time_point)
        self.addDensityAnalysisTab(time_point)
        self.addColocalizationTab(time_point)
        self.addBlobSizeTab(time_point)
        self.addIntensityAnalysisTab(time_point)
        self.add3DVisualizationTab(time_point)
        self.addStatsTab(time_point)

    def addDistanceAnalysisTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        all_distances, within_channel_distances, between_channel_distances = self.blob_analyzer.calculate_nearest_neighbor_distances(time_point)

        plt.figure(figsize=(10, 6))
        plt.hist(all_distances, bins=50, alpha=0.5, label='All Blobs')
        for ch, distances in within_channel_distances.items():
            plt.hist(distances, bins=50, alpha=0.5, label=f'Channel {ch}')
        plt.xlabel('Nearest Neighbor Distance')
        plt.ylabel('Frequency')
        plt.legend()
        plt.title(f'Nearest Neighbor Distance Distribution (Time: {time_point})')
        canvas = FigureCanvas(plt.gcf())
        layout.addWidget(canvas)

        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "Distance Analysis")

        plt.close()

    def addDensityAnalysisTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        overall_density, channel_densities = self.blob_analyzer.calculate_blob_density((30, 100, 100), time_point)

        textEdit = QTextEdit()
        textEdit.setReadOnly(True)
        textEdit.append(f"Time Point: {time_point}")
        textEdit.append(f"Overall Blob Density: {overall_density:.6f} blobs/unit^3")
        for ch, density in channel_densities.items():
            textEdit.append(f"Channel {ch} Density: {density:.6f} blobs/unit^3")

        layout.addWidget(textEdit)
        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "Density Analysis")

    def addBlobSizeTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        blob_sizes = self.blob_analyzer.calculate_blob_sizes(time_point)

        plt.figure(figsize=(10, 6))
        for ch, sizes in blob_sizes.items():
            plt.hist(sizes, bins=50, alpha=0.5, label=f'Channel {ch}')
        plt.xlabel('Blob Size')
        plt.ylabel('Frequency')
        plt.legend()
        plt.title(f'Blob Size Distribution (Time: {time_point})')
        canvas = FigureCanvas(plt.gcf())
        layout.addWidget(canvas)

        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "Blob Size Analysis")

        plt.close()

    def addStatsTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        textEdit = QTextEdit()
        textEdit.setReadOnly(True)

        all_distances, within_channel_distances, _ = self.blob_analyzer.calculate_nearest_neighbor_distances(time_point)
        overall_density, channel_densities = self.blob_analyzer.calculate_blob_density((30, 100, 100), time_point)
        blob_sizes = self.blob_analyzer.calculate_blob_sizes(time_point)

        textEdit.append(f"Statistics for Time Point: {time_point}")
        textEdit.append("\nOverall Statistics:")
        time_blobs = self.blob_analyzer.blobs[self.blob_analyzer.blobs[:, 5] == time_point]
        textEdit.append(f"Total number of blobs: {len(time_blobs)}")
        textEdit.append(f"Overall blob density: {overall_density:.6f} blobs/unit^3")
        if len(all_distances) > 0:
            textEdit.append(f"Mean nearest neighbor distance: {np.mean(all_distances):.2f}")
            textEdit.append(f"Median nearest neighbor distance: {np.median(all_distances):.2f}")
        else:
            textEdit.append("Not enough blobs to calculate nearest neighbor distances.")

        for ch in self.blob_analyzer.channels:
            textEdit.append(f"\nChannel {ch} Statistics:")
            channel_blobs = time_blobs[time_blobs[:, 4] == ch]
            textEdit.append(f"Number of blobs: {len(channel_blobs)}")
            textEdit.append(f"Blob density: {channel_densities[ch]:.6f} blobs/unit^3")
            if len(blob_sizes[ch]) > 0:
                textEdit.append(f"Mean blob size: {np.mean(blob_sizes[ch]):.2f}")
                textEdit.append(f"Median blob size: {np.median(blob_sizes[ch]):.2f}")
            else:
                textEdit.append("No blobs detected in this channel.")
            if ch in within_channel_distances and len(within_channel_distances[ch]) > 0:
                textEdit.append(f"Mean nearest neighbor distance: {np.mean(within_channel_distances[ch]):.2f}")
                textEdit.append(f"Median nearest neighbor distance: {np.median(within_channel_distances[ch]):.2f}")
            else:
                textEdit.append("Not enough blobs to calculate nearest neighbor distances.")

        layout.addWidget(textEdit)
        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "Statistics")

    def addIntensityAnalysisTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        intensities = self.blob_analyzer.calculate_blob_intensities(time_point)
        sizes = self.blob_analyzer.calculate_blob_sizes(time_point)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        for ch in self.blob_analyzer.channels:
            ax1.hist(intensities[ch], bins=50, alpha=0.5, label=f'Channel {ch}')
        ax1.set_xlabel('Intensity')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Blob Intensity Distribution')
        ax1.legend()

        for ch in self.blob_analyzer.channels:
            ax2.scatter(sizes[ch], intensities[ch], alpha=0.5, label=f'Channel {ch}')
        ax2.set_xlabel('Blob Size')
        ax2.set_ylabel('Intensity')
        ax2.set_title('Blob Size vs Intensity')
        ax2.legend()

        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "Intensity Analysis")

        plt.close()

    def add3DVisualizationTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        time_blobs = self.blob_analyzer.blobs[self.blob_analyzer.blobs[:, 5] == time_point]

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        for ch in self.blob_analyzer.channels:
            channel_blobs = time_blobs[time_blobs[:, 4] == ch]
            ax.scatter(channel_blobs[:, 0], channel_blobs[:, 1], channel_blobs[:, 2],
                       s=channel_blobs[:, 3]*10, alpha=0.5, label=f'Channel {ch}')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'3D Blob Positions (Time: {time_point})')
        ax.legend()

        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "3D Visualization")

        plt.close()

    def addColocalizationTab(self, time_point):
        tab = QWidget()
        layout = QVBoxLayout()

        basic_coloc = self.blob_analyzer.calculate_colocalization(distance_threshold=5, time_point=time_point)
        advanced_coloc = self.blob_analyzer.calculate_advanced_colocalization(time_point)

        textEdit = QTextEdit()
        textEdit.setReadOnly(True)
        textEdit.append(f"Time Point: {time_point}")
        textEdit.append("\nBasic Colocalization:")
        for (ch1, ch2), coloc in basic_coloc.items():
            textEdit.append(f"Channels {ch1} and {ch2}: {coloc:.2%}")

        textEdit.append("\nAdvanced Colocalization:")
        for (ch1, ch2), results in advanced_coloc.items():
            textEdit.append(f"Channels {ch1} and {ch2}:")
            pearson = results['pearson']
            textEdit.append(f"  Pearson's coefficient: {pearson:.4f}" if not np.isnan(pearson) else "  Pearson's coefficient: N/A")
            textEdit.append(f"  Manders' M1: {results['manders_m1']:.4f}")
            textEdit.append(f"  Manders' M2: {results['manders_m2']:.4f}")

        layout.addWidget(textEdit)
        tab.setLayout(layout)
        self.tabWidget.addTab(tab, "Colocalization")


class BlobResultsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Blob Detection Results")
        self.layout = QVBoxLayout(self)
        self.table = QTableWidget()
        self.layout.addWidget(self.table)

    def update_results(self, blobs):
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(["X", "Y", "Z", "Size", "Channel", "Time"])
        self.table.setRowCount(len(blobs))

        for i, blob in enumerate(blobs):
            y, x, z, r, channel, t = blob
            self.table.setItem(i, 0, QTableWidgetItem(f"{x:.2f}"))
            self.table.setItem(i, 1, QTableWidgetItem(f"{y:.2f}"))
            self.table.setItem(i, 2, QTableWidgetItem(f"{z:.2f}"))
            self.table.setItem(i, 3, QTableWidgetItem(f"{r:.2f}"))
            self.table.setItem(i, 4, QTableWidgetItem(f"{int(channel)}"))
            self.table.setItem(i, 5, QTableWidgetItem(f"{int(t)}"))

class ROI3D(gl.GLMeshItem):
    sigRegionChanged = pyqtSignal(object)

    def __init__(self, size=(10, 10, 10), color=(1, 1, 1, 0.3)):
        verts, faces = self.create_cube(size)
        super().__init__(vertexes=verts, faces=faces, smooth=False, drawEdges=True, edgeColor=color)
        self.size = size
        self.setColor(color)

    @staticmethod
    def create_cube(size):
        x, y, z = size
        verts = np.array([
            [0, 0, 0], [x, 0, 0], [x, y, 0], [0, y, 0],
            [0, 0, z], [x, 0, z], [x, y, z], [0, y, z]
        ])
        faces = np.array([
            [0, 1, 2], [0, 2, 3], [0, 1, 4], [1, 4, 5],
            [1, 2, 5], [2, 5, 6], [2, 3, 6], [3, 6, 7],
            [3, 0, 7], [0, 4, 7], [4, 5, 6], [4, 6, 7]
        ])
        return verts, faces

    def setPosition(self, pos):
        self.resetTransform()
        self.translate(*pos)
        self.sigRegionChanged.emit(self)

class TimeSeriesDialog(QDialog):
    def __init__(self, blob_analyzer, parent=None):
        super().__init__(parent)
        self.blob_analyzer = blob_analyzer
        self.setWindowTitle("Time Series Analysis")
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        plot_widget = pg.PlotWidget()
        layout.addWidget(plot_widget)

        time_points = np.unique(self.blob_analyzer.blobs[:, 5])
        channels = np.unique(self.blob_analyzer.blobs[:, 4])

        for channel in channels:
            blob_counts = [np.sum((self.blob_analyzer.blobs[:, 5] == t) & (self.blob_analyzer.blobs[:, 4] == channel))
                           for t in time_points]
            plot_widget.plot(time_points, blob_counts, pen=(int(channel), len(channels)), name=f'Channel {int(channel)}')

        plot_widget.setLabel('left', "Number of Blobs")
        plot_widget.setLabel('bottom', "Time Point")
        plot_widget.addLegend()

        self.setLayout(layout)

######################################################################################################
######################################################################################################
'''                                   MAIN LIGHTSHEETVIEWER CLASS                                  '''
######################################################################################################
######################################################################################################

class LightsheetViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initLogging()
        self.volume_processor = VolumeProcessor()
        self.data = None
        self.lastPos = None
        self.initDataGenerator()
        self.initSimulator()
        self.initColors()
        self.setupChannelControlsWidget()
        self.initUI()
        self.generateData()
        self.initTimer()
        self.blob_results_dialog = BlobResultsDialog(self)
        self.showBlobResultsButton.setVisible(False)
        self.biological_simulation_window = None
        self.toggleDownsamplingControls()  # Set initial state of downsampling controls


    def initColors(self):
        self.channel_colors = [
            (1, 0, 0, 1),    # Red
            (0, 1, 0, 1),    # Green
            (0, 0, 1, 1),    # Blue
            (1, 1, 0, 1),    # Yellow
            (1, 0, 1, 1),    # Magenta
            (0, 1, 1, 1),    # Cyan
            (0.5, 0.5, 0.5, 1),  # Gray
            (1, 0.5, 0, 1),  # Orange
            (0.5, 0, 0.5, 1),    # Purple
            (0, 0.5, 0.5, 1),    # Teal
        ]

    def initLogging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

    def initTimer(self):
        self.playbackTimer = QTimer(self)
        self.playbackTimer.timeout.connect(self.advanceTimePoint)

    def initUI(self):
        self.setWindowTitle('Lightsheet Microscopy Viewer')
        self.setGeometry(100, 100, 1600, 900)
        self.setCentralWidget(None)

        self.createDocks()
        self.organizeDocks()
        self.createMenuBar()

        self.connectViewEvents()
        self.check_view_state()

    def createDocks(self):
        self.createViewDocks()
        self.createDataGenerationDock()
        self.create3DVisualizationDock()
        self.createVisualizationControlDock()
        self.createPlaybackControlDock()
        self.createBlobVisualizationDock()
        self.createBlobDetectionDock()

    def organizeDocks(self):
        self.setDockOptions(QMainWindow.AllowNestedDocks | QMainWindow.AllowTabbedDocks)

        # Stack XY, XZ, and YZ views vertically on the left
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockXY)
        self.splitDockWidget(self.dockXY, self.dockXZ, Qt.Vertical)
        self.splitDockWidget(self.dockXZ, self.dockYZ, Qt.Vertical)

        # Add 3D view to the right of the 2D views
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock3D)

        # Add control docks to the far right
        self.addDockWidget(Qt.RightDockWidgetArea, self.dockDataGeneration)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dockVisualizationControl)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dockPlaybackControl)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dockBlobDetection)

        # Set the 3D view and control docks side by side
        self.splitDockWidget(self.dock3D, self.dockDataGeneration, Qt.Horizontal)
        self.tabifyDockWidget(self.dockDataGeneration, self.dockVisualizationControl)
        self.tabifyDockWidget(self.dockVisualizationControl, self.dockPlaybackControl)
        self.tabifyDockWidget(self.dockPlaybackControl, self.dockBlobDetection)

        # Add blob visualization dock below the 3D view
        self.splitDockWidget(self.dock3D, self.dockBlobVisualization, Qt.Vertical)

        # Adjust dock sizes
        self.resizeDocks([self.dockXY, self.dockXZ, self.dockYZ], [200, 200, 200], Qt.Vertical)
        self.resizeDocks([self.dock3D, self.dockDataGeneration], [800, 300], Qt.Horizontal)

    def initSimulator(self):
        self.biological_simulator = BiologicalSimulator(size=(30, 100, 100), num_time_points=10)

    def initDataGenerator(self):
        self.data_generator = DataGenerator()


    def runBiologicalSimulation(self, params):
        try:
            self.data = None  # Reset data at the start of simulation
            soma_center = tuple(s // 2 for s in self.biological_simulator.size)
            cell_type = params['cell_type']
            cell_shape, cell_interior, cell_membrane = self.generateCellShape(params, soma_center, cell_type)

            self.generateCellularStructures(params, cell_shape, cell_interior, cell_membrane, soma_center)
            self.simulateProteinDynamics(params)
            self.simulateCalciumSignal(params)

            self.updateUIForNewData()
            self.updateViews()
            self.create3DVisualization()

        except Exception as e:
            self.handleSimulationError(e)

    def generateCellShape(self, params, soma_center, cell_type):
        if cell_type == 'neuron':
            return self.biological_simulator.generate_cell_shape(
                cell_type,
                self.biological_simulator.size,
                params['pixel_size'],
                membrane_thickness=params['membrane_thickness'],
                cell_radius=params['cell_radius'],
                soma_radius=params['neuron']['soma_radius'],
                axon_length=params['neuron']['axon_length'],
                axon_width=params['neuron']['axon_width'],
                num_dendrites=params['neuron']['num_dendrites'],
                dendrite_length=params['neuron']['dendrite_length']
            )
        else:
            return self.biological_simulator.generate_cell_shape(
                cell_type,
                self.biological_simulator.size,
                params['pixel_size'],
                membrane_thickness=params['membrane_thickness'],
                cell_radius=params['cell_radius']
            )

    def generateCellularStructures(self, params, cell_shape, cell_interior, cell_membrane, soma_center):
        cellular_structures = params.get('cellular_structures', {})

        self.biological_simulator.cell_shape = cell_shape
        self.biological_simulator.size = cell_shape.shape

        if cellular_structures.get('cell_membrane', False):
            self.generateCellMembrane(cell_membrane)

        if cellular_structures.get('nucleus', False):
            self.generateNucleus(params, cell_interior, soma_center)

        if cellular_structures.get('er', False):
            self.generateER(params, cell_shape, soma_center)

        if cellular_structures.get('mitochondria', False):
            self.generateMitochondria(params)

        if cellular_structures.get('cytoskeleton', False):
            self.generateCytoskeleton(params)

    def generateCellMembrane(self, cell_membrane):
        membrane_timeseries = np.repeat(cell_membrane[np.newaxis, np.newaxis, :, :, :], self.biological_simulator.num_time_points, axis=0)
        self.addDataChannel(membrane_timeseries, "Cell Membrane")

    def generateNucleus(self, params, cell_interior, soma_center):
        nucleus_data, nucleus_center = self.biological_simulator.generate_nucleus(
            cell_interior,
            soma_center,
            params['nucleus_radius'],
            pixel_size=params['pixel_size']
        )
        nucleus_timeseries = np.repeat(nucleus_data[np.newaxis, np.newaxis, :, :, :], self.biological_simulator.num_time_points, axis=0)
        self.addDataChannel(nucleus_timeseries, "Nucleus")

    def generateER(self, params, cell_shape, soma_center):
        er_data = self.biological_simulator.generate_er(
            cell_shape,
            soma_center,
            params['nucleus_radius'],
            params['er_density'],
            params['pixel_size']
        )
        er_timeseries = np.repeat(er_data[np.newaxis, np.newaxis, :, :, :], self.biological_simulator.num_time_points, axis=0)
        self.addDataChannel(er_timeseries, "ER")

    def generateMitochondria(self, params):
        mito_count = params['mitochondria'].get('count', 50)
        mito_size_range = params['mitochondria'].get('size_range', (3, 8))
        mitochondria = self.biological_simulator.generate_mitochondria(mito_count, mito_size_range)
        mitochondria_timeseries = np.repeat(mitochondria[np.newaxis, np.newaxis, :, :, :], self.biological_simulator.num_time_points, axis=0)
        self.addDataChannel(mitochondria_timeseries, "Mitochondria")

    def generateCytoskeleton(self, params):
        actin_density = params['cytoskeleton'].get('actin_density', 0.05)
        microtubule_density = params['cytoskeleton'].get('microtubule_density', 0.02)
        actin, microtubules = self.biological_simulator.generate_cytoskeleton(actin_density, microtubule_density)
        actin_timeseries = np.repeat(actin[np.newaxis, np.newaxis, :, :, :], self.biological_simulator.num_time_points, axis=0)
        microtubules_timeseries = np.repeat(microtubules[np.newaxis, np.newaxis, :, :, :], self.biological_simulator.num_time_points, axis=0)
        self.addDataChannel(actin_timeseries, "Actin")
        self.addDataChannel(microtubules_timeseries, "Microtubules")

    def simulateProteinDynamics(self, params):
        if params['protein_diffusion']['enabled']:
            self.simulateProteinDiffusion(params)
        if params['active_transport']['enabled']:
            self.simulateActiveTransport(params)

    def simulateProteinDiffusion(self, params):
        diffusion_coefficient = params['protein_diffusion']['coefficient']
        initial_concentration = np.zeros(self.biological_simulator.size)
        initial_concentration[self.biological_simulator.size[0]//2,
                              self.biological_simulator.size[1]//2,
                              self.biological_simulator.size[2]//2] = 1.0  # Point source in the center
        diffusion_data = self.biological_simulator.simulate_protein_diffusion(
            diffusion_coefficient,
            initial_concentration
        )
        self.addDataChannel(diffusion_data[:, np.newaxis, :, :, :], "Protein Diffusion")

    def simulateActiveTransport(self, params):
        velocity = params['active_transport'].get('velocity', (1, 1, 1))
        use_microtubules = params['active_transport'].get('use_microtubules', True)
        initial_cargo = np.zeros(self.biological_simulator.size)
        initial_cargo[self.biological_simulator.size[0]//2, self.biological_simulator.size[1]//2, self.biological_simulator.size[2]//2] = 1.0
        transport_data = self.biological_simulator.simulate_active_transport(velocity, initial_cargo, use_microtubules)
        self.addDataChannel(transport_data[:, np.newaxis, :, :, :], "Active Transport")

    def simulateCalciumSignal(self, params):
        if params['calcium_signal']['type'] != 'None':
            calcium_signal_type = params['calcium_signal']['type']
            calcium_signal = self.biological_simulator.simulate_calcium_signal(calcium_signal_type.lower(), {})
            self.addDataChannel(calcium_signal[:, np.newaxis, :, :, :], "Calcium Signal")

    def addDataChannel(self, channel_data, channel_name):
        if self.data is None:
            self.data = channel_data
        else:
            self.data = np.concatenate((self.data, channel_data), axis=1)
        self.logger.info(f"Added {channel_name} data. New data shape: {self.data.shape}")

    def handleSimulationError(self, e):
        self.logger.error(f"Error in biological simulation: {str(e)}")
        self.logger.error(f"Full exception: {traceback.format_exc()}")
        QMessageBox.critical(self, "Simulation Error", f"An error occurred during simulation: {str(e)}")

    def createViewDocks(self):
        # XY View
        self.dockXY = QDockWidget("XY View", self)
        self.imageViewXY = pg.ImageView()
        self.imageViewXY.ui.roiBtn.hide()
        self.imageViewXY.ui.menuBtn.hide()
        self.imageViewXY.setPredefinedGradient('viridis')
        self.imageViewXY.timeLine.sigPositionChanged.connect(self.updateMarkersFromSliders)
        self.dockXY.setWidget(self.imageViewXY)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockXY)

        # XZ View
        self.dockXZ = QDockWidget("XZ View", self)
        self.imageViewXZ = pg.ImageView()
        self.imageViewXZ.ui.roiBtn.hide()
        self.imageViewXZ.ui.menuBtn.hide()
        self.imageViewXZ.setPredefinedGradient('viridis')
        self.imageViewXZ.timeLine.sigPositionChanged.connect(self.updateMarkersFromSliders)
        self.dockXZ.setWidget(self.imageViewXZ)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockXZ)

        # YZ View
        self.dockYZ = QDockWidget("YZ View", self)
        self.imageViewYZ = pg.ImageView()
        self.imageViewYZ.ui.roiBtn.hide()
        self.imageViewYZ.ui.menuBtn.hide()
        self.imageViewYZ.setPredefinedGradient('viridis')
        self.imageViewYZ.timeLine.sigPositionChanged.connect(self.updateMarkersFromSliders)
        self.dockYZ.setWidget(self.imageViewYZ)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dockYZ)


        # Set size policies for the image views
        for view in [self.imageViewXY, self.imageViewXZ, self.imageViewYZ]:
            view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            view.setMinimumSize(200, 200)

    def updateMarkersFromSliders(self):
        self.updateSliceMarkers()
        self.create3DVisualization()  # This might be heavy, consider optimizing if performance is an issue

    def updateSliceMarkers(self):
        if not hasattr(self, 'data') or self.data is None:
            return

        # Remove old slice marker items
        for item in self.slice_marker_items:
            try:
                self.blobGLView.removeItem(item)
            except:
                pass

        for item in self.main_slice_marker_items:
            try:
                self.glView.removeItem(item)
            except:
                pass

        self.slice_marker_items.clear()
        self.main_slice_marker_items.clear()

        if self.showSliceMarkersCheck.isChecked():
            z, y, x = self.data.shape[2:]  # z, y, x instead of depth, height, width

            z_slice = int(self.imageViewXY.currentIndex)
            y_slice = int(self.imageViewXZ.currentIndex)
            x_slice = int(self.imageViewYZ.currentIndex)

            # Create new markers
            x_marker = gl.GLLinePlotItem(pos=np.array([[x_slice, 0, 0], [x_slice, y, 0], [x_slice, y, z], [x_slice, 0, z]]),
                                         color=(1, 0, 0, 1), width=2, mode='line_strip')
            y_marker = gl.GLLinePlotItem(pos=np.array([[0, y_slice, 0], [x, y_slice, 0], [x, y_slice, z], [0, y_slice, z]]),
                                         color=(0, 1, 0, 1), width=2, mode='line_strip')
            z_marker = gl.GLLinePlotItem(pos=np.array([[0, 0, z_slice], [x, 0, z_slice], [x, y, z_slice], [0, y, z_slice]]),
                                         color=(0, 0, 1, 1), width=2, mode='line_strip')

            # Add new markers
            self.glView.addItem(x_marker)
            self.glView.addItem(y_marker)
            self.glView.addItem(z_marker)

            # Create and add new markers for blob visualization view
            x_marker_vis = gl.GLLinePlotItem(pos=np.array([[x_slice, 0, 0], [x_slice, y, 0], [x_slice, y, z], [x_slice, 0, z]]),
                                             color=(1, 0, 0, 1), width=2, mode='line_strip')
            y_marker_vis = gl.GLLinePlotItem(pos=np.array([[0, y_slice, 0], [x, y_slice, 0], [x, y_slice, z], [0, y_slice, z]]),
                                             color=(0, 1, 0, 1), width=2, mode='line_strip')
            z_marker_vis = gl.GLLinePlotItem(pos=np.array([[0, 0, z_slice], [x, 0, z_slice], [x, y, z_slice], [0, y, z_slice]]),
                                             color=(0, 0, 1, 1), width=2, mode='line_strip')

            self.blobGLView.addItem(x_marker_vis)
            self.blobGLView.addItem(y_marker_vis)
            self.blobGLView.addItem(z_marker_vis)

            self.glView.update()
            self.blobGLView.update()

            self.slice_marker_items.extend([x_marker_vis, y_marker_vis, z_marker_vis])
            self.main_slice_marker_items.extend([x_marker, y_marker, z_marker])

            self.logger.debug(f"Slice positions - X: {x_slice}, Y: {y_slice}, Z: {z_slice}")

    def createDataGenerationDock(self):
        self.dockDataGeneration = QDockWidget("Data Generation", self)
        dataGenWidget = QWidget()
        layout = QVBoxLayout(dataGenWidget)

        layout.addWidget(QLabel("Number of Volumes:"))
        self.numVolumesSpinBox = QSpinBox()
        self.numVolumesSpinBox.setRange(1, 100)
        self.numVolumesSpinBox.setValue(10)
        layout.addWidget(self.numVolumesSpinBox)

        layout.addWidget(QLabel("Number of Blobs:"))
        self.numBlobsSpinBox = QSpinBox()
        self.numBlobsSpinBox.setRange(1, 100)
        self.numBlobsSpinBox.setValue(30)
        layout.addWidget(self.numBlobsSpinBox)

        layout.addWidget(QLabel("Noise Level:"))
        self.noiseLevelSpinBox = QDoubleSpinBox()
        self.noiseLevelSpinBox.setRange(0, 1)
        self.noiseLevelSpinBox.setSingleStep(0.01)
        self.noiseLevelSpinBox.setValue(0.02)
        layout.addWidget(self.noiseLevelSpinBox)

        layout.addWidget(QLabel("Movement Speed:"))
        self.movementSpeedSpinBox = QDoubleSpinBox()
        self.movementSpeedSpinBox.setRange(0, 10)
        self.movementSpeedSpinBox.setSingleStep(0.1)
        self.movementSpeedSpinBox.setValue(0.5)
        layout.addWidget(self.movementSpeedSpinBox)

        # Add checkbox for structured data
        self.structuredDataCheck = QCheckBox("Generate Structured Data")
        self.structuredDataCheck.setChecked(False)
        layout.addWidget(self.structuredDataCheck)

        # Add options for channel ranges
        self.channelRangeWidgets = []
        for i in range(2):  # Assuming 2 channels for now
            channelLayout = QGridLayout()
            channelLayout.addWidget(QLabel(f"Channel {i+1} Range:"), 0, 0, 1, 6)

            channelLayout.addWidget(QLabel("X min"), 1, 0)
            xMinSpin = QSpinBox()
            channelLayout.addWidget(xMinSpin, 1, 1)

            channelLayout.addWidget(QLabel("X max"), 1, 2)
            xMaxSpin = QSpinBox()
            channelLayout.addWidget(xMaxSpin, 1, 3)

            channelLayout.addWidget(QLabel("Y min"), 2, 0)
            yMinSpin = QSpinBox()
            channelLayout.addWidget(yMinSpin, 2, 1)

            channelLayout.addWidget(QLabel("Y max"), 2, 2)
            yMaxSpin = QSpinBox()
            channelLayout.addWidget(yMaxSpin, 2, 3)

            channelLayout.addWidget(QLabel("Z min"), 3, 0)
            zMinSpin = QSpinBox()
            channelLayout.addWidget(zMinSpin, 3, 1)

            channelLayout.addWidget(QLabel("Z max"), 3, 2)
            zMaxSpin = QSpinBox()
            channelLayout.addWidget(zMaxSpin, 3, 3)

            for spin in [xMinSpin, xMaxSpin, yMinSpin, yMaxSpin, zMinSpin, zMaxSpin]:
                spin.setRange(0, 100)

            self.channelRangeWidgets.append((xMinSpin, xMaxSpin, yMinSpin, yMaxSpin, zMinSpin, zMaxSpin))
            layout.addLayout(channelLayout)

        self.generateButton = QPushButton("Generate New Data")
        self.generateButton.clicked.connect(self.generateData)
        layout.addWidget(self.generateButton)

        self.saveButton = QPushButton("Save Data")
        self.saveButton.clicked.connect(self.saveData)
        layout.addWidget(self.saveButton)

        self.loadButton = QPushButton("Load Data")
        self.loadButton.clicked.connect(self.loadData)
        layout.addWidget(self.loadButton)

        layout.addStretch(1)  # This pushes everything up
        self.dockDataGeneration.setWidget(dataGenWidget)

    def createBlobDetectionDock(self):
        self.dockBlobDetection = QDockWidget("Blob Detection", self)
        blobDetectionWidget = QWidget()
        layout = QVBoxLayout(blobDetectionWidget)

        # Blob Detection controls
        layout.addWidget(QLabel("Blob Detection:"))

        blobLayout = QGridLayout()

        blobLayout.addWidget(QLabel("Max Sigma:"), 0, 0)
        self.maxSigmaSpinBox = QDoubleSpinBox()
        self.maxSigmaSpinBox.setRange(1, 100)
        self.maxSigmaSpinBox.setValue(30)
        blobLayout.addWidget(self.maxSigmaSpinBox, 0, 1)

        blobLayout.addWidget(QLabel("Num Sigma:"), 1, 0)
        self.numSigmaSpinBox = QSpinBox()
        self.numSigmaSpinBox.setRange(1, 20)
        self.numSigmaSpinBox.setValue(10)
        blobLayout.addWidget(self.numSigmaSpinBox, 1, 1)

        blobLayout.addWidget(QLabel("Threshold:"), 2, 0)
        self.blobThresholdSpinBox = QDoubleSpinBox()
        self.blobThresholdSpinBox.setRange(0, 1)
        self.blobThresholdSpinBox.setSingleStep(0.01)
        self.blobThresholdSpinBox.setValue(0.5)
        blobLayout.addWidget(self.blobThresholdSpinBox, 2, 1)

        layout.addLayout(blobLayout)

        # Add checkbox for showing all blobs
        self.showAllBlobsCheck = QCheckBox("Show All Blobs")
        self.showAllBlobsCheck.setChecked(False)
        self.showAllBlobsCheck.stateChanged.connect(self.updateBlobVisualization)
        layout.addWidget(self.showAllBlobsCheck)

        # Add Blob Detection button
        self.blobDetectionButton = QPushButton("Detect Blobs")
        self.blobDetectionButton.clicked.connect(self.detect_blobs)
        layout.addWidget(self.blobDetectionButton)

        # Add button to show/hide blob results
        self.showBlobResultsButton = QPushButton("Show Blob Results")
        self.showBlobResultsButton.clicked.connect(self.toggleBlobResults)
        layout.addWidget(self.showBlobResultsButton)

        # Add Blob Analysis button
        self.blobAnalysisButton = QPushButton("Analyze Blobs")
        self.blobAnalysisButton.clicked.connect(self.analyzeBlobsasdkjfb)
        layout.addWidget(self.blobAnalysisButton)

        # Time Series Analysis button
        self.timeSeriesButton = QPushButton("Time Series Analysis")
        self.timeSeriesButton.clicked.connect(self.showTimeSeriesAnalysis)
        layout.addWidget(self.timeSeriesButton)

        layout.addStretch(1)  # This pushes everything up
        self.dockBlobDetection.setWidget(blobDetectionWidget)


    def createVisualizationControlDock(self):
        self.dockVisualizationControl = QDockWidget("Visualization Control", self)
        visControlWidget = QWidget()
        layout = QVBoxLayout(visControlWidget)

        # Set up channel controls widget
        self.setupChannelControlsWidget()
        layout.addWidget(self.channelControlsWidget)

        # Add checkbox for downsampling
        self.downsamplingCheckBox = QCheckBox("Enable Downsampling")
        self.downsamplingCheckBox.setChecked(False)  # Default is off
        self.downsamplingCheckBox.stateChanged.connect(self.toggleDownsamplingControls)
        self.downsamplingCheckBox.stateChanged.connect(self.updateVisualization)
        layout.addWidget(self.downsamplingCheckBox)

        # Downsampling control (now in a separate layout)
        downsamplingLayout = QHBoxLayout()
        downsamplingLayout.addWidget(QLabel("Max Points:"))
        self.downsamplingSpinBox = QSpinBox()
        self.downsamplingSpinBox.setRange(1000, 1000000)
        self.downsamplingSpinBox.setSingleStep(1000)
        self.downsamplingSpinBox.setValue(100000)
        self.downsamplingSpinBox.valueChanged.connect(self.updateVisualization)
        downsamplingLayout.addWidget(self.downsamplingSpinBox)
        layout.addLayout(downsamplingLayout)

        layout.addWidget(QLabel("Threshold:"))
        self.thresholdSpinBox = QDoubleSpinBox()
        self.thresholdSpinBox.setRange(0, 1)
        self.thresholdSpinBox.setSingleStep(0.1)
        self.thresholdSpinBox.setValue(0.2)
        self.thresholdSpinBox.valueChanged.connect(self.updateThreshold)
        layout.addWidget(self.thresholdSpinBox)

        # Point size control
        layout.addWidget(QLabel("Point Size:"))
        self.pointSizeSpinBox = QDoubleSpinBox()
        self.pointSizeSpinBox.setRange(0.1, 10)
        self.pointSizeSpinBox.setSingleStep(0.1)
        self.pointSizeSpinBox.setValue(2)
        self.pointSizeSpinBox.valueChanged.connect(self.updateVisualization)
        layout.addWidget(self.pointSizeSpinBox)

        layout.addWidget(QLabel("3D Rendering Mode:"))
        self.renderModeCombo = QComboBox()
        self.renderModeCombo.addItems(["Points", "Surface", "Wireframe"])
        self.renderModeCombo.currentTextChanged.connect(self.updateRenderMode)
        layout.addWidget(self.renderModeCombo)

        layout.addWidget(QLabel("Color Map:"))
        self.colorMapCombo = QComboBox()
        self.colorMapCombo.addItems(["Viridis", "Plasma", "Inferno", "Magma", "Grayscale"])
        self.colorMapCombo.currentTextChanged.connect(self.updateColorMap)
        layout.addWidget(self.colorMapCombo)

        self.showSliceMarkersCheck = QCheckBox("Show Slice Markers")
        self.showSliceMarkersCheck.stateChanged.connect(self.toggleSliceMarkers)
        layout.addWidget(self.showSliceMarkersCheck)

        layout.addWidget(QLabel("Clip Plane:"))
        self.clipSlider = QSlider(Qt.Horizontal)
        self.clipSlider.setMinimum(0)
        self.clipSlider.setMaximum(100)
        self.clipSlider.setValue(100)
        self.clipSlider.valueChanged.connect(self.updateClipPlane)
        layout.addWidget(self.clipSlider)

        # Add checkbox for synchronizing views
        self.syncViewsCheck = QCheckBox("Synchronize 3D Views")
        self.syncViewsCheck.setChecked(False)
        layout.addWidget(self.syncViewsCheck)

        # Add button for auto-scaling
        self.autoScaleButton = QPushButton("Auto Scale Views")
        self.autoScaleButton.clicked.connect(self.autoScaleViews)
        layout.addWidget(self.autoScaleButton)

        # Add buttons for orienting views
        self.topDownButton = QPushButton("Top-Down View")
        self.topDownButton.clicked.connect(self.setTopDownView)
        layout.addWidget(self.topDownButton)

        self.sideViewButton = QPushButton("Side View (XZ)")
        self.sideViewButton.clicked.connect(self.setSideView)
        layout.addWidget(self.sideViewButton)

        self.frontViewButton = QPushButton("Front View (YZ)")
        self.frontViewButton.clicked.connect(self.setFrontView)
        layout.addWidget(self.frontViewButton)

        layout.addStretch(1)  # This pushes everything up
        self.dockVisualizationControl.setWidget(visControlWidget)


    def createBlobVisualizationDock(self):
        self.dockBlobVisualization = QDockWidget("Blob Visualization", self)
        self.blobGLView = gl.GLViewWidget()
        self.dockBlobVisualization.setWidget(self.blobGLView)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dockBlobVisualization)

        # Add a grid to the view
        gx = gl.GLGridItem()
        gx.rotate(90, 0, 1, 0)
        self.blobGLView.addItem(gx)
        gy = gl.GLGridItem()
        gy.rotate(90, 1, 0, 0)
        self.blobGLView.addItem(gy)
        gz = gl.GLGridItem()
        self.blobGLView.addItem(gz)

        # Initialize empty lists to store blob and slice marker items
        self.blob_items = []
        self.slice_marker_items = []

        # Connect mouse events
        self.blobGLView.mousePressEvent = self.on3DViewMousePress
        self.blobGLView.mouseReleaseEvent = self.on3DViewMouseRelease
        self.blobGLView.mouseMoveEvent = self.on3DViewMouseMove
        self.blobGLView.wheelEvent = self.on3DViewWheel
        self.logger.debug(f"Blob visualization dock created. blobGLView: {self.blobGLView}")

    def create3DVisualizationDock(self):
        # 3D View
        self.dock3D = QDockWidget("3D View", self)
        self.glView = gl.GLViewWidget()
        self.dock3D.setWidget(self.glView)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock3D)

        #set camera positon
        self.glView.setCameraPosition(distance=50, elevation=30, azimuth=45)
        self.glView.opts['backgroundColor'] = pg.mkColor(20, 20, 20)  # Dark background


        # Add a grid to the view
        gx = gl.GLGridItem()
        gx.rotate(90, 0, 1, 0)
        self.glView.addItem(gx)
        gy = gl.GLGridItem()
        gy.rotate(90, 1, 0, 0)
        self.glView.addItem(gy)
        gz = gl.GLGridItem()
        self.glView.addItem(gz)

        # Initialize empty lists to store data and slice marker items
        self.data_items = []
        self.main_slice_marker_items = []

        self.glView.opts['fov'] = 60
        self.glView.opts['elevation'] = 30
        self.glView.opts['azimuth'] = 45

        # Connect mouse events
        self.glView.mousePressEvent = self.on3DViewMousePress
        self.glView.mouseReleaseEvent = self.on3DViewMouseRelease
        self.glView.mouseMoveEvent = self.on3DViewMouseMove
        self.glView.wheelEvent = self.on3DViewWheel
        self.logger.debug(f"3D visualization dock created. glView: {self.glView}")


    def visualize_blobs(self, blobs):
        # Remove old blob visualizations
        for item in self.blob_items:
            self.blobGLView.removeItem(item)
        self.blob_items.clear()

        current_time = self.timeSlider.value()

        # Define colors for each channel (you can adjust these)
        channel_colors = self.channel_colors

        # Add new blob visualizations
        for blob in blobs:
            y, x, z, r, channel, t, intensity = blob
            mesh = gl.MeshData.sphere(rows=10, cols=20, radius=r)

            # Get color based on channel
            base_color = channel_colors[int(channel) % len(channel_colors)]

            # Adjust color based on whether it's a current or past blob
            if t == current_time:
                color = base_color
            else:
                color = tuple(c * 0.5 for c in base_color[:3]) + (base_color[3] * 0.5,)  # Dimmed color for past blobs

            ## Optionally, adjust color based on intensity
            #alpha = min(1.0, intensity / 255.0)  # Assuming intensity is in 0-255 range
            #color = (*base_color[:3], alpha)


            # Add to blob visualization view
            blob_item_vis = gl.GLMeshItem(meshdata=mesh, smooth=True, color=color, shader='shaded')
            blob_item_vis.translate(z, x, y)  # Swapped y and z
            self.blobGLView.addItem(blob_item_vis)
            self.blob_items.append(blob_item_vis)

        self.blobGLView.update()

    def advanceTimePoint(self):
        current_time = self.timeSlider.value()
        if current_time < self.timeSlider.maximum():
            self.timeSlider.setValue(current_time + 1)
        elif self.loopCheckBox.isChecked():
            self.timeSlider.setValue(0)
        else:
            self.playbackTimer.stop()
            self.playPauseButton.setText("Play")

    def updateTimePoint(self, value):
        if self.data is not None:
            self.currentTimePoint = value
            self.updateViews()
            self.create3DVisualization()
            self.updateBlobVisualization()
        else:
            self.logger.warning("No data available to update time point")

    def createPlaybackControlDock(self):
        self.dockPlaybackControl = QDockWidget("Playback Control", self)
        playbackControlWidget = QWidget()
        layout = QVBoxLayout(playbackControlWidget)

        layout.addWidget(QLabel("Time:"))
        self.timeSlider = QSlider(Qt.Horizontal)
        self.timeSlider.setMinimum(0)
        #self.timeSlider.setMaximum(self.data.shape[0] - 1)  # Assuming first dimension is time
        self.timeSlider.valueChanged.connect(self.updateTimePoint)
        layout.addWidget(self.timeSlider)

        playbackLayout = QHBoxLayout()
        self.playPauseButton = QPushButton("Play")
        self.playPauseButton.clicked.connect(self.togglePlayback)
        playbackLayout.addWidget(self.playPauseButton)

        self.speedLabel = QLabel("Speed:")
        playbackLayout.addWidget(self.speedLabel)

        self.speedSpinBox = QDoubleSpinBox()
        self.speedSpinBox.setRange(0.1, 10)
        self.speedSpinBox.setSingleStep(0.1)
        self.speedSpinBox.setValue(1)
        self.speedSpinBox.valueChanged.connect(self.updatePlaybackSpeed)
        playbackLayout.addWidget(self.speedSpinBox)

        self.loopCheckBox = QCheckBox("Loop")
        playbackLayout.addWidget(self.loopCheckBox)

        layout.addLayout(playbackLayout)

        layout.addStretch(1)  # This pushes everything up
        self.dockPlaybackControl.setWidget(playbackControlWidget)


    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Adjust dock sizes to maintain aspect ratio if needed
        width = self.width()
        left_width = width // 3
        right_width = width - left_width

        # Resize the 2D view docks
        self.resizeDocks([self.dockXY, self.dockXZ, self.dockYZ], [left_width] * 3, Qt.Horizontal)

        # Resize the 3D view and control docks
        control_width = right_width // 4  # Allocate 1/4 of right side to controls
        self.resizeDocks([self.dock3D], [right_width - control_width], Qt.Horizontal)
        self.resizeDocks([self.dockDataGeneration, self.dockVisualizationControl, self.dockPlaybackControl],
                         [control_width] * 3, Qt.Horizontal)

    def createMenuBar(self):
        menuBar = self.menuBar()

        fileMenu = menuBar.addMenu('&File')

        loadAction = fileMenu.addAction('&Load Data')
        loadAction.triggered.connect(self.loadData)

        saveAction = fileMenu.addAction('&Save Data')
        saveAction.triggered.connect(self.saveData)

        importAction = fileMenu.addAction('&Import Microscope Data')
        importAction.triggered.connect(self.import_microscope_data)

        quitAction = fileMenu.addAction('&Quit')
        quitAction.triggered.connect(self.close)

        viewMenu = menuBar.addMenu('&View')

        for dock in [self.dockXY, self.dockXZ, self.dockYZ, self.dock3D,
                     self.dockDataGeneration, self.dockVisualizationControl, self.dockPlaybackControl]:
            viewMenu.addAction(dock.toggleViewAction())


        analysisMenu = menuBar.addMenu('&Analysis')
        timeSeriesAction = analysisMenu.addAction('Time Series Analysis')
        timeSeriesAction.triggered.connect(self.showTimeSeriesAnalysis)

        # Add a new menu item for the Biological Simulation window
        viewMenu = menuBar.addMenu('&View')
        self.showBioSimAction = QAction('Biological Simulation', self, checkable=True)
        self.showBioSimAction.triggered.connect(self.toggleBiologicalSimulationWindow)
        viewMenu.addAction(self.showBioSimAction)


    def toggleBiologicalSimulationWindow(self, checked):
        if checked:
            if self.biological_simulation_window is None:
                self.biological_simulation_window = QMainWindow(self)
                self.biological_simulation_window.setWindowTitle("Biological Simulation")
                simulation_widget = BiologicalSimulationWidget()
                simulation_widget.simulationRequested.connect(self.runBiologicalSimulation)
                self.biological_simulation_window.setCentralWidget(simulation_widget)
            self.biological_simulation_window.show()
        else:
            if self.biological_simulation_window:
                self.biological_simulation_window.hide()


    def generateData(self):
        try:
            params = self.getDataGenerationParams()
            if params['structured_data']:
                channel_ranges = [
                    ((ch['x'][0], ch['x'][1]), (ch['y'][0], ch['y'][1]), (ch['z'][0], ch['z'][1]))
                    for ch in params['channel_ranges']
                ]
                self.data = self.generateStructuredData(params, channel_ranges)
            else:
                self.data = self.generateUnstructuredData(params)

            self.logDataInfo()
            self.updateUIForNewData()
            self.updateViews()
            self.create3DVisualization()
            self.autoScaleViews()
            self.logger.info("Data generated and visualized successfully")

        except Exception as e:
            self.handleDataGenerationError(e)

    def getDataGenerationParams(self):
        return {
            'num_volumes': self.numVolumesSpinBox.value(),
            'num_blobs': self.numBlobsSpinBox.value(),
            'noise_level': self.noiseLevelSpinBox.value(),
            'movement_speed': self.movementSpeedSpinBox.value(),
            'structured_data': self.structuredDataCheck.isChecked(),
            'size': (30, 100, 100),  # (z, y, x)
            'num_channels': 2  # You can make this configurable if needed
        }

    def generateStructuredData(self, params):
        channel_ranges = self.getChannelRanges(params)
        return self.data_generator.generate_structured_multi_channel_time_series(
            num_volumes=params['num_volumes'],
            num_channels=params['num_channels'],
            size=params['size'],
            num_blobs=params['num_blobs'],
            intensity_range=(0.8, 1.0),
            sigma_range=(2, 6),
            noise_level=params['noise_level'],
            movement_speed=params['movement_speed'],
            channel_ranges=channel_ranges
        )

    def generateUnstructuredData(self, params):
        return self.data_generator.generate_multi_channel_time_series(
            num_volumes=params['num_volumes'],
            num_channels=params['num_channels'],
            size=params['size'],
            num_blobs=params['num_blobs'],
            intensity_range=(0.8, 1.0),
            sigma_range=(2, 6),
            noise_level=params['noise_level'],
            movement_speed=params['movement_speed']
        )

    def getChannelRanges(self, params):
        channel_ranges = []
        for widgets in self.channelRangeWidgets:
            xMin, xMax, yMin, yMax, zMin, zMax = [w.value() for w in widgets]
            channel_ranges.append(((xMin, xMax), (yMin, yMax), (zMin, zMax)))

        # If we have fewer channel ranges than channels, add default ranges
        while len(channel_ranges) < params['num_channels']:
            channel_ranges.append(((0, params['size'][2]), (0, params['size'][1]), (0, params['size'][0])))

        return channel_ranges

    def logDataInfo(self):
        self.logger.info(f"Generated data shape: {self.data.shape}")
        self.logger.info(f"Data min: {self.data.min()}, max: {self.data.max()}, mean: {self.data.mean()}")
        num_blobs = self.numBlobsSpinBox.value() * self.data.shape[1] * self.data.shape[0]
        self.logger.info(f"Generated {num_blobs} blobs")

    def handleDataGenerationError(self, e):
        self.logger.error(f"Error in data generation: {str(e)}")
        self.logger.error(f"Error type: {type(e).__name__}")
        self.logger.error(f"Error args: {e.args}")
        self.logger.error(f"Traceback: {traceback.format_exc()}")
        QMessageBox.critical(self, "Error", f"Failed to generate data: {str(e)}")

    def updateViews(self):
        if self.data is None:
            self.logger.warning("No data to update views")
            return
        t = self.timeSlider.value()
        threshold = self.thresholdSpinBox.value()

        self.logger.debug(f"Updating views for time point {t}")
        self.logger.debug(f"Data shape: {self.data.shape}")

        num_channels = self.data.shape[1]
        z, y, x = self.data.shape[2:]  # z, y, x instead of depth, height, width

        # Prepare 3D RGB images for each view
        combined_xy = np.zeros((z, y, x, 3))
        combined_xz = np.zeros((x, z, y, 3))  # Changed order to x, z, y
        combined_yz = np.zeros((y, z, x, 3))  # Changed order to y, z, x

        # Define colors for each channel (RGB)
        channel_colors = [
            (1, 0, 0), (0, 1, 0), (0, 0, 1),  # Red, Green, Blue
            (1, 1, 0), (1, 0, 1), (0, 1, 1),  # Yellow, Magenta, Cyan
            (0.5, 0.5, 0.5), (1, 0.5, 0),     # Gray, Orange
        ]

        for c in range(num_channels):
            if c < len(self.channelControls) and self.channelControls[c][0].isChecked():
                opacity = self.channelControls[c][1].value() / 100
                channel_data = self.data[t, c]

                # Normalize channel data
                channel_data = (channel_data - channel_data.min()) / (channel_data.max() - channel_data.min() + 1e-8)

                channel_data[channel_data < threshold] = 0

                # Apply color to the channel
                color = channel_colors[c % len(channel_colors)]
                colored_data = channel_data[:, :, :, np.newaxis] * color

                combined_xy += colored_data * opacity
                combined_xz += np.transpose(colored_data, (2, 0, 1, 3)) * opacity  # Corrected transpose
                combined_yz += np.transpose(colored_data, (1, 0, 2, 3)) * opacity  # Corrected transpose

        # Clip values to [0, 1] range
        combined_xy = np.clip(combined_xy, 0, 1)
        combined_xz = np.clip(combined_xz, 0, 1)
        combined_yz = np.clip(combined_yz, 0, 1)

        self.imageViewXY.setImage(combined_xy, autoLevels=False, levels=[0, 1])
        self.imageViewXZ.setImage(combined_xz, autoLevels=False, levels=[0, 1])
        self.imageViewYZ.setImage(combined_yz, autoLevels=False, levels=[0, 1])

        # Update slice ranges
        self.imageViewXY.view.setRange(xRange=(0, x), yRange=(0, y), padding=0)
        self.imageViewXZ.view.setRange(xRange=(0, x), yRange=(0, z), padding=0)
        self.imageViewYZ.view.setRange(xRange=(0, y), yRange=(0, z), padding=0)

        self.updateSliceMarkers()

    def create3DVisualization(self):
        try:
            self.clear3DVisualization()
            t = self.timeSlider.value()
            threshold = self.thresholdSpinBox.value()
            render_mode = self.getRenderMode()

            for c in range(self.data.shape[1]):
                if self.isChannelVisible(c):
                    try:
                        self.visualizeChannel(c, t, threshold, render_mode)
                    except Exception as e:
                        self.logger.error(f"Error visualizing channel {c}: {str(e)}")
                        self.logger.error(f"Traceback: {traceback.format_exc()}")

            # Set the correct camera position and center
            z, y, x = self.data.shape[2:]
            center = QVector3D(x/2, y/2, z/2)
            self.glView.setCameraPosition(pos=center, distance=max(x, y, z)*1.5, elevation=30, azimuth=45)
            self.glView.opts['center'] = center

            self.glView.update()

        except Exception as e:
            self.handle3DVisualization(e)

    def clear3DVisualization(self):
        for item in self.data_items:
            self.glView.removeItem(item)
        self.data_items.clear()

    def isChannelVisible(self, channel):
        return (channel < len(self.channelControls) and
                self.channelControls[channel][0].isChecked())

    def visualizeChannel(self, channel, time, threshold, render_mode):
        volume_data = self.data[time, channel]
        opacity = self.getChannelOpacity(channel)
        color = self.getChannelColor(channel)

        # Ensure color is in the correct format (r, g, b) without alpha
        color = color[:3]  # Remove alpha if present

        self.renderMesh(volume_data, threshold, color, opacity, render_mode)


    def renderPoints(self, volume_data, threshold, color, opacity):
        z, y, x = np.where(volume_data > threshold)
        pos = np.column_stack((x, y, z))

        if len(pos) > 0:
            colors = np.tile(color, (len(pos), 1))
            colors[:, 3] = opacity * (volume_data[z, y, x] - volume_data.min()) / (volume_data.max() - volume_data.min())

            scatter = gl.GLScatterPlotItem(pos=pos, color=colors, size=self.scatterPointSizeSpinBox.value())
            self.glView.addItem(scatter)
            self.data_items.append(scatter)

    def hasDegenerateTriangles(self, verts, faces):
        # Check if any triangle has zero area
        v0 = verts[faces[:, 0]]
        v1 = verts[faces[:, 1]]
        v2 = verts[faces[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)
        areas = 0.5 * np.sqrt((cross**2).sum(axis=1))
        return np.any(areas < 1e-10)

    def renderMesh(self, volume_data, threshold, color, opacity, render_mode):
        if render_mode == 'points':
            z, y, x = np.where(volume_data > threshold)
            pos = np.column_stack((x, y, z))

            if len(pos) > 0:
                original_points = len(pos)
                # Downsample if enabled and there are too many points
                if self.downsamplingCheckBox.isChecked():
                    max_points = self.downsamplingSpinBox.value()
                    if len(pos) > max_points:
                        indices = np.random.choice(len(pos), max_points, replace=False)
                        pos = pos[indices]
                        self.logger.info(f"Downsampled from {original_points} to {len(pos)} points")
                    else:
                        self.logger.info(f"Rendering {len(pos)} points (no downsampling needed)")
                else:
                    self.logger.info(f"Rendering all {len(pos)} points (downsampling disabled)")

                size = np.ones(len(pos)) * self.pointSizeSpinBox.value()
                color_array = np.tile(color + (opacity,), (len(pos), 1))
                scatter = gl.GLScatterPlotItem(pos=pos, size=size, color=color_array)
                self.glView.addItem(scatter)
                self.data_items.append(scatter)
            else:
                self.logger.info("No points above threshold for point rendering.")
        else:
            # For surface and wireframe, use marching cubes
            verts, faces = self.marchingCubes(volume_data, threshold)
            if len(verts) > 0 and len(faces) > 0:
                # Filter out degenerate triangles
                valid_faces = self.filterDegenerateTriangles(verts, faces)
                if len(valid_faces) == 0:
                    self.logger.warning("No valid triangles after filtering degenerate ones.")
                    return

                # Render as surface or wireframe
                mesh = gl.GLMeshItem(vertexes=verts, faces=valid_faces, smooth=True, drawEdges=render_mode=='wireframe')
                mesh.setColor(color + (opacity,))

                if render_mode == 'wireframe':
                    # For wireframe, set the face color to be more transparent
                    mesh.setColor(color + (opacity*0.1,))
                    # Set edge color (this is done during initialization, not with setEdgeColor)
                    mesh = gl.GLMeshItem(vertexes=verts, faces=valid_faces, smooth=True, drawEdges=True,
                                         edgeColor=color + (opacity,))

                self.glView.addItem(mesh)
                self.data_items.append(mesh)
            else:
                self.logger.info("No vertices or faces generated from marching cubes.")


    def calculateNormals(self, verts, faces):
        norm = np.zeros(verts.shape, dtype=verts.dtype)
        tris = verts[faces]
        n = np.cross(tris[::,1] - tris[::,0], tris[::,2] - tris[::,0])
        for i in range(3):
            norm[faces[:,i]] += n
        norm = norm / np.linalg.norm(norm, axis=1)[:, np.newaxis]
        return norm

    def filterDegenerateTriangles(self, verts, faces):
        v0 = verts[faces[:, 0]]
        v1 = verts[faces[:, 1]]
        v2 = verts[faces[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)
        areas = 0.5 * np.sqrt((cross**2).sum(axis=1))
        valid_faces = faces[areas > 1e-6]  # Increased threshold
        self.logger.info(f"Filtered out {len(faces) - len(valid_faces)} degenerate triangles out of {len(faces)}")
        return valid_faces

    def createEdgesFromFaces(self, faces):
        edges = set()
        for face in faces:
            for i in range(3):
                edge = (face[i], face[(i+1)%3])
                if edge[0] > edge[1]:
                    edge = (edge[1], edge[0])
                edges.add(edge)
        return np.array(list(edges))


    def marchingCubes(self, volume_data, threshold):
        try:
            verts, faces, _, _ = measure.marching_cubes(volume_data, threshold)
            self.logger.info(f"Marching cubes generated {len(verts)} vertices and {len(faces)} faces")
            return verts, faces
        except RuntimeError as e:
            self.logger.warning(f"Marching cubes failed: {str(e)}")
            return np.array([]), np.array([])

    def logDataStatistics(self, volume_data):
        self.logger.info(f"Data shape: {volume_data.shape}")
        self.logger.info(f"Data range: {volume_data.min()} to {volume_data.max()}")
        self.logger.info(f"Data mean: {volume_data.mean()}")
        self.logger.info(f"Data std: {volume_data.std()}")
        self.logger.info(f"Number of non-zero voxels: {np.count_nonzero(volume_data)}")


    def getChannelOpacity(self, channel):
        return self.channelControls[channel][1].value() / 100

    def getChannelColor(self, channel):
        return self.channel_colors[channel % len(self.channel_colors)]

    def getRenderMode(self):
        return self.renderModeCombo.currentText().lower()

    def handle3DVisualizationError(self, e):
        self.logger.error(f"Error in 3D visualization: {str(e)}")
        self.logger.error(f"Traceback: {traceback.format_exc()}")
        QMessageBox.critical(self, "Error", f"Failed to create 3D visualization: {str(e)}")

    def visualize_data_distribution(self):
        t = self.timeSlider.value()
        num_channels, depth, height, width = self.data.shape[1:]

        for c in range(min(num_channels, 3)):
            volume_data = self.data[t, c]
            plt.figure(figsize=(15, 5))
            plt.subplot(131)
            plt.imshow(np.max(volume_data, axis=0))
            plt.title(f'Channel {c} - XY Max Projection')
            plt.subplot(132)
            plt.imshow(np.max(volume_data, axis=1))
            plt.title(f'Channel {c} - XZ Max Projection')
            plt.subplot(133)
            plt.imshow(np.max(volume_data, axis=2))
            plt.title(f'Channel {c} - YZ Max Projection')
            plt.colorbar()
            plt.show()

        plt.figure()
        plt.hist(self.data[t].ravel(), bins=100)
        plt.title('Data Histogram')
        plt.xlabel('Intensity')
        plt.ylabel('Frequency')
        plt.show()

    def updateVisualization(self):
        self.create3DVisualization()

    def updateScatterPointSize(self, value):
        self.updateViews()
        self.create3DVisualization()

    def updateChannelVisibility(self):
        self.updateViews()
        self.create3DVisualization()

    def updateChannelOpacity(self):
        self.updateViews()
        self.create3DVisualization()

    def getColorMap(self):
        cmap_name = self.colorMapCombo.currentText().lower()
        if cmap_name == "grayscale":
            return pg.ColorMap(pos=[0.0, 1.0], color=[(0, 0, 0, 255), (255, 255, 255, 255)])
        else:
            return pg.colormap.get(cmap_name)

    def triangulate_points(self, points):
        from scipy.spatial import Delaunay
        tri = Delaunay(points)
        return tri.simplices

    def updateRenderMode(self):
        self.create3DVisualization()

    def updateColorMap(self):
        self.create3DVisualization()

    def toggleSliceMarkers(self, state):
        if state == Qt.Checked:
            self.updateSliceMarkers()
        else:
            for attr in ['x_marker', 'y_marker', 'z_marker']:
                if hasattr(self, attr):
                    self.glView.removeItem(getattr(self, attr))
                    delattr(self, attr)
        self.glView.update()

    def updateClipPlane(self, value):
        try:
            clip_pos = (value / 100) * 30
            mask = self.scatter.pos[:, 2] <= clip_pos
            self.scatter.setData(pos=self.scatter.pos[mask],
                                 color=self.scatter.color[mask])
            self.logger.info(f"Clip plane updated to position {clip_pos}")
        except Exception as e:
            self.logger.error(f"Error updating clip plane: {str(e)}")


    def updateThreshold(self, value):
        # if value >= 1:
        #     self.logger.warning("Array compute error")
        #     return

        if self.data is not None:
            self.updateViews()
            self.create3DVisualization()
        else:
            self.logger.warning("No data available to update display threshold")

    def updateBlobThreshold(self, value):
       self.filter_blobs()

    def saveData(self):
        try:
            filename, _ = QFileDialog.getSaveFileName(self, "Save Data", "", "TIFF Files (*.tiff);;NumPy Files (*.npy)")
            if filename:
                if filename.endswith('.tiff'):
                    self.data_generator.save_tiff(filename)
                elif filename.endswith('.npy'):
                    self.data_generator.save_numpy(filename)
                self.logger.info(f"Data saved to {filename}")
        except Exception as e:
            self.logger.error(f"Error saving data: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to save data: {str(e)}")

    def loadData(self):
        try:
            filename, _ = QFileDialog.getOpenFileName(self, "Load Data", "", "TIFF Files (*.tiff);;NumPy Files (*.npy)")
            if filename:
                if filename.endswith('.tiff'):
                    self.data = self.data_generator.load_tiff(filename)
                elif filename.endswith('.npy'):
                    self.data = self.data_generator.load_numpy(filename)

                self.timeSlider.setMaximum(self.data.shape[0] - 1)
                self.updateViews()
                self.create3DVisualization()
                self.autoScaleViews()  # Add this line
                self.logger.info(f"Data loaded from {filename}")

                # Stop playback when loading new data
                if self.playbackTimer.isActive():
                    self.playbackTimer.stop()
                    self.playPauseButton.setText("Play")

        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")
            QMessageBox.critical(self, "Error", f"Failed to load data: {str(e)}")

    def toggleDownsamplingControls(self):
        enabled = self.downsamplingCheckBox.isChecked()
        self.downsamplingSpinBox.setEnabled(enabled)

    def togglePlayback(self):
        if not hasattr(self, 'playbackTimer'):
            self.playbackTimer = QTimer(self)
            self.playbackTimer.timeout.connect(self.advanceTimePoint)

        if self.playbackTimer.isActive():
            self.playbackTimer.stop()
            self.playPauseButton.setText("Play")
        else:
            self.playbackTimer.start(int(1000 / self.speedSpinBox.value()))
            self.playPauseButton.setText("Pause")

    def updatePlaybackSpeed(self, value):
        if hasattr(self, 'playbackTimer') and self.playbackTimer.isActive():
            self.playbackTimer.setInterval(int(1000 / value))

    def updateUIForNewData(self):
        if self.data is not None:
            self.timeSlider.setMaximum(self.data.shape[0] - 1)

            # Update channel controls
            num_channels = self.data.shape[1]
            channel_names = ['Membrane', 'Nucleus', 'ER', 'Mitochondria', 'Actin', 'Microtubules', 'Calcium']

            # Clear existing channel controls
            for i in reversed(range(self.channelControlsLayout.count())):
                self.channelControlsLayout.itemAt(i).widget().setParent(None)
            self.channelControls.clear()

            for i in range(num_channels):
                channelLayout = QHBoxLayout()
                channelLayout.addWidget(QLabel(f"{channel_names[i] if i < len(channel_names) else f'Channel {i+1}'}:"))
                visibilityCheck = QCheckBox("Visible")
                visibilityCheck.setChecked(True)
                visibilityCheck.stateChanged.connect(self.updateChannelVisibility)
                channelLayout.addWidget(visibilityCheck)
                opacitySlider = QSlider(Qt.Horizontal)
                opacitySlider.setRange(0, 100)
                opacitySlider.setValue(100)
                opacitySlider.valueChanged.connect(self.updateChannelOpacity)
                channelLayout.addWidget(opacitySlider)
                self.channelControls.append((visibilityCheck, opacitySlider))

                channelWidget = QWidget()
                channelWidget.setLayout(channelLayout)
                self.channelControlsLayout.addWidget(channelWidget)
                self.logger.debug(f"Created control for channel {i}")

            self.updateViews()
            self.create3DVisualization()
        else:
            self.logger.warning("No data available to update UI")

    def closeEvent(self, event):
        # Stop the playback timer
        if hasattr(self, 'playbackTimer'):
            self.playbackTimer.stop()

        # Perform any other cleanup operations here
        # For example, you might want to save application settings

        # Log the application closure
        self.logger.info("Application closed")

        # Accept the event to allow the window to close
        event.accept()

        # Make sure to close the biological simulation window when the main window is closed
        if self.biological_simulation_window:
            self.biological_simulation_window.close()

        # Call the base class implementation
        super().closeEvent(event)

    def safeClose(self):
        self.close()  # This will trigger the closeEvent

    def exportProcessedData(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export Processed Data", "", "TIFF Files (*.tiff);;NumPy Files (*.npy)")
        if filename:
            if filename.endswith('.tiff'):
                tifffile.imwrite(filename, self.data)
            elif filename.endswith('.npy'):
                np.save(filename, self.data)

    def exportScreenshot(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Save Screenshot", "", "PNG Files (*.png)")
        if filename:
            screen = QApplication.primaryScreen()
            screenshot = screen.grabWindow(self.winId())
            screenshot.save(filename, 'png')

    def exportVideo(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export Video", "", "MP4 Files (*.mp4)")
        if filename:
            import imageio
            writer = imageio.get_writer(filename, fps=10)
            for t in range(self.data.shape[0]):
                self.timeSlider.setValue(t)
                QApplication.processEvents()
                screen = QApplication.primaryScreen()
                screenshot = screen.grabWindow(self.winId())
                writer.append_data(screenshot.toImage().convertToFormat(QImage.Format_RGB888).bits().asstring(screenshot.width() * screenshot.height() * 3))
            writer.close()

    def detect_blobs(self):
        if self.data is None:
            self.logger.warning("No data available for blob detection")
            return

        max_time =  self.timeSlider.maximum()+1
        num_channels = self.data.shape[1]

        # Get blob detection parameters from UI
        max_sigma = self.maxSigmaSpinBox.value()
        num_sigma = self.numSigmaSpinBox.value()
        #threshold = self.blobThresholdSpinBox.value()

        all_blobs = []
        for channel in range(num_channels):
            for t in range(max_time):
                # Get the current 3D volume for this channel
                volume = self.data[t, channel]

                # Detect blobs
                blobs = blob_log(volume, max_sigma=max_sigma, num_sigma=num_sigma)

                # Calculate intensity for each blob
                for blob in blobs:
                    y, x, z, r = blob
                    y, x, z = int(y), int(x), int(z)
                    r = int(r)

                    # Define a small region around the blob center
                    y_min, y_max = max(0, y-r), min(volume.shape[0], y+r+1)
                    x_min, x_max = max(0, x-r), min(volume.shape[1], x+r+1)
                    z_min, z_max = max(0, z-r), min(volume.shape[2], z+r+1)

                    # Extract the region
                    region = volume[y_min:y_max, x_min:x_max, z_min:z_max]

                    # Calculate the intensity (you can use different measures here)
                    intensity = np.mean(region)  # or np.max(region), np.sum(region), etc.

                    # Add blob information including intensity
                    all_blobs.append([y, x, z, r, channel, t, intensity])

        # Convert to numpy array
        all_blobs = np.array(all_blobs)

        # Store the original blob data
        self.original_blobs = all_blobs

        # Filter blobs based on threshold -> creates self.all_detected_blobs
        self.filter_blobs()

        self.logger.info(f"Detected {len(all_blobs)} blobs across all channels")

        # Store all detected blobs
        #self.all_detected_blobs = all_blobs

        # Display results
        self.display_blob_results(self.all_detected_blobs)

        # Visualize blobs
        self.updateBlobVisualization()

        # Show the blob results button and update its text
        self.showBlobResultsButton.setVisible(True)
        self.showBlobResultsButton.setText("Show Blob Results")

        return all_blobs


    def filter_blobs(self):
        if not hasattr(self, 'original_blobs'):
            return

        threshold = self.blobThresholdSpinBox.value()
        self.all_detected_blobs = self.original_blobs[self.original_blobs[:, 6] > threshold]

        # Update visualization
        self.updateBlobVisualization()


    def display_blob_results(self, blobs):
        result_dialog = QDialog(self)
        result_dialog.setWindowTitle("Blob Detection Results")
        layout = QVBoxLayout(result_dialog)

        table = QTableWidget()
        table.setColumnCount(7)
        table.setHorizontalHeaderLabels(["X", "Y", "Z", "Size", "Channel", "Time", "Intensity"])
        table.setRowCount(len(blobs))

        for i, blob in enumerate(blobs):
            y, x, z, r, channel, t, intensity = blob
            table.setItem(i, 0, QTableWidgetItem(f"{x:.2f}"))
            table.setItem(i, 1, QTableWidgetItem(f"{y:.2f}"))
            table.setItem(i, 2, QTableWidgetItem(f"{z:.2f}"))
            table.setItem(i, 3, QTableWidgetItem(f"{r:.2f}"))
            table.setItem(i, 4, QTableWidgetItem(f"{int(channel)}"))
            table.setItem(i, 5, QTableWidgetItem(f"{int(t)}"))
            table.setItem(i, 6, QTableWidgetItem(f"{intensity:.2f}"))

        layout.addWidget(table)

        close_button = QPushButton("Close")
        close_button.clicked.connect(result_dialog.close)
        layout.addWidget(close_button)

        result_dialog.exec_()

    def updateBlobVisualization(self):
        if not hasattr(self, 'all_detected_blobs') or self.all_detected_blobs is None:
            return

        # Clear previous visualizations
        for item in self.blob_items:
            self.blobGLView.removeItem(item)
        self.blob_items.clear()

        current_time = self.timeSlider.value()

        if self.showAllBlobsCheck.isChecked():
            blobs_to_show = self.all_detected_blobs
        else:
            blobs_to_show = self.all_detected_blobs[self.all_detected_blobs[:, 5] == current_time]

        self.visualize_blobs(blobs_to_show)

    def toggleBlobResults(self):
        if self.blob_results_dialog.isVisible():
            self.blob_results_dialog.hide()
            self.showBlobResultsButton.setText("Show Blob Results")
        else:
            if hasattr(self, 'all_detected_blobs'):
                self.blob_results_dialog.update_results(self.all_detected_blobs)
            self.blob_results_dialog.show()
            self.showBlobResultsButton.setText("Hide Blob Results")

    def clearDetectedBlobs(self):
        if hasattr(self, 'all_detected_blobs'):
            del self.all_detected_blobs
        self.updateBlobVisualization()
        self.blob_results_dialog.hide()
        self.showBlobResultsButton.setVisible(False)

    def analyzeBlobsasdkjfb(self):
        if hasattr(self, 'all_detected_blobs') and self.all_detected_blobs is not None:
            blob_analyzer = BlobAnalyzer(self.all_detected_blobs)
            analysis_dialog = BlobAnalysisDialog(blob_analyzer, self)
            analysis_dialog.setWindowTitle("Blob Analysis Results")
            analysis_dialog.setGeometry(100, 100, 800, 600)
            analysis_dialog.show()
        else:
            QMessageBox.warning(self, "No Blobs Detected", "Please detect blobs before running analysis.")

    def showTimeSeriesAnalysis(self):
        if hasattr(self, 'all_detected_blobs') and self.all_detected_blobs is not None:
            dialog = TimeSeriesDialog(BlobAnalyzer(self.all_detected_blobs), self)
            dialog.exec_()
        else:
            QMessageBox.warning(self, "No Data", "Please detect blobs first.")

    def autoScaleViews(self):
        if self.data is None:
            return

        # Get the bounds of the data
        z, y, x = self.data.shape[2:]  # Assuming shape is (t, c, z, y, x)
        center = QVector3D(x/2, y/2, z/2)

        # Calculate the diagonal of the bounding box
        diagonal = np.sqrt(x**2 + y**2 + z**2)

        # Set the camera position for both views
        for view in [self.glView, self.blobGLView]:
            view.setCameraPosition(pos=center, distance=diagonal*1.2, elevation=30, azimuth=45)
            view.opts['center'] = center

        self.glView.update()
        self.blobGLView.update()



    def connectViewEvents(self):
        for view in [self.glView, self.blobGLView]:
            if view is not None:
                view.installEventFilter(self)
        self.logger.debug("View events connected")

    def eventFilter(self, source, event):
        if source in [self.glView, self.blobGLView]:
            if event.type() == QEvent.MouseButtonPress:
                self.on3DViewMousePress(event, source)
                return True
            elif event.type() == QEvent.MouseButtonRelease:
                self.on3DViewMouseRelease(event, source)
                return True
            elif event.type() == QEvent.MouseMove:
                self.on3DViewMouseMove(event, source)
                return True
            elif event.type() == QEvent.Wheel:
                self.on3DViewWheel(event, source)
                return True
        return super().eventFilter(source, event)

    @pyqtSlot(QEvent)
    def on3DViewMousePress(self, event, source):
        self.lastPos = event.pos()
        self.logger.debug(f"Mouse press event on {source}")

    @pyqtSlot(QEvent)
    def on3DViewMouseRelease(self, event, source):
        self.lastPos = None
        self.logger.debug(f"Mouse release event on {source}")

    @pyqtSlot(QEvent)
    def on3DViewMouseMove(self, event, source):
        if self.lastPos is None:
            return

        diff = event.pos() - self.lastPos
        self.lastPos = event.pos()

        self.logger.debug(f"Mouse move event on {source}")

        if event.buttons() == Qt.LeftButton:
            self.rotate3DViews(diff.x(), diff.y(), source)
        elif event.buttons() == Qt.MidButton:
            self.pan3DViews(diff.x(), diff.y(), source)

    @pyqtSlot(QEvent)
    def on3DViewWheel(self, event, source):
        delta = event.angleDelta().y()
        self.logger.debug(f"Wheel event on {source}")
        self.zoom3DViews(delta, source)

    def rotate3DViews(self, dx, dy, active_view):
        views_to_update = [active_view]
        if self.syncViewsCheck.isChecked():
            views_to_update = [self.glView, self.blobGLView]

        for view in views_to_update:
            if view is not None and hasattr(view, 'opts'):
                view.opts['elevation'] -= dy * 0.5
                view.opts['azimuth'] += dx * 0.5
                view.update()
            else:
                self.logger.error(f"Invalid view object: {view}")

    def pan3DViews(self, dx, dy, active_view):
        views_to_update = [active_view]
        if self.syncViewsCheck.isChecked():
            views_to_update = [self.glView, self.blobGLView]

        for view in views_to_update:
            if view is not None and hasattr(view, 'pan'):
                view.pan(dx, dy, 0, relative='view')
            else:
                self.logger.error(f"Invalid view object for panning: {view}")

    def zoom3DViews(self, delta, active_view):
        views_to_update = [active_view]
        if self.syncViewsCheck.isChecked():
            views_to_update = [self.glView, self.blobGLView]

        for view in views_to_update:
            if view is not None and hasattr(view, 'opts'):
                view.opts['fov'] *= 0.999**delta
                view.update()
            else:
                self.logger.error(f"Invalid view object for zooming: {view}")

    def check_view_state(self):
        self.logger.debug(f"glView state: {self.glView}, has opts: {hasattr(self.glView, 'opts')}")
        self.logger.debug(f"blobGLView state: {self.blobGLView}, has opts: {hasattr(self.blobGLView, 'opts')}")
        self.logger.debug(f"Sync checked: {self.syncViewsCheck.isChecked()}")


    def setTopDownView(self):
        for view in [self.glView, self.blobGLView]:
            view.setCameraPosition(elevation=90, azimuth=0)
            view.update()

    def setSideView(self):
        for view in [self.glView, self.blobGLView]:
            view.setCameraPosition(elevation=0, azimuth=0)
            view.update()

    def setFrontView(self):
        for view in [self.glView, self.blobGLView]:
            view.setCameraPosition(elevation=0, azimuth=90)
            view.update()

    def setupChannelControlsWidget(self):
        self.channelControlsWidget = QWidget()
        self.channelControlsLayout = QVBoxLayout(self.channelControlsWidget)
        self.channelControls = []

    def import_microscope_data(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Microscope Data", "",
                                                   "All Files (*);;SIS Files (*.sis);;TIFF Files (*.tif *.tiff)")
        if file_path:
            try:
                importer = ImporterFactory.get_importer(file_path)
                metadata = importer.read_metadata()
                data = importer.read_data()

                if data is None:
                    raise ValueError("Failed to read data from the file.")

                # Update the viewer with the new data
                self.data = data
                self.updateUIForNewData()
                self.updateViews()
                self.create3DVisualization()

                # Display metadata
                self.display_metadata(metadata)

            except Exception as e:
                QMessageBox.critical(self, "Import Error", str(e))

    def display_metadata(self, metadata):
        if metadata is None:
            QMessageBox.warning(self, "Metadata", "No metadata available for this file.")
            return

        # Create a new dialog to display metadata
        dialog = QDialog(self)
        dialog.setWindowTitle("Metadata")
        layout = QVBoxLayout()
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)

        # Format metadata for display
        metadata_text = "\n".join([f"{key}: {value}" for key, value in metadata.items()])
        text_edit.setText(metadata_text)

        layout.addWidget(text_edit)
        dialog.setLayout(layout)
        dialog.exec_()

##############################################################################

def main():
    app = QApplication(sys.argv)
    viewer = LightsheetViewer()
    viewer.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
