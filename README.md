# LightSheet Microscopy Viewer and Biological Simulator

A comprehensive application for visualizing, analyzing, and simulating microscopy data with advanced biological modeling capabilities.

![LightSheet Viewer](https://via.placeholder.com/900x500)

## Overview

This project combines a powerful 3D visualization tool for microscopy data with sophisticated simulation capabilities for cellular structures and dynamics. It allows researchers to:

1. Visualize multi-dimensional microscopy data (3D volume, multiple channels, time series)
2. Generate synthetic cellular data with realistic biological properties
3. Simulate complex cellular structures and processes
4. Perform advanced analysis on microscopy data

## Features

### Visualization

- **Interactive 3D Rendering**: View volumetric data with customizable rendering modes (points, surface, wireframe)
- **Multi-View Interface**: Simultaneous XY, XZ, YZ, and 3D views for comprehensive visualization
- **Channel Management**: Independently control visibility, opacity, and color of each channel
- **Time Series Playback**: Animate time series data with adjustable playback speed

### Biological Simulation

- **Cellular Structures**: Generate realistic cell membranes, nuclei, endoplasmic reticulum, mitochondria, and cytoskeleton
- **Multiple Cell Types**: Simulate spherical, neuronal, epithelial, and muscle cells with appropriate morphologies
- **Protein Dynamics**: Model protein diffusion and active transport along cytoskeletal elements
- **Calcium Signaling**: Simulate various calcium signaling events (waves, puffs, blips)
- **Multi-Cell Interactions**: Model cell-cell adhesion, signaling, and collective behavior

### Data Analysis

- **Blob Detection**: Identify and analyze discrete structures within volumetric data
- **Colocalization Analysis**: Quantify spatial relationships between different channels
- **Intensity Profiling**: Measure intensity distribution along user-defined paths
- **Statistical Analysis**: Calculate and visualize statistical properties of detected features

### Data Management

- **File Import/Export**: Support for TIFF and other microscopy file formats
- **Synthetic Data Generation**: Create test data with controllable noise, density, and dynamics
- **Metadata Handling**: View and manage microscopy metadata

## Installation

### Prerequisites

- Python 3.7+
- PyQt5
- PyQtGraph
- NumPy
- SciPy
- scikit-image

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/lightsheet-viewer.git
   cd lightsheet-viewer
   ```

2. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv venv
   
   # On Windows
   venv\Scripts\activate
   
   # On macOS/Linux
   source venv/bin/activate
   ```

3. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Starting the Application

Run the main viewer application:

```bash
python lightsheetViewer.py
```

### Visualizing Data

1. **Loading Data**:
   - Use "File > Load Data" to open existing microscopy files
   - Use "File > Import Microscope Data" for specialized microscopy formats
   - Generate synthetic data using the Data Generation panel

2. **Viewing Controls**:
   - Use mouse to rotate, pan, and zoom in 3D views
   - Adjust channel visibility and appearance in Channel Controls panel
   - Set visualization parameters in Visualization Control panel

3. **Time Series Playback**:
   - Use the Playback Control panel to navigate and animate time series data

### Running Biological Simulations

1. Open the Biological Simulation window from the Simulation menu
2. Configure the desired cell type and structures
3. Set parameters for protein dynamics, organelles, and signaling
4. Click "Simulate" to run the simulation and visualize the results

### Analysis

1. **Blob Detection**:
   - Use the Blob Detection panel to find and analyze structures
   - Adjust detection parameters for sensitivity and specificity
   - View statistical analysis of detected features

2. **Intensity Analysis**:
   - Use the Intensity Profile Tool to measure intensity along a line
   - View intensity distributions with the Raw Data Viewer

## Architecture

The application follows a modular architecture:

- **UI Layer**: Manages the graphical interface and user interactions
- **Visualization Layer**: Handles 3D rendering and data display
- **Data Management Layer**: Controls data loading, generation, and processing
- **Simulation Layer**: Implements biological models and simulations
- **Analysis Layer**: Provides analytical tools and measurements

## Development and Extension

The codebase is designed to be extensible. Key areas for customization:

- **New Cell Types**: Extend the `BiologicalSimulator` class with additional cell morphologies
- **Custom Analysis**: Add specialized analysis algorithms to detect features of interest
- **Visualization Enhancements**: Implement additional rendering modes or views
- **File Format Support**: Add importers for specialized microscopy formats

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Built with [PyQt5](https://www.riverbankcomputing.com/software/pyqt/) and [PyQtGraph](http://pyqtgraph.org/)
- Uses [scikit-image](https://scikit-image.org/) for image processing
- Cell biology simulation models based on current scientific literature


## Contact

For questions, feedback, or contributions, please contact george.dickinson@gmail.com
