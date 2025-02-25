# LightSheet Microscopy Viewer and Biological Simulator
# User Manual

![Application Interface](https://via.placeholder.com/900x500)

## Table of Contents

1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [User Interface Overview](#user-interface-overview)
4. [Visualization Features](#visualization-features)
5. [Data Management](#data-management)
6. [Biological Simulation](#biological-simulation)
   - [Cell Types and Structures](#cell-types-and-structures)
   - [Protein Dynamics](#protein-dynamics)
   - [Organelle Simulation](#organelle-simulation)
   - [Calcium Signaling](#calcium-signaling)
   - [Multi-Cell Interactions](#multi-cell-interactions)
7. [Analysis Tools](#analysis-tools)
8. [Technical Reference](#technical-reference)
9. [Troubleshooting](#troubleshooting)
10. [Appendix: Biological Background](#appendix-biological-background)

## Introduction

The LightSheet Microscopy Viewer and Biological Simulator is an integrated platform designed for researchers working with microscopy data and cellular biology simulations. This software bridges the gap between experimental data visualization and theoretical modeling, allowing users to:

- Visualize and analyze multi-dimensional microscopy data
- Generate synthetic data for testing and validation
- Simulate cellular structures and dynamics based on biological principles
- Perform quantitative analysis on both real and simulated data

The simulator implements biologically accurate models of cellular components and processes, making it a valuable tool for hypothesis testing and experimental design in cell biology research.

## Getting Started

### System Requirements

- Operating System: Windows 10/11, macOS 10.14+, or Linux
- RAM: Minimum 8GB (16GB+ recommended for large datasets)
- GPU: Dedicated graphics card recommended for 3D rendering
- Storage: 1GB for installation, additional space for data

### Installation

1. **Install Python**: Ensure Python 3.7+ is installed on your system.

2. **Download the software**:
   ```bash
   git clone https://github.com/yourusername/lightsheet-viewer.git
   cd lightsheet-viewer
   ```

3. **Set up a virtual environment** (recommended):
   ```bash
   python -m venv venv
   
   # On Windows
   venv\Scripts\activate
   
   # On macOS/Linux
   source venv/bin/activate
   ```

4. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

5. **Launch the application**:
   ```bash
   python lightsheetViewer.py
   ```

### Quick Start Guide

1. **Generate sample data**:
   - Go to the "Data Generation" panel
   - Set the number of volumes, channels, and blobs
   - Click "Generate New Data"

2. **Explore the visualization**:
   - Use the mouse to rotate the 3D view
   - Toggle channel visibility in the Channel Controls
   - Adjust the threshold to highlight structures of interest

3. **Run a simple simulation**:
   - Select "Biological Simulation" from the Simulation menu
   - Choose "spherical" as the cell type
   - Enable the cell membrane and nucleus
   - Click "Simulate"

## User Interface Overview

The application interface is organized into several dockable panels:

### Main Views

- **XY View**: Top-down slice view of the volume
- **XZ View**: Side-view slice
- **YZ View**: Front-view slice
- **3D View**: Interactive 3D rendering of the volume
- **Blob Visualization**: 3D view of detected blob features

### Control Panels

- **Channel Controls**: Manage visibility, color, and opacity of data channels
- **Data Generation**: Create synthetic data with customizable parameters
- **Visualization Control**: Adjust rendering options, thresholds, and display modes
- **Playback Control**: Navigate through time series data
- **Blob Detection**: Configure and run feature detection algorithms

## Visualization Features

### 3D Rendering Modes

- **Points**: Render volume as discrete points in 3D space
- **Surface**: Create a solid surface using marching cubes algorithm
- **Wireframe**: Display the contour mesh as a wireframe

To change the rendering mode:
1. Go to the Visualization Control panel
2. Select your preferred mode from the "3D Rendering Mode" dropdown
3. Adjust the threshold to control the density/detail level

### Multi-Channel Visualization

The software supports visualization of multi-channel data with independent control of each channel:

- **Visibility**: Toggle checkbox for each channel
- **Color**: Click the color button to choose a custom color
- **Opacity**: Adjust the slider to control transparency

For optimal multi-channel visualization:
- Use contrasting colors for different channels
- Reduce opacity of outer structures to see internal features
- Apply appropriate thresholds to highlight specific signal ranges

### Cross-Sectional Views

The XY, XZ, and YZ views show 2D slices through the volume. These views are synchronized and interactive:

- **Navigation**: Use the time slider under each view to move through slices
- **Synchronization**: Slice positions are marked in the 3D view with colored planes
- **Intensity Profile**: Draw a line on any view to generate an intensity profile graph

## Data Management

### Loading Microscopy Data

The application supports several microscopy file formats:

1. Select "File > Import Microscope Data"
2. Navigate to your data file
3. Configure data properties (z-slices, channels, pixel size)
4. Data will load and display in all views

### Generating Synthetic Data

Create test data with controllable parameters:

1. Go to the Data Generation panel
2. Configure:
   - Number of volumes (time points)
   - Number of channels
   - Number of blobs (features)
   - Noise level and movement speed
3. For structured data (specific spatial distributions):
   - Enable "Generate Structured Data"
   - Set x, y, z ranges for each channel
4. Click "Generate New Data"

### Processing Options

Basic image processing tools are available under "Raw Data > Image Processing":

- **Gaussian Blur**: Apply smoothing to reduce noise
- **Median Filter**: Remove salt-and-pepper noise while preserving edges

## Biological Simulation

The biological simulation capabilities allow modeling of cellular structures and processes with realistic parameters based on cell biology principles.

### Cell Types and Structures

#### Available Cell Types

1. **Spherical Cell**
   - Simple spherical morphology
   - Uniform membrane and cytoplasm
   - Parameters: cell radius, membrane thickness

2. **Neuron**
   - Soma (cell body) with extending processes
   - Asymmetric morphology with axon and dendrites
   - Parameters: soma radius, axon length/width, dendrite number/length

3. **Epithelial Cell**
   - Polarized cell with distinct apical and basal surfaces
   - Structured lateral membranes
   - Parameters: height, membrane thickness

4. **Muscle Cell**
   - Elongated cylindrical morphology
   - Tapered ends
   - Parameters: length, diameter, tapering factor

#### Cellular Structures

The following subcellular components can be simulated:

1. **Cell Membrane**
   - Defines the cell boundary
   - Can be customized with varying thickness
   - Biologically accurate: Implemented as a phospholipid bilayer-like boundary

2. **Nucleus**
   - Central organelle containing genetic material
   - Typically located at cell center
   - Enclosed by nuclear envelope
   - Biologically accurate: Proper scaling relative to cell size (nucleocytoplasmic ratio)

3. **Endoplasmic Reticulum (ER)**
   - Network of membrane tubules and sheets
   - Higher density near nucleus
   - Connected to nuclear envelope
   - Biologically accurate: Forms a continuous network with appropriate distribution

4. **Mitochondria**
   - Multiple discrete organelles
   - Excluded from nucleus
   - Variable size and distribution
   - Biologically accurate: Appropriate size range and density throughout cytoplasm

5. **Cytoskeleton**
   - Actin filaments and microtubules
   - Organized network throughout cytoplasm
   - Microtubules radiate from cell center
   - Biologically accurate: Models proper network architecture and density

### Protein Dynamics

#### Diffusion Simulation

Simulate the movement of proteins through the cytoplasm via diffusion:

1. Enable "Protein Diffusion" in the Protein Dynamics tab
2. Set diffusion coefficient (typical range: 0.1-10 μm²/s)
   - Small proteins: 5-10 μm²/s
   - Medium proteins: 1-5 μm²/s
   - Large complexes: 0.1-1 μm²/s
3. The simulation applies a physically accurate diffusion equation using Gaussian kernels

#### Active Transport

Model directed movement of proteins along cytoskeletal elements:

1. Enable "Active Transport" in the Protein Dynamics tab
2. Set transport velocity (typical range: 0.1-5 μm/s)
   - Kinesin motors: ~0.5-1 μm/s
   - Dynein motors: ~0.5-1 μm/s
   - Myosin motors: ~0.1-5 μm/s depending on type
3. Choose between microtubule-based or actin-based transport
4. Transport follows the cytoskeletal network with realistic velocities

### Organelle Simulation

#### Mitochondria

Customize mitochondrial properties:

1. Set mitochondria count (typical range: 20-2000 depending on cell type)
2. Configure size range (typical values: 0.5-10 μm)
3. Distribution follows biologically accurate patterns:
   - Excluded from nucleus
   - Higher density in high-energy regions
   - Dynamic positioning

#### Cytoskeleton

Configure the cytoskeletal network:

1. Set actin density (typical range: 0.01-0.1)
2. Set microtubule density (typical range: 0.005-0.05)
3. The simulation generates:
   - Actin filaments as a branched network
   - Microtubules radiating from cell center
   - Appropriate connections to cell membrane and organelles

### Calcium Signaling

Simulate various types of calcium signaling events:

1. **Calcium Blips**
   - Localized, brief calcium release events
   - Represent single IP3R channel openings
   - Duration: 10-100 ms, amplitude: low

2. **Calcium Puffs**
   - Larger events involving multiple channels
   - Spatially concentrated
   - Duration: 100-500 ms, amplitude: medium

3. **Calcium Waves**
   - Propagating elevations in calcium
   - Travel across the cell
   - Duration: 1-10 s, amplitude: high, speed: 5-100 μm/s

Configure these parameters in the Calcium Signaling tab and observe spatiotemporal dynamics in the visualization.

### Multi-Cell Interactions

The Enhanced Biological Simulator allows modeling interactions between multiple cells:

1. **Cell-Cell Adhesion**
   - Simulates adhesion forces between cells
   - Controls cell aggregation and tissue formation
   - Strength parameter determines adhesion force

2. **Cell-Cell Signaling**
   - Models exchange of signaling molecules between cells
   - Diffusion through extracellular space
   - Receptor-mediated responses

To run a multi-cell simulation:
1. Select "Enhanced Biological Simulation" from the Simulation menu
2. Set the environment size and number of cells
3. Configure cell properties and protein dynamics
4. Run the simulation and observe emergent behaviors

## Analysis Tools

### Blob Detection

Detect and analyze discrete structures in volumetric data:

1. Configure detection parameters:
   - Max Sigma: Controls the maximum size of features to detect
   - Num Sigma: Number of intermediate scales to analyze
   - Threshold: Minimum intensity for feature detection

2. Click "Detect Blobs" to run the algorithm

3. Results include:
   - 3D visualization of detected blobs
   - Table of blob positions, sizes, and intensities
   - Statistical analysis of blob distributions

### Colocalization Analysis

Quantify spatial relationships between different channels:

1. Run blob detection on multi-channel data
2. Open "Blob Analysis" to see colocalization metrics:
   - Pearson's correlation coefficient
   - Manders' overlap coefficients
   - Distance-based colocalization measures
   - Channel-specific statistics

### Intensity Profile Analysis

Measure intensity distribution along a user-defined line:

1. Enable the "Intensity Profile Tool"
2. Draw a line on any 2D view
3. A graph will display showing intensity values along the line
4. Multi-channel data will show separate lines for each channel

## Technical Reference

### Rendering Performance Tips

For optimal performance with large datasets:

1. Enable downsampling for 3D rendering
2. Use points mode instead of surface rendering
3. Increase threshold to reduce the number of rendered points
4. Close unused panels to conserve memory

### Custom Data Import

For specialized microscopy formats:

1. Create a custom importer class that inherits from `MicroscopeDataImporter`
2. Implement `read_metadata()` and `read_data()` methods
3. Register your importer in `ImporterFactory.get_importer()`

### Simulation Parameters

Key parameters for biological accuracy:

| Structure | Parameter | Typical Range | Notes |
|-----------|-----------|---------------|-------|
| Cell | Radius | 5-50 μm | Depends on cell type |
| Nucleus | Radius | 2-15 μm | Usually 10-20% of cell volume |
| ER | Density | 0.05-0.2 | Higher near nucleus |
| Mitochondria | Count | 20-2000 | Cell-type dependent |
| Cytoskeleton | Actin density | 0.01-0.1 | Higher near membrane |
| Diffusion | Coefficient | 0.1-10 μm²/s | Molecule-size dependent |

## Troubleshooting

### Common Issues

1. **Application crashes when loading large files**
   - Try enabling downsampling
   - Increase system virtual memory
   - Use 64-bit Python installation

2. **3D rendering appears empty**
   - Adjust threshold values
   - Check channel visibility
   - Reset camera position with "Auto Scale Views"

3. **Simulation fails to initialize**
   - Check for valid parameter combinations
   - Ensure cell size fits within environment
   - Review log for specific error messages

4. **Poor rendering performance**
   - Reduce point size or switch to wireframe mode
   - Enable downsampling
   - Close unused application panels

## Appendix: Biological Background

### Cellular Structures and Functions

#### Cell Membrane
The plasma membrane forms a semi-permeable barrier around the cell, composed of a phospholipid bilayer embedded with proteins. It regulates the exchange of materials between the cell and its environment and maintains cell integrity. In the simulation, the membrane is modeled as a thin boundary with appropriate thickness (5-10 nm) and selective permeability to diffusing molecules.

#### Nucleus
The nucleus houses the cell's genetic material (DNA) and is surrounded by a double-membrane called the nuclear envelope. It controls gene expression and cellular activities through the synthesis of RNA. The simulator models the nucleus with proper scaling relative to cell size, maintains appropriate nucleocytoplasmic ratios (typically 0.1-0.3), and implements functional connections to the endoplasmic reticulum.

#### Endoplasmic Reticulum
The ER is a network of interconnected tubules and flattened sacs involved in protein synthesis, lipid metabolism, and calcium storage. The rough ER (studded with ribosomes) specializes in protein synthesis, while the smooth ER handles lipid synthesis and calcium regulation. The simulation creates a continuous ER network with higher density near the nucleus and connections to the nuclear envelope, as observed in actual cells.

#### Mitochondria
Mitochondria are double-membrane organelles responsible for cellular respiration and ATP production. They also play roles in calcium homeostasis, cell signaling, and apoptosis. The simulation models mitochondria as discrete organelles with appropriate size (0.5-10 μm) and distribution, typically excluding them from the nucleus and positioning them near areas with high energy demands.

#### Cytoskeleton
The cytoskeleton provides structural support, enables cell movement, and facilitates intracellular transport. It consists of three main components: microfilaments (actin), intermediate filaments, and microtubules. The simulator generates realistic cytoskeletal networks with appropriate architecture: actin forms a branched network concentrated near the membrane, while microtubules radiate from the microtubule organizing center near the nucleus.

### Protein Dynamics in Living Cells

#### Diffusion Principles
Diffusion is the random movement of molecules from areas of high concentration to low concentration. In cells, cytoplasmic diffusion is influenced by:

- Molecular size: Larger proteins diffuse more slowly
- Cytoplasmic viscosity: 3-5 times higher than water
- Molecular crowding: Reduces effective diffusion rates
- Temperature: Higher temperatures increase diffusion

The simulator implements physically accurate diffusion using the diffusion equation and appropriate coefficients based on molecular size.

#### Active Transport Mechanisms
Active transport involves the directed movement of molecules along cytoskeletal tracks, driven by motor proteins that consume ATP:

- Kinesins: Move cargo toward microtubule plus-ends (typically cell periphery)
- Dyneins: Move cargo toward microtubule minus-ends (typically cell center)
- Myosins: Move cargo along actin filaments

The simulation models realistic velocities, directionality, and path selection along the cytoskeletal network.

### Calcium Signaling
Calcium serves as a universal second messenger in cells, controlling processes ranging from muscle contraction to neurotransmitter release. Calcium signals vary in spatial spread, duration, and amplitude:

- **Calcium Blips**: Single-channel openings of IP3 receptors, typically lasting 10-100 ms
- **Calcium Puffs**: Coordinated opening of multiple channels, lasting 100-500 ms
- **Calcium Waves**: Regenerative signals that propagate across cells at 5-100 μm/s

The simulation implements these events with biologically accurate spatial and temporal parameters, modeling calcium release, diffusion, and re-uptake mechanisms.

### Cell-Cell Interactions in Tissues

#### Adhesion Mechanisms
Cells adhere to each other through specialized junction proteins:

- Tight junctions: Seal adjacent cells together
- Adherens junctions: Connect adjacent cells via cadherin proteins
- Desmosomes: Provide strong adhesion in tissues under mechanical stress
- Gap junctions: Allow direct communication between adjacent cells

The simulator models the combined effects of these mechanisms through parameterized adhesion forces.

#### Intercellular Signaling
Cells communicate through several mechanisms:

- Paracrine signaling: Secretion of factors that diffuse to nearby cells
- Juxtacrine signaling: Direct contact between adjacent cells
- Autocrine signaling: Cells responding to factors they themselves secrete

The enhanced biological simulator incorporates these signaling mechanisms with appropriate diffusion coefficients, receptor dynamics, and cellular responses.

---
