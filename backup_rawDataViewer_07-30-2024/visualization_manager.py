#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:37:34 2024

@author: george
"""

import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from PyQt5.QtCore import Qt, QLineF, QPointF
from skimage import measure
from PyQt5.QtWidgets import QApplication, QGraphicsLineItem
from PyQt5.QtGui import QPen, QColor
import logging

class VisualizationManager:
    def __init__(self, parent):
        self.parent = parent
        self.logger = parent.logger
        self.glView = None
        self.blobGLView = None
        self.data_items = []
        self.blob_items = []
        self.main_slice_marker_items = []
        self.slice_marker_items = []
        self.profile_roi = None  # Initialize profile_roi attribute
        self.profile_window = None  # Initialize profile_window attribute


    def create_3d_view(self):
        self.glView = gl.GLViewWidget()
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

        self.glView.opts['fov'] = 60
        self.glView.opts['elevation'] = 30
        self.glView.opts['azimuth'] = 45

        return self.glView

    def create_blob_view(self):
        self.blobGLView = gl.GLViewWidget()

        # Add a grid to the view
        gx = gl.GLGridItem()
        gx.rotate(90, 0, 1, 0)
        self.blobGLView.addItem(gx)
        gy = gl.GLGridItem()
        gy.rotate(90, 1, 0, 0)
        self.blobGLView.addItem(gy)
        gz = gl.GLGridItem()
        self.blobGLView.addItem(gz)

        return self.blobGLView

    def clear_3d_visualization(self):
        for item in self.data_items:
            self.glView.removeItem(item)
        self.data_items.clear()

    def update_3d_visualization(self, data, time_point, threshold, render_mode):
        self.clear_3d_visualization()

        for c in range(data.shape[0]):  # Iterate over channels
            if self.parent.isChannelVisible(c):
                try:
                    self.visualize_channel(data[c], c, time_point, threshold, render_mode)
                except Exception as e:
                    self.parent.logger.error(f"Error visualizing channel {c}: {str(e)}")

        self.glView.update()

    def visualize_channel(self, volume_data, channel, time, threshold, render_mode):
        opacity = self.parent.getChannelOpacity(channel)
        color = self.parent.getChannelColor(channel)
        color = color[:3]  # Remove alpha if present

        if render_mode == 'points':
            self.render_points(volume_data, threshold, color, opacity)
        else:
            self.render_mesh(volume_data, threshold, color, opacity, render_mode)


    def render_points(self, volume_data, threshold, color, opacity):
        z, y, x = np.where(volume_data > threshold)
        pos = np.column_stack((x, y, z))

        if len(pos) > 0:
            # Create colors array with 4 channels (RGBA)
            colors = np.zeros((len(pos), 4))
            colors[:, :3] = color  # Set RGB values

            # Calculate alpha based on intensity
            intensity = (volume_data[z, y, x] - volume_data.min()) / (volume_data.max() - volume_data.min())
            colors[:, 3] = opacity * intensity  # Set alpha values

            scatter = gl.GLScatterPlotItem(pos=pos, color=colors, size=self.parent.pointSizeSpinBox.value())
            self.glView.addItem(scatter)
            self.data_items.append(scatter)

    def render_mesh(self, volume_data, threshold, color, opacity, render_mode):
        verts, faces = self.marching_cubes(volume_data, threshold)
        if len(verts) > 0 and len(faces) > 0:
            # Reorder the vertices to match the orientation of the points rendering
            verts = verts[:, [2, 1, 0]]  # Change from (y, x, z) to (z, y, x)

            mesh = gl.GLMeshItem(vertexes=verts, faces=faces, smooth=True, drawEdges=render_mode=='wireframe')
            mesh.setColor(color + (opacity,))

            if render_mode == 'wireframe':
                mesh.setColor(color + (opacity*0.1,))
                mesh = gl.GLMeshItem(vertexes=verts, faces=faces, smooth=True, drawEdges=True,
                                     edgeColor=color + (opacity,))

            self.glView.addItem(mesh)
            self.data_items.append(mesh)

    def marching_cubes(self, volume_data, threshold):
        try:
            verts, faces, _, _ = measure.marching_cubes(volume_data, threshold)
            self.logger.info(f"Marching cubes generated {len(verts)} vertices and {len(faces)} faces")
            return verts, faces
        except RuntimeError as e:
            self.logger.warning(f"Marching cubes failed: {str(e)}")
            return np.array([]), np.array([])

    def clear_slice_markers(self):
        for item in self.slice_marker_items:
            try:
                self.blobGLView.removeItem(item)
            except ValueError:
                self.logger.warning(f"Item {item} not found in blobGLView")

        for item in self.main_slice_marker_items:
            try:
                self.glView.removeItem(item)
            except ValueError:
                self.logger.warning(f"Item {item} not found in glView")

        self.slice_marker_items.clear()
        self.main_slice_marker_items.clear()


    def update_views(self, data, time_point, threshold, channel_controls):
        if data is None:
            self.logger.warning("No data to update views")
            return None, None, None

        num_channels = data.shape[1]
        depth, height, width = data.shape[2:]

        # Prepare 3D RGB images for each view
        combined_xy = np.zeros((depth, height, width, 3))
        combined_xz = np.zeros((height, depth, width, 3))
        combined_yz = np.zeros((height, depth, width, 3))

        # # Define colors for each channel (RGB)
        # channel_colors = [
        #     (1, 0, 0), (0, 1, 0), (0, 0, 1),  # Red, Green, Blue
        #     (1, 1, 0), (1, 0, 1), (0, 1, 1),  # Yellow, Magenta, Cyan
        #     (0.5, 0.5, 0.5), (1, 0.5, 0),     # Gray, Orange
        # ]

        for c in range(num_channels):
            if c < len(channel_controls) and channel_controls[c][0].isChecked():
                opacity = channel_controls[c][1].value() / 100
                channel_data = data[time_point, c]

                # Normalize channel data
                channel_data = (channel_data - channel_data.min()) / (channel_data.max() - channel_data.min() + 1e-8)

                channel_data[channel_data < threshold] = 0

                # Apply color to the channel
                color = self.parent.getChannelColor(c)[:3]  # Get RGB values
                colored_data = channel_data[:, :, :, np.newaxis] * color

                combined_xy += colored_data * opacity
                combined_xz += np.transpose(colored_data, (1, 0, 2, 3)) * opacity
                combined_yz += np.transpose(colored_data, (2, 0, 1, 3)) * opacity


        return combined_xy, combined_xz, combined_yz



    def create_slice_markers(self, x_slice, y_slice, z_slice, width, height, depth):
        x_marker = gl.GLLinePlotItem(pos=np.array([[x_slice, 0, 0], [x_slice, height, 0], [x_slice, height, depth], [x_slice, 0, depth]]),
                                     color=(1, 0, 0, 1), width=2, mode='line_strip')
        y_marker = gl.GLLinePlotItem(pos=np.array([[0, y_slice, 0], [width, y_slice, 0], [width, y_slice, depth], [0, y_slice, depth]]),
                                     color=(0, 1, 0, 1), width=2, mode='line_strip')
        z_marker = gl.GLLinePlotItem(pos=np.array([[0, 0, z_slice], [width, 0, z_slice], [width, height, z_slice], [0, height, z_slice]]),
                                     color=(0, 0, 1, 1), width=2, mode='line_strip')

        self.glView.addItem(x_marker)
        self.glView.addItem(y_marker)
        self.glView.addItem(z_marker)

        x_marker_vis = gl.GLLinePlotItem(pos=np.array([[x_slice, 0, 0], [x_slice, height, 0], [x_slice, height, depth], [x_slice, 0, depth]]),
                                         color=(1, 0, 0, 1), width=2, mode='line_strip')
        y_marker_vis = gl.GLLinePlotItem(pos=np.array([[0, y_slice, 0], [width, y_slice, 0], [width, y_slice, depth], [0, y_slice, depth]]),
                                         color=(0, 1, 0, 1), width=2, mode='line_strip')
        z_marker_vis = gl.GLLinePlotItem(pos=np.array([[0, 0, z_slice], [width, 0, z_slice], [width, height, z_slice], [0, height, z_slice]]),
                                         color=(0, 0, 1, 1), width=2, mode='line_strip')

        self.blobGLView.addItem(x_marker_vis)
        self.blobGLView.addItem(y_marker_vis)
        self.blobGLView.addItem(z_marker_vis)

        self.slice_marker_items.extend([x_marker_vis, y_marker_vis, z_marker_vis])
        self.main_slice_marker_items.extend([x_marker, y_marker, z_marker])

    def visualize_blobs(self, blobs, current_time, channel_colors):
        # Remove old blob visualizations
        for item in self.blob_items:
            self.blobGLView.removeItem(item)
        self.blob_items.clear()

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

            # Add to blob visualization view
            blob_item_vis = gl.GLMeshItem(meshdata=mesh, smooth=True, color=color, shader='shaded')
            blob_item_vis.translate(z, x, y)  # Swapped y and z
            self.blobGLView.addItem(blob_item_vis)
            self.blob_items.append(blob_item_vis)

        self.blobGLView.update()



    def enable_intensity_profile_tool(self, image_view):
        if self.profile_roi is not None:
            image_view.removeItem(self.profile_roi)

        # Create a LineROI
        self.profile_roi = pg.LineROI([0, 0], [100, 100], width=2)
        self.profile_roi.setPen(pg.mkPen('r', width=2))
        image_view.addItem(self.profile_roi)

        # Connect the ROI's sigRegionChanged signal to update the profile
        self.profile_roi.sigRegionChanged.connect(lambda: self.update_profile(image_view))

    def update_profile(self, image_view):
        if self.profile_roi is None:
            self.logger.warning("Profile ROI is None")
            return

        # Get the line coordinates
        handle_positions = self.profile_roi.getSceneHandlePositions()
        self.logger.debug(f"Number of handle positions: {len(handle_positions)}")
        if len(handle_positions) < 2:
            self.logger.warning("Not enough handle positions for the ROI")
            return

        p1, p2 = handle_positions[0][1], handle_positions[-1][1]
        x1, y1 = image_view.getView().mapSceneToView(p1).x(), image_view.getView().mapSceneToView(p1).y()
        x2, y2 = image_view.getView().mapSceneToView(p2).x(), image_view.getView().mapSceneToView(p2).y()
        self.logger.debug(f"Line coordinates: ({x1}, {y1}) to ({x2}, {y2})")

        # Calculate the intensity profile
        x, y, intensity_values = self.calculate_intensity_profile((x1, y1), (x2, y2), image_view)

        # Show the intensity profile
        self.show_intensity_profile(x, y, intensity_values)

    def calculate_intensity_profile(self, start_point, end_point, image_view):
        image_data = image_view.imageItem.image

        self.logger.debug(f"Image data shape: {image_data.shape}")

        # Create a line between start and end points
        num_samples = 1000
        x = np.linspace(start_point[0], end_point[0], num_samples)
        y = np.linspace(start_point[1], end_point[1], num_samples)

        # Round to nearest integer and ensure within image bounds
        x = np.clip(np.round(x).astype(int), 0, image_data.shape[1] - 1)
        y = np.clip(np.round(y).astype(int), 0, image_data.shape[0] - 1)

        # Extract intensity values along the line
        intensity_values = image_data[y, x]

        self.logger.debug(f"Intensity values shape: {intensity_values.shape}")

        return x, y, intensity_values

    def show_intensity_profile(self, x, y, intensity_values):
        self.logger.debug(f"Intensity values shape in show_intensity_profile: {intensity_values.shape}")

        # Create a new window for the intensity profile if it doesn't exist
        if self.profile_window is None:
            self.profile_window = pg.PlotWidget(title="Intensity Profile")
            self.profile_window.setWindowTitle("Intensity Profile")
            self.profile_window.setLabel('left', 'Intensity')
            self.profile_window.setLabel('bottom', 'Position along line')

        # Clear the previous plot
        self.profile_window.clear()

        # Check if we have multiple channels
        if len(intensity_values.shape) > 1 and intensity_values.shape[1] > 1:
            for i in range(intensity_values.shape[1]):
                self.profile_window.plot(range(len(intensity_values)), intensity_values[:, i], pen=(i, intensity_values.shape[1]))
        else:
            # If it's a single channel or already flattened, plot as is
            self.profile_window.plot(range(len(intensity_values)), np.mean(intensity_values, axis=1) if len(intensity_values.shape) > 1 else intensity_values)

        # Make sure the window is visible
        self.profile_window.show()
        self.profile_window.raise_()

        # Force the GUI to update
        QApplication.processEvents()

        self.logger.debug("Intensity profile window should be visible now")
