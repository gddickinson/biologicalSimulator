#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 07:01:21 2024

@author: george
"""

import numpy as np
import logging
import tifffile

from scipy.ndimage import rotate
from typing import Tuple, List, Optional, Any, Dict


class DataGenerator:
    def __init__(self):
        self.data = None
        self.metadata = {}
        self.initLogging()

    def initLogging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

    def generate_multi_channel_volume(self, size=(100, 100, 30), num_channels=2, num_blobs=30,
                                      intensity_range=(0.5, 1.0), sigma_range=(2, 6), noise_level=0.02):
        volume = np.zeros((num_channels, *size))
        for c in range(num_channels):
            for _ in range(num_blobs):
                x, y, z = [np.random.randint(0, s) for s in size]
                sigma = np.random.uniform(*sigma_range)
                intensity = np.random.uniform(*intensity_range)

                x_grid, y_grid, z_grid = np.ogrid[-x:size[0]-x, -y:size[1]-y, -z:size[2]-z]
                blob = np.exp(-(x_grid*x_grid + y_grid*y_grid + z_grid*z_grid) / (2*sigma*sigma))
                volume[c] += intensity * blob

            volume[c] += np.random.normal(0, noise_level, size)
            volume[c] = np.clip(volume[c], 0, 1)
        return volume

    def _generate_single_volume(
        self,
        size: Tuple[int, int, int],
        blob_positions: np.ndarray,
        blob_velocities: np.ndarray,
        intensity_range: Tuple[float, float],
        sigma_range: Tuple[float, float],
        noise_level: float
    ) -> np.ndarray:
        z, y, x = size
        volume = np.zeros((z, y, x))

        for i, (bz, by, bx) in enumerate(blob_positions):
            sigma = np.random.uniform(*sigma_range)
            intensity = np.random.uniform(*intensity_range)

            zz, yy, xx = np.ogrid[
                max(0, int(bz-3*sigma)):min(z, int(bz+3*sigma)),
                max(0, int(by-3*sigma)):min(y, int(by+3*sigma)),
                max(0, int(bx-3*sigma)):min(x, int(bx+3*sigma))
            ]
            blob = np.exp(-((zz-bz)**2 + (yy-by)**2 + (xx-bx)**2) / (2*sigma*sigma))
            volume[zz, yy, xx] += intensity * blob

        volume += np.random.normal(0, noise_level, (z, y, x))
        return np.clip(volume, 0, 1)

    def generate_multi_channel_time_series(
        self,
        num_volumes: int,
        num_channels: int = 2,
        size: Tuple[int, int, int] = (30, 100, 100),
        num_blobs: int = 30,
        intensity_range: Tuple[float, float] = (0.5, 1.0),
        sigma_range: Tuple[float, float] = (2, 6),
        noise_level: float = 0.02,
        movement_speed: float = 1.0
    ) -> np.ndarray:
        z, y, x = size
        time_series = np.zeros((num_volumes, num_channels, z, y, x))
        blob_positions = np.random.rand(num_channels, num_blobs, 3) * np.array([z, y, x])
        blob_velocities = np.random.randn(num_channels, num_blobs, 3) * movement_speed

        for t in range(num_volumes):
            for c in range(num_channels):
                volume = self._generate_single_volume(
                    size, blob_positions[c], blob_velocities[c],
                    intensity_range, sigma_range, noise_level
                )
                time_series[t, c] = volume

                # Update blob positions
                blob_positions[c] += blob_velocities[c]
                blob_positions[c] %= [z, y, x]  # Wrap around the volume

        self.data = time_series
        self.logger.info(f"Generated time series with shape: {time_series.shape}")
        return time_series


    def generate_volume(self, size: Tuple[int, int, int] = (100, 100, 30),
                        num_blobs: int = 30,
                        intensity_range: Tuple[float, float] = (0.8, 1.0),
                        sigma_range: Tuple[float, float] = (2, 6),
                        noise_level: float = 0.02) -> np.ndarray:
        """Generate a single volume of data."""
        volume = np.zeros(size)
        for _ in range(num_blobs):
            x, y, z = [np.random.randint(0, s) for s in size]
            sigma = np.random.uniform(*sigma_range)
            intensity = np.random.uniform(*intensity_range)

            x_grid, y_grid, z_grid = np.ogrid[-x:size[0]-x, -y:size[1]-y, -z:size[2]-z]
            blob = np.exp(-(x_grid*x_grid + y_grid*y_grid + z_grid*z_grid) / (2*sigma*sigma))
            volume += intensity * blob

        volume += np.random.normal(0, noise_level, size)
        volume = np.clip(volume, 0, 1)
        return volume

    def generate_time_series(self, num_volumes: int,
                             size: Tuple[int, int, int] = (100, 100, 30),
                             num_blobs: int = 30,
                             intensity_range: Tuple[float, float] = (0.8, 1.0),
                             sigma_range: Tuple[float, float] = (2, 6),
                             noise_level: float = 0.02,
                             movement_speed: float = 1.0) -> np.ndarray:
        """Generate a time series of volumes with moving blobs."""
        time_series = np.zeros((num_volumes, *size))
        blob_positions = np.random.rand(num_blobs, 3) * np.array(size)
        blob_velocities = np.random.randn(num_blobs, 3) * movement_speed

        for t in range(num_volumes):
            volume = np.zeros(size)
            for i in range(num_blobs):
                x, y, z = blob_positions[i]
                sigma = np.random.uniform(*sigma_range)
                intensity = np.random.uniform(*intensity_range)

                x_grid, y_grid, z_grid = np.ogrid[-x:size[0]-x, -y:size[1]-y, -z:size[2]-z]
                blob = np.exp(-(x_grid*x_grid + y_grid*y_grid + z_grid*z_grid) / (2*sigma*sigma))
                volume += intensity * blob

                # Update blob position
                blob_positions[i] += blob_velocities[i]
                blob_positions[i] %= size  # Wrap around the volume

            volume += np.random.normal(0, noise_level, size)
            time_series[t] = np.clip(volume, 0, 1)

        self.data = time_series
        return time_series



    def generate_structured_multi_channel_time_series(self, num_volumes, num_channels=2, size=(30, 100, 100),
                                                      num_blobs=30, intensity_range=(0.5, 1.0), sigma_range=(2, 6),
                                                      noise_level=0.02, movement_speed=1.0, channel_ranges=None):
        z, y, x = size
        time_series = np.zeros((num_volumes, num_channels, z, y, x))

        if channel_ranges is None:
            channel_ranges = [((0, x), (0, y), (0, z)) for _ in range(num_channels)]

        blob_positions = []
        blob_velocities = []

        for c in range(num_channels):
            x_range, y_range, z_range = channel_ranges[c]
            channel_blobs = np.random.rand(num_blobs, 3)
            channel_blobs[:, 0] = channel_blobs[:, 0] * (x_range[1] - x_range[0]) + x_range[0]
            channel_blobs[:, 1] = channel_blobs[:, 1] * (y_range[1] - y_range[0]) + y_range[0]
            channel_blobs[:, 2] = channel_blobs[:, 2] * (z_range[1] - z_range[0]) + z_range[0]
            blob_positions.append(channel_blobs)
            blob_velocities.append(np.random.randn(num_blobs, 3) * movement_speed)

        for t in range(num_volumes):
            for c in range(num_channels):
                volume = np.zeros((z, y, x))
                x_range, y_range, z_range = channel_ranges[c]
                for i in range(num_blobs):
                    bx, by, bz = blob_positions[c][i]
                    sigma = np.random.uniform(*sigma_range)
                    intensity = np.random.uniform(*intensity_range)

                    xx, yy, zz = np.ogrid[max(0, int(bx-3*sigma)):min(x, int(bx+3*sigma)),
                                          max(0, int(by-3*sigma)):min(y, int(by+3*sigma)),
                                          max(0, int(bz-3*sigma)):min(z, int(bz+3*sigma))]
                    blob = np.exp(-((xx-bx)**2 + (yy-by)**2 + (zz-bz)**2) / (2*sigma*sigma))
                    volume[zz, yy, xx] += intensity * blob

                    # Update blob position
                    blob_positions[c][i] += blob_velocities[c][i]
                    blob_positions[c][i][0] = np.clip(blob_positions[c][i][0], x_range[0], x_range[1])
                    blob_positions[c][i][1] = np.clip(blob_positions[c][i][1], y_range[0], y_range[1])
                    blob_positions[c][i][2] = np.clip(blob_positions[c][i][2], z_range[0], z_range[1])

                volume += np.random.normal(0, noise_level, (z, y, x))
                time_series[t, c] = np.clip(volume, 0, 1)

        return time_series

    def simulate_angular_recording(self, angle: float) -> np.ndarray:
        """Simulate an angular recording by rotating the volume."""
        if self.data is None:
            raise ValueError("No data to rotate. Generate data first.")
        rotated_data = rotate(self.data, angle, axes=(1, 2), reshape=False, mode='constant', cval=0)
        return rotated_data

    def correct_angular_recording(self, angle: float) -> np.ndarray:
        """Correct an angular recording by rotating the volume back."""
        if self.data is None:
            raise ValueError("No data to correct. Generate data first.")
        corrected_data = rotate(self.data, -angle, axes=(1, 2), reshape=False, mode='constant', cval=0)
        return corrected_data

    def save_tiff(self, filename: str):
        """Save the data as a TIFF stack."""
        if self.data is None:
            raise ValueError("No data to save. Generate data first.")
        tifffile.imwrite(filename, self.data)

    def save_numpy(self, filename: str):
        """Save the data as a numpy array."""
        if self.data is None:
            raise ValueError("No data to save. Generate data first.")
        np.save(filename, self.data)

    def load_tiff(self, filename: str):
        """Load data from a TIFF stack."""
        self.data = tifffile.imread(filename)
        return self.data

    def load_numpy(self, filename: str):
        """Load data from a numpy array file."""
        self.data = np.load(filename)
        return self.data

    def apply_gaussian_filter(self, sigma: float):
        """Apply a Gaussian filter to the data."""
        from scipy.ndimage import gaussian_filter
        if self.data is None:
            raise ValueError("No data to filter. Generate or load data first.")
        self.data = gaussian_filter(self.data, sigma)
        return self.data

    def adjust_intensity(self, gamma: float):
        """Adjust the intensity of the data using gamma correction."""
        if self.data is None:
            raise ValueError("No data to adjust. Generate or load data first.")
        self.data = np.power(self.data, gamma)
        return self.data

    def get_metadata(self) -> dict:
        """Return metadata about the current dataset."""
        if self.data is None:
            return {}
        self.metadata.update({
            'shape': self.data.shape,
            'dtype': str(self.data.dtype),
            'min': float(np.min(self.data)),
            'max': float(np.max(self.data)),
            'mean': float(np.mean(self.data)),
            'std': float(np.std(self.data))
        })
        return self.metadata

    def set_metadata(self, key: str, value: Any):
        """Set a metadata value."""
        self.metadata[key] = value
