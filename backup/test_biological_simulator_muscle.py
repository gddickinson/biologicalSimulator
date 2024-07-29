import unittest
import numpy as np
from biological_simulation import BiologicalSimulator

class TestBiologicalSimulator(unittest.TestCase):
    def setUp(self):
        self.simulator = BiologicalSimulator(size=(30, 100, 100), num_time_points=10)

    def test_generate_muscle_cell_shape(self):
        cell_shape, cell_interior, cell_membrane = self.simulator.generate_cell_shape(
            cell_type='muscle',
            size=(30, 100, 100),
            pixel_size=(1,1,1),
            membrane_thickness=1,
            cell_radius=10
        )

        # Print shapes for debugging
        print(f"Shape of cell_shape: {cell_shape.shape}")
        print(f"Shape of cell_interior: {cell_interior.shape}")
        print(f"Shape of cell_membrane: {cell_membrane.shape}")

        # Check that the shapes are correct
        self.assertEqual(cell_shape.shape, (30, 100, 100), "cell_shape has incorrect shape")
        self.assertEqual(cell_interior.shape, (30, 100, 100), "cell_interior has incorrect shape")
        self.assertEqual(cell_membrane.shape, (30, 100, 100), "cell_membrane has incorrect shape")

        # Check that the cell_shape is approximately the sum of cell_interior and cell_membrane
        np.testing.assert_allclose(cell_shape, cell_interior + cell_membrane, atol=1e-6)

        # Check that the cell_interior and cell_membrane don't significantly overlap
        self.assertTrue(np.sum(cell_interior * cell_membrane) / np.sum(cell_shape) < 0.01)

        # Check that the cell is roughly cylindrical
        self.assertTrue(np.any(cell_shape[:, :, 0] > 0))  # Not empty at the start
        self.assertTrue(np.any(cell_shape[:, :, -1] > 0))  # Not empty at the end
        self.assertTrue(np.sum(cell_shape[:, :, 0] > 0) < np.prod(cell_shape.shape[:2]))  # Not filling the entire slice

        # Check that the cell extends along the x-axis and varies in shape
        self.assertFalse(np.all(cell_shape[:, :, 0] == cell_shape[:, :, -1]))
        self.assertFalse(np.all(cell_shape[:, :, 0] == cell_shape[:, :, cell_shape.shape[2]//2]))

        # Check that the cell tapers at the ends
        self.assertTrue(np.sum(cell_shape[:, :, 0]) < np.sum(cell_shape[:, :, cell_shape.shape[2]//2]))
        self.assertTrue(np.sum(cell_shape[:, :, -1]) < np.sum(cell_shape[:, :, cell_shape.shape[2]//2]))

        # Check that all values are between 0 and 1
        self.assertTrue(np.all(cell_shape >= 0) and np.all(cell_shape <= 1))
        self.assertTrue(np.all(cell_interior >= 0) and np.all(cell_interior <= 1))
        self.assertTrue(np.all(cell_membrane >= 0) and np.all(cell_membrane <= 1))

        # Print some additional information for debugging
        print(f"Sum of first slice: {np.sum(cell_shape[:, :, 0])}")
        print(f"Sum of middle slice: {np.sum(cell_shape[:, :, cell_shape.shape[2]//2])}")
        print(f"Sum of last slice: {np.sum(cell_shape[:, :, -1])}")
        print(f"Are first and last slices identical? {np.all(cell_shape[:, :, 0] == cell_shape[:, :, -1])}")
        print(f"Are first and middle slices identical? {np.all(cell_shape[:, :, 0] == cell_shape[:, :, cell_shape.shape[2]//2])}")

if __name__ == '__main__':
    unittest.main()
