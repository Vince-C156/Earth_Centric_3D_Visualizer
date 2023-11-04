import pyvista
from pyvista import examples
import numpy as np
from pyvistaqt import BackgroundPlotter
import os

class RenderRPO:

    def __init__(self):
        self.gmst_matrix = np.array([
            [0.540302, -0.841471, 0],
            [0.841471, 0.540302, 0],
            [0, 0, 1]
        ])
        # Assuming the inverse is simply the transpose since it's an orthogonal matrix
        self.gmst_matrix_inv = self.gmst_matrix.T

        self.chaser_stl_path, self.target_stl_path = 'models/DART_STL.stl', 'models/DART_STL.stl'
        self.pl, self.earth, self.chaser_mesh, self.target_mesh = self._BuildScene()
        print('Building scene...')
        print('Scene built.')
        print('Rendering scene...')
        #self.pl.show()

    def _BuildScene(self):
        # Light of the Sun.
        light = pyvista.Light()
        light.set_direction_angle(30, -20)
        earth, earth_texture = self._BuildEci()

        print('Building satellites...')
        print(f'Chaser: {self.chaser_stl_path}')
        print(f'Target: {self.target_stl_path}')
        chaser_mesh, target_mesh = self._BuildSatellites(self.chaser_stl_path, self.target_stl_path)

        # Add planets to Plotter.
        #pl = pyvista.Plotter()
        pl = BackgroundPlotter(lighting="none")
        cubemap = examples.download_cubemap_space_16k()
        _ = pl.add_actor(cubemap.to_skybox())
        pl.set_environment_texture(cubemap, True)
        pl.add_light(light)

        pl.add_mesh(earth, texture=earth_texture, smooth_shading=True)
        # Add the chaser and target to the scene with some color for differentiation
        pl.add_mesh(chaser_mesh, color='red')
        pl.add_mesh(target_mesh, color='blue')

        print('Returning plotter and objects.')
        return [pl, earth, chaser_mesh, target_mesh]
    
    def _BuildEci(self):
        earth = examples.planets.load_earth(radius=6371)  # Average radius of Earth in km
        earth_texture = examples.load_globe_texture()
        return earth, earth_texture

    def _BuildSatellites(self, chaser_stl_path, target_stl_path):
        # Load the STL files for chaser and target
        assert os.path.exists(self.chaser_stl_path), "Chaser STL file not found"
        assert os.path.exists(self.target_stl_path), "Target STL file not found"
        chaser_mesh = pyvista.read(chaser_stl_path)
        target_mesh = pyvista.read(target_stl_path)

        # Transformations (scaling, rotating, translating) can be applied to the meshes if needed

        return chaser_mesh, target_mesh

    def _ecef_to_eci(self, coord_ecef):
        coord_eci = self.gmst_matrix @ coord_ecef
        return coord_eci

    def _eci_to_ecef(self, coord_eci):
        coord_ecef = self.gmst_matrix_inv @ coord_eci
        return coord_ecef

    def set_chaser_position(self, position_eci):
        # Convert ECI to ECEF
        position_ecef = self._eci_to_ecef(position_eci)

        # Apply a transformation that moves the chaser to the new position
        chaser_current_position = self.chaser_mesh.center
        translation = position_ecef - chaser_current_position
        self.chaser_mesh.translate(translation, inplace=True)

        # After updating the position, focus the camera on the chaser
        self.focus_on_chaser()

    def focus_on_chaser(self):
        # Ensure that the mesh is valid
        if self.chaser_mesh.n_points == 0:
            print("Chaser mesh has no points, check the STL file.")
            return

        # Calculate the center of the chaser satellite
        chaser_center = self.chaser_mesh.center
        print(f"Chaser center: {chaser_center}")  # Debug info

        # Set a closer distance from the chaser for the camera
        distance_factor = 0.5
        distance = self.chaser_mesh.length * distance_factor

        # Define the camera position, focal point, and view up
        camera_position = [chaser_center[0], chaser_center[1] - distance, chaser_center[2]]
        focal_point = chaser_center
        view_up = [0, 0, 1]

        # Debug info
        print(f"Camera position: {camera_position}")
        print(f"Focal point: {focal_point}")
        print(f"View up: {view_up}")

        # Update the camera settings
        self.pl.camera_position = [camera_position, focal_point, view_up]
        self.pl.reset_camera()
        self.pl.update()

        # Debug: Add a sphere at the chaser center to check visibility
        self.pl.add_mesh(pyvista.Sphere(radius=1, center=chaser_center), color='white')






