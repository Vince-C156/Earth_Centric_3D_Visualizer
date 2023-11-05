import pyvista
from pyvista import examples
import numpy as np
from pyvistaqt import BackgroundPlotter
import os
import datetime
import math

class RenderRPO:

    def __init__(self):
        # Assuming the inverse is simply the transpose since it's an orthogonal matrix
        self.timer_count = 0

        self.chaser_stl_path, self.target_stl_path = 'models/DART_STL.stl', 'models/DART_STL.stl'
        self.pl, self.earth, self.chaser_mesh, self.target_mesh = self._BuildScene()
        self.light = self._set_light_position_by_date(datetime.datetime.now())
        self.pl.add_light(self.light)

        self.earthaxis = pyvista.AxesActor()
        self.earthaxis.origin = self.earth.center
        scale = 1000
        self.earthaxis.SetTotalLength(6371+scale, 6371+scale, 6371+scale)
        self.earthaxis.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor((1, 0, 0))  # Red for X
        self.earthaxis.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor((0, 1, 0))  # Green for Y
        self.earthaxis.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor((0, 0, 1)) 
        
        self.pl.add_axes()


        print('Building scene...')
        print('The units of the scene are in km.')
        print('Scene built.')
        print('Rendering scene...')
        #self.pl.show()

    def _BuildScene(self):
        # Light of the Sun.
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
        #pl.add_light(light)

        pl.add_mesh(earth, texture=earth_texture, smooth_shading=True)
        # Add the chaser and target to the scene with some color for differentiation
        pl.add_mesh(chaser_mesh, color='red')
        pl.add_mesh(target_mesh, color='blue')

        print('Returning plotter and objects.')
        return [pl, earth, chaser_mesh, target_mesh]
    
    def _BuildEci(self):
        earth = examples.planets.load_earth(radius=6371)  # Average radius of Earth 6,378.1370 km 
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

    def set_chaser_position(self, position_eci):
        # Convert ECI to ECEF
        #position_ecef = self._eci_to_ecef(position_eci)

        # Apply a transformation that moves the chaser to the new position
        chaser_current_position = self.chaser_mesh.center
        translation = position_eci - chaser_current_position
        self.chaser_mesh.translate(translation, inplace=True)

        # After updating the position, focus the camera on the chaser
        print(f'The chaser is {self._distanceFromEarthSurface(self.chaser_mesh.center)} km from the surface of the Earth.')

    def focus_on_chaser(self):
        # Ensure that the mesh is valid
        if self.chaser_mesh.n_points == 0:
            print("Chaser mesh has no points, check the STL file.")
            return

        # Calculate the center of the chaser satellite
        chaser_center = self.chaser_mesh.center
        print(f"Chaser center: {chaser_center}")  # Debug info

        # Calculate the direction vector from the Earth's center to the chaser
        earth_center = [0, 0, 0]  # Assuming Earth's center is at the origin
        direction_vector = [chaser_center[i] - earth_center[i] for i in range(3)]

        # Set a closer distance from the chaser for the camera
        distance_factor = 1.5  # You may need to adjust this to ensure the Earth is behind the camera
        distance = self.chaser_mesh.length * distance_factor

        # Set the camera position further back along the direction vector
        camera_position = [chaser_center[i] + distance * direction_vector[i] for i in range(3)]

        # The focal point is the chaser's center
        focal_point = chaser_center

        # Define the view up vector
        view_up = [0, 0, 1]  # Assuming the Z-axis is up

        # Debug info
        print(f"Camera position: {camera_position}")
        print(f"Focal point: {focal_point}")
        print(f"View up: {view_up}")

        # Update the camera settings
        self.pl.camera_position = [camera_position, focal_point, view_up]
        self.pl.reset_camera()
        self.pl.update()

        # Debug: Add a sphere at the chaser center to check visibility
        self.pl.add_mesh(
            pyvista.Sphere(radius=250, center=chaser_center),
            color='red',  # Red color
            opacity=0.5       # Semi-transparent
        )

    def animate(self, update_freq=30):
        """
        Animates the movement of the chaser and the target.
        :param update_freq: Frequency of updates in Hz
        """
        # Calculate the update interval from the frequency
        update_interval = round(1000 / update_freq) 

        # Start the animation using a timer
        timer_id = self.pl.app.startTimer(1000 * update_interval)  # interval in milliseconds
        self.pl.app.timerEvent = self.update_scene
        self.timer_count = 0

    def update_scene(self, event):
        """
        This method is called by the Qt timer at each interval.
        Here you can update the position of the chaser and target meshes.
        """
        # You would add logic here to update the position of the satellite meshes
        # For example, let's just move the chaser in the z direction for demonstration
        #chaser_position_eci = np.array([0, 0, self.timer_count * 10])  # Example increment
        #self.set_chaser_position(chaser_position_eci)

        # Redraw the plotter
        self.pl.update()
        self.pl.render()

        # Increment the counter or compute the new position based on time
        self.timer_count += 1

    def _distanceFromEarthSurface(self, position_eci):
        # Calculate the distance from the surface of the Earth to the satellite
        earth_center, position_eci = np.array(self.earth.center), np.array(position_eci)
        distance = np.linalg.norm(position_eci - earth_center) - 6371
        return distance

    def _set_light_position_by_date(self, date_time):
        """
        Set the position of the light source based on the given date and time.
        :param date_time: A datetime object specifying the date and time.
        """

        # Calculate the day of the year (1-365/366)
        day_of_year = date_time.timetuple().tm_yday
        
        # Assume a circular orbit: calculate the angle in radians the Earth would have moved on its orbit
        # Multiply by 2 * pi to get the angle in radians
        orbit_angle = 2 * math.pi * (day_of_year / 365.25)
        
        # Assume the distance from Earth to Sun is 1 AU (astronomical unit)
        # 1 AU is about 149.6 million km, but you can scale it down for the visualization
        # For example, let's say 1000 units in your scene are equivalent to 1 AU
        distance_scale = 1000

        # Calculate the Sun's position in the ECI frame (Earth Centered Inertial)
        # Assuming Earth is at the center of your scene, at (0, 0, 0)
        sun_x = distance_scale * math.cos(orbit_angle)
        sun_y = distance_scale * math.sin(orbit_angle)
        sun_position_eci = [sun_x, sun_y, 0]  # Assuming the orbit is in the XY plane
        
        # Set the light's direction to be coming from the Sun's position
        light = pyvista.Light(position=sun_position_eci, positional=True, cone_angle=180)
        light.set_direction_angle(0, 0)  # Direction angle is relative to the position

        # Add or update the light in the scene
        return light







