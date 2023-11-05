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
        self.chaser_sphere = pyvista.Sphere(radius=250, center=self.chaser_mesh.center)
        self.chaser_sphere_actor = self.pl.add_mesh(self.chaser_sphere, color='red')


        self.light = self._set_light_position_by_date(datetime.datetime.now())
        self.pl.add_light(self.light)

        self.earthaxis = pyvista.AxesActor()
        self.earthaxis.origin = self.earth.center
        scale = 1000
        self.earthaxis.SetTotalLength(6371+scale, 6371+scale, 6371+scale)
        self.earthaxis.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor((1, 0, 0))  # Red for X
        self.earthaxis.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor((0, 1, 0))  # Green for Y
        self.earthaxis.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor((0, 0, 1)) 
        
        #self.pl.add_axes()
        #self._initaxes()
        self.pl.add_axes(color='white')


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
        # Ensure that the mesh is valid and has points to focus on
        if self.chaser_mesh.n_points == 0:
            print("Chaser mesh has no points, check the STL file.")
            return

        # Calculate the center of the chaser mesh
        chaser_center = self.chaser_mesh.center

        # Calculate the appropriate distance from the chaser for the camera
        # This ensures that the chaser is within the camera's field of view
        distance_factor = 2.5
        distance_from_surface = self._distanceFromEarthSurface(chaser_center)
        if distance_from_surface < 0:
            print("Chaser is below the surface of the Earth. Adjusting distance to a positive value.")
            distance_from_surface = abs(distance_from_surface)

        # Add some altitude to the camera position to make sure it's above the chaser
        camera_altitude = (distance_from_surface + self.chaser_mesh.length) * distance_factor

        # Calculate the direction vector for camera placement
        # This vector should point from the chaser's center to the camera
        earth_center = np.array(self.earth.center)
        direction_vector = chaser_center - earth_center
        direction_vector = direction_vector / np.linalg.norm(direction_vector)  # Normalize the vector

        # Now we calculate the camera position
        camera_position = chaser_center + camera_altitude * direction_vector

        # The focal point is the chaser's center
        focal_point = chaser_center

        # Define the view up vector perpendicular to the direction_vector
        # This is a bit tricky as we need to ensure that the view up vector is indeed perpendicular
        # A simple cross product with another vector not aligned with direction_vector can give us a good approximation
        temp_vector = np.array([1, 0, 0]) if direction_vector[0] == 0 else np.array([0, 1, 0])
        view_up = np.cross(direction_vector, temp_vector)
        view_up = view_up / np.linalg.norm(view_up)  # Normalize the vector

        print(f'Camera Position: {camera_position}')
        print(f'Chaser Center: {chaser_center}')
        print(f'Focal Point: {focal_point}')
        print(f'View Up Vector: {view_up}')

        # Update the camera settings
        self.pl.camera_position = [camera_position, focal_point, view_up]
        self.pl.reset_camera()

        self.pl.remove_actor(self.chaser_sphere_actor)

        self.chaser_sphere = pyvista.Sphere(radius=250, center=self.chaser_mesh.center)
        self.chaser_sphere_actor = self.pl.add_mesh(self.chaser_sphere, color='red')
        self.pl.update()
        self.pl.render()


    def animate(self, update_freq=30):
        """
        Animates the movement of the chaser and the target.`
        :param update_freq: Frequency of updates in Hz
        """
        # Calculate the update interval from the frequency
        update_interval = round(1000 / update_freq) 

        # Start the animation using a timer
        timer_id = self.pl.app.startTimer(1000 * update_interval)  # interval in milliseconds
        self.pl.app.timerEvent = self.update_scene
        self.timer_count = 0

    def update_scene(self, event = None):
        """
        This method is called by the Qt timer at each interval.
        Here you can update the position of the chaser and target meshes.
        """
        # You would add logic here to update the position of the satellite meshes
        # For example, let's just move the chaser in the z direction for demonstration
        #chaser_position_eci = np.array([0, 0, self.timer_count * 10])  # Example increment
        #self.set_chaser_position(chaser_position_eci)

        # Redraw the plotter

        translation = np.array(self.chaser_mesh.center) - np.array(self.chaser_sphere.center)
        self.chaser_sphere.translate(translation)
        #self.pl.add_mesh(self.chaser_sphere, color='red')


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







