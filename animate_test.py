import sys
import numpy as np
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QTimer

from graphics import RenderRPO 

def generate_orbital_trajectory(semi_major_axis, eccentricity, num_points):
    """
    Generates a set of points representing an elliptical orbit.

    :param semi_major_axis: Semi-major axis of the ellipse (in km)
    :param eccentricity: Eccentricity of the ellipse
    :param num_points: Number of points to generate along the orbit
    :return: Array of points in the orbit
    """
    orbit_points = []
    for i in range(num_points):
        # True anomaly varies from 0 to 2*pi
        true_anomaly = 2 * np.pi * (i / num_points)
        
        # Compute the distance r using the formula r = a(1-e^2) / (1 + e*cos(v))
        # where a is the semi-major axis, e is the eccentricity, and v is the true anomaly.
        r = (semi_major_axis * (1 - eccentricity ** 2)) / (1 + eccentricity * np.cos(true_anomaly))
        
        # Convert polar coordinates (r, v) to Cartesian coordinates (x, y, z)
        x = r * np.cos(true_anomaly)
        y = r * np.sin(true_anomaly)
        z = 0  # Assuming a 2D ellipse in the XY-plane
        
        orbit_points.append([x, y, z])
    
    return np.array(orbit_points)

# Parameters for the orbit
semi_major_axis = 7000  # km, for example
eccentricity = 0.01  # nearly circular
num_points = 360  # number of points in the trajectory

# Generate the orbital trajectory points
trajectory = generate_orbital_trajectory(semi_major_axis, eccentricity, num_points)

# Create the application instance
app = QApplication(sys.argv)

# Create an instance of your RenderRPO class
render_rpo_instance = RenderRPO()
render_rpo_instance.pl.show()
# Setup a timer to update the position
timer = QTimer()
timer.setInterval(1000 // 30)  # Update frequency, here it's set for 30 Hz

index = 0  # To keep track of the current position in the trajectory
def update_position():
    global index
    render_rpo_instance.set_chaser_position(trajectory[index])
    render_rpo_instance.focus_on_chaser()
    render_rpo_instance.update_scene()  # Trigger a repaint which should call update_scene internally
    index = (index + 1) % num_points  # Loop over the trajectory points


# Connect the timer timeout signal to the update_position function
timer.timeout.connect(update_position)
timer.start()

# Start the application main loop
sys.exit(app.exec_())