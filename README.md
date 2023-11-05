# RenderRPO Class Documentation

## Overview
The `RenderRPO` class is designed to render a Relative Positioning Operation (RPO) scenario in space with the Earth, a chaser satellite, and a target satellite using the `pyvista` and `pyvistaqt` libraries.

## Dependencies
To use the `RenderRPO` class, the following libraries must be installed:
```python
import pyvista
from pyvista import examples
import numpy as np
from pyvistaqt import BackgroundPlotter
import os
import datetime
import math
```

## Class Definition
```python
class RenderRPO:

    def __init__(self):
        self.timer_count = 0

        self.chaser_stl_path, self.target_stl_path = 'models/DART_STL.stl', 'models/DART_STL.stl'
        self.pl, self.earth, self.chaser_mesh, self.target_mesh = self._BuildScene()
        self.light = self._set_light_position_by_date(datetime.datetime.now())
        self.pl.add_light(self.light)

        self.earthaxis = pyvista.AxesActor()
        self.earthaxis.origin = self.earth.center
        scale = 1000
        self.earthaxis.SetTotalLength(6371+scale, 6371+scale, 6371+scale)
        self.earthaxis.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor((1, 0, 0))
        self.earthaxis.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor((0, 1, 0))
        self.earthaxis.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor((0, 0, 1))
        
        self.pl.add_axes()

        print('Building scene...')
        print('The units of the scene are in km.')
        print('Scene built.')
        print('Rendering scene...')

    def _BuildScene(self):
        earth, earth_texture = self._BuildEci()

        print('Building satellites...')
        print(f'Chaser: {self.chaser_stl_path}')
        print(f'Target: {self.target_stl_path}')
        chaser_mesh, target_mesh = self._BuildSatellites(self.chaser_stl_path, self.target_stl_path)

        pl = BackgroundPlotter(lighting="none")
        cubemap = examples.download_cubemap_space_16k()
        _ = pl.add_actor(cubemap.to_skybox())
        pl.set_environment_texture(cubemap, True)

        pl.add_mesh(earth, texture=earth_texture, smooth_shading=True)
        pl.add_mesh(chaser_mesh, color='red')
        pl.add_mesh(target_mesh, color='blue')

        print('Returning plotter and objects.')
        return [pl, earth, chaser_mesh, target_mesh]

    # ... (Other methods continue in a similar fashion)
```

## Methods
- `__init__(self)`: Initializes the rendering plotter, loads the meshes, and sets up the scene.
- `_BuildScene(self)`: Constructs the visual scene with the Earth, chaser, and target meshes.
- `_BuildEci(self)`: Creates the Earth mesh with texture.
- `_BuildSatellites(self, chaser_stl_path, target_stl_path)`: Loads the STL models for the chaser and target satellites.
- `set_chaser_position(self, position_eci)`: Sets the chaser's position in Earth-Centered Inertial (ECI) coordinates.
- `focus_on_chaser(self)`: Adjusts the camera to focus on the chaser satellite.
- `animate(self, update_freq=30)`: Initiates the animation of the satellites.
- `update_scene(self, event)`: Updates the scene; meant to be called by a Qt timer.
- `_distanceFromEarthSurface(self, position_eci)`: Calculates the distance of an object from Earth's surface.
- `_set_light_position_by_date(self, date_time)`: Sets the lighting based on the current date and time to simulate sunlight.

## Usage
To use the `RenderRPO` class, instantiate it and call the appropriate methods for your simulation needs. Ensure that the `models/DART_STL.stl` paths are correctly set to the locations of your STL files.

## Notes
- The class is designed to work with the `pyvistaqt` BackgroundPlotter for interactive visualization.
- Error handling for file paths and STL file integrity should be considered when implementing in a production environment.
- The `_set_light_position_by_date` function is a simplification and may need more accurate astronomical calculations for precise simulations.

## Authors
- [Vincent Chen]
