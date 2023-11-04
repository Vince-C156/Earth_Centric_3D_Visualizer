from graphics import RenderRPO
import pyvista
import numpy as np

def main():
    # Usage
    rpo = RenderRPO()
    # Assuming you have a valid position in ECI coordinates for the chaser
    chaser_position_eci = np.array([7000, 5000, 0])  # Example position in km
    rpo.set_chaser_position(chaser_position_eci)
    rpo.pl.show()
    rpo.focus_on_chaser()
    input()

if __name__ == '__main__':
    main()