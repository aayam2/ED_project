# PHY204 Project
The completed notebook for PHY204 project: laser beam to the moon through the atmosphere

Techniques learnt and utilized:
- A 2D-spectral method to numerically solve a differential equation
- Implementing the split step method for additional (non-differential) terms in the equation
- Implementing phase screens to incorporate atmospheric turbulence
- Change of coordinates to a dynamic computational space to prevent the solution from overflowing the finite initial grid
- Implementation of a absorption coefficient around the borders of the grid to prevent the expanding noise from bouncing back and interfering with the wave.


Techniques learnt but not included in the final notebook in github:
- The Crank Nicolson scheme in 2D
- Creating a basic matplotlib animation
