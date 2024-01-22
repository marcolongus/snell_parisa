# Particle Simulation

This repository contains a simple particle simulation program with customizable parameters. The simulation consists of particles represented by the `agents.h` and `particle.h` classes, with configuration parameters defined in the `parameters.h` file. The simulation results are stored in the `data` folder.

## Instructions

1. Navigate to the `headers` folder, where you will find the following files:
   - `agents.h`
   - `particle.h`
   - `parameters.h`

2. Open the `parameters.h` file and locate the variables `velocity`, `alpha_left`, and `alpha_right`. These variables control the velocity and rotation rates of the particles. Modify these values according to your simulation requirements.

3. Once you have adjusted the parameters, you can compile the simulation by running the following command in the command prompt (cmd):
   ```
   python makefile.py
   ```

   This will compile the `agents.cpp` file and execute the simulation with the updated parameters.

4. The simulation results will be stored in the `data` folder.

## File Structure

- `headers/`
  - `agents.h`: Class definition for the agents in the particle simulation.
  - `particle.h`: Class definition for the particles in the simulation.
  - `parameters.h`: Configuration file where you can adjust simulation parameters.

- `data/`: Folder to store simulation results.

- `makefile.py`: Python script to compile the `agents.cpp` file.

Feel free to experiment with different parameter values to observe how they affect the behavior of the particle simulation. 