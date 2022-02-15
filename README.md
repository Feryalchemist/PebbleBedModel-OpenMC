# Pebble Bed Model-OpenMC

This Code will explain how to use the **Pebble_Model** module. This module was made to create Pebble Bed reactor core model for [**OpenMC Monte Carlo Simulation**](https://github.com/openmc-dev/openmc) with discretely modelled Pebble position. The data from DEM model is used as a reference of Pebble position.

Current feature of this code are:

1. Make HTR10 core model based on [reference](https://www-pub.iaea.org/MTCD/publications/PDF/te_1382_web/TE_1382_Part2.pdf)
2. Grouping Pebble based on its axial, radial position
3. Tracking the pebble position as its move
4. Model randomized TRISO and latticed TRISO
5. Mimic the OTTO recirculation scheme based on algorithm in [reference](https://iopscience.iop.org/article/10.1088/1742-6596/2048/1/012031)

Incoming feature of this code are:

1. Ability to control DEM simulation using [LIGGGHTS DEM](https://github.com/CFDEMproject/LIGGGHTS-PUBLIC)
2. Model other type of reactor
3. Post Processing

Refer and read more in this [link.](https://iopscience.iop.org/article/10.1088/1742-6596/2048/1/012031) for educational usage.
