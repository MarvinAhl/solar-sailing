# Solar Sailing Trajectory Design Tool
Author: Marvin Ahlborn</br>
Date: 2026-03-19

Implements heliocentric solar sailing only considering sun gravity and solar radiation pressure. The control is described by two angles and considered piece-wise constant. Monotonic Basin Hopping is used for optimizing solar sailing legs. The main goal is a 5-year Near-Earth Asteroid flyby tour. The combinatorial problem is solved with Beam Search.

An example tour starting from Earth in 2030-01-01 and visiting 16 asteroids before 2035-01-01.
<img src="https://github.com/MarvinAhl/solar-sailing/blob/main/output/nea_tour_3.png" alt="5-year solar sailing NEA tour" width="500"/>

Example of a direct Earth to Venus transfer in the heliocentric frame.
<img src="https://github.com/MarvinAhl/solar-sailing/blob/main/output/earth_venus.png" alt="Earth to Venus transfer" width="400"/>
