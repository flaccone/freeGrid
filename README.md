# freeGrid
freeGrid benchmark is meant to offer a common benchmark to test and compare different approaches to the conceptual design of form-resistant structures, from the classical man-based heuristic design, to the experimental-based design, to the computational based design.\
freeGrid sets three Design Baseline Gridshells (DBG) with their spring line partially not constrained (free-edge), subjected to symmetric and asymmetric load conditions.\
Participants are called to modify the DBG(s) and conceive Design Solution Gridshell(s) (DSG) in order to holistically improve their structural, buildability and sustainability  performances, weighted in a single, bulk quantitative metric.  


With the support of the [Italian Council for Steel Structures](https://www.collegiotecniciacciaio.it/)

Under the umbrella of the [International Association for Shell and Spatial Structures](https://iass-structures.org/)

With the industrial partner: [ArcelorMittal Steligence](https://steligence.arcelormittal.com/) 

For more and updated information, please visit the [website](https://sites.google.com/view/freegrid) or contact the [Steering Committee](mailto:freegrid@ctanet.it).


This repository includes two Python scripts:

freeGridConstraints.py performs an automatic checking for the fulfilment of the design geometrical constraints. \
The script needs as input data: \
•	the geometry of the design solution in .obj or .ply file format; \
•	the type DBG (0 for barrel, 1 for dome, 2 for hypar)

freeGridBuildability.py performs an automatic calculation of the buildability goal and performance metrics.\
The script needs as input data:\
•	the geometry of the design solution in .obj or .ply file format; \
•	the number of cross sections of the members.



Required Python >= 3.11 (64 bit). See requirements.txt for detailed package dependencies.

## Installation

Install the required Python packages using pip:

```sh
pip install -r requirements.txt
```
