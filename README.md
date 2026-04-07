# 3D ECI Interferometer Kinematics Simulator

This repository contains a high-precision, 3D Earth-Centered Inertial (ECI) kinematic simulation of classical interferometer experiments, specifically modeling the **Michelson-Morley** and **Michelson-Gale-Pearson** experiments.

By calculating the exact time-of-flight of light traveling between moving mirrors on the surface of a spherical, rotating Earth, this simulator achieves sub-nanometer, 50-decimal-place precision to calculate relativistic and kinematic fringe shifts.

## Core Concepts

* **Absolute Light Speed:** Light travels at the constant speed `c` relative to the Earth-Centered Inertial (ECI) reference frame, independently of the motion of source and receiver.
* **The ECI Frame:** The reference frame is centered on the Earth, aligned with its axis of rotation, and maintains a fixed orientation relative to the distant stars.
* **Dynamic Geometry:** The mirrors of the interferometer travel dynamically with the surface of the rotating Earth while the light is in transit from one mirror to the next.
* **Drift Testing:** The simulator allows you to introduce a hypothetical drift vector to accurately test and visualize the original Michelson-Morley hypothesis.

## Features

The model calculates the exact 3D distances and times of flight, naturally accounting for:
* **The Sagnac Effect:** Instantiated via Earth's physical rotation.
* **Geographical Flexibility:** Simulates experiments at arbitrary geographical latitudes.
* **Custom Apparatus Shapes:** Supports custom arm lengths and intersection angles (e.g., modeling the 60° arms of the Einstein Telescope).
* **Spherical Curvature:** Utilizes exact 3D spatial chord adjustments via an iterative Newton-Raphson approximation to ensure macroscopic flat distances accurately drape over Earth's curvature.
* **Automated Reporting:** Generates a compiled PDF report containing numerical results, top-down geometric setup schematics, and 360-degree rotation kinematic plots.

## Prerequisites & Installation

This project uses **Conda** to manage dependencies and ensure mathematical plotting libraries build correctly. 

1. Clone this repository to your local machine.
2. Navigate to the project directory in your terminal.
3. Create the isolated Conda environment from the provided configuration file:
   ```bash
   conda env create -f environment.yml
4. Activate the new environment:
   ```bash
    conda activate eci-simulator

## Usage
The simulator is driven by a YAML configuration file, allowing you to define multiple experimental scenarios without editing the Python codebase.

Ensure scenarios.yaml is present in the same directory as the script.

Run the simulator:
```bash
    python ECI_Interferometer_simulator_michelson_gale_morley.py
```

## Output
The script will output progress to the console and generate a file named Interferometer_Report.pdf in your directory. 

Depending on your YAML configuration, this report will include:
* Calculated Time Differences (ns)
* Calculated Distance Differences (nm)
* Total Fringe Shifts
* 2D Overhead Apparatus Schematics
* 360-Degree Kinematic Rotation Graphs