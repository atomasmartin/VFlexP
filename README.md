# VFlexP: A Vector-Based Flexible-Complexity Tool for Simulation and Small-Signal Analysis of Hybrid AC/DC Power Systems

This repository contains the VFlexP tool, a MATLAB/Simulink-based tool for the simulation and small-signal analysis of hybrid AC/DC power systems. The tool is designed to facilitate the modelling and analysis of complex power systems, with a focus on flexibility and ease of use.

## Key Components

1. **Functions**: A set of custom MATLAB functions for the initialisation and analysis of power systems.
2. **Modal Analysis Tool**: Simulink-based tool for power system modal analysis.
3. **Documentation**: Extensive HTML documentation.

## Features

- Modal analysis capabilities for power system stability studies
- Integration with MATLAB/Simulink for dynamic simulations
- Comprehensive documentation and examples

## Directory Structure

- `Documentation/`: User manuals, technical documentation, and figures
- `functions/`: Custom MATLAB functions for various analyses
- `bus_systems/`: Sample bus systems for testing
- `VFlexP.slx`: VFlexP simulink model

## Getting Started

Options:

a) Online testing

<a href="https://matlab.mathworks.com/open/github/v1?repo=atomasmartin/VFlexP" target="_blank">
  <img src="https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg" alt="Open in MATLAB Online">
</a>

b) Local installation
1. Clone this repository
2. Ensure you have MATLAB R2024a or later installed
3. Run 'install.m' in the root directory
4. Refer to the documentation in the `Documentation/` folder for usage instructions

Note: The MATPOWER-FUBM extension is required for the powerflow analysis implemented in the VFlexP tool. You can download it from [matpower-fubm GitHub
](https://github.com/AbrahamAlvarezB/matpower-fubm).

## Documentation

Comprehensive documentation is available in the `Documentation/` folder, including:

- Details on the modelling of power systems
- Instructions for using the custom functions
- Structure of the parameters of case files

## Contributing

Contributions to this project are welcome. Please refer to the contributing guidelines (TODO: add link) for more information.

## Acknowledgments

- Contributors to the MATPOWER-FUBM extension

## Citation

If you use this tool in your research, please cite the following paper:

Andrés Tomás-Martín, Carlos David Zuluaga-Ríos, Jorge Suárez-Porras, Behzad Kazemtabrizi, Javier García-Aguilar, Lukas Sigrist, Aurelio García-Cerrada, A Vector-Based Flexible-Complexity Tool for Simulation and Small-Signal Analysis of Hybrid AC/DC Power Systems, submitted to Sustainable Energy, Grid and Networks, 2025.



## Contact

For questions or support, please contact atomas@comillas.edu.
