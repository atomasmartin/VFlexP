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

a) Online testing of the tool

<a href="https://matlab.mathworks.com/open/github/v1?repo=atomasmartin/VFlexP" target="_blank">
  <img src="https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg" alt="Open in MATLAB Online">
</a>

b) Local installation
1. Clone this repository
2. Ensure you have MATLAB R2024a or later installed
3. Run 'install.m' in the root directory
4. Refer to the documentation in the `Documentation/` folder for usage instructions

## Documentation

Comprehensive documentation is available in the `Documentation/` folder, including:

- Details on the modelling of power systems
- Instructions for using the custom functions
- Structure of the parameters of case files

The documentation, available in HTML format for easy navigation, is also accessible online at [VFlexP Documentation](https://atomasmartin.github.io/VFlexP/).

## Contributing

This project is in early development, and bug reports and feature requests are welcome. If you encounter any issues or have suggestions for improvements, please open an issue on the GitHub repository.

## Acknowledgments

- Contributors to the MATPOWER-FUBM extension, required for the powerflow analysis implemented in the VFlexP tool. Full documentation of the MATPOWER-FUBM extension can be found in [MATPOWER-FUBM GitHub
](https://github.com/AbrahamAlvarezB/matpower-fubm).

- The development of this tool was possible thanks to the advice and support of the researchers within the Modelling, Analysis and Control group of the Institute for Research in Technology (IIT) at Comillas Pontifical University, with special mention to Luis Rouco, Fidel Fernández-Bernal and Ignacio Egido Cortes.

- This work has been partially financed through the research program S2018/EMT-4366 PROMINT-CAM on Smart Grids of Madrid Government, Spain, with 50\% support from the European Social Fund (ESF), and through Grant TED2021-130610B-C22 funded by MICIU/AEI/10.13039/501100011033 and by European Union NextGenerationEU/PRTR.

- The collaboration between Durham University and Comillas has been eased with a special grant to Andrés Tomás-Martín from Comillas Pontifical University for temporary stays of researchers in foreign research centres.  

- This work is also a collaboration with the project PID2021-125628OB-C21 funded by MICIU/AEI/10.13039/501100011033 and by ERDF/EU, and the Grant PRE2022-101606 funded by MICIU/AEI/10.13039/501100011033 and by ESF+.

## Citation

If you use this tool in your research, please cite the following paper:

Andrés Tomás-Martín, Carlos David Zuluaga-Ríos, Jorge Suárez-Porras, Behzad Kazemtabrizi, Javier García-Aguilar, Lukas Sigrist, Aurelio García-Cerrada, A Vector-Based Flexible-Complexity Tool for Simulation and Small-Signal Analysis of Hybrid AC/DC Power Systems, accepted for publication in IREP Special Issue of the journal Sustainable Energy, Grids and Networks, 2025.

If you use the static analysis included in this tool, implemented using the MATPOWER-FUBM extension, please cite:

Abraham Alvarez-Bustos, Behzad Kazemtabrizi, Mahmoud Shahbazi, Enrique Acha-Daza, "Universal branch model for the solution of optimal power flows in hybrid AC/DC grids", International Journal of Electrical Power & Energy Systems, Volume 126, Part A, 2021, 106543

R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems Research and Education," Power Systems, IEEE Transactions on, vol. 26, no. 1, pp. 12-19, Feb. 2011. doi: 10.1109/TPWRS.2010.2051168

## Contact

For questions or support, please contact atomas@comillas.edu.
