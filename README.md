# VFlexP: A Vector-Based Flexible-Complexity Tool for Simulation and Small-Signal Analysis of Hybrid AC/DC Power Systems

This repository contains **VFlexP**, a MATLAB/Simulink tool for the time-domain simulation, small-signal modal analysis, and state-relevance model reduction of hybrid AC/DC power systems. It is built on the same vector-based, flexible-complexity philosophy and shares its case file formats (`.xlsx`, `.dyr`, `.json`) and multi-tab workflow with its Python counterpart, **VFlexPy**.

## Key Features

* **Pure-`.m` Multi-Tab MATLAB App (`vflexp`) & Decoupled Scripting Core (`core/`).** Launch `vflexp` for an interactive 6-tab graphical interface (*Static Analysis*, *System Diagram*, *State Relevance*, *Complexity*, *Simulation*, and *Modal Analysis*) written entirely in plain-text `.m` classes (no binary `.mlapp`), or drive `VFlexPSession` and `core/vflexp_*.m` headlessly from MATLAB scripts and automated pipelines.
* **Vector-based modelling in Simulink.** Elements of the same type and level of detail are grouped into a single vectorised Simulink block whose parameters and signals are vectors/matrices.
* **Modular Simulink libraries with full `git diff` traceability.** All 31 primitive, generator, converter, network, and control subsystems live as standalone Simulink libraries under `Simulink/lib/`, linked (`LinkStatus = 'resolved'`) into `Simulink/VFlexP.slx` and parameterised via `ModelWorkspace` structures (`prm_*`, `op_*`). Loading cases, changing complexity, initialising, simulating, and linearising never mutate blocks or mark `.slx` files dirty. Every library is mirrored to a deterministic plain-text netlist (`Simulink/spec/*.spec.txt`) and configured with a `.gitattributes` `diff=slx` driver (`tools/slx_textconv.py`).
* **Flexible complexity & automatic reduction.** Every network element (AC/DC buses, lines, loads) and device (`MGES`, `SG`, `VSM`, `GFL`, `GFM`, `DC_PES`, `intconv_GFL`, `intconv_VSM`, `intconv_GFM`, `DFIG`, `IG`, `SC`, `MSG`, `DCBUCK`) selects its own level of detail (e.g. `order 8` / `EMT8` vs `order 4` / `RMS4` synchronous machines, `dynamic` vs `static` converter loops and network elements) via expression-driven Simulink Variant Subsystems, either manually or automatically from balanced-realisation state relevance thresholds.
* **Three-file case compatibility (`.xlsx` + `.dyr` + `.json`).** Load a MATPOWER-FUBM workbook (`<case>.xlsx`) and VFlexP automatically picks up its companion PSS/E-style dynamic data file (`<case>.dyr`, including `$name = value` symbolic variable definitions) and case settings (`<case>.json`).
* **Time-domain disturbances, exact linearisation, and parametric sweeps.** Simulate step/pulse disturbances and line switching events, linearise the initialised Simulink model analytically (`linearize`) to compute eigenvalues, damping ratios, participation factors, and state relevance, and run dynamic (`keep` OP) or static (`recompute` PF + OP) parametric sweeps.

## Documentation

| | |
|---|---|
| [`docs/index.html`](docs/index.html) | User documentation site: open it in a browser (figures, equations, block diagrams) |
| [`docs/agent/README.md`](docs/agent/README.md) | Plain-Markdown documentation for AI agents and scripts: workflow pipeline, architecture, compatibility |
| [`docs/agent/formats.md`](docs/agent/formats.md) | Specification of `<case>.xlsx`, `<case>.dyr`, `<case>.json`, `$name` variables, and simulation events |
| [`docs/agent/models.md`](docs/agent/models.md) | Equations, per-unit bases, reference frames, complexity vocabulary, and device models |
| [`docs/agent/reference.md`](docs/agent/reference.md) | Reference of `VFlexPSession`, `core/`, `ui/`, `functions/`, `tools/`, and `Simulink/lib/` libraries |

The HTML documentation is also accessible online at [VFlexP Documentation](https://atomasmartin.github.io/VFlexP/).

## Directory Structure

* `vflexp.m`: Main entry point to launch the 6-tab `VFlexPApp` (optionally preloading a case file)
* `install.m`: Configures MATLAB paths (`core/`, `ui/`, `functions/`, `Simulink/`, `tools/`, `MATPOWER-FUBM/`) and checks required toolboxes
* `core/`: Decoupled state container (`VFlexPSession.m`) and pure MATLAB actions (`vflexp_load_case.m`, `vflexp_run_pf.m`, `vflexp_run_opf.m`, `vflexp_build_and_init.m`, `vflexp_simulate.m`, `vflexp_linearise.m`, `vflexp_modal.m`, `vflexp_relevance.m`, `vflexp_apply_relevance.m`, `vflexp_sweep.m`)
* `ui/`: Programmatic `.m` graphical interface (`VFlexPApp.m`, 6 tabs, sweep component, and configuration/inspection dialogs)
* `Simulink/VFlexP.slx`: Top-level Simulink simulation and linearisation shell with resolved links to `Simulink/lib/`
* `Simulink/lib/`: 31 modular Simulink libraries (`lib_primitives.slx` + 1 standalone `.slx` per generator, converter, grid, and control model)
* `Simulink/spec/`: Deterministic plain-text netlist mirrors (`*.spec.txt`) for clean `git diff` reviews of Simulink models
* `functions/`: `.dyr` and `.json` parsers, `param_genrou.m`, MATPOWER converter (`fun_MATPOWER_to_Excel.m`), custom linearisation hooks, and pure `.m` Simulink mask initialisers (`functions/masks/`)
* `tools/`: Developer utilities (`vflexp_export_specs.m`, `slx_textconv.py`, `vflexp_build_libraries.m`, `vflexp_rebuild_top_model.m`)
* `bus_systems/`: Sample hybrid AC/DC and AC power system cases
* `docs/`: HTML user documentation (`index.html`, `doc_*.html`) and Markdown documentation for AI agents (`docs/agent/`)

## Getting Started

a) **Online testing of the tool**

<a href="https://matlab.mathworks.com/open/github/v1?repo=atomasmartin/VFlexP" target="_blank">
  <img src="https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg" alt="Open in MATLAB Online">
</a>

b) **Local installation & quick start**
1. Clone this repository.
2. Ensure you have MATLAB R2024a or later installed (with Simulink, Control System Toolbox, Simulink Control Design, and Optimization Toolbox).
3. Run `install` in the root directory to configure paths and verify MATPOWER-FUBM.
4. Launch the interactive app:
   ```matlab
   vflexp('bus_systems/4buses_hybrid.xlsx')
   ```
   Or run headlessly from a MATLAB script:
   ```matlab
   sess = vflexp_load_case('bus_systems/4buses_hybrid.xlsx');
   vflexp_run_pf(sess);
   vflexp_build_and_init(sess);
   vflexp_simulate(sess, false);
   vflexp_modal(sess);
   ```
5. Refer to [`docs/index.html`](docs/index.html) or [`docs/agent/README.md`](docs/agent/README.md) for detailed usage instructions.

## Contributing

Bug reports and feature requests are welcome. If you encounter any issues or have suggestions for improvements, please open an issue on the GitHub repository.

## Acknowledgments

- Contributors to the MATPOWER-FUBM extension, required for the powerflow analysis implemented in the VFlexP tool. Full documentation of the MATPOWER-FUBM extension can be found in [MATPOWER-FUBM GitHub](https://github.com/AbrahamAlvarezB/matpower-fubm).
- The development of this tool was possible thanks to the advice and support of the researchers within the Modelling, Analysis and Control group of the Institute for Research in Technology (IIT) at Comillas Pontifical University, with special mention to Luis Rouco, Fidel Fernández-Bernal and Ignacio Egido Cortes.
- This work has been partially financed through the research program S2018/EMT-4366 PROMINT-CAM on Smart Grids of Madrid Government, Spain, with 50\% support from the European Social Fund (ESF), and through Grant TED2021-130610B-C22 funded by MICIU/AEI/10.13039/501100011033 and by European Union NextGenerationEU/PRTR.
- The collaboration between Durham University and Comillas has been eased with a special grant to Andrés Tomás-Martín from Comillas Pontifical University for temporary stays of researchers in foreign research centres.
- This work is also a collaboration with the project PID2021-125628OB-C21 funded by MICIU/AEI/10.13039/501100011033 and by ERDF/EU, and the Grant PRE2022-101606 funded by MICIU/AEI/10.13039/501100011033 and by ESF+.

## Citation

If you use this tool in your research, please cite the following paper:

A. Tomás-Martín, C. D. Zuluaga-Ríos, J. Suárez-Porras, B. Kazemtabrizi, J. García-Aguilar, L. Sigrist, and A. García-Cerrada, “A vector-based flexible-complexity tool for simulation and small-signal analysis of hybrid AC/DC power systems,” Sustainable Energy, Grids and Networks, vol. 43, p. 101817, Sep. 2025, [doi: 10.1016/j.segan.2025.101817](https://doi.org/10.1016/j.segan.2025.101817)

If you use the static analysis included in this tool, implemented using the MATPOWER-FUBM extension, please cite:

Abraham Alvarez-Bustos, Behzad Kazemtabrizi, Mahmoud Shahbazi, Enrique Acha-Daza, "Universal branch model for the solution of optimal power flows in hybrid AC/DC grids", International Journal of Electrical Power & Energy Systems, Volume 126, Part A, 2021, 106543, [doi: 10.1016/j.ijepes.2020.106543](https://doi.org/10.1016/j.ijepes.2020.106543)

R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems Research and Education," Power Systems, IEEE Transactions on, vol. 26, no. 1, pp. 12-19, Feb. 2011. [doi: 10.1109/TPWRS.2010.2051168](https://doi.org/10.1109/TPWRS.2010.2051168)

## Contact

For questions or support, please contact atomas@comillas.edu.
