# VFlexP for AI agents

This folder is the plain-text documentation of **VFlexP**, written to be read by a language model, an automated MATLAB script, or a person in a terminal. The HTML site in `../` (`index.html`, `doc_Home.html`) covers the same ground for people with figures and rendered equations; nothing here depends on it.

| file | what it holds |
|---|---|
| `README.md` | this file: what VFlexP is, how to drive it from MATLAB/Simulink (`core/` and `vflexp`), the workflow pipeline, Simulink `git diff` traceability, and the repository map |
| `formats.md` | the files of a case: network workbook (`.xlsx`), PSS/E-style dynamic data (`.dyr`), case settings (`.json`), `$name` symbolic variables, and temporal disturbance events |
| `models.md` | the equations of every Simulink library block, the complexity vocabulary, reference frames, and balanced-realisation state relevance |
| `reference.md` | every class and function in `core/`, `ui/`, `functions/`, and `tools/`, plus the 31 modular Simulink libraries in `Simulink/lib/` |

## What VFlexP is

VFlexP is a MATLAB/Simulink tool for the time-domain simulation, small-signal modal analysis, and state-relevance complexity reduction of hybrid AC/DC power systems. It shares its case formats and design principles with its Python counterpart (`VFlexPy`):

* **Vector-based modelling in Simulink.** Elements of the same type and complexity level are grouped into a single vectorised Simulink block whose parameters and signals are arrays/vectors, avoiding one-block-per-device overhead and enabling fast compilation and analytical linearisation.
* **Flexible complexity.** Every network element (AC bus, AC line, AC load, DC bus, DC line, DC load) and every generating unit or converter (`MGES`, `SG`, `VSM`, `GFL`, `GFM`, `DC_PES`, `intconv_GFL`, `intconv_VSM`, `intconv_GFM`, `DFIG`, `IG`, `SC`, `MSG`, `DCBUCK`) selects its own level of detail through expression-driven Simulink Variant Subsystems and dynamic routing matrices.
* **Decoupled engine + pure-`.m` App + immutable Simulink libraries.** All session state lives in `VFlexPSession` (`core/`), the 6-tab GUI (`ui/`) is 100% plain-text `.m` code, and `Simulink/VFlexP.slx` links (`LinkStatus = 'resolved'`) to 31 modular `.slx` libraries in `Simulink/lib/` reading `prm_*` and `op_*` structures from `ModelWorkspace` without ever mutating blocks at runtime (`Dirty == 'off'`).

## How an agent or script runs VFlexP

### 1. Headless / Scripted Workflow (Recommended for AI agents and batch runs)

From the repository root in MATLAB:

```matlab
% 1. Add VFlexP (core/, ui/, functions/, Simulink/, tools/) and MATPOWER-FUBM to the path
run('install.m');

% 2. Load a case (.xlsx + optional sibling .dyr and .json, or a MATPOWER .m file)
sess = vflexp_load_case('bus_systems/4buses_hybrid.xlsx');

% 3. Solve AC/DC Power Flow (or Optimal Power Flow via vflexp_run_opf(sess))
pfRes = vflexp_run_pf(sess);

% 4. Build routing matrices, compute analytical / findop equilibrium, and populate ModelWorkspace
initRes = vflexp_build_and_init(sess);

% 5. (Optional) Add a time-domain disturbance event and run Simulink simulation
sess.eventList = [sess.eventList; {true, 0.2, "AC Bus", "102", "dP_load", 10.0, 0.0}];
simRes = vflexp_simulate(sess, true);    % true = refresh workspace events & init

% 6. Linearise and run small-signal modal analysis
modRes = vflexp_modal(sess);             % calls vflexp_linearise(sess) automatically if needed

% 7. Compute balanced-realisation state relevance and auto-simplify low-relevance states
relRes  = vflexp_relevance(sess);
applied = vflexp_apply_relevance(sess, 1e-3);  % dryRun=false applies & re-initialises

% 8. Run a parametric sweep ('keep' = dynamic only, 'recompute' = re-run PF + init per step)
swpSpec = struct('kind', 'modal', 'mode', 'keep', ...
                 'targetType', 'unitParam', 'paramName', 'GFM.mp', ...
                 'values', linspace(0.01, 0.05, 5));
swpRes = vflexp_sweep(sess, swpSpec);
```

### 2. Interactive 6-Tab MATLAB App

```matlab
run('install.m');
app = vflexp('bus_systems/4buses_hybrid.xlsx');   % or simply vflexp()
```

Double-clicking the `VFlexP/Grid` block inside `Simulink/VFlexP.slx` also executes `OpenFcn = 'vflexp;'`, bringing up the singleton `VFlexPApp` window.

## The pipeline

State during an analysis session has a strict single source of truth:
1. **`VFlexPSession` (`handle` object in MATLAB)**: holds `tables` (`Bus`, `Line`, `Gen`, `IntConv`, `SysParam`, `GenCost`, `ParamVar`, `SecCtrl`, `TerCtrl`, `VoltSag`), `rawTables` (preserving `$name` strings), `variables` (`$name` dictionary), `variableRanges`, `dyrData`, `unitOverrides`, `config` (`powerflow`, `sim`, `modal`, `relevance`), `eventList`, and cached results (`pfResult`, `initResult`, `simResult`, `linResult`, `modalResult`, `relevanceResult`, `sweepModalDynamic`, `sweepModalStatic`, `sweepRelDynamic`, `sweepRelStatic`).
2. **`VFlexP` Simulink `ModelWorkspace` (read-only cache populated by `vflexp_build_and_init`)**: booted from `functions/masks/vflexp_ws_defaults.m` (`DataSource = 'MATLAB Code'`) so `VFlexP.slx` always compiles cleanly on cold load, and updated in memory via `assignin(mdlWks, ...)` with `Grid_par`, `DC`, `prm_<Unit>`, `op_<Unit>`, `ev_bus_dP`, `ev_bus_dQ`, `ev_dc_dI`, `ev_gen_dref`, and `t_close_lines`. Zero `set_param` calls are ever issued on internal library blocks.

| step | MATLAB function | reads from `VFlexPSession` | writes to `VFlexPSession` / `ModelWorkspace` |
|---|---|---|---|
| 1. Load case | `vflexp_load_case(sess, path)` | `<case>.xlsx` (+ optional `<case>.dyr`, `<case>.json`) | `tables`, `rawTables`, `variables`, `dyrData`, `unitOverrides`, `config` |
| 2. Save case | `vflexp_save_case(sess, [xlsxPath])` | `rawTables` / `tables`, `variables`, `config` | `<case>.xlsx` (+ optional `<case>.json`) |
| 3. Power flow | `vflexp_run_pf(sess)` | `tables`, `variables`, `config.powerflow` | `pfResult`, updated `tables.Bus`, `Gen`, `Line`, `IntConv` |
| 4. Optimal power flow | `vflexp_run_opf(sess)` | `tables`, `GenCost`, `config.powerflow` | `pfResult`, optimal dispatch in `tables.Bus`, `Gen`, `IntConv` |
| 5. Build & initialise | `vflexp_build_and_init(sess)` | solved `pfResult`, `tables`, `unitOverrides`, `eventList` | `initResult` in session; `Grid_par`, `DC`, `prm_*`, `op_*`, `ev_*` in `ModelWorkspace` |
| 6. Simulate | `vflexp_simulate(sess, [reinit])` | `config.sim`, `eventList`, initialised `ModelWorkspace` | `simResult` (`t`, `signals`, `rawOut`) via `Simulink.SimulationInput` |
| 7. Linearise | `vflexp_linearise(sess, [forceReinit])` | initialised `VFlexP.slx`, `config.modal` | `linResult` (`linsys`, `A`, `B`, `C`, `D`, mapped `stateNames`, `stateMeta`) |
| 8. Modal analysis | `vflexp_modal(sess, [forceRelin], [exportPath])` | `linResult`, `config.modal` | `modalResult` (`modesTable`, `eigenvalues`, `freqHz`, `dampingPct`, `participation`, `V`, `W`) |
| 9. State relevance | `vflexp_relevance(sess, [forceRelin])` | `linResult`, `modalResult`, `config.relevance` | `relevanceResult` (`rankTable`, `stateRelevance`, `hankelSV`, `planTable`) |
| 10. Auto-simplify | `vflexp_apply_relevance(sess, [thr], ..., [dryRun])` | `relevanceResult`, `tables`, `unitOverrides` | `plan` table; updates complexity strings & re-initialises when `dryRun == false` |
| 11. Parameter sweep | `vflexp_sweep(sess, sweepSpec, [progressCb])` | `variables` / `unitOverrides` / `tables` | `sweepModalDynamic`, `sweepModalStatic`, `sweepRelDynamic`, or `sweepRelStatic` |

## Simulink Modular Libraries & `git diff` Traceability

To prevent opaque binary `.slx` diffs and broken library links:
1. **31 Modular Libraries (`Simulink/lib/`)**: `lib_primitives.slx` (`LPF_opt`, `PI_ctrl`, `LeadLag_opt`, `WOF_opt`, `RL_dq`, `RC_dq`, `Park_transforms`) plus 1 standalone `.slx` library per generator (`lib_AC_*.slx`, `lib_DC_*.slx`), converter (`lib_IC_*.slx`), network (`lib_AC_Grid.slx`, `lib_DC_Grid.slx`), and secondary/tertiary controller (`lib_SecCtrl_*.slx`, `lib_TerCtrl_*.slx`).
2. **External Mask Code (`functions/masks/*.m`)**: Mask initialisation logic lives in plain-text `.m` functions (`mask_init_*.m`, `vflexp_ws_defaults.m`).
3. **Deterministic Netlist Mirrors (`Simulink/spec/*.spec.txt`)**: After editing any `.slx` library or `VFlexP.slx` in Simulink, run `vflexp_export_specs()` (in `tools/`) to regenerate sorted, position-independent text specs listing every block, mask parameter, variant condition, and port connection.
4. **Automatic `git diff` Textconv (`.gitattributes` + `tools/slx_textconv.py`)**: Even before running `vflexp_export_specs()`, `git diff` renders structural changes inside `.slx` archives cleanly in plain text.

## Compatibility between VFlexP and VFlexPy

VFlexP and VFlexPy share the same case files (`.xlsx`, `.dyr`, `.json`) and model definitions:

* **Case files**: Both tools read the same MATPOWER-FUBM `.xlsx` sheets (`Bus data`, `Line data`, `Gen data`, `IntConv data`, `System parameters`), the same PSS/E-style `.dyr` records (`GENROU`, `GENSAL`, `ST4B`, `IEEEG1`, `PSS2B`, `GFL`, `VSM`, `GFM`, `DC_PES`, `INTCONV_GFL`, `INTCONV_VSM`, `INTCONV_GFM`, `DCBUS`, `DCLINE`), the same `$name = value` variable definitions, and the same `<case>.json` settings.
* **Per-unit base scaling**:
  * In `.dyr` files, device parameters are given on the device's own MVA base (`mBase` from `Gen data` for AC/DC generators, or `rateA` from `IntConv data` for interface converters, falling back to `Sb` when `0` or omitted).
  * `MGES` and `SG` blocks operate internally on `mBase` (converting currents and powers at their Simulink ports by `mBase./Sb` and `Sb./mBase`), whereas converter blocks operate on the system base `Sb`. `vflexp_load_case` / `vflexp_build_and_init` automatically convert `.dyr` converter impedances, admittances, inertias, and droop gains from `mBase` (or `rateA`) to `Sb`.
* **Synchronous machine conversion (`param_genrou.m`)**: Both tools support `CONVERSION=CLASSICAL` and `CONVERSION=EXACT` (IEEE Std 1110 two-circuit operational-inductance root matching), plus direct rotor circuit overrides (`Lfd`, `Rfd`, `Lkd1`, `Rkd1`, `Lkq1`, `Rkq1`, `Lkq2`, `Rkq2`).
* **Complexity aliases**: Both VFlexP (`order 8`, `order 4`, `dynamic`, `static omega`, `series RL static omega`, `DC dynamic`, `DC static`, `DC ConstR`, `DC ConstI`) and VFlexPy (`EMT8`, `RMS4`, `EMT`, `RMS`, `DC EMT`, `DC RMS`, `DC static omega`) complexity strings are accepted seamlessly.

## Code map

```
install.m                        configures MATLAB paths and checks toolboxes
vflexp.m                         launches or focuses the singleton 6-tab VFlexPApp
core/
  VFlexPSession.m                handle class holding all session state, tables, config, and results
  vflexp_load_case.m             loads .xlsx + sibling .dyr / .json (or MATPOWER .m) into VFlexPSession
  vflexp_save_case.m             writes session tables and settings to .xlsx (+ optional .json)
  vflexp_run_pf.m                solves hybrid AC/DC power flow via MATPOWER-FUBM
  vflexp_run_opf.m               solves hybrid AC/DC optimal power flow via MATPOWER-FUBM
  vflexp_build_events_ws.m       compiles session disturbance events into Simulink workspace vectors
  vflexp_build_and_init.m        builds Grid_par, DC, prm_*, and analytical/findop op_* in ModelWorkspace
  vflexp_simulate.m              executes Simulink time-domain simulation via SimulationInput
  vflexp_map_state_names.m       maps Simulink block state paths to human-readable device/bus names
  vflexp_linearise.m             linearises VFlexP.slx analytically via Simulink Control Design
  vflexp_modal.m                 computes eigenvalues, frequencies, damping, and participation factors
  vflexp_relevance.m             computes balanced-realisation state relevance (Tomás-Martín et al., 2023)
  vflexp_apply_relevance.m       simplifies device/network complexity below a state-relevance threshold
  vflexp_sweep.m                 runs 'modal' or 'relevance' sweeps in 'keep' or 'recompute' PF modes
ui/
  VFlexPApp.m                    main uifigure shell (toolbar, 6 tabs, singleton lifecycle)
  StaticTab.m                    Tab 1: editable network sheets, Run PF, Run OPF, summary
  DiagramTab.m                   Tab 2: interactive AC/DC topology graph and element inspector
  RelevanceTab.m                 Tab 3: single state-relevance chart/table, Set complexity, dynamic/static sweeps
  ComplexityTab.m                Tab 4: per-element complexity editor, presets, state-count breakdown
  SimulationTab.m                Tab 5: Play simulation, disturbance event table, Simulink launcher, scope grid
  ModalTab.m                     Tab 6: single linearisation (Damping-Freq / S-plane, participations, export) + sweeps
  SweepTab.m                     reusable parametric sweep sub-tab for Modal and State Relevance
  VariablesDialog.m              $name symbolic variable editor
  DeviceDialog.m                 per-unit dynamic parameter (.dyr / prm_*) editor
  LogDialog.m                    session log viewer
  StaticSettingsDialog.m         PF/OPF solver settings dialog
  DynamicSettingsDialog.m        Simulink horizon, solver, and tolerance settings dialog
  ModalSettingsDialog.m          I/O linearisation mode and modal filter settings dialog
  RelevanceSettingsDialog.m      state-relevance I/O mode and threshold settings dialog
Simulink/
  VFlexP.slx                     top-level simulation & linearisation shell (30 resolved library links)
  lib/                           31 modular .slx libraries (primitives, generators/, converters/, grid/, controls/)
  spec/                          32 deterministic .spec.txt netlists mirroring VFlexP.slx and Simulink/lib/
functions/
  fun_read_dyr.m                 parses PSS/E-style .dyr files and $name variables
  fun_read_case_json.m           reads <case>.json solver and power-flow settings
  param_genrou.m                 converts GENROU/GENSAL standard parameters to per-unit circuit values
  fun_MATPOWER_to_Excel.m        converts a MATPOWER case (.m or struct) to a VFlexP .xlsx workbook
  linearizedatan/cos/sin.m       exact custom linearisation callbacks for Simulink trigonometric blocks
  fun_MAS_tc_chart.m             multi-agent tertiary control chart helper
  masks/                         24 pure-.m mask initialisation functions and vflexp_ws_defaults.m
tools/
  vflexp_export_specs.m          exports deterministic .spec.txt netlists for all .slx models/libraries
  slx_textconv.py                Python git textconv driver for *.slx diffs (.gitattributes)
  vflexp_build_libraries.m       programmatic builder/verifier for Simulink/lib/*.slx
  vflexp_rebuild_top_model.m     programmatic builder/verifier for Simulink/VFlexP.slx
bus_systems/                     sample hybrid AC/DC and AC power system cases
docs/
  index.html, doc_*.html         HTML user documentation
  agent/                         Markdown documentation for AI agents and scripts
```
