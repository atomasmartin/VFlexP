# VFlexP Function, Class, and Library Reference

Complete reference of the entry points, `core/` session & action engine, `ui/` classes, `functions/` parsers & mask helpers, `tools/` spec exporters, and the 31 modular Simulink libraries in `Simulink/lib/`.

---

## 1. Entry Points (`install.m` and `vflexp.m`)

### `install()`
Adds `core/`, `ui/`, `functions/`, `functions/masks/`, `Simulink/`, `Simulink/lib/` (and subfolders), `tools/`, and `MATPOWER-FUBM/` to the MATLAB search path and verifies required toolboxes (`Simulink`, `Control System Toolbox`, `Simulink Control Design`, `Optimization Toolbox`).

### `appOut = vflexp([casePath])`
Launches or focuses the singleton `VFlexPApp` window.
* **Inputs**:
  * `casePath` (optional): path to an `.xlsx` workbook or MATPOWER `.m` case file. When supplied, loads the case into the session (and runs an initial power flow) whether launching a new window or reusing an already open `VFlexPApp` instance.
* **Outputs**: returns `app` (`VFlexPApp` handle) when `nargout > 0`.

---

## 2. Central Session State (`core/VFlexPSession.m`)

`VFlexPSession < handle` is the single in-memory source of truth for a study. Zero dependency on `gcb`, `bdroot`, or Simulink dialog masks.

### Key Properties
* `casePath`, `caseName`: active case file path and stem.
* `rawTables`: struct of raw string/cell tables (`Bus`, `Line`, `Gen`, `IntConv`, `SysParam`, `GenCost`, `ParamVar`, `SecCtrl`, `TerCtrl`, `VoltSag`) preserving symbolic `$name` expressions.
* `tables`: struct of numerically evaluated MATLAB tables after resolving `session.variables`.
* `variables`: struct mapping each symbolic variable name (without `$`) to its current scalar value.
* `variableRanges`: 5-column table (`Parameter`, `Value`, `From`, `To`, `N_steps`) for parameter sweeps.
* `dyrData`: parsed `<case>.dyr` structure from `fun_read_dyr` (or `[]`).
* `unitOverrides`: struct of interactive or programmatic overrides for device parameters and complexity flags (`unitOverrides.<UnitType>.<FieldName>`).
* `config`: nested settings struct:
  * `config.powerflow`: `tolerance`, `max_iterations`, `enforce_q_limits`, `write_branch`.
  * `config.sim`: `t_end`, `method` (`'ode15s'`, `'ode23tb'`, etc.), `rtol`, `atol`, `max_step`, `sc_structure`, `define_sc_in_matlab`, `tc_structure`.
  * `config.modal`: `input_channels`, `output_channels`, `sort_by`, `min_freq`.
  * `config.relevance`: `input_channels`, `output_channels`, `stability_margin`, `threshold`.
* `eventList`: 7-column disturbance table (`Enabled`, `Time`, `TargetType`, `TargetId`, `Channel`, `Magnitude`, `Duration`).
* Cached results: `pfResult`, `initResult`, `simResult`, `linResult`, `modalResult`, `relevanceResult`, `sweepModalDynamic`, `sweepModalStatic`, `sweepRelDynamic`, `sweepRelStatic`.
* `logs`: string array of timestamped session messages.

### Key Methods
* `obj = VFlexPSession([casePath])`: constructs a clean session (and calls `vflexp_load_case` if `casePath` is provided).
* `obj.resetDefaults()`: resets tables, config, events, and cached results.
* `obj.evalTables()`: resolves `$name` expressions in `rawTables` into numeric `tables`.
* `val = obj.evalScalar(expr, fallback)`: evaluates a numeric scalar or `$name` algebraic expression against `obj.variables`.
* `[TbusAC, TgenAC, TlineAC, TbusDC, TgenDC, TlineDC, Tintconv, TFb, Sb] = obj.splitDomains()`: splits hybrid network tables into AC and DC domains.
* `obj.log(fmt, ...)`: appends a formatted message to `obj.logs`.

---

## 3. Core Actions (`core/vflexp_*.m`)

### `session = vflexp_load_case(sessionOrPath, [casePath])`
Loads `<case>.xlsx` (plus optional sibling `<case>.dyr` and `<case>.json`) or a MATPOWER `<case>.m` file into `VFlexPSession`, populating `rawTables`, `tables`, `variables`, `variableRanges`, `dyrData`, and `config`.

### `savePath = vflexp_save_case(session, [savePath])`
Writes all session sheets (preserving `$name` formulas from `rawTables`) to `savePath` (`.xlsx`) and updates sibling `<stem>.json` solver/powerflow settings when present.

### `pfRes = vflexp_run_pf(session)`
Converts `session.tables` into a MATPOWER-FUBM `mpc` struct, ensures at least one slack bus per synchronous AC area, executes `runpf(mpc, mpopt)`, unpacks solved bus voltages, generator powers, branch flows, and converter `alpha1` losses back into `session.tables`, and stores `session.pfResult`.

### `opfRes = vflexp_run_opf(session)`
Builds the MATPOWER-FUBM `mpc` struct including polynomial/piecewise costs from `session.tables.GenCost` (generating default quadratic costs when omitted), solves `runopf(mpc, mpopt)`, updates optimal setpoints (`Vm`, `Va`, `Pg`, `Qg`, converter `PF`/`QF`/`VF_SET`/`VT_SET`), and stores `session.pfResult`.

### `evWs = vflexp_build_events_ws(session)`
Compiles enabled rows of `session.eventList` into `FromWorkspace` structures (`ev_bus_dP`, `ev_bus_dQ`, `ev_dc_dI`, `ev_gen_dref`) with sharp piecewise-constant step/pulse edges and AC line closing times (`t_close_lines`).

### `initRes = vflexp_build_and_init(session)`
Builds AC/DC routing and incidence structures (`Grid_par`, `DC`), constructs vectorised parameter structs (`prm_<Unit>`) scaled to `Sb` (or `mBase` for `MGES`/`SG`) with complexity codes, solves steady-state initial conditions analytically for network elements and converters and via isolated single-block `findop` for synchronous/induction machines, and writes `Grid_par`, `DC`, `prm_*`, `op_*`, `ev_*`, and `t_close_lines` into the `ModelWorkspace` of `Simulink/VFlexP.slx` without mutating any Simulink block (`Dirty == 'off'`, `30/30` library links `'resolved'`).

### `simRes = vflexp_simulate(session, [reinit])`
Runs a non-linear time-domain simulation of `Simulink/VFlexP.slx` using `Simulink.SimulationInput` configured with `session.config.sim` (`StopTime`, `Solver`, `RelTol`, `AbsTol`, `MaxStep`), unpacks logged AC/DC bus voltages, line currents, generator/converter powers, frequencies, and currents into `session.simResult`.

### `[labels, stateMeta] = vflexp_map_state_names(session, rawStateNames)`
Maps vectorised Simulink state paths and element indices back to physical device/bus/branch identifiers (`category`, `elementType`, `elementId`, `busId`, `subGroup`, `stateVar`, and human-readable `label`).

### `linRes = vflexp_linearise(session, [forceReinit])`
Temporarily zeroes exogenous event waveforms, evaluates exact block-by-block analytical linearisation (`linearize('VFlexP', io)`) at the initialised equilibrium, maps all state names via `vflexp_map_state_names`, and stores `(linsys, A, B, C, D, stateNames, stateMeta)` in `session.linResult`.

### `modalRes = vflexp_modal(session, [forceRelin], [exportPath])`
Computes right/left eigenvectors (`V`, `W`), eigenvalues (`eigenvalues`), natural frequencies (`freqHz`), damping ratios (`dampingPct`), normalised participation matrix (`participation`), and summary `modesTable`, storing `session.modalResult` and optionally exporting a multi-sheet `.xlsx` or `.mat` report to `exportPath`.

### `relRes = vflexp_relevance(session, [forceRelin])`
Computes the normalised input-output state relevance vector (`stateRelevance`, normalised so its maximum entry is `1`) over the selected `input_channels` and `output_channels` by combining balanced realisation Hankel singular values (`hankelSV` via `balreal` / Lyapunov Gramians, shifting unstable/integrator poles by `stability_margin` when required) with original and balanced participation factors (Tomás-Martín et al., *IJEPES*, 2023), returning `rankTable` and `planTable`.

### `out = vflexp_apply_relevance(session, [threshold], [stateRel], [stateMeta], [dryRun])`
Plans (`dryRun = true`) or applies (`dryRun = false`, default) automatic complexity reductions whenever **all** states belonging to an algebraic-capable sub-block have relevance strictly below `threshold` (e.g. `order 8 -> order 4 -> order 2`, `Dynamic VCC -> Static VCC`, `Dynamic CC/PLL -> Static CC/PLL`, dynamic lines/loads -> `static omega`, and dynamic AC/DC buses -> static when free of algebraic-loop blockers), then re-runs `vflexp_build_and_init(session)`.

### `sweepRes = vflexp_sweep(session, sweepSpec, [progressCb])`
Executes a parametric modal (`sweepSpec.kind = 'modal'`) or state-relevance (`sweepSpec.kind = 'relevance'`) sweep on a cloned session so the caller's base tables remain uncorrupted:
* `sweepSpec.mode`: `'keep'` (dynamic sweep keeping PF fixed) or `'recompute'` (re-running PF + unit initialisation at every step).
* `sweepSpec.targetType`: `'variable'` (`$name`), `'unitParam'` (e.g. `'GFM.mp'` or `'VSM.Ta(1)'`), `'loadScale'` (`'Pd_Qd'`), or `'tableCell'` (e.g. `'Gen.Pg(1)'`).
* `sweepSpec.values` (or `from`, `to`, `nSteps`): swept parameter vector, with Hungarian/greedy eigenvalue tracking across steps.

---

## 4. Parsers, Converters & Mask Helpers (`functions/`)

| file | purpose |
|---|---|
| `fun_read_dyr(dyr_file)` | Parses PSS/E-style `.dyr` files, `$name = value` declarations (slash-terminated or standalone line), positional `GENROU`/`GENSAL` records, named `KEY=VALUE` pairs, and `COMPLEXITY="..."` strings. |
| `fun_read_case_json(json_file)` | Reads `<case>.json` (`t_end`, `solver`, `powerflow`, and optional `variables` dictionary). |
| `param_genrou(...)` | Converts `GENROU`/`GENSAL` reactances and time constants into per-unit rotor/stator circuit parameters (`'classical'` or IEEE Std 1110 `'exact'`, plus explicit circuit overrides). |
| `fun_MATPOWER_to_Excel(mpc, [filename])` | Standalone converter from a MATPOWER/FUBM case struct or `.m` file into a VFlexP `.xlsx` workbook. |
| `linearizedatan.m`, `linearizedcos.m`, `linearizedsin.m` | Exact analytical Jacobian hooks attached to trigonometric blocks inside `Simulink/lib/`. |
| `fun_MAS_tc_chart.m` | Helper function invoked by `Simulink/lib/controls/lib_TerCtrl_MAS.slx`. |
| `functions/masks/vflexp_ws_defaults.m` | Default `ModelWorkspace` initialisation script (`DataSource = 'MATLAB Code'`) for `Simulink/VFlexP.slx`. |
| `functions/masks/mask_init_*.m` | 23 pure-`.m` Simulink mask initialisation functions backing the modular libraries in `Simulink/lib/`. |

---

## 5. Pure-`.m` Graphical Interface (`ui/`)

| class | role |
|---|---|
| `VFlexPApp` | Main `uifigure` window with top toolbar (*Load*, *Save*, *Variables*, *Settings*, *Open Simulink*, *Log*), 6-tab `uitabgroup`, and status bar. |
| `StaticTab` | Sub-tabs for `Bus`, `Line`, `Gen`, `IntConv`, `SysParam`, `GenCost`, `ParamVar`, *Run PF*, *Run OPF*, and voltage/loss summary. |
| `DiagramTab` | Interactive AC/DC network graph (`force`, `layered`, `circle`, coordinate layouts), colour-mapped by `Vm`, `Va`, or `Domain`, with click-to-inspect bus/branch details. |
| `RelevanceTab` | *Single analysis* (log bar chart & state table, threshold line, *Preview* & *Set complexity* button), *Sweep: dynamic model*, and *Sweep: static model*. |
| `ComplexityTab` | Per-element complexity tables (AC/DC buses, lines, loads, generators, interface converters), global `EMT` / `RMS` presets, state-count summary, and *Edit Params* launcher (`DeviceDialog`). |
| `SimulationTab` | *Initialise*, *Play* (`vflexp_simulate`), disturbance event editor (`session.eventList`), *Open Simulink*, channel selector, and multi-axes interactive scope grid. |
| `ModalTab` | *Single linearisation* (*Damping vs Frequency* & *Complex S-Plane* plots, sortable modes table, participation factors bar/table, `.xlsx`/`.mat` export), *Sweep: dynamic model*, and *Sweep: static model*. |
| `SweepTab` | Reusable sweep UI component driving `vflexp_sweep` for both `ModalTab` and `RelevanceTab`. |
| `VariablesDialog`, `DeviceDialog`, `LogDialog`, `StaticSettingsDialog`, `DynamicSettingsDialog`, `ModalSettingsDialog`, `RelevanceSettingsDialog` | Modal/non-modal helper dialogs in `ui/`. |

---

## 6. Modular Simulink Libraries (`Simulink/lib/`) & Git Spec Tools (`tools/`)

### Library Architecture (31 `.slx` files)
* **`Simulink/lib/lib_primitives.slx`**: reusable building blocks with parameterised initial state `x0` and expression-controlled bypass variants (`LPF_opt`, `PI_ctrl`, `LeadLag_opt`, `WOF_opt`, `RL_dq`, `RC_dq`, `Park_transforms`).
* **`Simulink/lib/generators/` (11 libraries)**: `lib_AC_MGES.slx`, `lib_AC_SG.slx`, `lib_AC_VSM.slx`, `lib_AC_GFL.slx`, `lib_AC_GFM.slx`, `lib_AC_DFIG.slx`, `lib_AC_IG.slx`, `lib_AC_SC.slx`, `lib_AC_MSG.slx`, `lib_DC_DC_PES.slx`, `lib_DC_DCBUCK.slx`.
* **`Simulink/lib/converters/` (6 libraries)**: `lib_IC_GFL.slx`, `lib_IC_VSM.slx`, `lib_IC_GFM.slx`, `lib_Conv_GFL.slx`, `lib_Conv_VSM.slx`, `lib_Conv_GFM.slx`.
* **`Simulink/lib/grid/` (2 libraries)**: `lib_AC_Grid.slx`, `lib_DC_Grid.slx`.
* **`Simulink/lib/controls/` (11 libraries)**: `lib_SecCtrl_NoSC.slx`, `lib_SecCtrl_Matrices.slx`, `lib_SecCtrl_Matrices1.slx`, `lib_SecCtrl_MatricesDelay.slx`, `lib_SecCtrl_MatricesDelayLim.slx`, `lib_SecCtrl_PilotNodes.slx`, `lib_SecCtrl_MAS.slx`, `lib_TerCtrl_NoTC.slx`, `lib_TerCtrl_Scheduled.slx`, `lib_TerCtrl_MAS.slx`, `lib_TerCtrl_MAS_edp.slx`.

Each active device slot in `Simulink/VFlexP.slx` maintains `LinkStatus = 'resolved'` to its library block and receives two `ModelWorkspace` structs on its mask: `prm_<Unit>` (physical parameters & numeric variant codes) and `op_<Unit>` (initial integrator/IC states).

### Spec & Diff Tooling (`tools/`)
* `outFiles = vflexp_export_specs([projRoot])`: exports `Simulink/VFlexP.slx` and all 31 `Simulink/lib/**/*.slx` libraries into deterministic, position-independent text netlists in `Simulink/spec/*.spec.txt`.
* `tools/slx_textconv.py`: Python script registered in `.gitattributes` (`*.slx diff=slx`) so `git diff` displays readable block/parameter/connection diffs directly from `.slx` ZIP archives.
