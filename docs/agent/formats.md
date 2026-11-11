# File Formats in VFlexP

A case in VFlexP consists of up to three files sharing the same stem:

| file | required | what it holds |
|---|---|---|
| `<case>.xlsx` | **yes** | network topology, steady-state dispatch, network complexity choices, and optional control/sweep tables |
| `<case>.dyr` | optional | PSS/E-style dynamic device records and `$name = value` symbolic variable definitions |
| `<case>.json` | optional | case settings: simulation horizon (`t_end`), Simulink solver tolerances/step, MATPOWER power-flow options, and optional `variables` |

When `vflexp_load_case('<case>.xlsx')` is called (or loaded via the *Load* button in `vflexp`), VFlexP automatically looks for `<case>.dyr` and `<case>.json` in the same folder and applies them into `VFlexPSession`. If no `.dyr` file is present, VFlexP uses the default unit parameter templates from `Simulink/lib/vflexp_lib_manifest.mat` (editable in the *Complexity* tab via `DeviceDialog` or `session.unitOverrides`). VFlexP can also load standard MATPOWER `.m` case files directly (`vflexp_load_case('case9.m')`).

---

## 1. The Network Workbook (`<case>.xlsx`)

An Excel workbook with sheets following the MATPOWER-FUBM convention plus VFlexP complexity columns.

### Required Sheets

#### `Bus data`
One row per AC or DC bus (`session.tables.Bus`):

| column | unit | meaning |
|---|---|---|
| `bus_i` | integer | unique bus number |
| `type` | `1`, `2`, `3`, `4` | `1` = PQ (or DC voltage-droop/power bus), `2` = PV (or DC voltage-controlled bus), `3` = Slack, `4` = isolated |
| `Pd`, `Qd` | MW, MVAr | active and reactive power demand at the bus (`Qd = 0` for DC buses) |
| `Gs`, `Bs` | MW, MVAr | shunt conductance and susceptance at `Vm = 1.0` pu (on DC buses, positive `Bs` specifies DC shunt capacitance when `.dyr` `DCBUS` is not provided) |
| `area` | integer | synchronous area number (`1, 2, ...`); DC grids usually carry their own area number |
| `Vm`, `Va` | pu, deg | voltage magnitude and angle (may be numbers or `$name` variable expressions) |
| `baseKV` | kV | nominal base voltage |
| `zone` | integer | loss zone |
| `Vmax`, `Vmin` | pu | upper and lower voltage limits for OPF |
| `Complexity_bus` | text | dynamic model of the bus shunt capacitance: `dynamic` (`EMT`), `static omega` (`RMS`), `DC dynamic` (`DC EMT`), or `DC static` (`DC RMS` / `DC static omega`) |
| `Complexity_load` | text | dynamic model of the load at the bus: `series RL dynamic` / `dynamic` (`EMT`), `series RL static omega` (`RMS`), `parallel RL dynamic`, `parallel RL static omega`, `dynamic PQ`, `PQ`, `DC DynamicP`, `DC ConstP`, `DC ConstR`, or `DC ConstI` |
| `Latitude`, `Longitude` | deg (optional) | geographic/layout coordinates used by `DiagramTab` when present |

#### `Line data`
One row per AC transmission line, transformer, or DC cable (`session.tables.Line`):

| column | unit | meaning |
|---|---|---|
| `fbus`, `tbus` | integer | from-bus and to-bus numbers |
| `r`, `x`, `b` | pu | series resistance, series reactance, and total charging susceptance on the area base `Sb` (`x = 0`, `b = 0` for DC lines in PF; positive `x` or `.dyr` `DCLINE` provides DC line inductance) |
| `rateA`, `rateB`, `rateC` | MVA | long-term, short-term, and emergency thermal ratings (`0` = unlimited in PF; enforced in OPF) |
| `ratio`, `angle` | pu, deg | off-nominal transformer tap ratio (`0` or `1` for lines) and phase-shift angle |
| `status` | `1` or `0` | `1` = in service, `0` = out of service |
| `angmin`, `angmax` | deg | minimum and maximum voltage angle difference |
| `Complexity` | text | `dynamic` (`EMT`), `static omega` (`RMS`), `static`, `DC dynamic` (`DC EMT`), or `DC static` (`DC RMS` / `DC static omega`) |

#### `Gen data`
One row per generating unit (`session.tables.Gen`), either AC or DC:

| column | unit | meaning |
|---|---|---|
| `bus` | integer | bus number to which the generator is connected |
| `type` | text | Simulink generator model: `MGES`, `SG`, `VSM`, `GFL`, `GFM`, `DFIG`, `IG`, `SC`, `MSG`, `DC_PES`, `DCBUCK`, or `inf` (infinite bus). If `.dyr` records exist for this bus, `vflexp_load_case` sets `type` automatically (`GENROU`/`GENSAL` -> `MGES`, etc.) |
| `Pg`, `Qg` | MW, MVAr | active and reactive power output (may be numbers or `$name` variable expressions) |
| `Qmax`, `Qmin` | MVAr | reactive power limits |
| `Vg` | pu | voltage magnitude setpoint |
| `mBase` | MVA | machine nominal apparent power base |
| `status` | `1` or `0` | `1` = in service, `0` = out of service |
| `Pmax`, `Pmin` | MW | active power limits |
| `Pc1..apf` | various | standard MATPOWER capability curve and ramp columns |

#### `System parameters`
One row per synchronous area (`session.tables.SysParam`):

| column | unit | meaning |
|---|---|---|
| `Area` | integer | area identifier matching `Bus data.area` |
| `Fb` | Hz | nominal base frequency of the area (e.g. `50` or `60`) |
| `Sb` | MVA | base apparent power of the area (e.g. `100`) |

### Optional Sheets

* **`IntConv data`** (`session.tables.IntConv`): MATPOWER-FUBM VSC interface-converter branches between a DC bus (`fbus`) and an AC bus (`tbus`). The first column `type` selects the Simulink converter model (`intconv_GFL`, `intconv_VSM`, or `intconv_GFM`), followed by the 36 FUBM branch columns (`fbus`, `tbus`, `r`, `x`, `b`, `rateA`, `rateB`, `rateC`, `ratio`, `angle`, `status`, `angmin`, `angmax`, `PF`, `QF`, `PT`, `QT`, `MU_SF`, `MU_ST`, `MU_ANGMIN`, `MU_ANGMAX`, `VF_SET`, `VT_SET`, `TAP_MAX`, `TAP_MIN`, `CONV`, `BEQ`, `K2`, `BEQ_MIN`, `BEQ_MAX`, `SH_MIN`, `SH_MAX`, `GSW`, `ALPHA1`, `ALPHA2`, `ALPHA3`).
* **`Gen cost data`** (`session.tables.GenCost`): MATPOWER generator cost polynomials (`model`, `startup`, `shutdown`, `n`, `c2`, `c1`, `c0`) used by `vflexp_run_opf`.
* **`Parameter variation`** (`session.tables.ParamVar` / `session.variableRanges`): 5-column table (`Parameter`, `Default_value`/`Value`, `From`, `To`, `N_steps`) defining symbolic variables and their sweep ranges.
* **`Secondary control`** (`session.tables.SecCtrl`) & **`Tertiary control`** (`session.tables.TerCtrl`): configuration of secondary and tertiary control structures.
* **`Voltage sag`** (`session.tables.VoltSag`): fault/sag sequence table (`Bus`, `Time (s)`, `Duration (s)`, `Depth (p.u.)`).

---

## 2. Dynamic Data File (`<case>.dyr`)

A plain-text PSS/E-style file holding variable definitions and one slash-terminated record (`... /`) per device or control module. Comments start with `//` or `@!`. Records may span multiple lines until `/`.

### Variable Definitions (`$name = value`)

At the top of the `.dyr` file (or between records), named symbolic parameters can be declared either with or without a trailing `/`:

```
$D_vsm = 0.01 /
$Pf_set = -40
```

* Every `$name = value` declaration registers `session.variables.name = value` (and populates `session.variableRanges` when not already specified in Excel).
* Any `.dyr` parameter value or workbook cell (`Vm`, `Va`, `Pg`, `Qg`, `r`, `x`, etc.) may reference `$name` (including arithmetic expressions such as `-1*$Pf_set`).
* Editing `session.variables` (via `VariablesDialog`, scripts, or `vflexp_sweep`) and calling `session.evalTables()` + `vflexp_build_and_init(session)` propagates updated values everywhere without touching `.slx` files.

### Record Syntax

Each record begins with:

```
<key> '<MODEL>' <id> <parameters...> /
```

* `<key>`:
  * Integer `<bus>` for generators, exciters, governors, stabilisers, and DC bus records (e.g. `1 'GENROU' 1 ... /`).
  * Hyphenated `<fbus>-<tbus>` for AC/DC interface converters and DC lines (e.g. `58-4 'INTCONV_GFL' 1 ... /` or `'58-4' 'INTCONV_GFL' 1 ... /`).
* `<MODEL>`: case-insensitive string in single quotes (`'GENROU'`, `'GENSAL'`, `'ST4B'`, `'IEEEG1'`, `'PSS2B'`, `'GFL'`, `'VSM'`, `'GFM'`, `'DC_PES'`, `'INTCONV_GFL'`, `'INTCONV_VSM'`, `'INTCONV_GFM'`, `'DCBUS'`, `'DCLINE'`).
* `<id>`: unit/circuit identifier (typically `1` or `'1'`).
* `<parameters...>`:
  * **Positional numeric values** (standard PSS/E order for `GENROU` and `GENSAL`), optionally followed by named `KEY=VALUE` pairs.
  * **Named `KEY=VALUE` pairs** for all converter models, control models (`ST4B`, `IEEEG1`, `PSS2B`), and DC network records (`DCBUS`, `DCLINE`), as well as optional overrides on `GENROU`/`GENSAL`.
  * **`COMPLEXITY="..."`**: semicolon-separated key-value string selecting the internal complexity variants of the device (e.g. `COMPLEXITY="MODEL=EMT8"`, `COMPLEXITY="MODEL=RMS4"`, `COMPLEXITY="CC=STATIC;PLL=DYNAMIC"`, `COMPLEXITY="VCC=DYNAMIC;REF=DC_ATM"`).

### Supported `.dyr` Models

#### Synchronous Machine (`GENROU` / `GENSAL`) -> Simulink `MGES`
Positional parameters (PSS/E order):
* `GENROU` (14 numbers): `Tpdo, Tppdo, Tpqo, Tppqo, H, D, Xd, Xq, Xpd, Xpq, Xppd, Xl, S10, S12`
* `GENSAL` (12 numbers): `Tpdo, Tppdo, Tppqo, H, D, Xd, Xq, Xpd, Xppd, Xl, S10, S12`

Optional named parameters on `GENROU`/`GENSAL`:
* `Ra` or `Rs`: stator armature resistance (pu on `mBase`, default `0.0025`).
* `CONVERSION`: `CLASSICAL` (default) or `EXACT` (IEEE Std 1110 operational-inductance root match via `param_genrou.m`).
* Direct rotor circuit overrides (pu on `mBase`): `Lfd`, `Rfd`, `Lkd1`, `Rkd1`, `Lkq1`, `Rkq1`, `Lkq2`, `Rkq2`.
* `cosn`: nominal power factor (default `0.85`).
* `COMPLEXITY`: `MODEL=EMT8` (or `order 8`, default) or `MODEL=RMS4` (or `order 4`).

#### Excitation System (`ST4B`) -> Simulink `MGES` AVR
* `T_r` (s, measurement filter time constant -> `T_meas`)
* `K_pr`, `K_ir` (pu, voltage regulator PI gains), `K_pm`, `K_im` (pu, inner bridge gains)
* `T_a` (s, bridge time constant -> `T_f`), `V_rmax`, `V_rmin` (pu -> `E_fdmax`, `E_fdmin`)

#### Steam/Hydro Turbine Governor (`IEEEG1`) -> Simulink `MGES` Governor
* `K` (pu, droop gain $= 1/R$), `T_1`, `T_2`, `T_3`, `T_4`, `T_5`, `T_6`, `T_7` (s), `K_1..K_8` (pu, shaft fractions), `U_o`, `U_c` (pu/s, valve rate limits), `P_max`, `P_min` (pu, turbine limits).

#### Dual-Input Power System Stabiliser (`PSS2B`) -> Simulink `MGES` PSS
* `K_s1`, `K_s2`, `K_s3`, `T_w1..T_w4`, `T_1..T_4`, `T_6..T_11`, `M`, `N`, `V_stmax`, `V_stmin`.

#### Grid-Following Generator (`GFL`) & Interface Converter (`INTCONV_GFL`)
Named parameters (on `mBase` / `rateA` in `.dyr`, scaled to `Sb` in `prm_*`):
* `Lf`, `Rf`, `Cf`, `kp_i`, `ki_i`, `kp_pll`, `ki_pll`, `Ta` (or `Tf`), `Rp`, `Kq`
* On `INTCONV_GFL`: `Kpvdc`, `Kivdc` (DC-voltage PI/droop gains -> `K_p`, `K_i`)
* `COMPLEXITY`: `CC=DYNAMIC|STATIC;PLL=DYNAMIC|STATIC`

#### Virtual Synchronous Machine (`VSM`) & Interface Converter (`INTCONV_VSM`)
Named parameters (on `mBase` / `rateA` in `.dyr`, scaled to `Sb` in `prm_*`):
* `Ta` ($= 2H$), `D`, `Rp` (`k_droop_f = 1/Rp`), `Kq` (`k_droop_u`), `Tf` (`iTf = 1/Tf`), `Rv`, `Lv`, `Rf`, `Lf`, `Cf`, `kp_v`, `ki_v`, `kp_i`, `ki_i`, `k_ffi`, `k_ffv`
* On `INTCONV_VSM`: `K_p_dc`, `k_p_atm`, `k_i_atm`, `P_atm_lim`, `T_f_dc`
* `COMPLEXITY`: `VCC=DYNAMIC|STATIC|DYNAMIC_LIMITING` (`DYNAMIC_LIMITING_P|DYNAMIC_LIMITING_DC;REF=CONSTANT|DC_ATM` on `INTCONV_VSM`)

#### Droop Grid-Forming Generator (`GFM`) & Interface Converter (`INTCONV_GFM`)
Named parameters (on `mBase` / `rateA` in `.dyr`, scaled to `Sb` in `prm_*`):
* `Tf` (power measurement LPF time constant), `mp`, `nq` (active/reactive droop coefficients), `VI_R`, `VI_L` (virtual impedance), `WOF_Tf` (washout damping filter), `Rf`, `Lf`, `Cf`, `kp_v`, `ki_v`, `kp_i`, `ki_i`
* On `INTCONV_GFM`: DC-voltage droop/PI parameters (`K_p_dc`, `K_i_dc`)
* `COMPLEXITY`: `VCC=DYNAMIC|STATIC`

#### DC Primary Energy Source (`DC_PES`) & DC Network (`DCBUS`, `DCLINE`)
* `DC_PES`: `Tau` (s), `Kp` (pu), `Ki` (pu/s).
* `<dcbus> 'DCBUS' 1 C=<pu> /`: overrides `DC.CbusDC` for that DC bus (default `0.02` pu).
* `<fbus>-<tbus> 'DCLINE' 1 L=<pu> /`: overrides `DC.L_Lines` for that DC line (default `0.01` pu).

---

## 3. Case Settings File (`<case>.json`)

An optional JSON file beside `<case>.xlsx` configuring simulation, solver, power-flow defaults, and optional symbolic variables:

```json
{
  "variables": {
    "Pf_set": -40.0
  },
  "pf.run": {
    "solver": "newton",
    "tolerance": 1e-8,
    "max_iterations": 20,
    "enforce_q_limits": false
  },
  "dyn.simulate": {
    "t_end": 15.0,
    "method": "ode15s",
    "rtol": 1e-5,
    "atol": 1e-7,
    "max_step": 0.01
  }
}
```

Both nested (`VFlexPy`-style `"dyn.simulate"` / `"pf.run"`) and flat keys are accepted by `fun_read_case_json.m` and stored in `session.config`.

---

## 4. Time-Domain Disturbance Events (`session.eventList`)

In `SimulationTab` (or headlessly via `session.eventList`), temporal disturbances are represented as rows of a 7-column MATLAB `table`:

| column | type | meaning |
|---|---|---|
| `Enabled` | `logical` | whether the event is active in `vflexp_simulate` |
| `Time` | `double` (s) | onset instant $t_{\text{ev}} \ge 0$ |
| `TargetType` | `string` | `"AC Bus"`, `"DC Bus"`, `"AC Generator"`, `"DC Generator"`, `"IntConv"`, or `"AC Line"` |
| `TargetId` | `string` | bus number (e.g. `"102"`), generator bus/index, or `"fbus-tbus"` line key |
| `Channel` | `string` | `"dP_load"` (MW), `"dQ_load"` (MVAr), `"dP"` / `"dQ"` (pu), `"dI"` (pu DC), `"dPref"` / `"dVref"` / `"dQref"`, or `"close"` (for `"AC Line"`) |
| `Magnitude` | `double` | step/pulse amplitude (in MW/MVAr for `*_load`, or pu for `dP`/`dQ`/`dI`/`d*ref`) |
| `Duration` | `double` (s) | `0` (or `Inf`) for a permanent step; `> 0` for a rectangular pulse of width `Duration` |
