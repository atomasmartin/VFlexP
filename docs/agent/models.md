# Models and Equations in VFlexP

This document describes the mathematical conventions, reference frames, complexity levels, and modular Simulink block models implemented in `Simulink/lib/`, `functions/masks/`, and `core/`.

---

## 1. Per-Unit Conventions, Reference Frames, and Modular Parameterisation

### Per-Unit Bases
* Every synchronous AC area $a$ has a nominal frequency $F_{b,a}$ (Hz), angular base frequency $\omega_{b,a} = 2\pi F_{b,a}$ (rad/s), and apparent power base $S_{b,a}$ (MVA), defined in the `System parameters` sheet (`session.tables.SysParam`).
* Network elements (lines, shunt capacitors, loads) and converter blocks (`VSM`, `GFL`, `GFM`, `DC_PES`, `intconv_GFL`, `intconv_VSM`, `intconv_GFM`) operate in per unit on the **system base $S_b$**.
* Synchronous machine blocks (`MGES`, `SG`) operate internally in per unit on their **machine base `mBase`** (from `Gen data`), converting injected stator currents to the system base via $\mathbf{i}_{\text{sys}} = (\text{mBase}/S_b)\,\mathbf{i}_{\text{mach}}$ and terminal power to the machine base via $P_{\text{mach}} = (S_b/\text{mBase})\,P_{\text{sys}}$.
* When `.dyr` records are built by `vflexp_build_and_init.m`, converter parameters specified on `mBase` (or `rateA` for interface converters) are automatically scaled to $S_b$ using $\text{scale} = S_b / S_{\text{unit}}$ (impedances $\times \text{scale}$, admittances/capacitances/inertias $/ \text{scale}$, droop resistances $\times \text{scale}$).

### Park ($dq$) Transformation and Power-Invariant Convention
VFlexP uses a rotating $dq$ reference frame where active and reactive power are:

$$P = v_d i_d + v_q i_q, \qquad Q = v_q i_d - v_d i_q$$

At the steady-state power-flow solution $(V_m, \theta)$ of each bus:

$$v_{d0} = V_m \cos\theta, \qquad v_{q0} = V_m \sin\theta$$

Each synchronous area has its own rotating reference frame whose frequency $\omega_{\text{area}}$ is given by the first generator in that area (`Grid_par.FirstGen_each_area`). A unit with local rotor or virtual angle $\delta_i$ transforms terminal voltages from the area frame $(d, q)$ to its local frame $(d_i, q_i)$ by:

$$\begin{bmatrix} v_{di} \\ v_{qi} \end{bmatrix} = \begin{bmatrix} \cos\delta_i & \sin\delta_i \\ -\sin\delta_i & \cos\delta_i \end{bmatrix} \begin{bmatrix} v_d \\ v_q \end{bmatrix}, \qquad \frac{\mathrm{d}\delta_i}{\mathrm{d}t} = \omega_b (\omega_i - \omega_{\text{area}})$$

Because the first generator of each area defines the area reference angle ($\delta_1 \equiv 0$), its angle integrator state is eliminated in `MGES`/`SG` via the selector matrix `Eliminate_delta1` (`T_elim_delta1`).

### Immutable Library Links & Expression-Driven Variants (`Simulink/lib/`)
* **Primitive building blocks (`lib_primitives.slx`)**: `PI_ctrl`, `LPF_opt`, `LeadLag_opt`, `WOF_opt`, `RL_dq`, and `RC_dq` declare an explicit initial-condition parameter `x0` on their mask (so internal `Integrator` blocks permanently hold `InitialCondition = 'x0'`) and switch between `Active` (`any(T > 0)`) and `Bypass` (`all(T == 0)`) using Simulink expression-mode `Variant Subsystem` blocks.
* **Zero block mutation**: `vflexp_build_and_init.m` writes `prm_<Unit>` (physical arrays + numeric complexity variant codes + `n_units`) and `op_<Unit>` (steady-state `x0` vectors) into the `ModelWorkspace` of `Simulink/VFlexP.slx`. Unused device groups (`n_units == 0`) activate an internal `Disabled` variant outputting zero current, keeping all 30 library links resolved (`LinkStatus = 'resolved'`) and `Dirty = 'off'`.

---

## 2. Network Models (`lib_AC_Grid.slx` and `lib_DC_Grid.slx`)

Transmission lines use a $\pi$-model whose half-line charging susceptances $b/2$ are lumped into the parallel $RC$ shunt admittance of their terminal buses.

### Inductive Branch ($R$–$L$ Line or Series Load) in $dq$
* **Dynamic (`dynamic` / `EMT` / `series RL dynamic`)**:
  $$\frac{L}{\omega_b}\frac{\mathrm{d}i_d}{\mathrm{d}t} = v_d - R i_d + \omega L i_q, \qquad \frac{L}{\omega_b}\frac{\mathrm{d}i_q}{\mathrm{d}t} = v_q - R i_q - \omega L i_d$$
  *(Purely resistive AC loads with $Q_d \le 0$ on dynamic capacitor buses automatically use the exact algebraic branch to avoid $1/L$ division by zero.)*
* **Quasi-stationary with frequency coupling (`static omega` / `RMS` / `series RL static omega`)**:
  $$\begin{bmatrix} i_d \\ i_q \end{bmatrix} = \frac{1}{R^2 + (\omega L)^2} \begin{bmatrix} R & \omega L \\ -\omega L & R \end{bmatrix} \begin{bmatrix} v_d \\ v_q \end{bmatrix}$$
* **Static at nominal frequency (`static`)**: same algebraic matrix with $\omega = 1\text{ pu}$.

### Capacitive Bus Shunt ($G$–$C$ Bus) in $dq$
* **Dynamic (`dynamic` / `EMT`)**:
  $$\frac{C}{\omega_b}\frac{\mathrm{d}v_d}{\mathrm{d}t} = i_d - G v_d + \omega C v_q, \qquad \frac{C}{\omega_b}\frac{\mathrm{d}v_q}{\mathrm{d}t} = i_q - G v_q - \omega C v_d$$
* **Algebraic (`static omega` / `RMS`)**:
  $$\begin{bmatrix} v_d \\ v_q \end{bmatrix} = \frac{1}{G^2 + (\omega C)^2} \begin{bmatrix} G & \omega C \\ -\omega C & G \end{bmatrix} \begin{bmatrix} i_d \\ i_q \end{bmatrix}$$

### DC Grid (`lib_DC_Grid.slx`)
* **Dynamic DC bus (`DC dynamic` / `DC EMT`)**: $\dfrac{C_{\text{dc}}}{\omega_b}\dfrac{\mathrm{d}v_{\text{dc}}}{\mathrm{d}t} = i_{\text{dc}} - G_{\text{dc}} v_{\text{dc}}$, with per-bus capacitance `DC.CbusDC` (from `.dyr` `DCBUS`, positive `Bs`, or default $0.02\text{ pu}$).
* **Dynamic DC line (`DC dynamic` / `DC EMT`)**: $\dfrac{L_{\text{dc}}}{\omega_b}\dfrac{\mathrm{d}i_{\text{dc}}}{\mathrm{d}t} = \Delta v_{\text{dc}} - R_{\text{dc}} i_{\text{dc}}$, with per-line inductance `DC.L_Lines` (from `.dyr` `DCLINE`, positive `x`, or default $0.01\text{ pu}$).
* **Static DC network (`DC static` / `DC RMS` / `DC static omega`)**: algebraic conductance matrix $\mathbf{v}_{\text{dc}} = \mathbf{G}_{\text{dc}}^{-1}\mathbf{i}_{\text{dc}}$.
* **DC Load models**: `DC ConstR` (constant resistance $R_{\text{dc}} = V_{\text{dc},0}^2 / P_{\text{dc},0}$), `DC ConstI` (constant current), `DC ConstP` (algebraic constant power), and `DC DynamicP` (filtered constant-power load).

---

## 3. Synchronous Machine with Controls (`lib_AC_MGES.slx` and `lib_AC_SG.slx`)

### Parameter Conversion (`param_genrou.m`)
Standard generator reactances and time constants $(X_d, X_q, X'_d, X'_q, X''_d, X_l, T'_{d0}, T''_{d0}, T'_{q0}, T''_{q0})$ are converted to per-unit circuit parameters $(L_{md}, L_{mq}, L_{ffd}, R_{fd}, L_{11d}, R_{1d}, L_{11q}, R_{1q}, L_{22q}, R_{2q})$ by `param_genrou.m`:
* **`classical` (default)**: explicit series-approximation expressions used in classical PSS/E and textbook models.
* **`exact`**: exact two-circuit operational inductance root matching (IEEE Std 1110) solving the coupled quadratic relations so $L_d(s)$ and $L_q(s)$ reproduce $(X'_d, X''_d, T'_{d0}, T''_{d0})$ and $(X'_q, X''_q, T'_{q0}, T''_{q0})$ without the $T'' \ll T'$ approximation.
* **Explicit overrides**: if `.dyr` supplies `Lfd, Rfd, Lkd1, Rkd1, Lkq1, Rkq1, Lkq2, Rkq2`, they override the converted values directly.

### Flux-Linkage Inductance Matrix ($\mathbf{L}$)
The 6 electrical windings $(d, q, fd, kd_1, kq_1, kq_2)$ satisfy $\boldsymbol{\psi} = \mathbf{L}\mathbf{i}$ with:

$$\mathbf{L} = \begin{bmatrix}
-(L_{md}+L_l) & 0 & L_{md} & L_{md} & 0 & 0 \\
0 & -(L_{mq}+L_l) & 0 & 0 & L_{mq} & L_{mq} \\
-L_{md} & 0 & L_{ffd} & L_{md} & 0 & 0 \\
-L_{md} & 0 & L_{md} & L_{11d} & 0 & 0 \\
0 & -L_{mq} & 0 & 0 & L_{11q} & L_{mq} \\
0 & -L_{mq} & 0 & 0 & L_{mq} & L_{22q}
\end{bmatrix}$$

### Machine Complexity Variants (`MGES` and `SG`)
* **`order 8` (`EMT8`)**: All 6 flux linkages $(\psi_{sd}, \psi_{sq}, \psi_{fd}, \psi_{kd1}, \psi_{kq1}, \psi_{kq2})$ are dynamic states integrated by $1/\omega_b$, plus rotor speed $\omega_r$ and rotor angle $\delta_r$:
  $$\frac{1}{\omega_b}\frac{\mathrm{d}\psi_{sd}}{\mathrm{d}t} = v_d + R_s i_d + \omega_r \psi_{sq}, \qquad \frac{1}{\omega_b}\frac{\mathrm{d}\psi_{sq}}{\mathrm{d}t} = v_q + R_s i_q - \omega_r \psi_{sd}$$
  $$\frac{1}{\omega_b}\frac{\mathrm{d}\psi_{fd}}{\mathrm{d}t} = \frac{R_{fd}}{L_{md}} E_{fd} - R_{fd} i_{fd}, \qquad \frac{1}{\omega_b}\frac{\mathrm{d}\psi_{k}}{\mathrm{d}t} = -R_k i_k \quad (k \in \{kd_1, kq_1, kq_2\})$$
* **`order 4` (`RMS4`)**: Quasi-stationary stator ($\mathrm{d}\psi_{sd}/\mathrm{d}t = \mathrm{d}\psi_{sq}/\mathrm{d}t = 0$) with 4 rotor flux states $(\psi_{fd}, \psi_{kd1}, \psi_{kq1}, \psi_{kq2})$ plus $(\omega_r, \delta_r)$. Eliminating $(\psi_{sd}, \psi_{sq})$ from $\boldsymbol{\psi} = \mathbf{L}\mathbf{i}$ yields the Schur-complement rotor relation $\boldsymbol{\psi}_r = \mathbf{L}_{\text{rot}}\mathbf{i}_r + \mathbf{M}_{r,s}\mathbf{i}_s$ and internal flux $\boldsymbol{\psi}_{\text{int}} = \mathbf{M}_{s,r}\mathbf{i}_r$, solving the algebraic stator voltage equation with transient saliency at each step:
  $$\begin{bmatrix} R_s & -\omega_r L''_q \\ \omega_r L''_d & R_s \end{bmatrix} \begin{bmatrix} i_d \\ i_q \end{bmatrix} = -\begin{bmatrix} v_d + \omega_r \psi_{\text{int},q} \\ v_q - \omega_r \psi_{\text{int},d} \end{bmatrix}$$
* **`order 2` (`SG` only)**: Classical swing-equation model with constant internal voltage behind transient reactance.

### Controls in `MGES`
1. **Terminal Voltage Transducer & AVR (`ST4B`-compatible)**: measurement lag $1/(1 + s T_{\text{meas}})$ (`T_meas`), PI + lead-lag regulator (`K_pr`, `K_ir`, `K_d`, `T_d`, `T_a`, `T_b`, `T_e`, `T_f`) with anti-windup limits `[E_fdmin, E_fdmax]`.
2. **Turbine Governor (`IEEEG1`-compatible)**: droop $K = 1/R$, lead-lag $(1 + s T_2)/(1 + s T_1)$, rate- and position-limited servomotor, and four steam stages $(T_4..T_7)$ weighted by shaft fractions $(K_1+K_2, K_3+K_4, K_5+K_6, K_7+K_8)$.
3. **Power System Stabiliser (`PSS2B`)**: dual-input speed/power washouts ($T_{w1}..T_{w4}$), transducers ($T_6, T_7$), ramp-tracking filter $\left(\frac{1 + s T_8}{1 + s T_9}\right)^M$, and three lead-lag stages $(T_1/T_2, T_3/T_4, T_{10}/T_{11})$ with gain $K_{s1}$ and limits `[V_stmin, V_stmax]`.

---

## 4. Converter-Interfaced Generators and AC/DC Interface Converters

### Virtual Synchronous Machine (`lib_AC_VSM.slx` and `lib_IC_VSM.slx`)
* **Power Controller (`PC`)**:
  $$T_a \frac{\mathrm{d}\omega}{\mathrm{d}t} = P_{\text{ref}} - P_{\text{meas}} - (D + k_{\text{droop},f})(\omega - \omega_{\text{ref}}), \qquad v_{\text{ref}} = V_0 + k_{\text{droop},u}(Q_{\text{ref}} - Q_{\text{meas}})$$
  with power measurement filter `iTf = 1/Tf` and virtual impedance (`Rv`, `Lv`).
* **Voltage and Current Controller (`VCC`) Variants**: `dynamic` (cascaded PI voltage + PI current + dynamic $LC$ filter), `static` (algebraic inner loop + static filter), and `dynamic_limiting` / `dynamic_limiting_P` / `dynamic_limiting_DC`.
* **DC-Side Coupling & ATM Control (`intconv_VSM`)**: includes FUBM converter loss polynomial (`alpha1`, `alpha2`, `alpha3`) and optional `DC_ATM` outer DC-voltage / active-power fast-tracking control (`K_p_dc`, `k_p_atm`, `k_i_atm`, `P_atm_lim`, `T_f_dc`).

### Droop Grid-Forming Converter (`lib_AC_GFM.slx` and `lib_IC_GFM.slx`)
* **Active/Reactive Droop & Virtual Impedance**: low-pass filtered active and reactive power measurements (`Tf`), frequency droop `mp`, voltage droop `nq`, virtual impedance (`VI_R`, `VI_L`), washout damping filter (`WOF_Tf`), and `dynamic` vs `static` `VCC` inner loops.

### Grid-Following Converter (`lib_AC_GFL.slx` and `lib_IC_GFL.slx`)
* **Phase-Locked Loop (`PLL`)**:
  * `dynamic`: PI controller (`Kp_pll`, `Ki_pll`) driving $v_{oq}^{\text{local}} \to 0$ and integrating $\theta_{\text{pll}}$.
  * `static`: exact algebraic angle tracking $\theta_{\text{pll}} = \text{atan2}(v_q, v_d)$ (initialised cleanly at $t = 0$ via `InitialCondition` block `op_*.PLL_dynamic_frame` with zero transport delay).
* **Power & DC-Voltage Controller (`PC`)**: filtered power references (`iTf = 1/Ta`), frequency/voltage droops (`1/Rp`, `Kq`), plus DC-voltage PI/droop (`K_p`, `K_i`) on `intconv_GFL`.
* **Current Controller (`CC`)**: `dynamic` (PI current loop + dynamic filter) or `static` (algebraic current injection).

### DC Primary Energy Source (`lib_DC_DC_PES.slx`)
* Dispatchable DC source/storage behind first-order lag $\tau$ (`Tau`) and DC voltage droop + integral controller (`Kp`, `Ki`):
  $$\tau \frac{\mathrm{d}P_{\text{dc}}}{\mathrm{d}t} = P_{\text{dc,0}} + \left(K_p + \frac{K_i}{s}\right)(v_{\text{dc,0}} - v_{\text{dc}}) - P_{\text{dc}}, \qquad i_{\text{dc}} = \frac{P_{\text{dc}}}{v_{\text{dc}}}$$

---

## 5. Exact Linearisation, Modal Analysis, and State Relevance

* **Linearisation (`vflexp_linearise.m`)**: Uses Simulink Control Design's `linearize('VFlexP', io)` (plus custom analytical trigonometric Jacobians `linearizedsin/cos/atan.m`) at the initialised equilibrium point to obtain $(\mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D})$ and physical state metadata (`vflexp_map_state_names.m`).
* **Eigenvalues & Participation Factors (`vflexp_modal.m`)**: Computes right eigenvectors $\mathbf{V}$ ($\mathbf{A}\mathbf{V} = \mathbf{V}\boldsymbol{\Lambda}$), left eigenvectors $\mathbf{W} = \mathbf{V}^{-1}$, natural frequencies $f_i = |\text{Im}(\lambda_i)|/(2\pi)$, damping ratios $\zeta_i = -\text{Re}(\lambda_i)/|\lambda_i|$, and normalised participation factors $P_{ki} = |V_{ki} W_{ik}| / \sum_j |V_{ji} W_{ij}|$.
* **State Relevance & Automatic Complexity Reduction (`vflexp_relevance.m`, `vflexp_apply_relevance.m`)**: Computes the input-output relevance of every physical state by combining the Hankel singular values of the balanced realisation (`balreal`) with the participation factors of the original and balanced state-space models (Tomás-Martín et al., *IJEPES*, 2023), normalised to a maximum of 1, and automatically downgrades sub-blocks whose states all fall below a user-selected threshold.
