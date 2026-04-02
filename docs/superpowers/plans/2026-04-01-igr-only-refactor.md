# IGR-Only Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove all non-IGR code paths from MFC, making IGR the only solver. Remove incompatible features (bubbles, IBM, hypoelasticity, hyperelasticity, surface tension, acoustic sources, phase change, MHD, cylindrical coords, alt_soundspeed, characteristic BCs). Remove the `igr` toggle itself. Clean up examples, tests, and Python toolchain.

**Architecture:** Layered removal in 5 stages, each producing a compilable codebase. Stage 1 removes non-IGR solver infrastructure. Stage 2 removes incompatible feature modules. Stage 3 inlines IGR guards. Stage 4 cleans the Python toolchain. Stage 5 cleans examples and tests.

**Tech Stack:** Fortran 2008+ with Fypp preprocessing, Python toolchain, CMake build (glob-based, so file deletions auto-propagate).

**Build system note:** CMake uses `file(GLOB ...)` to discover `.fpp` files. Deleting a `.fpp` file automatically removes it from the build. The only CMakeLists.txt edits needed are for special-case exclusions and compiler flag references to specific files.

---

## Stage 1: Remove Non-IGR Solver Infrastructure

### Task 1: Delete non-IGR solver module files

These modules are only initialized/called when `igr=false`. Deleting them is safe because all call sites are inside `.not. igr` guards that will be removed in subsequent tasks.

**Files to delete:**
- `src/simulation/m_weno.fpp` (1,573 lines — WENO/WENO-Z/TENO reconstruction)
- `src/simulation/m_muscl.fpp` (365 lines — MUSCL reconstruction)
- `src/simulation/m_riemann_solvers.fpp` (4,657 lines — HLL/HLLC/HLLD/LF Riemann solvers)
- `src/simulation/m_cbc.fpp` (1,456 lines — characteristic boundary conditions)
- `src/simulation/m_compute_cbc.fpp` (368 lines — CBC computation kernels)
- `src/simulation/m_viscous.fpp` (1,363 lines — non-IGR viscous stress module)
- `src/simulation/include/inline_riemann.fpp` (116 lines — inline Riemann macros)

**~9,900 lines removed.**

- [ ] **Step 1:** Delete the 7 files listed above using `git rm`.

```bash
cd /storage/scratch1/6/sbryngelson3/mfc-igr
git rm src/simulation/m_weno.fpp \
       src/simulation/m_muscl.fpp \
       src/simulation/m_riemann_solvers.fpp \
       src/simulation/m_cbc.fpp \
       src/simulation/m_compute_cbc.fpp \
       src/simulation/m_viscous.fpp \
       src/simulation/include/inline_riemann.fpp
```

- [ ] **Step 2:** Commit.

```bash
git commit -m "remove non-IGR solver modules: WENO, MUSCL, Riemann, CBC, viscous"
```

---

### Task 2: Strip non-IGR code from m_rhs.fpp

This is the main branching point. Remove all `if (.not. igr)` blocks (allocation, computation, deallocation) and keep only the IGR path. Also remove `use` statements for deleted modules.

**Files:**
- Modify: `src/simulation/m_rhs.fpp`

Key regions to modify (line numbers from current file):

**Initialization (`s_initialize_rhs_module`):**
- Remove `use m_weno`, `use m_muscl`, `use m_riemann_solvers`, `use m_cbc`, `use m_viscous` imports
- Lines 141-150: Remove `if (.not. igr)` block allocating `q_cons_qp%vf` and `q_prim_qp%vf`
- Lines 166-186: Remove `if (.not. igr)` block for `q_prim_qp` pointer setup and `ACC_SETUP_VFs`
- Lines 200-271: Remove `if (.not. igr)` block allocating `flux_n`, `flux_src_n`, `flux_gsrc_n`
- Lines 273-491: Remove `if ((.not. igr) .or. dummy)` block allocating all WENO reconstruction data (`qL_rsx_vf`, `qR_rsx_vf`, `dqL_prim_*`, `tau_Re_vf`, etc.)
- Lines 516-519: Remove `if (.not. igr)` block allocating `gm_alphaL_n`, `gm_alphaR_n`

**Computation (`s_compute_rhs`):**
- Lines 564-600: Remove `if (.not. igr .or. dummy)` block copying `q_cons_vf` to `q_cons_qp`
- Lines 602-606: Keep IGR buffer population but remove `if (igr .or. dummy)` guard (make unconditional)
- Lines 607-615: Remove `if (.not. igr .or. dummy)` block (conservative-to-primitive conversion + buffer population)
- Lines 630-636: Remove `if ((viscous .and. .not. igr) .or. dummy)` block (non-IGR viscous call)
- Lines 646-674: Keep IGR computation block, remove `if (igr .or. dummy)` guard (make unconditional)
- Lines 675-862: Remove entire `if ((.not. igr) .or. dummy)` block (~190 lines: WENO reconstruction, Riemann solver, advection source, hypoelasticity, chemistry diffusion, additional physics, bubble dynamics, all non-IGR source terms)
- Lines 874-887: Remove `if (.not. igr .or. dummy)` block copying `q_prim_qp` back to `q_prim_vf`

**Finalization (`s_finalize_rhs_module`):**
- Lines 1750-1769: Remove `if (.not. igr)` block deallocating `q_cons_qp`/`q_prim_qp`
- Lines 1773-1855: Remove `if (.not. igr)` block deallocating WENO reconstruction data
- Lines 1862-1903: Remove `if (.not. igr)` block deallocating flux arrays

Also remove declarations of variables that are now unused: `flux_n`, `flux_src_n`, `flux_gsrc_n`, `q_cons_qp`, `q_prim_qp`, `qL_rsx_vf`, `qR_rsx_vf`, `qL_rsy_vf`, `qR_rsy_vf`, `qL_rsz_vf`, `qR_rsz_vf`, `dqL_prim_*`, `dqR_prim_*`, `tau_Re_vf`, `gm_alphaL_n`, `gm_alphaR_n`, and any procedure pointers to deleted modules (`s_riemann_solver`, `s_reconstruct_cell_boundary_values`, etc.).

Also remove subroutines that are only called from the non-IGR path:
- `s_compute_advection_source_term` (calls CBC, only from non-IGR block)
- `s_compute_additional_physics_rhs` (viscous/surface_tension/chem_diffusion, only from non-IGR block)
- `s_get_viscous` (WENO viscous reconstruction, only from non-IGR block)
- Any other subroutines in m_rhs.fpp that are only reachable from the deleted non-IGR block

- [ ] **Step 1:** Read m_rhs.fpp fully to understand current structure.
- [ ] **Step 2:** Remove `use` statements for deleted modules.
- [ ] **Step 3:** Remove non-IGR initialization blocks and unused variable declarations.
- [ ] **Step 4:** Simplify computation section: remove non-IGR blocks, make IGR blocks unconditional.
- [ ] **Step 5:** Remove non-IGR finalization blocks.
- [ ] **Step 6:** Remove now-dead subroutines (`s_compute_advection_source_term`, `s_get_viscous`, `s_compute_additional_physics_rhs`, etc.).

---

### Task 3: Strip non-IGR code from m_start_up.fpp (simulation)

Remove initialization and finalization of deleted modules.

**Files:**
- Modify: `src/simulation/m_start_up.fpp`

Key changes:
- Remove `use m_weno`, `use m_muscl`, `use m_riemann_solvers`, `use m_cbc`, `use m_viscous` imports
- Lines 920-921: Remove `if (igr .or. dummy)` guard around `s_initialize_igr_module()` — make unconditional
- Lines 923-928: Remove `if (.not. igr .or. dummy)` block that initializes WENO or MUSCL based on `recon_type`
- Line 929: Remove `if (.not. igr .or. dummy)` block initializing CBC module
- Line 930: Remove `if (.not. igr .or. dummy)` block initializing Riemann solvers module
- Lines 864-866: Remove `if (viscous .and. (.not. igr))` block initializing viscous module
- Lines 1090-1091: Remove `if (igr)` guard around `s_finalize_igr_module()` — make unconditional
- Lines 1093-1098: Remove `else` block finalizing CBC, Riemann solvers, WENO/MUSCL
- Lines 1108-1109: Remove `if (viscous .and. (.not. igr))` block finalizing viscous module

- [ ] **Step 1:** Read m_start_up.fpp and apply changes.
- [ ] **Step 2:** Remove dead `use` statements.
- [ ] **Step 3:** Make IGR init/finalize unconditional, remove non-IGR init/finalize blocks.

---

### Task 4: Strip non-IGR checks from m_checker.fpp (simulation)

Remove WENO/MUSCL/Riemann solver validation that no longer applies.

**Files:**
- Modify: `src/simulation/m_checker.fpp`

- Remove all `@:PROHIBIT` checks related to `weno_order`, `riemann_solver`, `recon_type`, `muscl_order`, `muscl_lim`, `wave_speeds`, `avg_state`.
- Keep IGR-specific checks (`nv_uvm_igr_temps_on_gpu`, `igr_iter_solver`, etc.).

- [ ] **Step 1:** Read m_checker.fpp and remove non-IGR validation checks.

---

### Task 5: Update CMakeLists.txt special cases

Remove references to deleted files in compiler-specific workarounds.

**Files:**
- Modify: `CMakeLists.txt`

- Lines 502-506: Remove reference to `m_cbc` in cross-file inlining exclusion
- Lines 505: Remove mention of `m_riemann_solvers`, `m_viscous`, `m_weno` from comment
- Lines 721-729: Remove Cray compiler workarounds for `m_bubbles_EL` and `m_phase_change` (these files will be deleted in Stage 2, but the CMake references should be cleaned up)

- [ ] **Step 1:** Read relevant CMakeLists.txt sections and clean up references.

---

### Task 6: Verify Stage 1 builds

- [ ] **Step 1:** Format and precheck.

```bash
./figr.sh format -j 8
./figr.sh precheck -j 8
```

- [ ] **Step 2:** Build all three targets.

```bash
./figr.sh build -j 8
```

- [ ] **Step 3:** Fix any compilation errors (missing symbols, dangling references).

- [ ] **Step 4:** Commit.

```bash
git add -A
git commit -m "strip non-IGR solver code from m_rhs, m_start_up, m_checker, CMakeLists"
```

---

## Stage 2: Remove Incompatible Feature Modules

### Task 7: Delete incompatible feature module files

**Files to delete from `src/simulation/`:**
- `m_bubbles.fpp` (596 lines — shared bubble dynamics)
- `m_bubbles_EE.fpp` (326 lines — Euler-Euler bubble source terms)
- `m_bubbles_EL.fpp` (1,643 lines — Lagrangian bubble tracking)
- `m_bubbles_EL_kernels.fpp` (386 lines — Lagrangian bubble kernels)
- `m_ibm.fpp` (1,233 lines — immersed boundary method)
- `m_ib_patches.fpp` (1,148 lines — IBM patch geometries) — also in `src/common/`
- `m_compute_levelset.fpp` (659 lines — level-set for IBM) — also in `src/common/`
- `m_hypoelastic.fpp` (423 lines — hypoelastic stress)
- `m_hyperelastic.fpp` (268 lines — hyperelastic deformation)
- `m_surface_tension.fpp` (390 lines — capillary source fluxes)
- `m_acoustic_src.fpp` (687 lines — acoustic source injection)
- `m_qbmm.fpp` (1,080 lines — QBMM for polydisperse bubbles)
- `include/inline_capillary.fpp` (21 lines — capillary force macros)

**Files to delete from `src/common/`:**
- `m_phase_change.fpp` (~1,100 lines — phase transition relaxation)

**Files to delete from `src/pre_process/`:**
- `m_check_ib_patches.fpp` (227 lines — IBM patch validation)

**~10,200 lines removed.**

- [ ] **Step 1:** Delete all listed files with `git rm`.

```bash
git rm src/simulation/m_bubbles.fpp \
       src/simulation/m_bubbles_EE.fpp \
       src/simulation/m_bubbles_EL.fpp \
       src/simulation/m_bubbles_EL_kernels.fpp \
       src/simulation/m_ibm.fpp \
       src/simulation/m_hypoelastic.fpp \
       src/simulation/m_hyperelastic.fpp \
       src/simulation/m_surface_tension.fpp \
       src/simulation/m_acoustic_src.fpp \
       src/simulation/m_qbmm.fpp \
       src/simulation/include/inline_capillary.fpp \
       src/common/m_phase_change.fpp \
       src/pre_process/m_check_ib_patches.fpp
```

Note: `m_ib_patches.fpp` and `m_compute_levelset.fpp` exist in `src/common/`. Check if they are used by any remaining code before deleting. They are currently excluded from post_process builds (CMakeLists.txt lines 410-413), suggesting they are only used for IBM. Delete them:

```bash
git rm src/common/m_ib_patches.fpp \
       src/common/m_compute_levelset.fpp
```

- [ ] **Step 2:** Commit.

```bash
git commit -m "remove incompatible feature modules: bubbles, IBM, elasticity, surface tension, acoustic, phase change, QBMM"
```

---

### Task 8: Strip feature references from simulation modules

Remove `use` statements, initialization/finalization calls, and feature-gated code blocks for deleted features.

**Files to modify:**

**`src/simulation/m_start_up.fpp`:**
- Remove `use` statements for all deleted modules
- Remove initialization calls: `s_initialize_bubbles_EE_module`, `s_initialize_bubbles_EL_module`, `s_initialize_ibm_module`, `s_initialize_hypoelastic_module`, `s_initialize_hyperelastic_module`, `s_initialize_surface_tension_module`, `s_initialize_acoustic_src_module`, `s_initialize_qbmm_module`, `s_initialize_phase_change_module`
- Remove finalization calls for all the above
- Remove all conditional blocks gated on: `bubbles_euler`, `bubbles_lagrange`, `ib`, `hypoelasticity`, `hyperelasticity`, `surface_tension`, `acoustic_source`, `qbmm`, `relax`, `mhd`

**`src/simulation/m_rhs.fpp` (remaining references after Task 2):**
- Remove any remaining `use` statements for deleted modules
- Remove any remaining feature-gated blocks for deleted features (bubble RHS, IBM forcing, hypoelastic RHS, surface tension flux, acoustic source injection, phase change relaxation)
- Remove bubble-related variables and allocations
- Remove IBM-related code paths

**`src/simulation/m_global_parameters.fpp`:**
- Remove variable declarations for: `bubbles_euler`, `bubbles_lagrange`, `bubble_model`, `polytropic`, `polydisperse`, `qbmm`, `nb`, `R0ref`, `ib`, `num_ibs`, `hypoelasticity`, `hyperelasticity`, `surface_tension`, `sigma`, `acoustic_source`, `num_source`, `relax`, `relax_model`, `mhd`, `hyper_cleaning`, `cyl_coord`, `alt_soundspeed`, `probe_wrt`
- Remove all bubble index variables (`bub_idx`, `bub_fld`, etc.)
- Remove IBM-related variables
- Remove stress index variables for hypoelasticity
- Remove MHD index variables (`B_idx`, etc.)
- Remove subroutines that compute index ranges for deleted features
- Keep: `igr_order`, `igr_iter_solver`, `num_igr_iters`, `num_igr_warm_start_iters`, `alf_factor`, `igr_pres_lim`, `nv_uvm_igr_temps_on_gpu`

**`src/simulation/m_checker.fpp`:**
- Remove all `@:PROHIBIT` checks for deleted features

**`src/simulation/m_time_steppers.fpp`:**
- Remove any bubble/IBM/feature-specific time stepping logic

**`src/simulation/m_data_output.fpp`:**
- Remove probe write code (probe_wrt incompatible with IGR)
- Remove bubble/IBM/MHD-specific output logic

**`src/simulation/m_derived_variables.fpp`:**
- Remove any feature-specific derived variable computation for deleted features

**`src/simulation/m_sim_helpers.fpp`:**
- Remove bubble/IBM/feature-specific helper code

- [ ] **Step 1:** Read each file, identify all references to deleted features.
- [ ] **Step 2:** Remove `use` statements, variable declarations, initialization/finalization calls.
- [ ] **Step 3:** Remove feature-gated code blocks.
- [ ] **Step 4:** Remove now-unused subroutines.

---

### Task 9: Strip feature references from common modules

**Files to modify:**

**`src/common/m_global_parameters.fpp`:**
- Remove variable declarations for all deleted features (same list as simulation target)
- Remove feature-specific index computation subroutines
- Remove bubble/IBM/MHD/elasticity parameter groups

**`src/common/m_variables_conversion.fpp`:**
- Remove bubble-specific conversion code
- Remove MHD-specific conversion code
- Remove hypoelastic/hyperelastic stress variable conversion
- Remove surface tension variable handling

**`src/common/m_boundary_common.fpp`:**
- Keep `s_populate_F_igr_buffers` (IGR-specific, ~lines 1368-1532)
- Remove IBM boundary handling
- Remove bubble boundary handling
- Remove characteristic BC handling (those BC types are incompatible with IGR)

**`src/common/m_mpi_common.fpp`:**
- Remove bubble-specific MPI packing/unpacking
- Remove IBM-specific MPI communication
- Remove QBMM-specific MPI handling

**`src/common/m_derived_types.fpp`:**
- Remove types: `bub_bounds_info`, `bubble_mv`, `subgrid_bubble_physical_parameters`, `bubbles_lagrange_parameters`
- Remove IBM-related types
- Remove any types only used by deleted features

**`src/common/m_helper_basic.fpp`:**
- Remove bubble/IBM helper functions

**`src/common/m_checker_common.fpp`:**
- Remove validation for deleted features

**`src/common/m_model.fpp`:**
- Check for IBM-specific model handling and remove if present

- [ ] **Step 1:** Read each common file, identify references to deleted features.
- [ ] **Step 2:** Remove types, variables, subroutines for deleted features.
- [ ] **Step 3:** Keep IGR-specific code intact (especially `s_populate_F_igr_buffers`).

---

### Task 10: Strip feature references from pre_process

**Files to modify:**

**`src/pre_process/m_global_parameters.fpp`:**
- Remove declarations for all deleted features

**`src/pre_process/m_start_up.fpp`:**
- Remove initialization of deleted feature modules
- Remove `use` statements for deleted modules

**`src/pre_process/m_icpp_patches.fpp`:**
- Remove bubble-specific initial condition patches
- Remove IBM-related patch setup

**`src/pre_process/m_assign_variables.fpp`:**
- Remove bubble variable assignment
- Remove stress/elastic variable assignment
- Remove MHD variable assignment

**`src/pre_process/m_check_patches.fpp`:**
- Remove validation for bubble/IBM/elastic/MHD patches

**`src/pre_process/m_data_output.fpp`:**
- Remove bubble/IBM/elastic/MHD output

- [ ] **Step 1:** Read each pre_process file, strip deleted feature references.

---

### Task 11: Strip feature references from post_process

**Files to modify:**

**`src/post_process/m_global_parameters.fpp`:**
- Remove declarations for all deleted features

**`src/post_process/m_start_up.fpp`:**
- Remove initialization of deleted feature modules
- Remove `use` statements for deleted modules

**`src/post_process/m_data_output.fpp`:**
- Remove bubble/IBM/elastic/MHD/surface_tension output logic

**`src/post_process/m_derived_variables.fpp` (if exists):**
- Remove feature-specific derived variable computation

- [ ] **Step 2:** Remove CMakeLists.txt exclusions for `m_compute_levelset.fpp` and `m_ib_patches.fpp` (lines 410-413) since those files are now deleted.

- [ ] **Step 1:** Read each post_process file, strip deleted feature references.

---

### Task 12: Verify Stage 2 builds

- [ ] **Step 1:** Format and precheck.

```bash
./figr.sh format -j 8
./figr.sh precheck -j 8
```

- [ ] **Step 2:** Build all three targets.

```bash
./figr.sh build -j 8
```

- [ ] **Step 3:** Fix any compilation errors.

- [ ] **Step 4:** Commit.

```bash
git add -A
git commit -m "strip all references to deleted feature modules from simulation, common, pre_process, post_process"
```

---

## Stage 3: Inline IGR Guards and Remove `igr` Parameter

### Task 13: Remove `igr` parameter and inline all guards

Now that non-IGR code is gone, the `igr` boolean is always effectively true. Remove it and clean up all guards.

**Files to modify (all targets):**

**`src/simulation/m_global_parameters.fpp`:**
- Remove `igr` variable declaration (keep `igr_order`, `igr_iter_solver`, etc.)
- Remove `igr` from namelist in `m_start_up.fpp`

**`src/simulation/m_start_up.fpp`:**
- Remove `igr` from namelist
- Remove any remaining `if (igr)` guards — make the guarded code unconditional
- Remove any remaining `if (.not. igr)` guards — delete the guarded code entirely

**`src/simulation/m_rhs.fpp`:**
- Remove any remaining `if (igr .or. dummy)` guards — make unconditional
- Remove any remaining `if (.not. igr .or. dummy)` guards — delete

**`src/simulation/m_igr.fpp`:**
- Check for any self-references to `igr` parameter and remove

**`src/simulation/m_data_output.fpp`:**
- Remove any `igr` guards

**`src/common/m_boundary_common.fpp`:**
- Remove `igr` guards around `s_populate_F_igr_buffers` calls

**`src/common/m_helper_basic.fpp`:**
- Remove `igr` references

**`src/common/m_global_parameters.fpp`:**
- Remove `igr` declaration

**`src/pre_process/m_global_parameters.fpp`:**
- Remove `igr` declaration

**`src/post_process/m_global_parameters.fpp`:**
- Remove `igr` declaration

**All targets `m_start_up.fpp`:**
- Remove `igr` from all namelists

**Search strategy:** `grep -rn 'igr' src/` to find ALL remaining references. Each one must be either:
- Removed (if it's the `igr` boolean or a guard on it)
- Kept (if it's `igr_order`, `igr_iter_solver`, `igr_pres_lim`, etc.)

- [ ] **Step 1:** Grep for all `igr` references across `src/`.
- [ ] **Step 2:** Remove `igr` boolean declarations from all `m_global_parameters.fpp`.
- [ ] **Step 3:** Remove `igr` from all namelists.
- [ ] **Step 4:** Inline or remove all `if (igr)` / `if (.not. igr)` guards.
- [ ] **Step 5:** Rename `s_populate_F_igr_buffers` and similar to just reflect they're the standard path (optional, cosmetic).

---

### Task 14: Verify Stage 3 builds

- [ ] **Step 1:** Format and precheck.

```bash
./figr.sh format -j 8
./figr.sh precheck -j 8
```

- [ ] **Step 2:** Build.

```bash
./figr.sh build -j 8
```

- [ ] **Step 3:** Commit.

```bash
git add -A
git commit -m "remove igr parameter toggle, inline all IGR guards as unconditional"
```

---

## Stage 4: Clean Python Toolchain

### Task 15: Strip dead parameters from definitions.py

**File:** `toolchain/mfc/params/definitions.py`

Remove parameter registrations, constraints, and dependencies for:

**WENO/MUSCL/Riemann (solver infrastructure):**
- `weno_order`, `weno_eps`, `mapped_weno`, `wenoz`, `teno`, `mp_weno`, `teno_CT`, `wenoz_q`, `weno_avg`, `null_weights`, `weno_Re_flux`
- `recon_type`, `muscl_order`, `muscl_lim`
- `riemann_solver`, `wave_speeds`, `avg_state`, `low_Mach`

**Incompatible features:**
- `bubbles_euler`, `bubbles_lagrange`, `polytropic`, `polydisperse`, `qbmm`, `nb`, `R0ref`, `bubble_model`
- `surface_tension`, `sigma`
- `hypoelasticity`, `hyperelasticity`
- `ib`, `num_ibs`, `ib_state_wrt`
- `acoustic_source`, `num_source`
- `mhd`, `hyper_cleaning`, `hyper_cleaning_speed`, `hyper_cleaning_tau`, `Bx0`
- `relax`, `relax_model`, `palpha_eps`, `ptgalpha_eps`
- `alt_soundspeed`
- `cyl_coord`
- `probe_wrt`

**The `igr` toggle itself:**
- `igr` (remove — always on now)

**Keep:**
- `igr_order`, `igr_iter_solver`, `num_igr_iters`, `num_igr_warm_start_iters`, `alf_factor`, `igr_pres_lim`, `nv_uvm_igr_temps_on_gpu`

Also update `CASE_OPT_PARAMS` set to remove deleted parameters and `igr`.

Also remove entries from `CONSTRAINTS` and `DEPENDENCIES` dicts for deleted parameters.

Also clean up any indexed parameter families:
- `patch_ib(i)%*` — all IBM patch parameters
- `acoustic(i)%*` — all acoustic source parameters
- `bub_pp%*` — bubble physics parameters
- `fluid_pp(i)%Re(j)` — keep (viscous is still supported)
- Remove bubble-related fluid_pp attributes

- [ ] **Step 1:** Read definitions.py fully.
- [ ] **Step 2:** Remove parameter registrations for dead parameters.
- [ ] **Step 3:** Remove constraints for dead parameters.
- [ ] **Step 4:** Remove dependencies for dead parameters.
- [ ] **Step 5:** Update `CASE_OPT_PARAMS`.

---

### Task 16: Strip dead validation from case_validator.py

**File:** `toolchain/mfc/case_validator.py`

**Functions to remove entirely:**
- `check_weno()` — WENO validation
- `check_muscl()` — MUSCL validation
- `check_hypoelasticity()` — hypoelastic validation
- `check_phase_change()` — phase change validation
- `check_ibm()` — IBM validation
- `check_surface_tension()` — surface tension validation
- `check_mhd()` — MHD validation
- `check_weno_simulation()` — WENO simulation validation
- `check_muscl_simulation()` — MUSCL simulation validation
- `check_bubbles_euler_simulation()` — bubble simulation validation
- `check_mhd_simulation()` — MHD simulation validation
- `check_acoustic_source()` — acoustic source validation
- `check_alt_soundspeed()` — alt soundspeed validation
- `check_bubbles_lagrange()` — Lagrangian bubble validation
- `check_hyperelasticity()` — hyperelastic validation
- `check_surface_tension_post()` — post-process surface tension validation
- `check_bubbles_euler()` — bubble Euler validation
- `check_riemann_solver()` — Riemann solver validation
- `check_low_Mach()` — low Mach validation

**Functions to modify:**
- `check_igr()` — remove the `igr` parameter check (it's always on), keep grid resolution checks for `igr_order`
- `check_igr_simulation()` — remove checks for features that no longer exist. Keep IGR parameter validation (`num_igr_iters`, `igr_iter_solver`, `alf_factor` bounds).
- `check_boundary_conditions()` — remove cylindrical coordinate checks, remove characteristic BC checks (those BCs no longer exist)
- `check_parallel_io()` — remove `down_sample` requiring `igr` (igr is always on)
- Main validation function — remove calls to deleted check functions

- [ ] **Step 1:** Read case_validator.py fully.
- [ ] **Step 2:** Delete dead validation functions.
- [ ] **Step 3:** Simplify remaining validation functions.
- [ ] **Step 4:** Remove calls to deleted functions from the main validation entry point.

---

### Task 17: Strip dead descriptions from descriptions.py

**File:** `toolchain/mfc/params/descriptions.py`

Remove description entries for all parameters deleted in Task 15. Remove pattern-based descriptions for:
- `patch_ib()` patterns
- `acoustic()` patterns
- `bub_pp%` patterns
- `lag_*_wrt` patterns
- All WENO/MUSCL/Riemann descriptions
- All bubble/IBM/elastic/MHD descriptions

- [ ] **Step 1:** Read descriptions.py and remove dead entries.

---

### Task 18: Clean up other toolchain files

Check and clean:

**`toolchain/mfc/run/input.py`:**
- Remove any special handling for deleted parameters in namelist generation

**`toolchain/mfc/case.py`:**
- Remove references to deleted parameters

**`toolchain/mfc/state.py` or similar:**
- Remove feature-specific state handling

**`toolchain/mfc/packer/packer.py` (if exists):**
- Remove feature-specific packing

Search with: `grep -rn 'bubbles_euler\|bubbles_lagrange\|hypoelasticity\|hyperelasticity\|surface_tension\|acoustic_source\|mhd\|weno_order\|riemann_solver\|recon_type\|cyl_coord\|alt_soundspeed\|probe_wrt' toolchain/`

- [ ] **Step 1:** Grep toolchain/ for all dead parameter references and clean up.

---

### Task 19: Verify Stage 4

- [ ] **Step 1:** Run Python linting.

```bash
./figr.sh lint
```

- [ ] **Step 2:** Validate an IGR example case.

```bash
./figr.sh validate examples/3D_IGR_TaylorGreenVortex/case.py
```

- [ ] **Step 3:** Commit.

```bash
git add -A
git commit -m "strip dead parameters, validation, and descriptions from Python toolchain"
```

---

## Stage 5: Clean Examples, Tests, and Golden Files

### Task 20: Delete non-IGR examples and benchmarks

Delete all 131 non-IGR example directories and 4 non-IGR benchmark directories.

Keep 7 IGR examples:
- `examples/2D_IGR_2fluid/`
- `examples/2D_IGR_triple_point/`
- `examples/3D_IGR_33jet/`
- `examples/3D_IGR_jet/`
- `examples/3D_IGR_jet_1fluid/`
- `examples/3D_IGR_TaylorGreenVortex/`
- `examples/3D_IGR_TaylorGreenVortex_nvidia/`

Keep 1 IGR benchmark:
- `benchmarks/igr/`

- [ ] **Step 1:** Generate deletion list by finding all non-IGR example dirs.

```bash
# List all example dirs, exclude IGR ones
ls -d examples/*/ | grep -v IGR | xargs git rm -r
```

- [ ] **Step 2:** Delete non-IGR benchmarks.

```bash
git rm -r benchmarks/5eq_rk3_weno3_hllc benchmarks/hypo_hll benchmarks/ibm benchmarks/viscous_weno5_sgb_acoustic
```

- [ ] **Step 3:** Update IGR example case.py files to remove the now-unnecessary `"igr": "T"` parameter (since it no longer exists). Keep all other IGR-specific parameters (`igr_order`, etc.).

- [ ] **Step 4:** Remove any `"alt_soundspeed": "F"` or other removed parameters from IGR case files.

- [ ] **Step 5:** Commit.

```bash
git add -A
git commit -m "remove non-IGR examples and benchmarks"
```

---

### Task 21: Rewrite test generator for IGR-only

**File:** `toolchain/mfc/test/cases.py`

This is the most complex part of Stage 5. The test generator must be rewritten to only produce IGR-compatible test cases.

**Remove test generation functions:**
- `alter_weno()` — WENO test variants
- `alter_muscl()` — MUSCL test variants
- `alter_riemann_solvers()` — Riemann solver test variants
- `alter_low_Mach_correction()` — low Mach test variants
- `alter_ib()` — IBM test variants
- `alter_bubbles()` — bubble test variants
- `alter_lag_bubbles()` — Lagrangian bubble test variants
- `alter_hypoelasticity()` — hypoelastic test variants
- `alter_phasechange()` — phase change test variants
- `alter_acoustic_src()` — acoustic source test variants
- `alter_2d()` — remove axisymmetric/cylindrical parts (keep Cartesian 2D)
- `alter_3d()` — remove cylindrical parts
- `alter_grcbc()` — GRCBC test variants (characteristic BCs)
- `mhd_cases()` — entire MHD test suite
- `chemistry_cases()` — check if chemistry is IGR-compatible; if not, remove

**Modify test generation functions:**
- `alter_igr()` — becomes the primary solver variant generator. Expand to cover more IGR configurations.
- `alter_bcs()` — remove characteristic BC types (bc values -12 to -5)
- `alter_num_fluids()` — keep (IGR supports 1 and 2 fluids)
- `alter_viscosity()` — keep (IGR supports viscous flows)
- `alter_elliptic_smoothing()` — keep (IGR feature)
- `alter_body_forces()` — keep (IGR-compatible)
- `alter_mixlayer_perturb()` — keep if IGR-compatible
- `alter_bc_patches()` — keep, remove incompatible BC types

**Modify BASE_CFG:**
- Remove: `recon_type`, `weno_order`, `weno_eps`, `mapped_weno`, `wenoz`, `teno`, `mp_weno`, `null_weights`, `weno_Re_flux`, `weno_avg`, `riemann_solver`, `wave_speeds`, `avg_state`, `alt_soundspeed`
- Add: `igr_order`, `igr_iter_solver`, `num_igr_iters`, `num_igr_warm_start_iters`, `alf_factor`
- Remove: `bubbles_euler`, `ib`, etc. (no longer valid parameters)

**Modify `foreach_dimension()`:**
- Remove calls to deleted `alter_*` functions
- Ensure `alter_igr()` is called appropriately
- Remove 1D test generation (IGR requires 2D+) — or keep 1D if IGR supports it (check `case_validator.py` — the grid size check `m + 1 >= igr_order` allows 1D if m is large enough, but the `alter_igr` currently only runs for `len(dimInfo[0]) > 1`)

**Modify `foreach_example()`:**
- Update skip list to only skip examples that are too large for testing
- Remove non-IGR examples from the scan (they'll be deleted)

- [ ] **Step 1:** Read cases.py fully.
- [ ] **Step 2:** Update BASE_CFG for IGR-only defaults.
- [ ] **Step 3:** Remove dead `alter_*` functions and `mhd_cases()`.
- [ ] **Step 4:** Update `foreach_dimension()` to remove calls to dead functions.
- [ ] **Step 5:** Update `foreach_example()` skip list.
- [ ] **Step 6:** Verify test list generates without errors: `./figr.sh test -l`.

---

### Task 22: Delete all old golden files and regenerate

All existing golden files are for non-IGR tests. They need to be regenerated.

- [ ] **Step 1:** Delete all existing golden files.

```bash
git rm -r tests/
```

- [ ] **Step 2:** Regenerate golden files for the new IGR-only test suite.

```bash
./figr.sh test --generate -j 8
```

- [ ] **Step 3:** Verify the new tests pass.

```bash
./figr.sh test -j 8
```

- [ ] **Step 4:** Commit.

```bash
git add tests/ toolchain/mfc/test/cases.py
git commit -m "rewrite test suite for IGR-only, regenerate golden files"
```

---

### Task 23: Final verification and cleanup

- [ ] **Step 1:** Full precheck.

```bash
./figr.sh format -j 8
./figr.sh precheck -j 8
```

- [ ] **Step 2:** Full build.

```bash
./figr.sh build -j 8
```

- [ ] **Step 3:** Run full test suite.

```bash
./figr.sh test -j 8
```

- [ ] **Step 4:** Grep for any remaining references to deleted features that were missed.

```bash
grep -rn 'bubbles_euler\|bubbles_lagrange\|hypoelasticity\|hyperelasticity\|surface_tension\|acoustic_source\|\.not\. igr\|riemann_solver\|weno_order\|recon_type\|m_weno\|m_muscl\|m_riemann_solvers\|m_cbc\|m_viscous\|m_ibm\|m_bubbles\|m_qbmm\|m_phase_change\|m_acoustic_src\|m_hypoelastic\|m_hyperelastic\|m_surface_tension\|m_compute_cbc' src/ toolchain/
```

- [ ] **Step 5:** Fix any remaining references found.

- [ ] **Step 6:** Final commit.

```bash
git add -A
git commit -m "final cleanup: remove remaining dead references"
```

- [ ] **Step 7:** Push.

```bash
git push origin IGR-only
```
