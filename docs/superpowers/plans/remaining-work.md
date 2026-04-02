# Remaining IGR-Only Refactor Work

## Current state
- Stages 1-3 DONE: All non-IGR Fortran deleted, all feature modules removed, igr boolean removed
- All three targets (pre_process, simulation, post_process) compile clean
- 9 commits on IGR-only branch

## Stage 4: Python toolchain (4 steps)

### Step 4a: Clean definitions.py
File: `toolchain/mfc/params/definitions.py`
- Remove parameter registrations for deleted features
- Remove CONSTRAINTS entries for deleted params  
- Remove DEPENDENCIES entries for deleted params
- Remove from CASE_OPT_PARAMS set
- Remove indexed families: patch_ib, acoustic, bub_pp
- Remove igr boolean (keep igr_order etc.)

### Step 4b: Clean case_validator.py
File: `toolchain/mfc/case_validator.py`
- Delete validation functions for removed features
- Simplify check_igr (remove igr toggle check)
- Simplify check_igr_simulation (remove checks for non-existent features)
- Remove calls to deleted functions from main validator

### Step 4c: Clean descriptions.py
File: `toolchain/mfc/params/descriptions.py`
- Remove descriptions for all deleted parameters

### Step 4d: Clean remaining toolchain files
- Grep toolchain/ for dead refs and fix
- Run `./mfc.sh lint` to verify
- Run `./mfc.sh validate examples/3D_IGR_TaylorGreenVortex/case.py`
- Commit

## Stage 5: Examples, tests, golden files (3 steps)

### Step 5a: Delete non-IGR examples and benchmarks
- Delete 131 non-IGR example dirs
- Delete 4 non-IGR benchmark dirs
- Clean IGR example case.py files (remove igr="T" param, remove alt_soundspeed etc.)
- Commit

### Step 5b: Rewrite test generator
File: `toolchain/mfc/test/cases.py`
- Update BASE_CFG for IGR defaults
- Remove dead alter_* functions
- Remove mhd_cases(), chemistry_cases() if incompatible
- Update foreach_dimension() and foreach_example()
- Verify: `./mfc.sh test -l`
- Commit

### Step 5c: Regenerate golden files
- Delete tests/ directory
- `./mfc.sh test --generate -j 8`
- `./mfc.sh test -j 8`
- Final commit and push
