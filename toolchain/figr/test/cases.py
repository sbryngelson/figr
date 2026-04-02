"""
IGR test case generator.

Uses the upstream Sod shock tube base configuration (3-patch, ghost
extrapolation BCs) adapted for the IGR-only mini-app.
"""

from .case import CaseGeneratorStack, define_case_d


def list_cases():
    stack = CaseGeneratorStack()
    cases = []

    # 2D base: Sod shock tube with ghost-extrapolation BCs
    stack.push(
        "2D",
        {
            "m": 49,
            "n": 39,
            "p": 0,
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "num_patches": 3,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.05,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 0.1,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%y_centroid": 0.45,
            "patch_icpp(2)%length_x": 1.0,
            "patch_icpp(2)%length_y": 0.7,
            "patch_icpp(2)%alpha_rho(1)": 0.5,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 0.5,
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%x_centroid": 0.5,
            "patch_icpp(3)%y_centroid": 0.9,
            "patch_icpp(3)%length_x": 1.0,
            "patch_icpp(3)%length_y": 0.2,
            "patch_icpp(3)%alpha_rho(1)": 0.125,
            "patch_icpp(3)%alpha(1)": 1.0,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": 0.1,
        },
    )

    # 1-fluid IGR tests
    stack.push(
        "1 Fluid(s)",
        {
            "num_fluids": 1,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
        },
    )

    igr_base = {
        "alf_factor": 10,
        "num_igr_iters": 10,
        "num_igr_warm_start_iters": 10,
    }

    # IGR order 3, Jacobi
    cases.append(
        define_case_d(
            stack,
            ["IGR", "igr_order=3", "Jacobi"],
            {**igr_base, "igr_order": 3, "igr_iter_solver": 1},
        )
    )

    # IGR order 5, Jacobi
    cases.append(
        define_case_d(
            stack,
            ["IGR", "igr_order=5", "Jacobi"],
            {**igr_base, "igr_order": 5, "igr_iter_solver": 1},
        )
    )

    # IGR order 5, Gauss-Seidel
    cases.append(
        define_case_d(
            stack,
            ["IGR", "igr_order=5", "Gauss Seidel"],
            {**igr_base, "igr_order": 5, "igr_iter_solver": 2},
        )
    )

    # Viscous IGR
    stack.push("Viscous", {"viscous": "T", "fluid_pp(1)%Re(1)": 1e4})

    for order in [3, 5]:
        cases.append(
            define_case_d(
                stack,
                ["IGR", f"igr_order={order}", "Jacobi"],
                {**igr_base, "igr_order": order, "igr_iter_solver": 1},
            )
        )

    stack.pop()  # Viscous

    # 2-fluid IGR
    stack.pop()  # 1 Fluid(s)
    stack.push(
        "2 Fluid(s)",
        {
            "num_fluids": 2,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.5,
            "fluid_pp(2)%pi_inf": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 0.8,
            "patch_icpp(1)%alpha_rho(2)": 0.2,
            "patch_icpp(1)%alpha(1)": 0.8,
            "patch_icpp(1)%alpha(2)": 0.2,
            "patch_icpp(2)%alpha_rho(1)": 0.4,
            "patch_icpp(2)%alpha_rho(2)": 0.1,
            "patch_icpp(2)%alpha(1)": 0.8,
            "patch_icpp(2)%alpha(2)": 0.2,
            "patch_icpp(3)%alpha_rho(1)": 0.1,
            "patch_icpp(3)%alpha_rho(2)": 0.025,
            "patch_icpp(3)%alpha(1)": 0.8,
            "patch_icpp(3)%alpha(2)": 0.2,
        },
    )

    for order in [3, 5]:
        cases.append(
            define_case_d(
                stack,
                ["IGR", f"igr_order={order}", "Jacobi"],
                {**igr_base, "igr_order": order, "igr_iter_solver": 1},
            )
        )

    stack.pop()  # 2 Fluid(s)
    stack.pop()  # 2D

    # 3D IGR tests
    stack.push(
        "3D",
        {
            "m": 29,
            "n": 29,
            "p": 29,
            "dt": 1e-4,
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "z_domain%beg": 0.0,
            "z_domain%end": 1.0,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "num_patches": 3,
            "num_fluids": 1,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%z_centroid": 0.125,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%length_z": 0.25,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%x_centroid": 0.5,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%z_centroid": 0.5,
            "patch_icpp(2)%length_x": 1.0,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%length_z": 0.5,
            "patch_icpp(2)%alpha_rho(1)": 0.5,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 0.5,
            "patch_icpp(3)%geometry": 9,
            "patch_icpp(3)%x_centroid": 0.5,
            "patch_icpp(3)%y_centroid": 0.5,
            "patch_icpp(3)%z_centroid": 0.875,
            "patch_icpp(3)%length_x": 1.0,
            "patch_icpp(3)%length_y": 1.0,
            "patch_icpp(3)%length_z": 0.25,
            "patch_icpp(3)%alpha_rho(1)": 0.125,
            "patch_icpp(3)%alpha(1)": 1.0,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%vel(3)": 0.0,
            "patch_icpp(3)%pres": 0.1,
        },
    )

    for order in [3, 5]:
        cases.append(
            define_case_d(
                stack,
                ["1 Fluid(s)", "IGR", f"igr_order={order}", "Jacobi"],
                {**igr_base, "igr_order": order, "igr_iter_solver": 1},
            )
        )

    stack.pop()  # 3D

    return cases
