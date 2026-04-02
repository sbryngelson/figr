"""
Minimal IGR test case generator.
"""

from .case import CaseGeneratorStack, define_case_d


def list_cases():
    stack = CaseGeneratorStack()
    cases = []

    # 2D IGR tests (periodic)
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
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "num_patches": 2,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%alpha_rho(1)": 0.125,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 0.1,
            "patch_icpp(2)%alter_patch(1)": "T",
        },
    )

    # IGR order 3, Jacobi
    stack.push(
        "1 Fluid",
        {
            "num_fluids": 1,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
        },
    )

    cases.append(
        define_case_d(
            stack,
            ["IGR order=3 Jacobi"],
            {
                "igr_order": 3,
                "igr_iter_solver": 1,
                "alf_factor": 10,
                "num_igr_iters": 10,
                "num_igr_warm_start_iters": 10,
                "elliptic_smoothing": "T",
                "elliptic_smoothing_iters": 10,
            },
        )
    )

    # IGR order 5, Jacobi
    cases.append(
        define_case_d(
            stack,
            ["IGR order=5 Jacobi"],
            {
                "igr_order": 5,
                "igr_iter_solver": 1,
                "alf_factor": 10,
                "num_igr_iters": 10,
                "num_igr_warm_start_iters": 10,
                "elliptic_smoothing": "T",
                "elliptic_smoothing_iters": 10,
            },
        )
    )

    # IGR order 5, Gauss-Seidel
    cases.append(
        define_case_d(
            stack,
            ["IGR order=5 GS"],
            {
                "igr_order": 5,
                "igr_iter_solver": 2,
                "alf_factor": 10,
                "num_igr_iters": 10,
                "num_igr_warm_start_iters": 10,
                "elliptic_smoothing": "T",
                "elliptic_smoothing_iters": 10,
            },
        )
    )

    # 2D viscous IGR
    stack.push(
        "Viscous",
        {
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1e4,
        },
    )
    cases.append(
        define_case_d(
            stack,
            ["IGR viscous order=5"],
            {
                "igr_order": 5,
                "igr_iter_solver": 1,
                "alf_factor": 10,
                "num_igr_iters": 10,
                "num_igr_warm_start_iters": 10,
                "elliptic_smoothing": "T",
                "elliptic_smoothing_iters": 10,
            },
        )
    )
    stack.pop()  # Viscous

    # 2D two-fluid IGR
    stack.pop()  # 1 Fluid
    stack.push(
        "2 Fluids",
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
            "patch_icpp(2)%alpha_rho(1)": 0.1,
            "patch_icpp(2)%alpha_rho(2)": 0.025,
            "patch_icpp(2)%alpha(1)": 0.8,
            "patch_icpp(2)%alpha(2)": 0.2,
        },
    )
    cases.append(
        define_case_d(
            stack,
            ["IGR 2fluid order=5"],
            {
                "igr_order": 5,
                "igr_iter_solver": 1,
                "alf_factor": 10,
                "num_igr_iters": 10,
                "num_igr_warm_start_iters": 10,
                "elliptic_smoothing": "T",
                "elliptic_smoothing_iters": 10,
            },
        )
    )
    stack.pop()  # 2 Fluids

    stack.pop()  # 2D

    # 3D IGR test (periodic)
    stack.push(
        "3D",
        {
            "m": 24,
            "n": 24,
            "p": 24,
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "z_domain%beg": 0.0,
            "z_domain%end": 1.0,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "num_patches": 2,
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%z_centroid": 0.5,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%length_z": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%z_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%length_z": 1.0,
            "patch_icpp(2)%alpha_rho(1)": 0.125,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 0.1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "num_fluids": 1,
            "fluid_pp(1)%gamma": 2.5,
            "fluid_pp(1)%pi_inf": 0.0,
        },
    )

    cases.append(
        define_case_d(
            stack,
            ["IGR 3D order=5"],
            {
                "dt": 1e-4,
                "igr_order": 5,
                "igr_iter_solver": 1,
                "alf_factor": 10,
                "num_igr_iters": 10,
                "num_igr_warm_start_iters": 10,
            },
        )
    )

    stack.pop()  # 3D

    return cases
