# figr

**figr** is a GPU-accelerated finite volume solver for compressible multi-component flows using Information Geometric Regularization (IGR).
It's a miniapp based on [MFC](https://github.com/mflowcode/MFC).
IGR replaces traditional viscous or numerical shock-capturing methods with an elliptic regularization of the Euler equations, enabling high-order accurate representations of shock-laden flows.

## Method

figr discretizes the compressible Euler equations on structured Cartesian grids using a finite volume method.
Instead of solving Riemann problems at cell interfaces, IGR introduces a regularization parameter that smooths discontinuities via an elliptic PDE solve at each time step.
The method supports IGR orders 3 and 5, Jacobi and Gauss-Seidel iterative solvers, multi-component flows, viscous flows, 1D/2D/3D problems, MPI parallelism, and GPU acceleration via OpenACC and OpenMP target offload.

## References

- R. Cao and F. Schäfer, "Information geometric regularization of the barotropic Euler equation," *SIAM Journal on Scientific Computing*, 2026. [arXiv:2308.14127](https://arxiv.org/abs/2308.14127)
- R. Cao and F. Schäfer, "Information geometric regularization of unidimensional pressureless Euler equations yields global strong solutions," 2024. [arXiv:2411.15121](https://arxiv.org/abs/2411.15121)
- B. Wilfong, A. Radhakrishnan, H. Le Berre, D. J. Vickers, T. Prathi, N. Tselepidis, B. Dorschner, R. Budiardja, B. Cornille, S. Abbott, F. Schäfer, and S. H. Bryngelson, "Simulating many-engine spacecraft: Exceeding 1 quadrillion degrees of freedom via information geometric regularization," *Proceedings of SC '25: The International Conference for High Performance Computing, Networking, Storage and Analysis*, 14–24, 2025. [arXiv:2505.07392](https://arxiv.org/abs/2505.07392)

## Quick start

```bash
# Build (CPU)
./figr.sh build -j $(nproc)

# Build (GPU, OpenACC)
. ./figr.sh load -c <cluster> -m g
./figr.sh build -j $(nproc) --gpu acc

# Run an example
./figr.sh run examples/2D_IGR_double_mach/case.py -n 4

# Run tests
./figr.sh test -j 8
```

## License

MIT License
