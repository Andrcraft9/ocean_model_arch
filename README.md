# Example of ocean model architecture. PSyKAl and Block Load-Balancing approach.
# Using modern Fortran (2003 at least).

# Shallow Water only example.

Structure:

```
control   - Algorithm Layer       : work with grid_type and ocean_type

interface - Parallel System Layer : transfer data to kernel

kernel    - Kernel Layer          : work with 1D, 2D and 3D fortran arrays, per one block

core      - Data types
```

```
service   - Specific utills (grid and mask construction, etc)

shared    - Common utills (MPI, configs, etc)

legacy    - Code from INMOM
```

I/O, Time Manager, Configs - from ocean model INMOM
