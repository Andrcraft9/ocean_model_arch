### Example of ocean model architecture. 
### PSyKAl and Block Load-Balancing approach.
### Using modern Fortran (2003 at least).

### Shallow Water only example.

Structure:

```
control         - Algorithm Layer       : work with grid_type and ocean_type

core, interface - Parallel System Layer : data types and transfer data to kernel (core/kernel_interface.f90 - main envokes for kernels)

kernel          - Kernel Layer          : work with 1D, 2D and 3D fortran arrays, per one block



```

```
service   - Specific utills (grid and mask construction, etc)

shared    - Common utills (MPI, etc)

configs   - Common configs

tools     - Shell for legacy

legacy    - Original code from INMOM
```

I/O, Time Manager, Configs - from ocean model INMOM
