# TDEP.jl
Temperature Dependent Effective Potential Bindings in Julia

This project is still under development. At the moment only the `extract_forceconstants` command is verified to work. More to come soon!

```julia
using TDEP

efc = ExtractForceConstants(secondorder_cutoff = 2.0)
rundir = "/home/emeitz/tests/TDEP/tdep_jl"
execute(efc, rundir)
```


TO-DO:
--------
- Support MPI!!!
- Re-implment checks for required folders to get better errors
- Finish wrapping all executables
- Support local builds of TDEP/MPI
- Do something cool with libolle
