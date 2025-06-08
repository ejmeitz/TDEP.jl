# TDEP.jl
Temperature Dependent Effective Potential Bindings in Julia

To install this package (and TDEP!) run in the Julia REPL:
```julia
using Pkg
Pkg.add(url = "https://github.com/ejmeitz/TDEP.jl.git")
```

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
- Re-implment checks for required input files to get useful errors without opening logfile
- Finish wrapping all executables
- Support local builds of TDEP/MPI
- Integrate with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl) ecosystem
- Possible integrations with [Brillouin.jl](https://github.com/thchr/Brillouin.jl) or [Spglib.jl](https://github.com/singularitti/Spglib.jl) or [SimpleCrystals.jl](https://github.com/ejmeitz/SimpleCrystals.jl)
- Do something cool with libolle (implement sTDEP loop etc.)
