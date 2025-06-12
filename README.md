# TDEP.jl
Temperature Dependent Effective Potential Bindings in Julia
> :warning: **This project is still under development** More to come soon!

----------------------------------------------
To install this package (and TDEP!) run in the Julia REPL:
```julia
using Pkg
Pkg.add(url = "https://github.com/ejmeitz/TDEP.jl.git")
```

```julia
using TDEP

efc = ExtractForceConstants(secondorder_cutoff = 2.0)
rundir = "/home/emeitz/tests/TDEP/tdep_jl"
execute(efc, rundir)
```

MPI is also supported. Unless otherwise specified, a comptabile MPI is automaticalled downloaded. To use it simply specify the number of cores when calling `execute`
```julia
ncores = 10
execute(efc, rundir, ncores)
```
It is possible to use a local MPI and/or a locally built TDEP. I will write docs on how to do that soon!

TO-DO:
--------
- Re-implment checks for required input files to get useful errors without opening logfile
- Finish wrapping all executables
- Support local builds of TDEP/MPI
- Integrate with [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl) ecosystem
- Possible integrations with [Brillouin.jl](https://github.com/thchr/Brillouin.jl) or [Spglib.jl](https://github.com/singularitti/Spglib.jl) or [SimpleCrystals.jl](https://github.com/ejmeitz/SimpleCrystals.jl)
- Do something cool with libolle (implement sTDEP loop etc.)
