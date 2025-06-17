module DFTExt

using TDEP
using Unitful
using DFTK

# Needed for sTDEP with DFT as force calculator
# Example params:
 # model_kwargs = (; functionals = LDA()),
 # basis_kwargs = (; Ecut=10, kgrid=(2, 2, 2)),
 # scf_kwargs = (; tol = 1e-4)
function TDEP.single_point_forces(
        sys::AbstractSystem{3}, 
        model_kwargs,
        basis_kwargs,
        scf_kwargs,
    )

    calc = DFTKCalculator(
        model_kwargs=model_kwargs,
        basis_kwargs=basis_kwargs,
        scf_kwargs=scf_kwargs
    )

    scfres = DFTK.compute_scf(sys, calc, nothing)
    forces = compute_forces_cart(scfres) * u"hartree/bohr"
    energy = scfres.energies.total * u"hartree" # just potential part since single point?

    return uconvert.(u"eV", energy), uconvert.(u"eV / Å", forces)
                           
end

end