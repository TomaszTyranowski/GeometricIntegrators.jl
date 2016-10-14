__precompile__()

module Integrators

    using ..Equations
    using ..Solvers
    using ..Utils

    import ..Solvers.function_stages!

    include("utils/macro_utils.jl")

    export Tableau, TableauRK, TableauERK, TableauDIRK, TableauFIRK, TableauSIRK,
           TableauPRK, TableauSARK, TableauSPARK, TableauGLM,
           showTableau, writeTableauToFile, readTableauERKFromFile

    include("integrators/tableaus.jl")

    export Solution, SolutionODE, SolutionPODE, SolutionDAE, SolutionPDAE,
           reset!, set_initial_conditions!,
           createHDF5, writeSolutionToHDF5

    include("integrators/solutions.jl")

    export Integrator, IntegratorERK, IntegratorDIRK, IntegratorFIRK,
           IntegratorPRK, IntegratorSARK, IntegratorSPARK,
           integrate, integrate!, function_stages!

    include("integrators/integrators.jl")
    include("integrators/integrators_erk.jl")
    include("integrators/integrators_dirk.jl")
    include("integrators/integrators_firk.jl")
    include("integrators/integrators_prk.jl")
    include("integrators/integrators_sark.jl")
    include("integrators/integrators_spark.jl")

end