module ContextualEBM

import SymbolicIndexingInterface
using SciCompDSL, DiffEqBase

inlcude("VSS.jl")
export VSSSystem, VSSSolution, restart!, getRestartedIntegrators, VSSProblem, VariableStructureSystemProblem

include("CVSS.jl")
export CVSSSystem, CVSSProblem

include("RVSS.jl")
export RVSSSystem, RVSSProblem

export solve

end
