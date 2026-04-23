module ContextualEBM

import SymbolicIndexingInterface
using SciCompDSL, DiffEqBase, ModelingToolkitBase, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Contexts
using Satisfiability

include("VSS.jl")
export VSSSystem, VSSSolution, restart!, getRestartedIntegrators, VSSProblem, VariableStructureSystemProblem

include("CVSS.jl")
include("CVSSModel.jl")
export CVSSSystem, CVSSProblem, @cvssmodel, @cvsscompile, CVSSModel

include("RVSS.jl")
include("RVSSModel.jl")
export RVSSSystem, RVSSProblem, @rvssmodel, @cvsscompile, RVSSModel

export solve

end
