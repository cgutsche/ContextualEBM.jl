import SymbolicIndexingInterface

struct VariableStructureSystemProblem end

function restartDef()
    restartedIntegrators = []
    function restart!(integ::SciMLBase.DEIntegrator)
        push!(restartedIntegrators, integ)
        terminate!(integ)
    end
    function getRestartedIntegrators()
        return copy(restartedIntegrators)
    end
    return restart!, getRestartedIntegrators
end

restart!, getRestartedIntegrators = restartDef()


mutable struct VSSProblem
    f::Any
    m::ModelingToolkit.Model
    u0::Vector
    tspan::Tuple
    p::Any
    kwargs::Any
    problem_type::Any
end

function VSSProblem(f, m, u0::Vector, tspan::Tuple; kwargs...)
    VSSProblem(f, m, u0, tspan, SciMLBase.NullParameters(), kwargs, VariableStructureSystemProblem())
end

mutable struct CVSSProblem
    f::Dict{Union{Contexts.Context, Contexts.AbstractContextRule}, ODESystem}
    u0::Dict{Num, Float64}
    tspan::Tuple{Float64, Float64}
    p::Any
    kwargs::Any
    problem_type::Any
end

function CVSSProblem(f::Dict, u0::Vector, tspan::Tuple; kwargs...)
    CVSSProblem(f, Dict(u0), tspan, SciMLBase.NullParameters(), kwargs, VariableStructureSystemProblem())
end

struct VSSSolution{T, N, uType, uType2, DType, tType, rateType, P, A, IType, S, AC <: Vector{<:Union{Nothing, Vector{Int}}}, R, O, Sols}  <: SciMLBase.AbstractODESolution{T, N, uType}
    u::uType
    u_analytic::uType2
    errors::DType
    t::tType
    k::rateType
    prob::P
    alg::A
    interp::IType
    dense::Bool
    tslocation::Int
    stats::S
    statsVector::Vector{S}
    alg_choice::AC
    retcode::Vector{SciMLBase.ReturnCode.T}
    resid::R
    original::O
    singleSolutions::Sols
end

function solve(problem::CVSSProblem, alg::DiffEqBase.AbstractODEAlgorithm; kwargs...)
    t_max::Float64 = problem.tspan[2]
    t_end::Float64 = 0.0
    solutions::Vector{ODESolution} = []
    contexts::Vector{Union{Contexts.Context, Contexts.AbstractContextRule}} = collect(keys(problem.f))
    modelProbMap::Dict = Dict()
    while t_end < t_max
        activeContexts::Vector{Union{Contexts.Context, Contexts.AbstractContextRule}} = Context[filter(x -> isActive(x), contexts)...]
        if length(activeContexts) == 1
            nothing
        elseif length(activeContexts) == 0
            @warn "No active Contexts. Simulation stopped before time span ended."
            break
        else
            @warn "Active contexts do not specify exactly 1 active model. Only first possible model is simulated. You might want to define the context-to-model assignment disjunctly."
        end
        @info "Simulating model for context < $(activeContexts[1]) > at t = $(t_end)." 
        activeModel::ODESystem = problem.f[activeContexts[1]]
        u0 = solutions == [] ? problem.u0 : Dict(filter(x -> x !== nothing, [toexpr(v) in toexpr.(unknowns((solutions[end]).prob.f.sys)) ? v => solutions[end][v.metadata[Symbolics.VariableSource][2]][end] : nothing for v in unknowns(activeModel)]))    
        if activeModel in keys(modelProbMap)
            #prob = ODEProblem(activeModel, u0, (t_end, t_max))
            prob = remake(modelProbMap[activeModel], u0 = u0, tspan = (t_end, t_max))
        else
            prob = ODEProblem(activeModel, u0, (t_end, t_max))
            modelProbMap[activeModel] = prob
        end
        integ = init(prob, alg; kwargs...)
        push!(solutions,  ModelingToolkit.solve!(integ))
        t_end = solutions[end][t, end]
        if solutions[end].retcode == ReturnCode.Terminated
            restartedIntegrators = getRestartedIntegrators()
            if integ in restartedIntegrators
                filter!(x -> x !== integ, restartedIntegrators)
            else
                break
            end
        else
            break
        end
    end

    stats::SciMLBase.DEStats = solutions[1].stats
    t_sol::Vector{Float64} = solutions[1].t
    u::Vector{Vector{Float64}} = solutions[1].u
    symbolics_idx::Dict = Dict([toexpr(v) => SymbolicIndexingInterface.variable_index(solutions[1].prob.f.sys, v) for v in unknowns(solutions[1].prob.f.sys)])
    """println(symbolics_idx)
    symbolics = collect(keys(symbolics_idx))
    u_new = nothing
    for sol in solutions[2:end]
        stats = SciMLBase.merge(stats, sol.stats)
        t_sol = [t_sol; sol.t]
        u_new = copy(sol.u)
        symbolics_idx_new = Dict([toexpr(v) => SymbolicIndexingInterface.variable_index(sol.prob.f.sys, v) for v in unknowns(sol.prob.f.sys)])
        println("A: ", symbolics_idx_new)
        new_variables = filter(v -> !(toexpr(v) in toexpr.(symbolics)),  unknowns(sol.prob.f.sys))
        missing_variables = filter(v -> !(toexpr(v) in toexpr.(unknowns(sol.prob.f.sys))), symbolics)
        for v in new_variables
            println()
            symbolics_idx[toexpr(v)] = max(collect(values(symbolics_idx)) + 1)
            push!(u_new[max(collect(values(symbolics_idx)) + 1)], fill(NaN, length(sol.u[1]))...)
        end
        for (s, j) in symbolics_idx
            if toexpr(s) in toexpr.(missing_variables)
                println(j)
                push!(u_new[j], fill(NaN, length(sol.u[1]))...)
            else
                println(j)
                push!(u_new[j], sol.u[symbolics_idx_new[toexpr(s)]]...)
            end
        end
        u = [u; u_new]
    end"""

    probs = [sol.prob for sol in solutions]
    alg = [sol.alg for sol in solutions]
    interp = [sol.interp for sol in solutions]
    dense  = all([sol.dense for sol in solutions])
    tslocation = 0
    statsVector = [sol.stats for sol in solutions]
    alg_choice = [sol.alg_choice for sol in solutions]
    retcode = [sol.retcode for sol in solutions]
    resid = [sol.resid for sol in solutions]
    original = [sol.original for sol in solutions]

    T = eltype(eltype(u))
    N = length((size(u[1])..., length(u)))


    VSSSolution{T, N, typeof(u), Nothing, Nothing, typeof(t_sol), typeof([]), typeof(probs), 
                typeof(alg), typeof(interp), typeof(stats), typeof(alg_choice), typeof(resid), 
                typeof(original), typeof(solutions)}(u, nothing, nothing, t_sol, [], probs, alg, interp, dense,
                                  tslocation, stats, statsVector, alg_choice, retcode, resid, original, solutions)

end

function solve(problem::VSSProblem, solver=Rosenbrock23(); kwargs...)
    t_max = problem.tspan[2]
    
    sol_temp = nothing
    t_end = problem.tspan[1]
    solutions::Vector{ODESolution} = []
    @mtkbuild system = (problem.m)()
    prob = ODEProblem(system, problem.u0, (t_end, t_max))
    sol_temp = ModelingToolkit.solve(prob, solver; kwargs...)
    t_end = sol_temp[t, end]
    push!(solutions, sol_temp)
    stats = solutions[1].stats
    t_sol = solutions[1].t
    u = solutions[1].u
    symbolics = unknowns(solutions[1].prob.f.sys)
    u_new = nothing
    while true
        if t_end == t_max
            break
        end
        @mtkbuild system = (problem.m)()
        global system = system
        u0 = Dict(filter(x -> x !== nothing, [toexpr(v) in toexpr.(unknowns(sol_temp.prob.f.sys)) ? v => sol_temp[v.metadata[Symbolics.VariableSource][2]][end] : nothing for v in unknowns(system)]))
        prob = remake(ODEProblem(system, u0, (t_end, t_max)))
        sol_temp = ModelingToolkit.solve(prob, solver; kwargs...)
        t_end = sol_temp[t, end]
        push!(solutions, sol_temp)
        
        stats = SciMLBase.merge(stats, sol_temp.stats)
        t_sol = [t_sol; sol_temp.t]
        u_new = [sol_temp.u[i] for i in eachindex(sol_temp.u)]
        
        u = [u; u_new]
    end

    prob = [sol.prob for sol in solutions]
    alg = [sol.alg for sol in solutions]
    interp = [sol.interp for sol in solutions]
    dense  = all([sol.dense for sol in solutions])
    tslocation = 0
    stats = stats
    statsVector = [sol.stats for sol in solutions]
    alg_choice = [sol.alg_choice for sol in solutions]
    retcode = [sol.retcode for sol in solutions]
    resid = [sol.resid for sol in solutions]
    original = [sol.original for sol in solutions]

    T = eltype(eltype(u))
    N = length((size(u[1])..., length(u)))

    VSSSolution{T, N, typeof(u), Nothing, Nothing, typeof(t_sol), typeof([]), typeof(prob), 
                typeof(alg), typeof(interp), typeof(stats), typeof(alg_choice), typeof(resid), 
                typeof(original), typeof(solutions)}(u, nothing, nothing, t_sol, [], prob, alg, interp, dense,
                                  tslocation, stats, statsVector, alg_choice, retcode, resid, original, solutions)

end

function Base.getindex(sol::VSSSolution, sym::Num)
    v = []

    for s in sol.singleSolutions
        try
            push!(v, s[sym])
        catch e
            if (SymbolicIndexingInterface.is_independent_variable(s.prob.f.sys, sym) | 
                       SymbolicIndexingInterface.is_observed(s.prob.f.sys, sym) | 
                       SymbolicIndexingInterface.is_variable(s.prob.f.sys, sym))
                push!(v, zeros(length(s.u[1])).+ NaN)
            else
                throw("ArgumentError: $sym is neither an observed nor an unknown variable of one of the contextual models.")
            end
        end
    end
    vcat(v...)
end


function VSSSystem(m, u0::Vector, tspan::Tuple; kwargs...)
    @mtkbuild system = m()
    (system, VSSProblem(system, m, u0, tspan, SciMLBase.NullParameters(), kwargs, VariableStructureSystemProblem()))
end
