
struct RVSSSystem
    f::Dict{Vector, System}
    r::Vector{Pair{System, Any}}
    m::SciCompDSL.Model
end


function Base.getproperty(obj::RVSSSystem, symbol::Symbol)
    if hasproperty(obj, symbol)
        return getfield(obj, symbol)
    else
        f = getfield(obj, :f)
        for subsystem in values(f)
            if hasproperty(subsystem, symbol)
                try 
                    return getproperty(subsystem, symbol)
                catch e
                    nothing
                end
            end
        end
        error("Property $symbol not found in RVSSSystem.")
    end
end

mutable struct RVSSProblem
    f::RVSSSystem
    u0::Dict{Num, Float64}
    tspan::Tuple{Float64, Float64}
    p::Any
    kwargs::Any
    problem_type::Any
end

function RVSSProblem(f::RVSSSystem, u0::Vector, tspan::Tuple; kwargs...)
    RVSSProblem(f, Dict(u0), tspan, SciMLBase.NullParameters(), kwargs, VariableStructureSystemProblem())
end

function find_representative(sys)
    results = Vector()
    sp = ModelingToolkit.get_ps(sys)
    for p in sp
        if string(nameof(p)) == "_representative"
            val = ModelingToolkit.getdefault(p)
            push!(results, sys => val.representative)
            break 
        end
    end

    for subsys in ModelingToolkit.get_systems(sys)
        append!(results, find_representative(subsys))
    end
    return results
end

function recursiveRoleTypeFinder(dict::Dict)
    results::Vector{Type{<:Role}} = Vector{Type{<:Role}}()
    for (k, v) in dict
        if v isa Role
            push!(results, typeof(v))
        elseif v isa Dict
            push!(results, recursiveRoleTypeFinder(v)...)
        end
    end
    results
end

function recursiveRoleTypeFinder(dict::Nothing)
    results::Vector{Type{<:Role}} = Vector{Type{<:Role}}()
end

struct RVSSRepresentative
    representative::Any
end

function solve(problem::RVSSProblem, alg::DiffEqBase.AbstractODEAlgorithm; kwargs...)
    t_max::Float64 = problem.tspan[2]
    t_end::Float64 = 0.0
    solutions::Vector{ODESolution} = []
    roleDist::Vector = []
    representatives::Vector = problem.f.r
    prob::Union{ODEProblem, Nothing} = nothing
    activeModel::Union{System, Nothing} = nothing
    prevModel::Union{System, Nothing} = nothing
    modelProbMap::Dict{System, ODEProblem} = Dict()
    u0_varsMap::Dict{System, Any} = Dict()
    u0_vars::Any = nothing
    allVars = Dict()
    u0::Dict{Num, Any} = Dict()
    first::Bool = true
    u0_init::Dict{Num, Any} = Dict(problem.u0)
    initialKwargs = Dict([problem.kwargs...])
    while t_end < t_max
        prevModel = activeModel
        roleDist = []
        for p in representatives
            rs = getRoles(p[2])
            push!(roleDist, p[2] => recursiveRoleTypeFinder(rs))
        end
        if !haskey(problem.f.f, roleDist)
            @mtkcompile temp = problem.f.m()
            activeModel = get!(problem.f.f, roleDist, temp)
            u0_vars = unknowns(structural_simplify(ModelingToolkit.generate_initializesystem(activeModel); fully_determined = false))
            u0_varsMap[activeModel] = u0_vars
        else
            activeModel = problem.f.f[roleDist]
            u0_vars = get!(u0_varsMap, activeModel, unknowns(structural_simplify(ModelingToolkit.generate_initializesystem(activeModel); fully_determined = false)))
        end
        if first
            @info "Starting role-based VSS simulation."  
        elseif activeModel != prevModel
            @info "Detected role change at t = $(t_end). Restarting simulation." 
        end
        if solutions == []
            u0 = u0_init
        else
            u0 = Dict()
            for v in u0_vars
                if haskey(allVars, toexpr(v))
                    u0[v] = allVars[toexpr(v)]
                else
                    @warn "no initial value found for $v."
                end
            end
        end
        if activeModel in keys(modelProbMap)
            prob = modelProbMap[activeModel]
            prob = remake(modelProbMap[activeModel]; u0 = u0, tspan = (t_end, t_max), build_initializeprob=false)
        else
            g::Dict{Num, Any} = guesses(activeModel)
            for k_u0 in keys(u0)
                if haskey(g, k_u0) && !haskey(u0_init, k_u0)
                    delete!(u0, k_u0)
                end
            end
            if first && haskey(initialKwargs, :guesses)
                prob = ODEProblem(activeModel, u0, (t_end, t_max), guesses = [initialKwargs[:guesses]...]) #ODEProblem(activeModel, u0, (t_end, t_max), problem.kwargs...) # , build_initializeprob
            else
                prob = ODEProblem(activeModel, u0, (t_end, t_max), build_initializeprob=false)
            end
            modelProbMap[activeModel] = prob
        end
        first = false  
        integ = init(prob, alg; kwargs...)
        push!(solutions, ModelingToolkit.solve!(integ))
        
        if true
            for v in unknowns(activeModel)
                allVars[toexpr(v)] = solutions[end][v.metadata[Symbolics.VariableSource][2]][end]
            end
            for v in observables(activeModel)
                allVars[toexpr(v)] = SymbolicIndexingInterface.getsym(prob, v)(solutions[end])[end]#solutions[end][v.metadata[Symbolics.VariableSource][2]][end]
            end
        end
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
    stats::Union{Nothing, SciMLBase.DEStats} = solutions[1].stats
    t_sol::Vector{Float64} = solutions[1].t
    u::Vector{Vector{Float64}} = solutions[1].u
    symbolics_idx::Dict = Dict([toexpr(v) => SymbolicIndexingInterface.variable_index(solutions[1].prob.f.sys, v) for v in unknowns(solutions[1].prob.f.sys)])
    ###
    ### TO DO: Correct concatenation of u (works without by saving single solutions)
    ###
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