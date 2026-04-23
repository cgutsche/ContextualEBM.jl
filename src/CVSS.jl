
struct CVSSSystem
    f::Dict{Contexts.AbstractContext, System}
end

function Base.getproperty(obj::CVSSSystem, symbol::Symbol)
    if hasproperty(obj, symbol)
        return getfield(obj, symbol)
    else
        f = getfield(obj, :f)
        for subsystem in values(f)
            if hasproperty(subsystem, symbol)
                return getproperty(subsystem, symbol)
            end
        end
        error("Property $symbol not found in CVSSSystem.")
    end
end

mutable struct CVSSProblem
    f::Dict{Union{Contexts.AbstractContext}, ODEProblem}
    u0::Dict{Num, Float64}
    tspan::Tuple{Float64, Float64}
    p::Any
    kwargs::Any
    problem_type::Any
    variableMap::Dict{System, Set{Expr}}
    initVariables::Dict{System, Vector{Num}}
end

function CVSSProblem(f::CVSSSystem, u0::Vector, tspan::Tuple; kwargs...)
	problems = Dict{Contexts.AbstractContext, ODEProblem}()
	stateMigrationVars = Dict{System, Set{Expr}}()
	initVars = Dict{System, Vector{Num}}()
	u0 = Dict(u0)
	for (context, model) in f.f
		u0_vars = unknowns(model)
		initVars[model] = u0_vars
		u0_prob = [var => var in keys(u0) ? u0[var] : 0.0 for var in u0_vars]
		problems[context] = ODEProblem(model, u0_prob, tspan; build_initializeprob= false)
		stateMigrationVars[model] = Set([toexpr.(unknowns(model)); toexpr.(map(eq -> eq.lhs, observed(model)))])
	end
	CVSSProblem(problems, u0, tspan, SciMLBase.NullParameters(), kwargs, VariableStructureSystemProblem(), stateMigrationVars, initVars)
end

function solve(problem::CVSSProblem, alg::DiffEqBase.AbstractODEAlgorithm; printLog = true, kwargs...)
    t_max::Float64 = problem.tspan[2]
    t_end::Float64 = 0.0
    solutions::Vector{ODESolution} = []
    contexts::Vector{Contexts.AbstractContext} = collect(keys(problem.f))
    modelProbMap::Dict{System, ODEProblem} = Dict()
    allVars::Dict{Expr, Float64} = Dict()
    warn = true
    u0_init::Dict{Num, Any} = Dict(problem.u0)
    initialKwargs = Dict([problem.kwargs...])
    first = true
    variablesToSave::Vector{Num} = union(collect(values(problem.initVariables))...)
    sizehint!(allVars, length(problem.u0))
    stateMigrationVars = problem.variableMap
    while t_end < t_max
        activeContexts::Vector{Contexts.AbstractContext} = Vector{Contexts.AbstractContext}([filter(x -> isActive(x), contexts)...])
        if length(activeContexts) == 1
            nothing
        elseif length(activeContexts) == 0
            @warn "No active Contexts. Simulation stopped before time span ended."
            break
        else
            @warn "Active contexts do not specify exactly 1 active model. Only first possible model is simulated. You might want to define the context-to-model assignment disjunctly."
        end
        if printLog
            @info "Simulating model for context < $(activeContexts[1]) > at t = $(t_end)." 
        end
        activeModel::System = problem.f[activeContexts[1]].f.sys
        if solutions == []
            u0 = u0_init
            g::Dict{Num, Any} = guesses(activeModel)
            for k_u0 in keys(u0)
                if haskey(g, k_u0) && !haskey(u0_init, k_u0)
                    delete!(u0, k_u0)
                end
            end
            prob = ODEProblem(activeModel, u0, (t_end, t_max), initialKwargs...)
        else
            u0 = Dict()
            for v in problem.initVariables[activeModel]
                if haskey(allVars, toexpr(v))
                    u0[v] = allVars[toexpr(v)]
                end
            end
            prob = remake(problem.f[activeContexts[1]]; u0 = u0, tspan = (t_end, t_max), build_initializeprob=false)

        end
        integ = init(prob, alg; kwargs...)
        push!(solutions,  ModelingToolkit.solve!(integ))
        
        sol = solutions[end]
        for var in variablesToSave
            varExpr = toexpr(var)
            if varExpr in stateMigrationVars[activeModel]
                allVars[varExpr] = sol[var][end]
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
    stats::SciMLBase.DEStats = solutions[1].stats
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