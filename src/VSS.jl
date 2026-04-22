struct VariableStructureSystemProblem end

mutable struct VSSProblem
    f::Any
    m::SciCompDSL.Model
    u0::Vector
    tspan::Tuple
    p::Any
    kwargs::Any
    problem_type::Any
end

function VSSProblem(f, m, u0::Vector, tspan::Tuple; kwargs...)
    VSSProblem(f, m, u0, tspan, SciMLBase.NullParameters(), kwargs, VariableStructureSystemProblem())
end

function restartDef()
    restartedIntegrators = []
    function restart!(integ::SciMLBase.DEIntegrator)
        push!(restartedIntegrators, integ)
        terminate!(integ)
    end
    function getRestartedIntegrators()
        return restartedIntegrators
    end
    return restart!, getRestartedIntegrators
end

restart!, getRestartedIntegrators = restartDef()

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
                push!(v, zeros(size(s.u)[1]).+ NaN)
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