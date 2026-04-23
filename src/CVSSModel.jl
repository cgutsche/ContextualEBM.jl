

macro cvssmodel(name, expr)
	esc(_cvss_model_macro(__module__, name, expr))
end

struct CVSSModel
    f::Dict{Contexts.AbstractContext, SciCompDSL.Model}
end
(m::CVSSModel)() = m.f

function extractArgs(expr::Vector{Union{Symbol, Expr}})
	function contextList(e::Symbol)
		return [e]  # Return a single-element array containing the symbol
	end
	
	function contextList(e::Expr)
		list = contextList.(filter(x -> (x != :(!)) & (x != :(&)) & (x != :(|)), e.args))
		if list == []
			return [e]  # If no arguments, return the expression itself
		end
		list = reduce(vcat, list)
		return list
	end
	# Extracts all arguments from a vector of expressions
	args = Set{Symbol}()
	for e in expr
		symbols = contextList(e)
		push!(args, symbols...)
	end
	return collect(args)
end

function get_atomic_models(exprs::Vector{Union{Symbol, Expr}})
	function alternative(exprs::BoolExpr...)
		alts = []
		for expr in exprs
			push!(alts, and(expr, and([not.(filter(x -> repr(x) != repr(expr), exprs))...])))
		end
		or(alts)
	end
	function exclusion(exprs::BoolExpr...)
		exprs = collect(exprs)
		excls = []
		for i in 1:(length(exprs)-1)
			for j in (i+1):length(exprs)
				push!(excls, or(not(exprs[i]), not(exprs[j])))
			end
		end
		if length(excls) == 1
			return excls[1]
		end
		and(excls)
	end
	function free(exprs::BoolExpr...)
		or([not(exprs[1]), exprs...])
	end
	function requirement(expr1::BoolExpr, expr2::BoolExpr)
		or(not(expr1), expr2)
	end
	function strongInclusion(expr1::BoolExpr, expr2::BoolExpr)
		or(not(expr1), expr2)
	end
	function getContextsFromConstraint(constraint::Contexts.Constraint)
		contextsInConstraints = Meta.parse.([String.(s...) for s in split.(repr.(constraint.contexts), "ContextType()"; keepempty=false)])
	end
	function addConstraint(constraint::Contexts.Exclusion, contextList::Vector{Symbol})
		exclusion([contexts[var_map[context]] for context in contextList]...)
	end
	function addConstraint(constraint::Contexts.Alternative, contextList::Vector{Symbol})
		alternative([contexts[var_map[context]] for context in contextList]...)
	end
	function addConstraint(constraint::Contexts.Requirement, contextList::Vector{Symbol})
		requirement(contexts[var_map[contextList[1]]], contexts[var_map[contextList[2]]])
	end
	function addConstraint(constraint::Contexts.Inclusion, contextList::Vector{Symbol})
		strongInclusion(contexts[var_map[contextList[1]]], contexts[var_map[contextList[2]]])
		
	end
	function convertExpr(expr::Union{Symbol, Expr})
		if expr isa Symbol
			return contexts[var_map[expr]]
		elseif expr isa Expr
			if expr.head == :(call)
				if expr.args[1] == :(!)
					return not(convertExpr(expr.args[2]))
				elseif expr.args[1] == :(&)
					return and(convertExpr.(expr.args[2:end])...)
				elseif expr.args[1] == :(|)
					return or(convertExpr.(expr.args[2:end])...)
				else
					error("Unsupported expression type: $(typeof(expr))")
				end
			else
				error("Unsupported expression type: $(typeof(expr))")
			end
		end
	end
	function generateExpr(assignments::Vector{Dict{Symbol, Bool}})
		if length(assignments) < 1
			return :()
		end
		expr = nothing
		for (var, value) in assignments[1]
			if isnothing(expr)
				if value
					expr = var
				else
					expr = Expr(:call, :(!), var)
				end
			else
				if value
					expr = Expr(:call, :(&), var, expr)
				else
					subexpr = Expr(:call, :(!), var)
					expr = Expr(:call, :(&), subexpr, expr)
				end
			end
		end
		if length(assignments) > 1
			for assignment in assignments[2:end]
				local_expr = nothing
				for (var, value) in assignment
					if isnothing(local_expr)
						if value
							local_expr = var
						else
							local_expr = Expr(:call, :(!), var)
						end
					else
						if value
							local_expr = Expr(:call, :(&), var, local_expr)
						else
							subexpr = Expr(:call, :(!), var)
							local_expr = Expr(:call, :(&), subexpr, local_expr)
						end
					end
				end
				expr = Expr(:call, :(|), expr, local_expr)
			end
		end
		return expr	
	end

	variables = extractArgs(exprs)

	# translate Contexts.jl-constraints to boolean expressions

	constraintList = Contexts.getConstraints()
	contextsInConstraints = nothing
	for constraint in constraintList
		contextsInConstraints = getContextsFromConstraint(constraint)
		if any(x -> x in variables, contextsInConstraints)
			if !issubset(contextsInConstraints, variables)
				push!(variables, filter(x -> !(x in variables), contextsInConstraints)...)
			end
		end
	end

	n = size(variables, 1)
	var_map = Dict{Symbol, Int}(variables[i] => i for i in 1:n)
	@satvariable(contexts[1:n], Bool)

	constraints = []
	containedContexts = Set{Symbol}()
	for constraint in constraintList
		contextsInConstraints = getContextsFromConstraint(constraint)
		if all(x -> x in variables, contextsInConstraints)
			newConstraint = addConstraint(constraint, contextsInConstraints)
			push!(constraints, newConstraint)
			push!(containedContexts, contextsInConstraints...)
		end
	end
	
	# Add all contexts that are not in the constraints
	missingContexts = filter(x -> !(x in containedContexts), variables)
	if length(missingContexts) != 0
		push!(constraints, free([contexts[var_map[c]] for c in missingContexts]...))
	end
	

	# Convert expressions to boolean variables
	conditions = Dict([expr => convertExpr(expr) for expr in exprs])
	solutions::Dict{Dict, Vector{Union{Expr, Symbol}}} = Dict()
	for (expr, cond) in conditions
		local_constraints = copy(constraints)
		# Add the condition to the constraints
		push!(local_constraints, cond)
		open(Z3()) do interactive_solver 
			assert!(interactive_solver, local_constraints...)
			i = 1
			status, assignment = sat!(interactive_solver)
			while status == :SAT
				# Try to solve the problem
				assign!(contexts, assignment)
				if haskey(solutions, Dict([k => value(contexts[v]) for (k, v) in var_map]))
					push!(solutions[Dict([k => value(contexts[v]) for (k, v) in var_map])], expr)
				else
					solutions[Dict([k => value(contexts[v]) for (k, v) in var_map])] = [expr]
				end
				# Use assert! to exclude the solution we just found. 
				assert!(interactive_solver, not(and(contexts .== value(contexts))))
				status, assignment = sat!(interactive_solver)
				i += 1
			end
		end
	end
	
	# Check if the "all false" case (skeleton only) is valid
	open(Z3()) do interactive_solver
		assert!(interactive_solver, constraints...)
		# Check if it's satisfiable with all contexts false
		push_assertion = and(contexts .== false)
		assert!(interactive_solver, push_assertion)
		status, assignment = sat!(interactive_solver)
		if status == :SAT
			# All contexts false is a valid solution
			all_false_dict = Dict([k => false for (k, v) in var_map])
			if !haskey(solutions, all_false_dict)
				# Mark this as "no context active" - use an empty list of active contexts
				solutions[all_false_dict] = Symbol[]
			end
		end
	end
	
	model_conditions = Dict()
	for (c, d) in solutions
		if haskey(model_conditions, d)
			push!(model_conditions[d], c)
		else
			model_conditions[d] = [c]
		end
	end
	return Dict([d => generateExpr(c) for (d, c) in model_conditions])
end


function _cvss_model_macro(mod, name, expr)
	Base.remove_linenums!(expr)
	skeletonExpr::Expr = :()
	contextMap = Dict()
	returnExpr = quote end
	for arg in expr.args
		if arg.head == :macrocall
			Base.remove_linenums!(arg.args)
			if arg.args[1] == Symbol("@skeleton")
				if skeletonExpr != :()
					error("More than one Skeleton has been specified.")
				end
				modExpr = arg.args[end]
				Base.remove_linenums!(modExpr)
				skeletonExpr = modExpr
			elseif arg.args[1] == Symbol("@context")
				modExpr = arg.args
				Base.remove_linenums!.(modExpr[4].args)
				contextMap[modExpr[3]] = modExpr[4]
			else
				error("Blocks must begin with @skeleton or @context.")
			end
		end
	end

	exprs::Vector{Union{Symbol, Expr}} = collect(keys(contextMap))
	atomic_models = get_atomic_models(exprs)
	dictExpr = Expr(:call, :Dict)
	i = 1
	for (models, contexts) in atomic_models
		viewName = Symbol(name, "_mode_$i")
		i += 1
		modelContent = copy(skeletonExpr)
		for model in models 
			contextContent = copy(contextMap[model])
			contentData = Dict()
			for arg in contextContent.args
				contentData[arg.args[1]] = arg
			end
			for (contentType, content) in contentData
				flag = false
				for arg in modelContent.args
					if arg.args[1] == contentType
						push!(arg.args[3].args, (content.args[3].args)...)
						flag = true
					end
				end
				if !flag
					push!(modelContent.args, content)
				end
			end
		end
		push!(dictExpr.args, :($contexts => $viewName))
		push!(returnExpr.args, SciCompDSL._model_macro(mod, viewName, modelContent, false))
	end
	if keys(atomic_models) == []
		push!(returnExpr.args, SciCompDSL._model_macro(mod, name, skeletonExpr, false))
	end
	dictReturnExpr = Expr(:(=), :($name), :(CVSSModel($dictExpr)))
	push!(returnExpr.args, dictReturnExpr)
	
	returnExpr
end

macro cvsscompile(exprs...)
    expr = exprs[1]
    name = expr.args[1]
    name2 = QuoteNode(Symbol(expr.args[1], "_"))
    kwargs = Base.tail(exprs)
    kwargs = map(kwargs) do ex
        @assert ex.head == :(=)
        Expr(:kw, ex.args[1], ex.args[2])
    end
    if isempty(kwargs)
        kwargs = ()
    else
        kwargs = (Expr(:parameters, kwargs...),)
    end
    esc(quote
        $name = Dict()
        for (i, (context, model)) in enumerate($(expr.args[2]))
            temp = model(; name=Symbol($name2, i))
            $name[context] = eval(Expr(:call, mtkcompile, $kwargs..., temp))
        end
        $name = CVSSSystem($name)
    end)
end