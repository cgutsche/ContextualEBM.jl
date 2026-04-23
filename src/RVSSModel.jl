macro rvssmodel(name, expr)
	esc(_rvss_model_macro(__module__, name, expr))
end

function argList(expr::Symbol)
	[expr]
end
function argList(expr::Expr)
	if expr.args[1] == :(!)
		return [expr]
	else
		return reduce(vcat, [filter(x -> (x != :(&)) & (x != :(|)), argList.(expr.args))...])
	end
end


function _rvss_model_macro(mod, name, expr)
	Base.remove_linenums!(expr)
	skeletonExpr::Expr = :()
	skeletonCodeMap = Dict()
	roleTeamMap = Dict()
	roleMap = Dict()
	returnExpr = quote end
	for arg in expr.args
		if arg.head == :macrocall
			Base.remove_linenums!(arg.args)
			if arg.args[1] == Symbol("@base")
				if skeletonExpr != :()
					error("More than one base has been specified.")
				end
				modExpr = arg.args[end]
				Base.remove_linenums!(modExpr)
				skeletonExpr = modExpr
			elseif arg.args[1] == Symbol("@role")
				modExpr = arg.args
                if modExpr[3] == :not
                    if modExpr[4] isa Expr
                        if !(modExpr[4].args[1] == :in)
                            error("When adding Team types to role block definitions, use \"<role> in <team>\" syntax.")
                        end
                        roleTeamMap[modExpr[4].args[2]] = modExpr[4].args[3]
                        Base.remove_linenums!.(modExpr[5].args)
                        roleMap[modExpr[4].args[2] => false] = modExpr[5]
                    else
                        Base.remove_linenums!.(modExpr[5].args)
                        roleMap[modExpr[4] => false] = modExpr[5]
                    end
                else
                    if modExpr[3] isa Expr
                        if !(modExpr[3].args[1] == :in)
                            error("When adding Team types to role block definitions, use \"<role> in <team>\" syntax.")
                        end
                        roleTeamMap[modExpr[3].args[2]] = modExpr[3].args[3]
                        Base.remove_linenums!.(modExpr[4].args)
                        roleMap[modExpr[3].args[2] => true] = modExpr[4]
                    else
                        Base.remove_linenums!.(modExpr[4].args)
                        roleMap[modExpr[3] => true] = modExpr[4]
                    end
                end
			else
				error("Blocks must begin with @base, @role, or @representative.")
			end
		end
	end
	for arg in skeletonExpr.args
		skeletonCodeMap[arg.args[1]] = arg
	end
    hasStrucParams = haskey(skeletonCodeMap, Symbol("@structural_parameters"))
    push!(get!(skeletonCodeMap, Symbol("@structural_parameters"), :(@structural_parameters begin end)).args[end].args, :(representative = nothing))
    if isempty(roleTeamMap)
        push!(skeletonCodeMap[Symbol("@structural_parameters")].args[end].args, :(team = nothing))
    end
    if !hasStrucParams
        push!(skeletonExpr.args, skeletonCodeMap[Symbol("@structural_parameters")])
    end
    hasParams = haskey(skeletonCodeMap, Symbol("@constants"))
    push!(get!(skeletonCodeMap, Symbol("@constants"), :(@constants begin end)).args[end].args, :(_representative::RVSSRepresentative = RVSSRepresentative(representative)))
    if isempty(roleTeamMap)
        push!(skeletonCodeMap[Symbol("@constants")].args[end].args, :(_team::Contexts.AbstractTeam = team))
    end
    if !hasParams
        push!(skeletonExpr.args, skeletonCodeMap[Symbol("@constants")])
    end
    for (rolePair, content) in roleMap
        role = rolePair[1]
        playsRole = rolePair[2]
		for arg in content.args
			flag = false
            for (skeletonContentType, skeletonContent) in skeletonCodeMap
                if arg.args[1] == skeletonContentType
                    if isempty(roleTeamMap)
                        if playsRole
                            expr = quote
                                    if hasRole(representative, $role, team)
                                    end
                                end
                        else
                            expr = quote
                                    if !hasRole(representative, $role, team)
                                    end
                                   end
                        end
                    else
                        if haskey(roleTeamMap, role)
                            teamtype = roleTeamMap[role]
                        else
                            error("Team type specification missing for role $role.")
                        end
                        if playsRole
                            expr = quote
                                    if hasRole(representative, $role, $teamtype)
                                    end
                                end
                        else
                            expr = quote
                                if !hasRole(representative, $role, $teamtype)
                                end
                            end
                        end
                    end
                    if expr.args[2].args[2] isa Expr
                        push!(expr.args[2].args[2].args, arg.args[3].args...)
                    else
                        expr.args[2].args[2] = arg.args[3]
                    end
                    push!(skeletonContent.args[3].args, expr.args[2])
                    flag = true
				end
			end
            if !flag
                if isempty(roleTeamMap)
                    if playsRole
                        expr = quote
                                if hasRole(representative, $role, team)
                                end
                            end
                    else
                        expr = quote
                                if !hasRole(representative, $role, team)
                                end
                                end
                    end
                else
                    if haskey(roleTeamMap, role)
                        teamtype = roleTeamMap[role]
                    else
                        error("Team type specification missing for role $role.")
                    end
                    if playsRole
                        expr = quote
                                if hasRole(representative, $role, $teamtype)
                                end
                            end
                    else
                        expr = quote
                                if !hasRole(representative, $role, $teamtype)
                                end
                                end
                    end
                end
                push!(expr.args[2].args[2].args, arg.args[3].args...)
                newContent = deepcopy(collect(values(skeletonCodeMap))[1])
                newContent.args = newContent.args[1:3]
                newContent.args[1] = arg.args[1]
                newContent.args[3] = expr
                push!(skeletonExpr.args, newContent)
                skeletonCodeMap[arg.args[1]] = newContent
            end
		end
	end
    println(skeletonExpr)
	push!(returnExpr.args, SciCompDSL._model_macro(mod, name, skeletonExpr, false))
	returnExpr
end

function _rvvs_compile_helper(name, sym, kwargs)
    map = Dict()
    temp = sym(; name=Symbol("$name"))
    r::Vector{Pair{System, Any}} = find_representative(temp)
    roleDist = []
    for p in r
        rs = getRoles(p[2])
        push!(roleDist, p[2] => recursiveRoleTypeFinder(rs))
    end
    map[roleDist] = eval(quote eval(Expr(:call, mtkcompile, $kwargs..., $temp)) end)
    return map, r, sym
end

macro rvsscompile(exprs...)
    expr = exprs[1]
    name = expr.args[1]
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
    sym = (Symbol(string(split(repr(expr.args[2]), "()")[1][3:end])))
    namesym = Meta.parse(":$name")
    esc(quote
        $name = RVSSSystem(_rvvs_compile_helper($namesym, $sym, $kwargs)...)
    end)
end