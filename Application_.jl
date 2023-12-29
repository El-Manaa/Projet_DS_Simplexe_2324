using LinearAlgebra:I,dot
using Crayons:@crayon_str,COLORS,CrayonWrapper
using Random:shuffle!
using PrettyTables

function type_problème()::Bool
    local s::String;
    while true
        s = Base.prompt("Choose the problem type ('max' or 'min')");
        if s in ("max","min")
            break
        end
    end
    return s == "max"
end

function type_inégalité(nc::Int)::Vector{Symbol}
    local s::String;
    local V = Vector{Union{Symbol,Nothing}}(nothing,nc);
    println("INFO: Each constraint written as (∑aᵢxᵢ 𝑅 b) must have-\nits non-strict inequality relation symbol 𝑅 written as (>= or <= or ==).")
    for i in 1:nc
        while true
            s = Base.prompt("Enter the inequality type of the constraint $i");
            try
                if s in ("<=",">=","==")
                    V[i] = Symbol(s);
                    break;
                else
                    throw(DomainError(s));
                end
            catch e
                if isa(e,DomainError)
                    println("Invalid type of inequality. Please try again.")
                end
            end
        end
    end
    return V |> Vector{Symbol}
end

function positivité_variables(nv::Int)::Vector{Symbol}
    local s::String;
    local V = Vector{Union{Symbol,Nothing}}(nothing,nv);
    println("INFO: Each variable's positivity symbol is (>= or <=).")
    for i in 1:nv
        while true
            s = Base.prompt("Enter the symbol of the variable $i");
            try
                if s in ("<=",">=")
                    V[i] = Symbol(s);
                    break
                else
                    throw(DomainError(s));
                end
            catch e
                if isa(e,DomainError)
                    println("Invalid symbol. Please try again.")
                end
            end
        end
    end
    return V |> Vector{Symbol}
end

function nbr_vars_cont()::Tuple{Int,Int}
    local nv::Int,nc::Int;
    while true
        s = Base.prompt("Enter the total number of your variables");
        try
            nv = parse(Int,s);
            if nv <= 0
                throw(DomainError(nv))
            end
            break
        catch e
            if isa(e,DomainError)
                println("Invalid number of variables. Please try again.")
            end
        end
    end
    while true
        s = Base.prompt("Enter the total number of your constraints");
        try
            nc = parse(Int,s);
            if nc <= 0
                throw(DomainError(nc))
            end
            break
        catch e
            if isa(e,DomainError)
                println("Invalid number of constraints. Please try again.")
            end
        end
    end
    return (nv,nc)
end

function saisirMatrice(nv::Int,nc::Int,est_max::Bool)::Matrix{Real}
    local M::Matrix{Real} = zeros(nc+1,nv+1);
    local s::String,ligne::Matrix{Real};
    println("INFO: A constraint { a₁x₁ + ... + aₙxₙ 𝑅 b , where 𝑅 = inequality relation} must be written as a vector [a₁ ... aₙ b].")
    for i in 1:nc
        while true
            s = Base.prompt("Enter the value of the constraints $i");
            try
                ligne = Meta.parse(s) |> eval |> Matrix{Real};
                if size(ligne,2) != nv+1
                    throw(ArgumentError(ligne))
                else
                    break
                end
            catch e
                if isa(e,ArgumentError)
                    println("The constraint $i is invalid. Please try again.")
                end
            end
        end
        M[i,1:nv+1] = ligne[1:nv+1];
    end

    println("INFO: An objective function {c₁x₁ + ... + cₙxₙ = $(est_max ? "Z(max)" : "W(min)")} must be written as a vector [c₁ ... cₙ 0]");
    local fo::Matrix{Real};
    while true
        s = Base.prompt("Enter the objective function");
        try
            fo = Meta.parse(s) |> eval |> Matrix{Real};
            if (size(fo,2) != nv+1) || (fo[end] != 0)
                throw(ArgumentError(s))
            else
                break
            end
        catch e
            if isa(e,ArgumentError)
                println("The objective function is invalid. Please try again.")
            end
        end
    end
    M[end,1:nv+1] = fo[1:nv+1]
    return M
end

function pivot(M::Matrix{Real},f::Function)
    local j = indexin(f(M[end,1:end-1]),M[end,1:end-1]);
    local rt = M[1:end-1,end] ./ M[1:end-1,j];
    local i = indexin(minimum(filter(x -> x >= 0,rt)),rt);
    return (i,j);
end

function est_vec_mat_id(V::Vector{Real})::Bool
    return (count(==(0),V) == length(V) - 1) && (count(==(1),V) == 1)
end

function standardiserSimplexe(M::Matrix{Real},est_max::Bool,v_var::Vector{Symbol},v_cont::Vector{Symbol})::Matrix{Real}
    A = copy(M)
    if est_max
        for i in 1:size(v_var,1)
            if v_var[i] == :<=
                A[:,i] .*= -1
            end
        end
        for i in 1:size(v_cont,1)
            if v_cont[i] == :(==)
                A = hcat(A[:,1:end-1],[Int(j == i) |> Real for j in 1:size(M,1)],[-Int(j == i) |> Real for j in 1:size(M,1)],A[:,end])
            elseif v_cont[i] == :>=
                A = hcat(A[:,1:end-1],[-Int(j == i) |> Real for j in 1:size(M,1)],A[:,end])
            else
                A = hcat(A[:,1:end-1],[Int(j == i) |> Real for j in 1:size(M,1)],A[:,end])
            end        
        end
    else
        A[end,:] .*= -1
        A = standardiserSimplexe(A,true,v_var,v_cont)
    end
    return A
end

function afficherSimplexe(M::Matrix{Real},est_max::Bool)
    local nopt = any((est_max ? (.>) : (.<))(M[end,:],0))
    local D = M |> x -> round.(x,digits=2) |> x -> string.(x);
    D[end,end] = (est_max ? "Z" : "W") * ("+" ^ (M[end,end] > 0)) * D[end,end] ^ !in(D[end,end],("0","0.0"));
    if nopt
        ij = pivot(M,(est_max ? maximum : minimum))
        pivot_l = Highlighter((data,i,j) -> (i == ij[1][1]) && (j == size(data)[2]),crayon"blue")
        pivot_c = Highlighter((data,i,j) -> (j == ij[2][1]) && (i == size(data)[1]),crayon"yellow")
        pivot_lc = Highlighter((data,i,j) -> (i == ij[1][1]) && (j == ij[2][1]),crayon"green")
        pretty_table(
            D;
            body_hlines=[size(M,1)-1,size(M,1)],
            linebreaks=true,
            header=append!(["x$(Char(0x2080) + i)" for i in 1:size(M[1,1:end-1])[1]],["b"]),
            highlighters = (pivot_l,pivot_c,pivot_lc)
        );
    else
        local Colors = filter(x -> x ∉ [
            :black,:white,:default,
            :dark_gray,:light_gray],keys(Crayons.COLORS)) |> collect
        shuffle!(Colors)
        local X :: Vector{Union{String,CrayonWrapper}} = ["x$(Char(0x2080) + i)" for i in 1:size(S[1,1:end-1])[1]];
        local hlt :: Vector{Highlighter} = [];
        local X_colors :: Vector{Crayon} = [];
        for j in 1:size(M)[2]-1
            if M[end,j] == 0 && est_vec_mat_id(M[1:end-1,j])
                i = findall(==(1),M[1:end-1,j])[1][1]
                push!(X_colors,Crayon(foreground=Colors[j],bold=true))
                #Colors = filter(x -> x != Colors[j],Colors)
                push!(hlt,Highlighter((data,x,y) -> (x == i) && y ∈ (j,size(M)[2]),X_colors[end])) 
            else
                push!(X_colors,Crayon(foreground=:default,bold=true))
            end
        end
        push!(X_colors,Crayon(foreground=:default))
        pretty_table(
            D;
            body_hlines=[size(M,1)-1,size(M,1)],
            linebreaks=true,
            header=append!(X,["b"]),
            header_crayon=X_colors,
            highlighters = hlt |> Tuple
        );
    end
end

function formaterBase(M::Matrix{Real})
    local V :: Vector{Int} = [];
    for j in 1:size(M)[2]-1
        if M[end,j] == 0 && est_vec_mat_id(M[1:end-1,j])
            push!(V,j)
        end
    end
    return "{$(join(V,", "))}"
end


function resoudre!(A::Matrix{Real},est_max::Bool)
    local k = 0
    println("iteration i = $k, basis J = $(formaterBase(A))")
    afficherSimplexe(A,est_max)
    while any((est_max ? (.>) : (.<))(A[end,1:end-1],0))
        ij = pivot(A,(est_max ? (maximum) : (minimum)));
        if A[ij...][1] == Real(0)
            break
        else
            A[ij[1]...,:] ./= A[ij...]
            for i in 1:size(A)[1]
                if i == ij[1]...
                    continue
                else
                    A[i,:] .-= A[ij[1]...,:] .* A[i,ij[2]...]
                end
            end
            k += 1;
            println()
            println("iteration i = $k, base J = $(formaterBase(A))")
            afficherSimplexe(A,est_max)
        end
    end
end

function est_réalisable(x::Vector{Real},M::Matrix{Real},v_vars::Vector{Symbol},v_cont::Vector{Symbol})::Bool
    return ((M[1:end-1,1:end-1] * x[1:end-1]) .|> (M[1:end-1,end] .|> eval.(v_cont)) |> all) && (x[1:end-1] .|> ([0] .|> eval.([v_vars;repeat([:>=],inner=length(v_cont))])) |> all)
end


function est_optimale(M::Matrix{Real},est_max::Bool)::Bool
    return all((est_max ? (.<=) : (.>=))(M[end,1:end-1],0))
end



function resultats(M::Matrix,N::Matrix,est_max::Bool,v_var::Vector{Symbol},v_cont::Vector{Symbol})
    local V = Vector{Real}(zeros(Real,size(M,2)))
    for j in 1:size(M,2)-1
        if M[end,j] == 0
            for i in 1:size(M)[1]
                if M[i,j] == 1
                    V[j] = round(M[i,end];digits=2);
                    break
                end
            end
        else
            V[j] = Real(0);
        end
    end
    
    for i in 1:size(V,1)-1
        print("x$i = $(V[i]); ");
    end
    V[end] = round(M[end,end] |> abs ;digits=2)
    println("\n$(est_max ? "Z" : "W")* = $(V[end])")
    if est_réalisable(V,N,v_var,v_cont)
        println("This solution is feasible.")
        if est_optimale(M,est_max)
            println("\n$(est_optimale(M,est_max) ? "And" : "But") the basis J = $(formaterBase(M)) $(est_optimale(M,est_max) ? "is" : "is not") an optimal basis.")
        end
    else
        println("No feasible solution.")
    end
end
