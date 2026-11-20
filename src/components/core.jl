#=
Lightweight component layer (feature/component-core).

A component declares
- its internal unknowns: `unknowns(c)`, a NamedTuple name => VariableTransform;
- its ports: `ports(c)`, a NamedTuple name => VariableTransform (values shared with other
  components through nodes);
- its residuals: `residuals!(r, c, x, p)` on physical values (`x`: internal unknowns, `p`: ports,
  both NamedTuples), unscaled, and a nominal magnitude for each one, `residual_scales(c)`;
- a physical initial guess for its internal unknowns, `guess(c; kwargs...)`.

A `ComponentSystem` joins components through nodes: each node is one unknown shared by the ports
it connects, or a fixed value. It maps the unknowns from ℝ with the transforms, calls the
components and divides each residual by its nominal magnitude.
=#

"""
    AbstractComponent

Supertype of the components assembled by [`ComponentSystem`](@ref). A component implements
[`unknowns`](@ref), [`ports`](@ref), [`residual_scales`](@ref), [`residuals!`](@ref) and, if it
has internal unknowns, [`guess`](@ref).
"""
abstract type AbstractComponent end

"""
    unknowns(c::AbstractComponent) -> NamedTuple

Internal unknowns of the component, name => [`VariableTransform`](@ref), in residual order.
"""
unknowns(::AbstractComponent) = (;)

"""
    ports(c::AbstractComponent) -> NamedTuple

Ports of the component, name => [`VariableTransform`](@ref). A port takes the value of the
node it is connected to.
"""
ports(::AbstractComponent) = (;)

"""
    residual_scales(c::AbstractComponent) -> Tuple

Nominal magnitude of each residual of the component (same units as the residual). The system
solves `residual / scale = 0`.
"""
function residual_scales end

"""
    residuals!(r, c::AbstractComponent, x::NamedTuple, p::NamedTuple)

Write the (unscaled) residuals of the component into `r`, from the physical values of its
internal unknowns `x` and of its ports `p`.
"""
function residuals! end

"""
    guess(c::AbstractComponent; kwargs...) -> NamedTuple

Physical values of the internal unknowns of the component, for an initial guess.
"""
guess(c::AbstractComponent; kwargs...) = (;)

nresiduals(c::AbstractComponent) = length(residual_scales(c))

# one component inside a system: where its unknowns, residuals and port nodes are
struct ComponentBlock{C <: AbstractComponent, P <: Tuple, S <: Tuple}
    component::C
    uoffset::Int            # internal unknowns are u[uoffset .+ (1:n)]
    roffset::Int            # residuals are du[roffset .+ (1:m)]
    portnodes::P            # node index of each port
    scales::S
end

"""
    ComponentSystem(components::NamedTuple; nodes::NamedTuple, fixed::NamedTuple = (;))

Square nonlinear system from `components` (name => component). `nodes` connects the ports:
node name => tuple of `component_name => port_name` pairs; every port belongs to exactly one
node. A node listed in `fixed` (node name => physical value) is a known value, the others are
unknowns.

The unknowns, in ℝ, are a `Vector{Float64}`: the internal unknowns of each component (in the
order of `components`), then the free nodes. Build it with [`initial_guess`](@ref), evaluate with
[`residual!`](@ref), and read the physical values by name with [`physical`](@ref). The type of the
system is inferred from the types of the components, so code that builds a system inside a
function stays concretely typed (node transforms: all ports of a given transform type, e.g.
`Free()`, give a concrete type).

Example, evaporator coupled to a borehole:

    sys = ComponentSystem((evaporator = evap, borehole = bhx);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin),
                 T_warm = (:evaporator => :Ts_in,  :borehole => :Tout),
                 T_ground = (:borehole => :Tb,)),
        fixed = (T_ground = 281.15,))         # undisturbed ground, e.g. ≈ 8 °C (Stockholm)

Fixed values are changed with [`set_fixed`](@ref).
"""
struct ComponentSystem{C <: NamedTuple, B <: Tuple, NN <: Val, NT <: Tuple, F <: Tuple, I <: Tuple}
    components::C
    blocks::B
    nodenames::NN           # Val((names...)): node names in the type, so that results keyed by them are inferred
    nodetransforms::NT
    nodefixed::F            # `nothing` (free) or the fixed value, per node
    nodeindex::I            # index in u of each free node (0 if fixed)
    n::Int
end

function ComponentSystem(components::NamedTuple; nodes::NamedTuple, fixed::NamedTuple = (;))
    nodenames = keys(nodes)
    for k in keys(fixed)
        k in nodenames || throw(ArgumentError("fixed node $k is not in `nodes`"))
    end
    # each port in exactly one node
    owner = Dict{Tuple{Symbol, Symbol}, Int}()
    for (j, k) in enumerate(nodenames), pr in nodes[k]
        cname, pname = pr
        haskey(components, cname) || throw(ArgumentError("node $k: no component $cname"))
        haskey(ports(components[cname]), pname) || throw(ArgumentError("node $k: component $cname has no port $pname"))
        haskey(owner, (cname, pname)) && throw(ArgumentError("port $cname.$pname is in two nodes"))
        owner[(cname, pname)] = j
    end
    for cname in keys(components), pname in keys(ports(components[cname]))
        haskey(owner, (cname, pname)) || throw(ArgumentError("port $cname.$pname is not connected to a node"))
    end
    # node transform: the one of its first port; all ports of a node must agree
    nodetransforms = map(nodenames) do k
        ts = [ports(components[c])[p] for (c, p) in nodes[k]]
        all(==(ts[1]), ts) || throw(ArgumentError("node $k joins ports with different transforms: $ts"))
        ts[1]
    end
    nodefixed = map(k -> haskey(fixed, k) ? Float64(fixed[k]) : nothing, nodenames)

    # layout of u: components with internal unknowns, then the free nodes. Built by recursion
    # over the tuples (no counters captured in closures), so that the system type is inferred.
    blocks = _blocks(keys(components), values(components), owner, 0, 0)
    nu, nr = _block_sizes(blocks)
    nodeindex = _node_index(nodefixed, nu + 1)
    n = nu + count(v -> v === nothing, nodefixed)
    n == nr || throw(ArgumentError("system is not square: $n unknowns, $nr residuals"))
    return ComponentSystem(components, blocks, Val(nodenames), nodetransforms, nodefixed, nodeindex, n)
end

_blocks(::Tuple{}, ::Tuple{}, owner, uoff, roff) = ()
function _blocks(names::Tuple, comps::Tuple, owner, uoff, roff)
    c, cname = first(comps), first(names)
    b = ComponentBlock(c, uoff, roff, map(p -> owner[(cname, p)], keys(ports(c))), map(float, residual_scales(c)))
    return (b, _blocks(Base.tail(names), Base.tail(comps), owner, uoff + length(unknowns(c)), roff + nresiduals(c))...)
end

# number of internal unknowns and of residuals of all blocks
_block_sizes(blocks) = (sum(b -> length(unknowns(b.component)), blocks; init = 0),
                        sum(b -> length(b.scales), blocks; init = 0))

# index in u of each free node (0 for a fixed node), from `next` on
_node_index(::Tuple{}, next) = ()
_node_index(fixed::Tuple, next) = first(fixed) === nothing ?
    (next, _node_index(Base.tail(fixed), next + 1)...) : (0, _node_index(Base.tail(fixed), next)...)

Base.length(sys::ComponentSystem) = sys.n

"""
    set_fixed(sys::ComponentSystem; node = value, ...) -> ComponentSystem

Copy of `sys` with new values of fixed nodes (e.g. the ground temperature at each time step of a
dynamic simulation). The copy has the same type as `sys`, so nothing is compiled again, and the
unknowns keep their layout: the previous solution is a warm start.
"""
function set_fixed(sys::ComponentSystem; values...)
    names = _nodenames(sys)
    for k in keys(values)
        i = findfirst(==(k), names)
        i === nothing && throw(ArgumentError("set_fixed: no node $k"))
        sys.nodefixed[i] === nothing && throw(ArgumentError("set_fixed: node $k is free, not fixed"))
    end
    nodefixed = map(names, sys.nodefixed) do k, v
        v === nothing ? v : (haskey(values, k) ? Float64(values[k]) : v)
    end
    return ComponentSystem(sys.components, sys.blocks, sys.nodenames, sys.nodetransforms, nodefixed, sys.nodeindex, sys.n)
end

function Base.show(io::IO, sys::ComponentSystem)
    print(io, "ComponentSystem(", join(keys(sys.components), ", "), "; ", sys.n, " unknowns, nodes: ",
          join((v === nothing ? string(k) : "$k = $v" for (k, v) in zip(_nodenames(sys), sys.nodefixed)), ", "), ")")
end

# physical value of every node, in the element type of u
_node_values(sys::ComponentSystem, u) =
    map(sys.nodetransforms, sys.nodefixed, sys.nodeindex) do t, v, i
        v === nothing ? from_R(t, u[i]) : convert(eltype(u), v)
    end

_internal_values(b::ComponentBlock, u) =
    _internal_values(unknowns(b.component), b.uoffset, u)
_internal_values(ts::NamedTuple{K}, off, u) where {K} =
    NamedTuple{K}(map((t, i) -> from_R(t, u[off + i]), values(ts), ntuple(identity, length(ts))))

_port_values(b::ComponentBlock, nodevals) =
    NamedTuple{keys(ports(b.component))}(map(j -> nodevals[j], b.portnodes))

"""
    residual!(du, u, sys::ComponentSystem)

Scaled residuals of the system at the unknowns `u` (in ℝ, layout as in
[`ComponentSystem`](@ref)). Signature compatible with `NonlinearProblem(f!, u0, sys)`.
"""
function residual!(du, u, sys::ComponentSystem)
    nodevals = _node_values(sys, u)
    map(sys.blocks) do b
        m = length(b.scales)
        r = view(du, b.roffset+1:b.roffset+m)
        residuals!(r, b.component, _internal_values(b, u), _port_values(b, nodevals))
        for i in 1:m
            r[i] /= b.scales[i]
        end
        nothing
    end
    return nothing
end

"""
    initial_guess(sys::ComponentSystem; nodes = (;), <component name> = (; kwargs...), ...)

Unknowns in ℝ (a `Vector{Float64}`) from physical values: for each component with internal
unknowns, the keyword arguments of its [`guess`](@ref) method; `nodes`, the values of the free
nodes.
"""
function initial_guess(sys::ComponentSystem; nodes::NamedTuple = (;), kwargs...)
    u = zeros(sys.n)
    # blocks are in the order of the components (no identity search: `===` on a large inline
    # immutable block was seen to return false for the same value on Julia 1.13). Recursion over
    # the tuples, not a loop: a loop over more than three component types widens the block to
    # an abstract type, and `guess` (property calls, Clapeyron) was inferred on `Any` (cycle of
    # four components: 560 s of compilation).
    _guess_blocks!(u, keys(sys.components), sys.blocks, _component_kwargs(sys.components, values(kwargs)))
    for (k, t, v, i) in zip(_nodenames(sys), sys.nodetransforms, sys.nodefixed, sys.nodeindex)
        v === nothing || continue
        haskey(nodes, k) || throw(ArgumentError("no initial value for node $k"))
        u[i] = to_R(t, nodes[k])
    end
    return u
end

# keyword arguments of `guess` for each component (an empty NamedTuple when not given), in the
# order of the components; generated so that each element has a concrete type
@generated function _component_kwargs(::NamedTuple{K}, kw::NamedTuple{KW}) where {K, KW}
    return Expr(:tuple, (k in KW ? :(kw.$k) : :((;)) for k in K)...)
end

_guess_blocks!(u, ::Tuple{}, ::Tuple{}, ::Tuple{}) = u
function _guess_blocks!(u, names::Tuple, blocks::Tuple, kws::Tuple)
    _guess_block!(u, first(names), first(blocks), first(kws))
    return _guess_blocks!(u, Base.tail(names), Base.tail(blocks), Base.tail(kws))
end

function _guess_block!(u, cname, b::ComponentBlock, kw::NamedTuple)
    ts = unknowns(b.component)
    isempty(ts) && return nothing
    x = guess(b.component; kw...)
    keys(x) == keys(ts) || throw(ArgumentError("guess for $cname has keys $(keys(x)), expected $(keys(ts))"))
    v = map(to_R, values(ts), values(x))               # tuple map: each transform concretely typed
    for i in eachindex(v)
        u[b.uoffset + i] = v[i]
    end
    return nothing
end

"""
    physical(sys::ComponentSystem, u) -> NamedTuple

Physical values of the unknowns `u`: one NamedTuple per component with internal unknowns, and
`nodes` with the values of all nodes (free and fixed).
"""
function physical(sys::ComponentSystem, u)
    nodevals = _node_values(sys, u)
    comps = _drop_empty(NamedTuple{keys(sys.components)}(map(b -> _internal_values(b, u), sys.blocks)))
    return merge(comps, (nodes = NamedTuple{_nodenames(sys)}(nodevals),))
end

_nodenames(::ComponentSystem{C, B, Val{NN}}) where {C, B, NN} = NN

# drop the entries without internal unknowns (empty NamedTuples), at the type level so that the
# result is inferred
_drop_empty(nt::NamedTuple{()}) = nt
function _drop_empty(nt::NamedTuple{K}) where {K}
    rest = _drop_empty(Base.tail(nt))
    return first(nt) isa NamedTuple{()} ? rest : merge(NamedTuple{(first(K),)}((first(nt),)), rest)
end

"""
    port_values(sys::ComponentSystem, u, name::Symbol) -> NamedTuple

Physical values of the ports of component `name` at the unknowns `u`.
"""
function port_values(sys::ComponentSystem, u, name::Symbol)
    b = sys.blocks[findfirst(==(name), keys(sys.components))]
    return _port_values(b, _node_values(sys, u))
end

########################
# validity             #
########################

"""
    validity(c::AbstractComponent, x, p) -> Vector

Physical-validity issues of component `c` at the physical values `x` (internal unknowns) and `p`
(ports), as NamedTuples `(quantity, value, limit, message)`. Empty if none. Reported only: the
solver is not constrained by them.
"""
validity(::AbstractComponent, x, p) = NamedTuple[]

"""
    validity_issues(sys::ComponentSystem, u) -> Vector

Validity issues of all components at the unknowns `u` (e.g. a solution), as NamedTuples
`(component, quantity, value, limit, message)`; empty if none. See [`validity`](@ref).
"""
function validity_issues(sys::ComponentSystem, u)
    nodevals = _node_values(sys, u)
    issues = NamedTuple[]
    for (name, b) in zip(keys(sys.components), sys.blocks)
        for i in validity(b.component, _internal_values(b, u), _port_values(b, nodevals))
            push!(issues, (component = name, i...))
        end
    end
    return issues
end

_issue(quantity, value, limit, message) = (quantity = quantity, value = value, limit = limit, message = message)

# secondary-fluid temperatures: freezing point and tabulated range
function _secondary_issues!(issues, sf::SecondaryFluid, temps::NamedTuple)
    for (k, T) in pairs(temps)
        if isfinite(sf.Tfreeze) && T < sf.Tfreeze
            push!(issues, _issue(k, T, sf.Tfreeze, "secondary fluid below its freezing point"))
        end
        if T < sf.Tmin || T > sf.Tmax
            push!(issues, _issue(k, T, T < sf.Tmin ? sf.Tmin : sf.Tmax,
                                 "secondary fluid outside its property table (properties extrapolated flat)"))
        end
    end
    return issues
end

########################
# continuation         #
########################

"""
    continuation(build, λs, u0, alg; secant = true, max_halvings = 6, kwargs...) -> Vector

Solve the systems `build(λ)` (a [`ComponentSystem`](@ref) for each parameter value `λ`, e.g. a
flow factor or a source temperature) for `λ` in `λs`, in order, starting from the guess `u0` for
`build(first(λs))`. Each solve starts from the previous solution, extrapolated linearly in `λ`
from the last two solutions (`secant = true`). When a step fails, it is halved (at most
`max_halvings` times) and the intermediate values are solved first.

Returns one NamedTuple `(λ, u, ok, iterations)` per value of `λs` solved (the continuation
stops after the first value that fails; its `u` is then the last converged solution, at an
intermediate `λ`). `kwargs` go to
`solve`. Needs NonlinearSolve (package extension).
"""
function continuation end
