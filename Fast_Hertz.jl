module Fast_Hertz

using JuLIP
using JuLIP.Potentials: pairs
import JuLIP: energy, forces, cutoff
using NeighbourLists
using StaticArrays

const SVecG{T} = SVector{3, T}
const SMatG{T} = SMatrix{3, 3, T, 9}

export Fast_Hertz_Potential

@pot struct Fast_HP <: AbstractCalculator
    rcut::Float64    # cutoff when calculating neighbourlist
    rtol2::Float64   # if max_displacement > rtol, recalculate the neighbourlist rtol2 = rtol^2
    rsum::Vector{Float64} # rsum = [2*rmin; rmin + rmax; 2*rmax]
    Lhalf::Float64
end

function Fast_Hertz_Potential(L, r0::Float64, r1::Float64; rcutfact = 2.5)
    rcut = rcutfact * r0 # rcut is rcutfact*2*R_{max}, used when searching possible pairs
    rtol = 0.5 * rcut    # rtol is the neighbourlist recomputing tolerance
    rtol2 = rtol * rtol
    rsum = [r1; (r0 + r1)/2; r0]
    Lhalf = L/2.0
    return Fast_HP(rcut, rtol2, rsum, Lhalf)
end

function energy(V::Fast_HP, at::AbstractAtoms)
    nlist_i, nlist_j, nlist_Z = update_neighbourlist!(V, at)
    bins = at.cell'
    inv_bins = inv(bins)
    return energy_inner(nlist_i, nlist_j, nlist_Z, at.X, V, bins, inv_bins)
end
# compute the maximum displacement of all particles
@inline max_dr2(X::Vector{JVecF}, X_Old::Vector{JVecF}) = 
                    maximum( [ sum( (X[i] - X_Old[i]) .* (X[i] - X_Old[i]) ) for i = 1:length(X) ] )

@inline PairType(i::Vector{Int}, j::Vector{Int}, Z::Vector{Int}) = [ Z[i[k]] + Z[j[k]] - 1 for k = 1:length(i) ]
                    

function update_neighbourlist!(V::Fast_HP, at::AbstractAtoms)
    if has_data(at, :nlist_X)
        X_Old = get_data(at, :nlist_X)::Vector{JVecF}
        
        X = at.X::Vector{JVecF}
        max_d = max_dr2(X, X_Old)
        if max_d < V.rtol2
            nlist_i = get_data(at, :nlist_i)::Vector{Int}
            nlist_j = get_data(at, :nlist_j)::Vector{Int}
            nlist_Z = get_data(at, :nlist_Z)::Vector{Int}
            return nlist_i, nlist_j, nlist_Z
        end
    end
    nlist = neighbourlist(at, V.rcut)
    nlist_Z = PairType(nlist.i, nlist.j, at.Z)
    set_data!(at, :nlist_X, copy(at.X))
    set_data!(at, :nlist_i, nlist.i)
    set_data!(at, :nlist_j, nlist.j)
    set_data!(at, :nlist_Z, nlist_Z)
    return nlist.i, nlist.j, nlist_Z
end

@inline wrap_to_unit(inv_bins::SMatG{T}, x::SVecG{T}) where {T <: Real} = mod.(inv_bins * x, 1)
@inline adjust_dx(bins::SMatG{T}, dx::SVecG{T}) where {T <: Real} = bins * (dx - (dx .> 0.5))
@inline dx2dx(dx::SVecG{T}, inv_bins::SMatG{T}, bins::SMatG{T}) where {T <: Real} = 
                                    adjust_dx(bins, wrap_to_unit(inv_bins, dx))                


function energy_inner(nlist_i, nlist_j, nlist_Z, X, V::Fast_HP, bins, inv_bins)
    E = 0.0
    for k = 1:length(nlist_i)
        dx = X[nlist_j[k]] - X[nlist_i[k]]
        r = norm(dx)
        if r > V.Lhalf
            dx = dx2dx(dx, inv_bins, bins)
            r = norm(dx)
        end
        nlist_s = r - V.rsum[nlist_Z[k]]
        E += 1.0 * nlist_s * nlist_s * (nlist_s < 0)
    end
    return 0.5*E
end


end
