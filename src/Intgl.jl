
using QCBase
using SpecialFunctions
using Combinatorics: doublefactorial
using LinearAlgebra: norm, eigen
using StaticArrays
const ang2bohr = 1.8897261246257702



function cartesian_dimesion_all(L::Integer)
    ((L + 1) * (L + 2) * (L + 3)) ÷ 6
end

function cartesian_dimension_oneShell(L::Integer)
    ((L + 1) * (L + 2)) ÷ 2
end


@inline function E(i::Int64, j::Int64, t::Int64, Qx::Float64, a::Float64, b::Float64, n::Int64=0, Ax::Float64=0.0)::Float64
    p = a + b
    u = a*b/p
    if n == 0
        if (t < 0) || (t > (i + j))
            return 0.0
        elseif i == j == t == 0
            return exp(-u*Qx*Qx)
        elseif j == 0
            return (1/(2*p))*E(i-1,j,t-1,Qx,a,b) - (u*Qx/a)*E(i-1,j,t,Qx,a,b) + \
                   (t+1)*E(i-1,j,t+1,Qx,a,b)
        else
            return (1/(2*p))*E(i,j-1,t-1,Qx,a,b) + (u*Qx/b)*E(i,j-1,t,Qx,a,b) + \
                   (t+1)*E(i,j-1,t+1,Qx,a,b)
        end
    else
        return E(i+1,j,t,Qx,a,b,n-1,Ax) + Ax*E(i,j,t,Qx,a,b,n-1,Ax)
    end
end


function boys(m::Float64, T::Float64)::Float64
    return hyp1f1(m+0.5, m+1.5, -T)/(2.0*m+1.0)
end


function gaussian_product_center(a::Float64, A, b::Float64, B)
    return (a.*A .+ b.*B) ./ (a+b)
end



@inline function R(t::Int64, u::Int64, v::Int64, n::Int64, p::Float64, PCx::Float64, PCy::Float64, PCz::Float64, RPC::Float64)::Float64
    T = p * RPC * RPC
    val = 0.0
    if t == 0 && u == 0 && v == 0
        val += (-2 * p)^n * boys(n, T)
    elseif t == u == 0
        if v > 1
            val += (v - 1) * R(t, u, v - 2, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCz * R(t, u, v - 1, n + 1, p, PCx, PCy, PCz, RPC)
    elseif t == 0
        if u > 1
            val += (u - 1) * R(t, u - 2, v, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCy * R(t, u - 1, v, n + 1, p, PCx, PCy, PCz, RPC)
    elseif t > 0
        if t > 1
            val += (t - 1) * R(t - 2, u, v, n + 1, p, PCx, PCy, PCz, RPC)
        end
        val += PCx * R(t - 1, u, v, n + 1, p, PCx, PCy, PCz, RPC)
    end
    return val
end



