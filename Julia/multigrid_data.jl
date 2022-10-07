
const max_level = 16
const base_grid = 3
const pre_relax = 4
const post_relax = 4
const n_vcycle = 20

const xmin = 0.
const xmax = 1.
const ρ_c = 1.28e-3
const N = 1.
const K = 1.e2
const alpha2 = (N+1)*K/(4π)
const alpha = sqrt(alpha2)
const r_s = π*alpha



###

mutable struct mg_data_type
    nx :: Int64 
    dx :: Float64 
    dx2 :: Float64 
    r_dx :: Float64 
    r_dx2 :: Float64 
    x :: Vector{Float64}
    r_x :: Vector{Float64}
    q :: Vector{Float64}
    src :: Vector{Float64}
    res :: Vector{Float64}
    err :: Vector{Float64}
    exact :: Vector{Float64}
end


###

function allocate_data(level :: Int64)

    n = (base_grid-1)*2^(max_level-level)+1

    nx = n
    dx = (xmax-xmin)/(n-1)
    dx2 = dx^2
    r_dx = 1/dx
    r_dx2 = r_dx^2

    x = zeros(n)
    for i = 1:n
        x[i] = xmin + dx*(i-1)
    end

    r_x = zeros(n)
    r_x[1] = 0.
    for i = 2:n
        r_x[i] = 1/x[i]
    end

    q = zeros(n)
    src = zeros(n)
    res = zeros(n)
    err = zeros(n)
    exact = zeros(n)

    return mg_data_type(nx,dx,dx2,r_dx,r_dx2,x,r_x,q,src,res,err,exact)
end


###

function grid()

    u = []
    for l = 1:max_level
        push!(u,allocate_data(l))
    end

    return u
end


###

function initialize_finest_grid(u1::mg_data_type)

    n = u1.nx
    rs2 = r_s^2

    for i = 1:n
        u1.q[i] = 0.
        u1.exact[i] = 0.
    end

    i = 1
    u1.src[i] = 4π*ρ_c*rs2
    u1.exact[i] = -8π*ρ_c*alpha2 
    for i = 2:n-1
        x = u1.x[i]
        r_x = u1.r_x[i]
        temp = x/(1-x)
        if x > 0.5 
            u1.src[i] = 0. 
            u1.exact[i] = -4π*alpha2*ρ_c/temp 
        else 
            u1.src[i] = 4π*ρ_c*rs2*(sin(π*temp)/(π*temp))/(1-x)^4
            u1.exact[i] = -4π*ρ_c*alpha2*(1+sin(pi*temp)/(pi*temp))
        end
    end

    return u1
end
