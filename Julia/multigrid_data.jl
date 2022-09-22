
const max_level = 16
const base_grid = 3
const pre_relax = 4
const post_relax = 4
const n_vcycle = 20

xmin = 0.
xmax = 1.
ρ_c = 1.28e-3
r_s = 8.


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

    return mg_data_type(nx,dx,dx2,r_dx,r_dx2,x,r_x,q,src,res,err)
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

function initialize_finest_grid(u1)

    n = u1.nx
    rs2 = r_s^2

    for i = 1:n
        u1.q[i] = 0.
    end

    for i = 1:n-1
        x = u1.x[i]
        if x > 0.5
            u1.src[i] = 0.
        else
            u1.src[i] = 4π*ρ_c*rs2*(1-(x/(1-x))^2)/(1-x)^4
        end
    end

    return u1
end
