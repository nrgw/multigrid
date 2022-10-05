
include("multigrid_data.jl")


###

function relaxation(u :: mg_data_type)
    nx = u.nx

    i = 1
    u.q[i] = u.q[i+1] - u.dx2/6*u.src[i]
    for j = 3:-1:2
        for i = j:2:nx-1
            u.q[i] = 0.5*( u.q[i+1]+u.q[i-1] 
                        + u.dx*u.r_x[i]*(u.q[i+1]-u.q[i-1]) 
                        - u.dx2*u.src[i] )
        end
    end

    return u
end    


###

function residual(u :: mg_data_type)
    nx = u.nx
    u.res[nx] = 0

    i = 1
    u.res[i] = u.src[i] - 6*(u.q[i+1]-u.q[i])*u.r_dx2
    for i = 2:nx-1
        u.res[i] = u.src[i] - (u.q[i+1]-2*u.q[i]+u.q[i-1])*u.r_dx2 - u.r_x[i]*(u.q[i+1]-u.q[i-1])*u.r_dx
    end

    return u
end


###

function residual_norm(u :: mg_data_type)
    nx = u.nx

    i = 1
    temp = abs(u.src[i]-6*(u.q[i+1]-u.q[i])*u.r_dx2)
    L_1 = temp
    L_inf = temp

    for i = 2:nx-1
        temp = abs( u.src[i]-(u.q[i+1]-2*u.q[i]+u.q[i-1])*u.r_dx2
                    - u.r_x[i]*(u.q[i+1]-u.q[i-1])*u.r_dx )
        L_1 = L_1 + temp
        L_inf = max(L_inf,temp)
    end
    L_1 = L_1/nx

    return L_1, L_inf
end


###

function relaxation_maxlevel(u :: mg_data_type)
    nx = u.nx

    u.q[1] = u.q[3] - u.dx2/(1+u.dx*u.r_x[2])*(u.src[1]/3+u.src[2])
    u.q[2] = u.q[1] + u.dx2*u.src[1]/6

    return u
end


###

function restriction(n_coarse::Int64,n_fine::Int64,q_fine::Vector{Float64})
    #q_coarse[n_coarse] = 0
    q_coarse = zeros(n_coarse)

    q_coarse[1] = 0.5*(q_fine[1]+q_fine[2])     ## reflective bc
    for i = 2:(n_coarse-1)
        i_fine = 2*i - 1
        q_coarse[i] = 0.25*(q_fine[i_fine-1] + q_fine[i_fine+1] + 2*q_fine[i_fine])
    end

    return q_coarse
end


###

function prolongation(n_coarse::Int64,n_fine::Int64,q_coarse::Vector{Float64})
    q_fine = zeros(n_fine)

    for i = 1:n_coarse
        i_fine = 2*i - 1
        q_fine[i_fine] = q_coarse[i]
    end

    for i = 2:2:n_fine
        q_fine[i] = 0.5*(q_fine[i-1]+q_fine[i+1])
    end

    return q_fine
end



###

function rms_error(u :: mg_data_type)
    rms_exact = 0.

    i = 1
    temp1 = u.src[i]-6*(u.q[i+1]-u.q[i])*u.r_dx2 
    rms_res = temp1^2 

    for i = 2:u.nx-1
        temp1 = (u.src[i]-(u.q[i+1]-2*u.q[i]+u.q[i-1])*u.r_dx2-u.r_x[i]*(u.q[i+1]-u.q[i-1])*u.r_dx)
        rms_res = rms_res + temp1^2
    end 

    for i = 1:u.nx 
        rms_exact = rms_exact + (u.q[i]-u.exact[i])^2
    end 

    temp2 = 1. /u.nx
    rms_res = sqrt(rms_res*temp2)
    rms_exact = sqrt(rms_exact*temp2)
    return rms_res, rms_exact 
end