
include("multigrid_data.jl")
include("multigrid_functions.jl")


mg_data = grid()
mg_data[1] = initialize_finest_grid(mg_data[1])

time_taken = []

for vcycle = 1:n_vcycle
    elapsed_time = @elapsed begin

    for level = 1:max_level-1

        for i = 1:pre_relax
            mg_data[level] = relaxation(mg_data[level])
        end

        mg_data[level] = residual(mg_data[level])
        mg_data[level+1].src = restriction(mg_data[level+1].nx,mg_data[level].nx,mg_data[level].res)
        mg_data[level+1].q .= 0

    end


    ###

    level = max_level
    mg_data[level] = relaxation_maxlevel(mg_data[level])

    ###

    for level = max_level-1:-1:1
        mg_data[level].err = prolongation(mg_data[level+1].nx,mg_data[level].nx,mg_data[level+1].q)
        mg_data[level].q = mg_data[level].q .+ mg_data[level].err

        for i = 1:post_relax
            mg_data[level] = relaxation(mg_data[level])
        end
    end

    level = 1

    #L_1, L_inf = residual_norm(mg_data[level])
    #println("$vcycle, $L_1, $L_inf")
    rms_res, rms_exact = rms_error(mg_data[level])
    println("$vcycle, $rms_res, $rms_exact")


    end
    @show elapsed_time
    push!(time_taken,elapsed_time)

end


println("Time average: ",sum(time_taken[2:n_vcycle])/(n_vcycle-1))
