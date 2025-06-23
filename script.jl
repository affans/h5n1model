# make sure salloc is used, see /home/slurm-connect for details
# run with: salloc -p <partition> -N <nodes> -n <tasks> bash
# for ODIN: salloc -N 5 --ntasks-per-node=10 bash
# some global packages that need not be in project

using Distributed
using SlurmClusterManager # needs to be installed in the global environment
using DelimitedFiles

addprocs(SlurmManager(), 
         exeflags="--project=$(ENV["SLURM_SUBMIT_DIR"])")

println("Hello world - $(gethostname()) on process $(myid())")
@everywhere println("Hello World from workers on $(gethostname()) with process $(myid())")

@everywhere using Revise
@everywhere includet("model.jl")

function run_calibration()
    @info "Total number of processors: $(nprocs())"
    @info "Running calibration from process $(myid())"
    betavals = [0.015, 0.02, 0.025, 0.04, 0.05]
    cd = pmap(1:500) do x
        # @info "Starting simulation $x on host $(gethostname()), id: $(myid())"
        # flush(stdout)
        total_infect = calibrate(0.02)
    end
    @info "Calibration completed on process $(myid())"
    @info "Total infected individuals: $(sum(cd))"
    @info "Average infected individuals: $(mean(cd))"
    return cd
end

function run_sims() 
    @info "Total number of processors: $(nprocs())"
    @info "Running sims from process $(myid())"
    cd = pmap(1:500) do x
        ihouse, ifarm, icomm = time_loop(beta=0.02)
        itotal = ihouse + ifarm + icomm
        return itotal
    end
    _mean = reduce(hcat, cd)
    return mean(_mean, dims=2)
end
