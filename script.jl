# make sure salloc is used to request resources, see /home/slurm-connect for details
# run with: salloc -p <partition> -N <nodes> -n <tasks> bash
# for ODIN: salloc -N 5 --ntasks-per-node=10 bash
# some global packages that need not be in project
using Distributed
using SlurmClusterManager # needs to be installed in the global environment
using DelimitedFiles
using Gnuplot
using CSV
Gnuplot.options.term = "sixelgd truecolor enhanced size 800,400 scroll font arial 10"

addprocs(SlurmManager(), 
         exeflags="--project=$(ENV["SLURM_SUBMIT_DIR"])")

println("salloc ENV directory: $(ENV["SLURM_SUBMIT_DIR"])")         
println("Hello world - $(gethostname()) on process $(myid())")
@everywhere println("Hello World from workers on $(gethostname()) with process $(myid())")

@everywhere using Revise
@everywhere includet("model.jl")

function run_calibration(beta, isoday, isoprop)
    @info "Total number of processors: $(nprocs())"
    @info "Running calibration from process $(myid())"
    #betavals = [0.015, 0.02, 0.025, 0.04, 0.05]
    cd = pmap(1:5000) do x
        # @info "Starting simulation $x on host $(gethostname()), id: $(myid())"
        # flush(stdout)
        total_infect = calibrate(beta; 
            iso_day = isoday, 
            iso_prop = isoprop
        ) # run the calibration
    end
    @info "Calibration completed on process $(myid())"
    @info "Total infected individuals: $(sum(cd))"
    @info "Average infected individuals: $(mean(cd))"
    @info "variance of infected individuals: $(var(cd))\n"
    return cd
end

function run_calibration_scenarios() 
    @info "Running calibration scenarios from process $(myid())"
    # run the calibration for different scenarios
    cfgs = [(;r=12, beta=0.26), (;r=15, beta=0.295), (;r=18, beta=0.329)]
    props = 0.2:0.1:0.8
    for cfg in cfgs
        data_collection = zeros(Int64, 5000, length(props) + 1) # to store the results
        data_collection[:, 1] = run_calibration(cfg.beta, -1, 0.0)
        _header = ["baseline"]
        for (i, isoprop) in enumerate(0.2:0.1:0.8)
            @info "column: $(i + 1): $(isoprop)"
            push!(_header, string(isoprop*100))
            data_collection[:, i + 1] = run_calibration(cfg.beta, 3, isoprop)
        end
        fname = "r$(cfg.r).csv"
        CSV.write(fname, CSV.Tables.table(data_collection); header=_header)
    end
end

function plot_figure1(sim_data) 
    # need to bootstrap the calibration results to look at R with isolation
    
    @gp "reset"
    # border
    @gp :- "set style line 101 lc rgb '#808080' lt 1 lw 1"
    @gp :- "set border 3 front ls 101"
    @gp :- "set tics nomirror out scale 0.75"

    # grid
    @gp :- "set style line 102 lc rgb '#808080' lt 0 lw 1"
    @gp :- "set grid back ls 102"

    # box plots
    # @gp :- "set style data boxplot" :- # better to use with boxplot
    @gp :- "set style boxplot medianlinewidth 2.5"
    @gp :- "set style boxplot outliers" :-
    @gp :- "set boxwidth 0.5 absolute" :-
    @gp :- "set style fill solid 0.80 border lt -1" :-
    @gp :- "set pointsize 0.5" :-
    
    # labels and keys
    @gp :- "unset key" :-
    @gp :- "set xtics   ('Baseline' 1, '20%%' 2, '30%%' 3, '40%%' 4, '50%%' 5, '60%%' 6, '70%%' 7, '80%%' 8)" :-
    #@gp :- "set format x '%.0f%%'"  :- 
    for (i, x) in enumerate(eachcol(sim_data))
        #println(mean(x))
        bs1 = bootstrap(mean, x, BasicSampling(2000))
        bsdata = bs1.t1[1]
        @gp :- [i for _ = 1:2000] bsdata "with boxplot lt 51 lc variable" :-
    end
    @gp 
    
end

function run_sims(_beta,
                  _iso_day,
                  _iso_prop,
                  _vac_scn1,
                  _vac_scn2) 
    @info "Running simulations with params: beta:$(_beta), iso day:$(_iso_day), prop: $(_iso_prop), vac1: $_vac_scn1, vac2: $_vac_scn2"
    cd = pmap(1:50) do x
        ihouse, ifarm, icomm, a1= time_loop(;
            simid = x, 
            beta=_beta, # 0.324:1.8, 0.255 = 1.2
            iso_day = _iso_day,  # 0 means on-set of symptomatic
            iso_prop = _iso_prop,
            vac_scn1 = _vac_scn1,
            vac_scn2 = _vac_scn2,
            vac_cov = 0.80
        )
        itotal = ihouse + ifarm + icomm
        return itotal, a1
    end

    return cd
end

function get_incidence(sims) 
    m = reduce(hcat, [x[1] for x in sims])
    mean(m, dims=2)
end

function get_reffective(sims)
    _allsims = map(sims) do y
        sortedarr = sort(y[2], by = x -> x[2]) # sort by time
        count_inf = zeros(Int64, 365) # create time vectors 
        total_inf = zeros(Int64, 365)
        for (inf, time) in sortedarr
            total_inf[time] += inf
            count_inf[time] += 1
        end
        withnans = total_inf ./ count_inf
    end
    allsims = reduce(hcat, _allsims)
    map(eachrow(allsims)) do x
        mean(filter(!isnan, x))
    end
end

function run_incidence_scenarios()
    # scenarios to run:
    cfgs = [(;r=12, beta=0.26), (;r=15, beta=0.295), (;r=18, beta=0.329)]

    # scenario combinations
    iso_props = 0.5:0.1:0.8
    vaxscen = (FARMONLY, FARM_AND_HH)
    vaxscen2 = (A1, A2, A3)   

    for cfg in cfgs
        @info "running simulation: $cfg"
        res_inf = zeros(Float64, 365, 1+4+24) # 1 baseline, 4 isolation (50% to 80%), 32 vaccine scenarios
        _header=["baseline"]

        # baseline
        sim_data = run_sims(cfg.beta, -1, 0.0, NONE, A0)
        res_inf[:, 1] = get_incidence(sim_data)

        # run baseline with isolation, no vaccine 
        for (i, isoprop) in enumerate(iso_props)
            push!(_header, string(isoprop*100))
            sim_data = run_sims(cfg.beta, 3, isoprop, NONE, A0)
            res_inf[:, i+1]= get_incidence(sim_data)
        end

        for (i, x) in enumerate(Base.Iterators.product(iso_props, vaxscen, vaxscen2))
            _isoprop = x[1]
            _vaxscen = x[2]
            _vaxscen2 = x[3]
            push!(_header, "$(string(_isoprop*100))_$(string(_vaxscen))_$(string(_vaxscen2))")
            sim_data = run_sims(cfg.beta, 3, _isoprop, _vaxscen, _vaxscen2)
            res_inf[:, i+5]= get_incidence(sim_data)
        end

        fname = "r$(cfg.r)-incidence.csv"
        CSV.write(fname, CSV.Tables.table(res_inf); header=_header)
    end
end

# @eval SlurmClusterManager function Distributed.launch(manager::SlurmManager, params::Dict, instances_arr::Array, c::Condition)
#     try
#         warn_if_unexpected_params(params)

#         exehome = params[:dir]
#         exename = params[:exename]
#         exeflags = params[:exeflags]

#         _srun_cmd_without_env = `srun -D $exehome $exename --threads 1 $exeflags --worker`

#         @static if Base.VERSION >= v"1.6.0"
#           # Pass the key-value pairs from `params[:env]` to the `srun` command:
#           env2 = _new_environment_additions(Dict{String,String}(params[:env]))
#           srun_cmd_with_env = addenv(_srun_cmd_without_env, env2)
#         else
#           # See discussion above for why we don't support this functionality on Julia 1.5 and earlier.
#           if haskey(params, :env)
#             @warn "SlurmClusterManager.jl does not support params[:env] on Julia 1.5 and earlier" Base.VERSION
#           end
#           srun_cmd_with_env = _srun_cmd_without_env
#         end

#         # Pass cookie as stdin to srun; srun forwards stdin to process
#         # This way the cookie won't be visible in ps, top, etc on the compute node
#         @debug "srun command: $srun_cmd_with_env"
#         manager.srun_proc = open(srun_cmd_with_env, write=true, read=true)
#         write(manager.srun_proc, cluster_cookie())
#         write(manager.srun_proc, "\n")

#         t = @async for i in 1:manager.ntasks
#           @debug "connecting to worker $i out of $(manager.ntasks)"

#           line = readline(manager.srun_proc)
#           @debug "Worker $i output: $line"
#           m = match(r"(\d*)#(.*)", line)
#           m === nothing && error("could not parse $line")
#           @debug "Worker $i matched: $m"
#           config = WorkerConfig()
#           config.port = parse(Int, m[1])
#           config.host = strip(m[2])

#           @debug "Worker $i ready on host $(config.host), port $(config.port)"

#           push!(instances_arr, config)
#           notify(c)
#         end

#         # workers must be launched before timeout otherwise interrupt
#         status = timedwait(() -> istaskdone(t), manager.launch_timeout)
#         if status !== :ok
#           @async Base.throwto(t, ErrorException("launch_timeout exceeded"))
#         end
#         wait(t)

#         # redirect output
#         @async while !eof(manager.srun_proc)
#           line = readline(manager.srun_proc)
#           println(line)
#         end

#         # wait to make sure that srun_proc exits before main program to avoid slurm complaining
#         # avoids "Job step aborted: Waiting up to 32 seconds for job step to finish" message
#         finalizer(manager) do manager
#           wait(manager.srun_proc)
#           # need to sleep briefly here to make sure that srun exit is recorded by slurm daemons
#           # TODO find a way to wait on the condition directly instead of just sleeping
#           sleep(manager.srun_post_exit_sleep)
#         end

#     catch ex
#         @error "Error launching Slurm job" exception=ex
#         rethrow()
#     end
# end
