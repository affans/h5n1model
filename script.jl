#!/usr/bin/env julia
#SBATCH --job-name=abm-h5n1          # Job name
#SBATCH --nodes=7                    # Request n nodes
##SBATCH --ntasks-per-node=90         # 10 tasks per node
#SBATCH --ntasks=90                  # Total number of tasks
#SBATCH --output=soutput.log         # Log output file output_%j.log

# how to run this script using either sbatch or salloc:
# SALLOC: For interactive use, you can use salloc to request resources and then run this script.
# on ODIN (in the correct working directory)
# > salloc  --nodes 5 --ntasks=50
# > julia --project=. 
# > using Revise; includet("script.jl")  
# and then run your functions ()

# sbatch can be used to submit the script to the SLURM scheduler (and run in the background)
# on ODIN, (in the correct working directory)
# > sbatch script.jl 
# note that this launches `julia script.jl` so the environment loaded in julia will be the global 
# use !isinteractive() to check if the script is run in a non-interactive mode and switch directory

using Pkg
if !isinteractive()
    # the script is run in a non-interactive mode, e.g. via julia scripts.jl or sbatch script.jl
    # if sbatch script.jl is used, the environment will be the global one (since #!/usr/bin/env julia on top of line)
    # and so we need to switch to the project directory 
    @info("starting script, non interactive use - sbatch")
    @info("@__DIR__: $(@__DIR__)") # this is where the script is copied and executed from
    @info("current working directory: $(pwd())") # this is the working directory (slurm sets this but also provides )
    
    # Pkg.activate(@__DIR__) # does not work since sbatch is "run" from project: /var/lib/slurm/slurmd/job00577/Project.toml
    # Pkg.activate("$(ENV["SLURM_SUBMIT_DIR"])") # works - activate the project in the SLURM_SUBMIT_DIR
    Pkg.activate(pwd()) # activate the project in the SLURM_SUBMIT_DIR
    @info("activated project: $(Pkg.project().path)")
end

using SlurmClusterManager # global env
using Distributed
using DelimitedFiles
using CSV
using Bootstrap
using DataFrames
using Gnuplot  # global env
Gnuplot.options.term = "sixelgd truecolor enhanced size 1000,400 scroll font arial 10"

# slurm needs to have the correct ENV VARs set. 
addprocs(SlurmManager(),
    exeflags="--project=$(ENV["SLURM_SUBMIT_DIR"])")

@info("Hello world - $(gethostname()) on process $(myid())")
@everywhere @info("Hello World from workers on $(gethostname()) with process $(myid())")
@everywhere using Revise
@everywhere includet("$(pwd())/model.jl")

function run_calibration(nsims, beta, isoday, isoprop)
    # for calibration, set isoday=-1 and isoprop=0.0 
    # we give these options to see what happens to R if we isolate at 
    # a certain day with proportion (assuming baseline beta) 
    # that's what the run_calibration_scenarios() does 
    cd = pmap(1:nsims) do x
        # @info "Starting simulation $x on host $(gethostname()), id: $(myid())"
        # flush(stdout)
        total_infect, tih, tif, tic  = calibrate(beta;
            iso_day=isoday,
            iso_prop=isoprop
        ) # run the calibration
        return total_infect, tih, tif, tic
    end
    
    total_infected = [x[1] for x in cd] # get the total number of infected
    @info "Average: $(round(mean(total_infected), digits=2)), Variance: $(round(var(total_infected), digits=3))"
    
    # print proportions of infected in households, farms, communities
    for c in (2, 3, 4)  # 2 is households, 3 is farms, 4 is communities
        _ti = [x[c] for x in cd] # get the total number of infected in households, farms, communities
        prop = mean(filter(!isnan, _ti ./ total_infected))
        @info "Mean proportion, id $c: $(prop)"
    end
    return total_infected
end

# runs a set of simulations with the given parameters
function run_sims(nsims, 
                _beta,    
                _init_inf,
                _iso_day,
                _iso_prop,
                _vac_scn1,
                _vac_scn2,
                _vac_time)
    @info """starting $nsims simulations with params: 
              beta => $(_beta), 
              iso day => $(_iso_day), 
              iso prop => $(_iso_prop), 
              vac scenario => $_vac_scn1, 
              vac type => $_vac_scn2, 
              vac start => $_vac_time,
              init infections => $(_init_inf)
          """
    cd = pmap(1:nsims) do x
        ihouse, ifarm, icomm = time_loop(;
            simid=x,
            beta = _beta, # 0.324:1.8, 0.255 = 1.2
            init_inf = _init_inf,
            iso_day = _iso_day,  # 0 means on-set of symptomatic
            iso_prop = _iso_prop,
            vac_scn1 = _vac_scn1,
            vac_scn2 = _vac_scn2,
            vac_cov = 0.80, 
            vac_time = _vac_time
        )
        itotal = ihouse + ifarm + icomm
        return itotal
    end
    return cd
end

function run_calibration_scenarios(nsims)
    # calculates the reproduction number for different scenarios
    # it loops over all R and Iso proportion values

    @info "Running calibration scenarios from process $(myid())"
    # run the calibration for different scenarios
    cfgs = [(; r=12, beta=0.265), (; r=15, beta=0.3), (; r=18, beta=0.325)]
    props = 0.2:0.1:0.8

    fld_name = "./output/secondary_infections"
    if !isdir(fld_name)
        @info "Creating folder: $fld_name"
        mkpath(fld_name)
    end
   
    # iso_day = 2
    for cfg in cfgs
        @info "running calibration configuration: r=$(cfg.r), beta=$(cfg.beta)"
        datares = zeros(Int64, nsims, length(props) + 1) # to store the results
        datares[:, 1] = run_calibration(nsims, cfg.beta, -1, 0.0)
        _header = ["baseline"]
        for (i, isoprop) in enumerate(0.2:0.1:0.8)
            @info "   R = $(cfg.r), iso prop: $(isoprop), column: $(i + 1)"
            push!(_header, string(isoprop * 100))
            datares[:, i+1] = run_calibration(nsims, cfg.beta, 2, isoprop)
        end
        fname = "$fld_name/r$(cfg.r)_isoday2_secondaryinfections.csv"
        CSV.write(fname, CSV.Tables.table(datares); header=_header)
    end

    # iso_day = 3
    for cfg in cfgs
        @info "running calibration configuration: r=$(cfg.r), beta=$(cfg.beta) for day 3"
        datares = zeros(Int64, nsims, length(props) + 1) # to store the results
        datares[:, 1] = run_calibration(nsims, cfg.beta, -1, 0.0)
        _header = ["baseline"]
        for (i, isoprop) in enumerate(0.2:0.1:0.8)
            @info "   R = $(cfg.r), iso prop: $(isoprop), column: $(i + 1)"
            push!(_header, string(isoprop * 100))
            datares[:, i+1] = run_calibration(nsims, cfg.beta, 3, isoprop)
        end
        fname = "$fld_name/r$(cfg.r)_isoday3_secondaryinfections.csv"
        CSV.write(fname, CSV.Tables.table(datares); header=_header)
    end
end

function get_hcat_incidence(sims)
    # concatenates the simulation results in a matrix
    m = reduce(hcat, sims)
end

function get_vcat_incidence(sims)
    # concatenates the simulation results in a single vector
    m = reduce(vcat, sims)
end

function run_incidence_scenarios(nsims)
    # scenario combinations
    beta_values = ((12, 0.265), (15, 0.3), (18, 0.325))
    iso_props = 0.5:0.1:0.8
    vaxscen = (FARMONLY, FARM_AND_HH)
    vaxtypes = (A1, A2, A3)
    #init_infections = (1, 2)
    init_infections = (1, )
    iso_day = (2, 3) # isolation day, 2 means on-set of symptomatic
    vac_time = (1, 42) # vaccination time, 42 days after the start of the simulation

    # create a vector of configurations
    file_configs = vec(collect(Base.Iterators.product(beta_values, vaxtypes, init_infections, iso_day, vac_time)))
   
    # loop through the configurations: each cfg is a tuple of (beta_values, vaxtypes, init_infections)
    # use the configuration to create a filename
    # columns per file: baseline, isolation 50 - 80 (4 columns), farm only vaccine (50 to 80) (4 columns), farm_HH vaccine (50 to 80) (4 columns)
    # total of 1 + 4 + 4 + 4 = 13 columns
    # the way this for loop works is that it repeats some configurations
    ctr = 0
    for cfg in file_configs   
        rval = cfg[1][1] 
        beta = cfg[1][2]
        vaxtype = cfg[2]
        init_inf = cfg[3] # initial infections
        iso_day = cfg[4] # time of isolation
        vt = cfg[5] # vaccination time
        
        # folder/filename setup
        _fprefix = vt == 1 ? "preemptive" : "reactive"
        fld_name = "./output/$(_fprefix)_isoday$(iso_day)"
        if !isdir(fld_name)
            @info "Creating folder: $fld_name"
            mkpath(fld_name)
        end
        fname = "r$(rval)_$(vaxtype)_i$(init_inf)_vt$(vt)_isoday$(iso_day)_incidence.csv"
        @info("Generating file: $fname")
        ctr += 1
        
        # create data structure for the results
        datares = zeros(Int64, 365 * nsims, 13)

        # run baseline scenario
        @info("   ctr: $ctr col 1: (running) baseline scenario")
        sim_data = run_sims(nsims, beta, init_inf, -1, 0.0, NONE, A0, 0) 
        datares[:, 1] = get_vcat_incidence(sim_data)
        
        # run isolation scenarios
        for (i, _isoprop) in enumerate(iso_props)
            ctr += 1
            @info("   ctr: $ctr col $(i + 1): (running) isoprop: $(_isoprop), no vax")
            sim_data = run_sims(nsims, beta, init_inf, iso_day, _isoprop, NONE, A0, 0)
            datares[:, i+1] = get_vcat_incidence(sim_data)            
        end

        # run vaccination scenarios
        for (i, x) in enumerate(Base.Iterators.product(iso_props, vaxscen))
            ctr += 1
            @info("   ctr: $ctr col $(i + 5): (running) iso, vax")
            _isoprop = x[1]
            _vaxscen = x[2]
            sim_data = run_sims(nsims, beta, init_inf, iso_day, _isoprop, _vaxscen, vaxtype, vt)
            datares[:, i+5] = get_vcat_incidence(sim_data)
        end
        CSV.write("$fld_name/$fname", CSV.Tables.table(datares); header=["baseline", 
            "isolation_50", "isolation_60", "isolation_70", "isolation_80",
            "farmonly_50", "farmonly_60", "farmonly_70", "farmonly_80",
            "farmhh_50", "farmhh_60", "farmhh_70", "farmhh_80"])
    end
end

function plot_r_figure() 
    # read the data in output folder and plot the R values using the secondary_infections files 
    # WIP 
    isoday = 3
    r12data = CSV.read("./output/secondary_infections/r12_isoday$(isoday)_secondaryinfections.csv", DataFrame)
    r15data = CSV.read("./output/secondary_infections/r15_isoday$(isoday)_secondaryinfections.csv", DataFrame)
    r18data = CSV.read("./output/secondary_infections/r18_isoday$(isoday)_secondaryinfections.csv", DataFrame)

    @gp "reset"
    @gp :- "set multiplot layout 1, 3" :- # 1 row, 3 columns
    @gp :- "set style line 101 lc rgb '#808080' lt 1 lw 1" :- # border
    @gp :- "set border 3 front ls 101" :-
    @gp :- "set tics nomirror out scale 0.75" :-
    @gp :- "set style line 102 lc rgb '#808080' lt 0 lw 1" :- # grid
    @gp :- "set grid back ls 102" :-
    # @gp :- "set style data boxplot" :- # better to use with boxplot
    @gp :- "set style boxplot medianlinewidth 2.5" :- # box plots
    @gp :- "set style boxplot nooutliers" :-
    #@gp :- "set boxwidth 0.5 absolute" :- # set by the column spec passed in 
    @gp :- "set style fill solid 0.80 border lt -1" :-
    @gp :- "set pointsize 0.5" :-

    @gp :- "unset key" :-
    @gp :- "set xtics   ('BL' 1, '20%%' 2, '30%%' 3, '40%%' 4, '50%%' 5, '60%%' 6, '70%%' 7, '80%%' 8)" :-
    for (pn, simdata) in enumerate([r12data, r15data, r18data])
        #@gp :- "set ylabel 'R'" :-
        #@gp :- "set xlabel 'Isolation Proportion'" :-
        #@gp :- "set format x '%.0f%%'"  :- 
        for colnum in 1:8
            
            bs1 = bootstrap(mean, simdata[:, colnum], BasicSampling(1000)).t1[1]
            xval = repeat([colnum], length(bs1)) # x value, in array form
            #@gp :- pn xval bs1  "with boxplot lt 51 lc variable" :-
            @gp :- pn xval bs1  "with boxplot lt 51 lc 'black'" :-
        end
    end
    @gp
end

function plot_incidence1()
    vt = "reactive"
    vt1 = 42
    iso = 3
    fldname = "./output/$(vt)_isoday$(iso)"
    r12data = CSV.read("$fldname/r12_A1_i1_vt$(vt1)_isoday$(iso)_incidence.csv", DataFrame)
    r15data = CSV.read("$fldname/r15_A1_i1_vt$(vt1)_isoday$(iso)_incidence.csv", DataFrame)
    r18data = CSV.read("$fldname/r18_A1_i1_vt$(vt1)_isoday$(iso)_incidence.csv", DataFrame)


    @gp "reset" 
    @gp :- "set multiplot layout 1, 4" :- # 1 row, 3 columns
    @gp :- "set style line 101 lc rgb '#808080' lt 1 lw 1" :- # border
    @gp :- "set border 3 front ls 101" :-
    @gp :- "set tics nomirror out scale 0.75" :-
    @gp :- "set style line 102 lc rgb '#808080' lt 0 lw 1" :- # grid
    @gp :- "set grid back ls 102" :-
    @gp :- "set termoption dashed" :- 
    #@gp :- "set xtics ('1' 1, '30' 2, '30%%' 3, '40%%' 4, '50%%' 5, '60%%' 6, '70%%' 7, '80%%' 8)" :-
    @gp :- "set xtics 30" :- 
    # plot baseline for all r values 
    _rdata = [r12data, r15data, r18data]
    _pdata = ["R 1.2", "R 1.5", "R 1.8"]
    _ltype = ["lc rgb '#969696'", "lc rgb '#636363'", "lc rgb '#252525'"]
    for (rd, pd, lt) in zip(_rdata, _pdata, _ltype)
        _data = rd[:, :baseline]
        h = reduce(hcat, Iterators.partition(_data, 365))
        daily_mean = cumsum(vec(mean(h, dims=2))) # mean over the columns
        @gp :- 1 1:150 daily_mean[1:150] "with line $lt title '$pd'" :-
    end

    
    figcolms = [:isolation_50, :isolation_60, :isolation_70, :isolation_80] 
    figcolrs = ["lc rgb '#377eb8'", "lc rgb '#4daf4a'", "lc rgb '#984ea3'", "lc rgb '#ff7f00'"]
    fittitle = ["50%", "60%", "70%", "80%"]
    for (i, rd) in enumerate([r12data, r15data, r18data])
        for (colm, clr, ft) in zip(figcolms, figcolrs, fittitle)
            _data = rd[:, colm]
            h = reduce(hcat, Iterators.partition(_data, 365))
            daily_mean = cumsum(vec(mean(h, dims=2))) # mean over the columns
            @gp :- (i + 1) 1:180 daily_mean[1:180] "with line $clr title '$ft'" :-
        end
    end
    
    @gp 
    
    # manual bootstrap
    #mbt = sample(1:1000, 1000, replace=true)
    #mbt_mean = vec(mean(h[:, mbt], dims=2))
    
    #@gp :- 1:365 daily_high "with lines lw 1 lc 'red' title 'R12 95% CI'" :-
    #@gp :- 1:365 daily_low "with lines lw 1 lc 'blue' title 'R12 5% CI'"
    return r12data
end


# run the functions if launched through sbatch
if !isinteractive()
    # run the calibration scenarios
    @info ("Running Calibration Scenarios")
    run_calibration_scenarios(5000)

    @info "Running Incidence Scenarios"
    run_incidence_scenarios(1000)
    exit(0)
end