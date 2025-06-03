using StatsBase, Random, Distributions, DelimitedFiles
using StaticArrays

# make sure the order doesn't change
@enum INFTYPE UNDEF=0 SUS=1 CARC=2 CARW=3 CARY=4 IMD=5 REC=6

# define an agent and all agent properties
Base.@kwdef mutable struct Human
    idx::Int64 = 0
    age::Int64 = 0 # age in weeks
    agegroup::Int8 = 0 # store the age group for easier access to other properties
    hid::Int64 = 0 # household id 
    inf::INFTYPE = SUS # infection status
    swap::INFTYPE = UNDEF # swap state
    tis::Int16 = 0 # time in state
    st::Int16 = typemax(Int16) # max time in state 
    vac::Int16 = 0 # dose count
    eff::Float64 = 0.0 # vaccine efficacy
    exp::Int64 = 0 # waning time of vaccine
end
Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

## system parameters
Base.@kwdef mutable struct ModelParameters
    beta::Vector{Float64} = [0.5, 0.5, 0.5] # for Carriage 
end

# constant variables
const POPSIZE = 153525
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters


include("helpers.jl") # include the helper functions
include("debug.jl")
init_state() = init_state(ModelParameters())
function init_state(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    # ip.popsize == 0 && error("no population size given")
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
    # resize the human array to change population size
    resize!(humans, POPSIZE)
    @debug "Initialization state complete..."
    return
end

init_agents() = init_agents(0)
function init_agents(simid) 
    _POPDISTR = @SVector [1530, 1545, 1565, 1670, 1700, 1820, 1755, 1790, 1820, 1935, 1875, 1925, 
                1925, 1910, 1890, 1800, 1895, 1900, 1820, 1855, 1995, 2135, 2040, 2060, 
                2030, 2120, 1985, 2020, 2125, 2140, 2120, 2230, 2240, 2065, 2105, 2190, 
                2160, 2195, 2225, 2210, 2190, 2025, 2005, 1965, 1845, 1870, 1920, 1840, 
                1835, 1855, 1880, 1895, 1870, 1915, 1645, 1820, 1885, 1970, 1875, 1755, 
                1850, 1945, 1845, 1850, 1790, 1720, 1810, 1735, 1595, 1520, 1460, 1515, 
                1430, 1500, 1410, 1135, 1060, 1050, 925, 860, 820, 735, 730, 755, 570, 530, 
                445, 395, 440, 370, 405, 315, 185, 170, 140, 145, 75, 40, 40, 35, 45]
    
    @info "Initialzing new zero population"
    # there is no shuffling going on here, so at index 32250 adults start 
    # this is hardcoded in the model, so don't change it 
    _age_vec = inverse_rle(0:100, _POPDISTR)
    for i in 1:POPSIZE
        humans[i] = Human()
        x = humans[i]
        x.idx = i
        x.age = _age_vec[i]
        x.agegroup = 0 # age group is not set yet, update if needed
    end 
    return 
end

function create_houseids() 
    # encoding scheme -- the house id encodes the household size in the last digit
    # houseid = baseid * 10 + size 
    #   base_id is a unique identifier for the house (e.g., sequential number).
    #   size_code is a number from 1 to 5 indicating household size.
    # to find household size from a house_id: size = house_id % 10
    # to get the unique baseid: base_id = house_id // 10

    # the houses here can accomodate 142285 individuals, we have a fixed size of POPSIZE = 153525, 
    # so we have less houses than popsize - but the last house size is 5+ so throw the rest in there
    household_sizes = @SVector [13180, 16335, 7265, 7660, 8800] #sum 53240
    sizes = inverse_rle(1:5, household_sizes) # size of each house
    bases = 1:sum(household_sizes) # unique ID of each house 
    houseids = bases .* 10 .+ sizes
    return houseids
end

function init_households() 
    house_ids = create_houseids()

    # we have these IDs separated because the initialized population is not randomized
    children_ids = shuffle!(collect(1:32250)) # children are 0-17 years old, adults are 18+ years old
    adult_ids = shuffle!(collect(32250:POPSIZE))
    children_idx = 1
    adult_idx = 1
    # configuration options for each house of different sizes. Each tuple is (num_adults, num_children)
    comp_options = @SVector [[(1, 0)], 
                    [(1, 1), (2, 0)],
                    [(2, 1), (1, 2), (3, 0)],
                    [(2, 2), (3, 1), (4, 0)],
                    [(2, 3), (3, 2), (4, 1), (5, 0)]]

    for house_id in house_ids
        size_code = house_id % 10 # get the size of the house using our encoding scheme
        
        # see which configuration works best. for each config, we need to check if we have enough adults and children 
        all_options = comp_options[size_code] # get the options for this size codecomp_options
        viable = false 
        local o
        while !viable
            o = rand(all_options)
            if length(adult_ids) >= o[1] && length(children_ids) >= o[2]
                viable = true
            end
        end
        num_adults, num_children = o # splat the tuple
        selected_adults = poplast!(adult_ids, num_adults) 
        selected_children = poplast!(children_ids, num_children)
     
        for h in selected_children
            humans[h].hid = house_id
        end
        for h in selected_adults
            humans[h].hid = house_id
        end
    end

    # Assign adults without a house to a 5+ house randomly. 
    leftover_adults = filter(x -> x.hid == 0 && x.age >= 18, humans)
    large_house_ids = filter(hid -> hid % 10 >= 5, house_ids) # get the house IDs with size 5+
    for h in leftover_adults
        h.hid = rand(large_house_ids) # assign to the last large house
    end
        
    # run error checks 
    #household_errorchecks()
    return 
end 