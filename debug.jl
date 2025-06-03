# Summary and Debug Functions 
# remove packages if not needed for final version
using UnicodePlots 

function dbg_household() 
    house_ids = create_houseids()
    # all houseIDs should be used (i.e., every house should have at least one person)
    house_ids_assigned = Set(humans[i].hid for i in 1:POPSIZE)
    for hid in house_ids
        if !(hid in house_ids_assigned)
            @error "House ID $(hid) is not assigned to any human!"
        end
    end

    # no one should be without a house 
    leftover_children = filter(h -> h.hid == 0, humans)
    length(leftover_children) > 0 && @error "There are still $(length(leftover_children)) children left without a house!"
    @info "Household error checks passed."
    
    # Average distribution of adults and children in households
    count_adults = [0, 0, 0, 0, 0]
    count_children = [0, 0, 0, 0, 0]
    houses_sizes = [13180, 16335, 7265, 7660, 8800] #sum 53240
    for h in humans 
        lives_in = h.hid % 10 
        if h.age < 18 
            count_children[lives_in] += 1
        else
            count_adults[lives_in] += 1
        end
    end
    average_adults = round.(count_adults ./ houses_sizes; digits=4)
    average_children = round.(count_children ./ houses_sizes; digits=4)
    @info "Average adults per house size: $(average_adults)"
    @info "Average children per house size: $(average_children)"   
    
    dict = Dict(k => Vector{Int64}() for k in house_ids)
    for h in humans 
        if h.age >= 18
            push!(dict[h.hid], h.age)
        end
    end
    result = map(mean, values(dict))
    histogram(result, nbins=100, vertical=true)
    #return dict
    #for hid in house_ids
    #     housecomp = findall(x -> x.hid == hid && x.age >= 18, humans)
    #     houseage = [humans[x].age for x in housecomp]
    # end
end
