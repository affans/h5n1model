# Summary and Debug Functions 
# remove packages if not needed for final version
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
    for h in humans
        if h.hid == 0
            @error "Human $(h.idx) has no house assigned!"
        end
    end
    @info "Household error checks passed. All individuals are assigned to a house"

    # print how many houses of each size 
    dict_a = Dict(k => Vector{Int64}() for k in house_ids)
    shid = [h.hid for h in humans]
    count_freq = countmap(values(countmap(shid)))
    # Convert to vector and sort by the number of elements (the key)
    sorted_counts = sort(collect(count_freq), by=x -> x[1])

    # Print
    overall_avg_size = 0 
    for (count, howmany) in sorted_counts
        overall_avg_size += count * howmany
        println("Houses of size $count: $howmany")
    end
    println("Overall average house size: $(round(overall_avg_size / length(house_ids), digits=4))")

    
    
    # Average distribution of adults and children in households
    count_adults = [0, 0, 0, 0, 0]
    count_children = [0, 0, 0, 0, 0]
    count_all = [0, 0, 0, 0, 0]
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
    #histogram(result, nbins=100, vertical=true)
    #return dict
    #for hid in house_ids
    #     housecomp = findall(x -> x.hid == hid && x.age >= 18, humans)
    #     houseage = [humans[x].age for x in housecomp]
    # end
end

function dbg_farms() 
    farmers = findall(x -> x.fid > 0, humans)
    total_farmers = length(farmers)
    println("Total farmers: $total_farmers")

    # unique number of farms

    # we should check if the age of the farmers follows the distribution
    ages = [humans[i].age for i in farmers]
    ranges = [(20, 24), (25, 34), (35, 44), (45, 54), (55, 64), (65, 100)]
    counts = [count(x -> r[1] <= x <= r[2], ages) for r in ranges]
    for (r, c) in zip(ranges, counts)
        println("Range $(r): $c elements, proportion: $(round(c / total_farmers * 100))")
    end

    NUM_FARMS = 259 # from Google worksheet 
    over259 = findall(x -> x.fid > 259, humans) 
    length(over259) > 0 && @error "There are humans with fid > 259! This is not expected. Check the data."
    
    dict = Dict(k => 0 for k in 1:259)
    for h in humans 
        if h.fid > 0
            dict[h.fid] += 1
        end
    end
    
    # go through dict and see which farm is over 50 
    for (k, v) in dict
        if v > 50
            println("Farm $k has $v employees")
        end
    end
end