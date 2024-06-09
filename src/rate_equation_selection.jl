using Dates, CSV, DataFrames, Distributed, HypothesisTests
include("rate_equation_fitting.jl")


function prepare_data(data::DataFrame, metab_names)

    # Check if the column source exists and add it if it doesn't
    if !hasproperty(data, :source)
           #Add source column that uniquely identifies a figure from publication
        data.source .= data.Article .* "_" .* data.Fig
    end

    # Remove Na's
    data = data[.!isnan.(data.Rate), [:Rate, metab_names..., :source]]

    #Only include Rate > 0 because otherwise log_ratio_predict_vs_data() will have to divide by 0
    filter!(row -> row.Rate != 0, data)

    # Check if all values in metab_names are columns in the data
    missing_columns = setdiff(metab_names, Symbol.(names(data)))
    @assert isempty(missing_columns) "The following metab columns are missing from the data: $(join(missing_columns, ", "))"
    
    return data
end



function data_driven_rate_equation_selection(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    range_number_params::Tuple{Int,Int},
    forward_model_selection::Bool;
    n_reps_opt::Int = 20,
    maxiter_opt::Int = 50_000,
    model_selection_method = "denis",
    p_val_threshold = .4,
)
    
    data = prepare_data(data, metab_names)
    
    #generate param_removal_code_names by converting each mirror parameter for a and i into one name
    #(e.g. K_a_Metabolite1 and K_i_Metabolite1 into K_Metabolite1)
    param_removal_code_names = (
        [
            Symbol(replace(string(param_name), "_a_" => "_allo_")) for
            param_name in param_names if
            !contains(string(param_name), "_i") && param_name != :Vmax
        ]...,
    )

    #generate all possible combination of parameter removal codes
    param_subsets_per_n_params = calculate_all_parameter_removal_codes(param_names, range_number_params)

    if model_selection_method == "denis"
        # results = fit_rate_equation_selection_denis(
        #     general_rate_equation,
        #     data,
        #     metab_names,
        #     param_names,
        #     param_removal_code_names, 
        #     range_number_params,
        #     forward_model_selection,
        #     n_reps_opt,
        #     maxiter_opt,
        #     param_subsets_per_n_params,
        #     )
        println("finish stage 1!!")
        test_res = CSV.read("/home/ec2-user/code/DataDrivenEnzymeRateEqs.jl/test/Data_for_tests/pkm2_test_results_df.csv", DataFrame)
        train_res = CSV.read("/home/ec2-user/code/DataDrivenEnzymeRateEqs.jl/test/Data_for_tests/pkm2_train_results_df.csv", DataFrame)
        results = (train_results = train_res, test_results = test_res)
        
        best_n_params, best_subset = find_optimal_n_params(results.test_results, p_val_threshold)
        println("Best subset")
        println(best_subset)

        # find best_subset row in train_results
        best_subset_row = filter(row -> row.nt_param_removal_codes == best_subset, results.train_results)
        println(best_subset_row)


    elseif model_selection_method == "cv_denis"
        figs = unique(data.source) 
        results_figs_df = pmap(
            dropped_fig -> fit_rate_equation_selection_per_fig(
                general_rate_equation,
                data,
                metab_names,
                param_names,
                param_removal_code_names,
                range_number_params,
                forward_model_selection,
                n_reps_opt,
                maxiter_opt,
                param_subsets_per_n_params,
                all_param_removal_codes,
                dropped_fig
                ), 
            figs
        )
        results = vcat(results_figs_df...)

        best_n_params = find_best_n_params(results)

        # TODO: add train and choose best subset out of all subsets with best_n_params using all data


    elseif model_selection_method == "cv_all_subsets"
        results =  fit_rate_equation_selection_all_subsets(
            general_rate_equation,
            data,
            meta_names,
            param_names,
            param_removal_code_names,
            n_reps_opt, 
            maxiter_opt
            )

        # TODO: for each n params: keep the best model in terms of train loss
        # TODO: choose best num of params 
        # TODO: accordingly, choose best subset 

    end
    # TODO: decide how to choose best n params -> one sample differences wilcoxon test, need to choose threshold (p=.36?)
    return (results = results, best_n_params = best_n_params, best_subset_row = best_subset_row)
end

function get_nt_subset(df, num)
    # Filter the DataFrame where n_params equals num
    filtered_df = filter(row -> row.num_params == num, df)

    return filtered_df.nt_param_removal_codes[1]

end

# function find_best_n_params(df_results::DataFrame, print_res = true)
#     println("find best n params")
#     # Calculate average test loss for each n_params
#     avg_values = combine(groupby(df_results, :num_params), :test_loss_leftout_fig => mean => :avg_test_loss)

#     min_row = argmin(avg_values.avg_test_loss)
#     best_n_params =  avg_values[min_row, :].num_params
#     println("Best n params")
#     println(best_n_params)

#     best_subset = get_nt_subset(df_results, best_n_params)

#     if print_res == true
#         println("Avg CV error for each n removed params:")
#         println(sort(avg_values, :avg_test_loss))
#     end
#     return (best_n_params = best_n_params, best_subset = best_subset)
# end

function find_optimal_n_params(df_results::DataFrame, p_value_threshold::Float64) :: Int
    # Group by number of parameters and calculate average test loss
    grouped = groupby(df_results, :num_params)
    avg_losses = combine(grouped, :test_loss => mean => :avg_test_loss)
    # Sort by number of parameters
    sort!(avg_losses, :num_params)
    println("Avg CV error for each n params:")
    println(avg_losses)
    # Find the row with the minimum average test loss
    idx_min_loss = argmin(avg_losses.avg_test_loss)
    n_param_minimal_loss = avg_losses[idx_min_loss, :num_params]
    losses_minimal_loss = filter(row -> row.num_params == n_param_minimal_loss, df_results).test_loss

    current_n_params = n_param_minimal_loss
    # Start checking from the model just below the minimal average loss model downwards
    for i in idx_min_loss-1:-1:1        
        current_n_params = avg_losses[i, :num_params]
        # Perform Wilcoxon signed-rank test on test losses
        losses_current = filter(row -> row.num_params == current_n_params, df_results).test_loss
        # compare with best n params: 
        test_result = SignedRankTest(losses_current, losses_minimal_loss)
        pval = pvalue(test_result)

        # If the difference is not significant, continue; else, stop and return last non-significant model's params
        if pval <= p_value_threshold
            current_n_params = avg_losses[i+1, :num_params]
            break  # Stop if a significant difference is found
        end
    end
    
    best_n_params = current_n_params
    best_subset = get_nt_subset(df_results, best_n_params)
    
    return (best_n_params = best_n_params, best_subset = best_subset)

    return current_n_params
end

function train_and_choose_best_subset(data,param_subsets_per_n_params,  best_n_params; n_repetiotions_opt = 20, maxiter_opt = 50_000, print_res = false)
    nt_param_removal_codes = param_subsets_per_n_params[best_n_params]

    results_array = pmap(
        nt_param_removal_code -> train_rate_equation(
            general_rate_equation,
            data,
            metab_names,
            param_names;
            n_iter = n_repetiotions_opt,
            maxiter_opt = maxiter_opt,
            nt_param_removal_code = nt_param_removal_code,
        ),
        nt_param_removal_codes,
    )

    df_results = DataFrame(results_array)
    df_results.num_params = fill(best_n_params, nrow(df_results))
    df_results.nt_param_removal_codes = nt_param_removal_codes
    # cols: n_params, param_subset, train_loss, params
    println(first(df_results, 5)) 

    best_param_subset = DataFrame(results_df[argmin(results_df.train_loss),:])
    println("Best subset: $(best_param_subset.param_subset)")

    return best_param_subset
end


function fit_rate_equation_selection_denis(
        general_rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}},
        param_removal_code_names, 
        range_number_params::Tuple{Int,Int},
        forward_model_selection::Bool,
        n_repetiotions_opt::Int,
        maxiter_opt::Int,
        param_subsets_per_n_params,
        )

        if forward_model_selection
            num_param_range = (range_number_params[2]):-1:range_number_params[1]
        elseif !forward_model_selection
            num_param_range = (range_number_params[1]):1:range_number_params[2]
        end
        starting_param_removal_codes = param_subsets_per_n_params[num_param_range[1]]
    
        previous_param_removal_codes = starting_param_removal_codes
        println("About to start loop with num_params: $num_param_range")
        
        df_train_results = DataFrame()
        df_test_results = DataFrame()
        for num_params in num_param_range
            println("Running loop with num_params: $num_params")
    
            #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
            if forward_model_selection
                nt_param_removal_codes = forward_selection_next_param_removal_codes(
                    param_subsets_per_n_params,
                    previous_param_removal_codes,
                    num_params,
                    param_names,
                    param_removal_code_names,
                )
            elseif !forward_model_selection
                nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                    param_subsets_per_n_params,
                    previous_param_removal_codes,
                    num_params,
                    param_names,
                    param_removal_code_names,
                )
            end
            println("nt_param_removel_codes", length(nt_param_removal_codes))
            # TODO: change to pmap after debugging
            #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
            results_array = pmap(
                nt_param_removal_code -> train_rate_equation(
                    general_rate_equation,
                    data,
                    metab_names,
                    param_names;
                    n_iter = n_repetiotions_opt,
                    maxiter_opt = maxiter_opt,
                    nt_param_removal_code = nt_param_removal_code,
                ),
                nt_param_removal_codes,
            )
    
            #convert results_array to DataFrame
            df_results = DataFrame(results_array)
            df_results.num_params = fill(num_params, nrow(df_results))
            df_results.nt_param_removal_codes = nt_param_removal_codes
            df_train_results = vcat(df_train_results, df_results)
    
            # Optinally consider saving results to csv file for long running calculation of cluster
            # CSV.write(
            #     "$(Dates.format(now(),"mmddyy"))_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
            #     df_results,
            # )

            #if all train_loss are Inf, then skip to next loop
            if all(df_results.train_loss .== Inf)
                previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
                continue
            end

            #store top 10% for next loop as `previous_param_removal_codes`
            filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)
            previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
    
            #calculate loocv test loss for top subset for each `num_params`
            best_nt_param_removal_code =
                df_results.nt_param_removal_codes[argmin(df_results.train_loss)]
            
            # TODO: move test_results out from the loop
            # test_results = pmap(
            #     removed_fig -> loocv_rate_equation(
            #         removed_fig,
            #         general_rate_equation,
            #         data,
            #         metab_names,
            #         param_names;
            #         n_iter = n_repetiotions_opt,
            #         maxiter_opt = maxiter_opt,
            #         nt_param_removal_code = best_nt_param_removal_code,
            #     ),
            #     unique(data.source),
            # )

            # df_results = DataFrame(test_results)
            # df_results.num_params = fill(num_params, nrow(df_results))
            # df_results.nt_param_removal_codes =
            #     fill(best_nt_param_removal_code, nrow(df_results))

            df_results = DataFrame(:num_params => [num_params], :nt_param_removal_codes => [best_nt_param_removal_code])
            df_test_results = vcat(df_test_results, df_results)
        end
        println("there are ", size(df_test_results)[1], "best models")
        # calculate loocv test loss for top subsets:
        # Prepare the data for pmap
        subsets_to_fit = [(row.nt_param_removal_codes, removed_fig, row.num_params) for row in eachrow(df_test_results) for removed_fig in unique(data.source)]

        results = pmap(
            subset -> loocv_rate_equation(
                subset[2], #removed_fig
                general_rate_equation,
                data,
                metab_names,
                param_names;
                n_iter = n_repetiotions_opt,
                maxiter_opt = maxiter_opt,
                nt_param_removal_code = subset[1],
            ), 
        subsets_to_fit
        )
        
        result_dfs = DataFrame[]
        for (res, subset) in zip(results, subsets_to_fit)
            res_df = DataFrame([res])
            res_df[!, :nt_param_removal_codes] = [subset[1]]
            res_df[!, :num_params] = [subset[3]]
            push!(result_dfs, res_df)
        end

        df_test_results = vcat(result_dfs...)
        println("size of df_test_results: ", size(df_test_results))

        return (train_results = df_train_results, test_results = df_test_results)
    

end


function fit_rate_equation_selection_per_fig(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    range_number_params::Tuple{Int,Int},
    forward_model_selection::Bool,
    n_repetiotions_opt::Int,
    maxiter_opt::Int,
    param_subsets_per_n_params,
    all_param_removal_codes,
    test_fig
    )

    train_data = data[data.source.!=test_fig, :]
    test_data = data[data.source.==test_fig, :]

    if forward_model_selection
        num_param_range = (range_number_params[2]):-1:range_number_params[1]
        starting_param_removal_codes = param_subsets_per_n_params[range_number_params[2]]
    elseif !forward_model_selection
        num_param_range = (range_number_params[1]):1:range_number_params[2]
        starting_param_removal_codes = param_subsets_per_n_params[range_number_params[1]]
    end

    previous_param_removal_codes = starting_param_removal_codes
    println("About to start loop with num_params: $num_param_range")
    
    df_train_results = DataFrame()
    df_test_results = DataFrame()
    for num_params in num_param_range
        println("Running loop with num_params: $num_params")

        #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
        if forward_model_selection
            nt_param_removal_codes = forward_selection_next_param_removal_codes(
                param_subsets_per_n_params,
                previous_param_removal_codes,
                num_params,
                param_names,
                param_removal_code_names,
            )
        elseif !forward_model_selection
            nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                param_subsets_per_n_params,
                previous_param_removal_codes,
                num_params,
                param_names,
                param_removal_code_names,
            )
        end

        #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
        results_array = map(
            nt_param_removal_code -> train_rate_equation(
                general_rate_equation,
                train_data,
                metab_names,
                param_names;
                n_iter = n_repetiotions_opt,
                maxiter_opt = maxiter_opt,
                nt_param_removal_code = nt_param_removal_code,
            ),
            nt_param_removal_codes,
        )

        #convert results_array to DataFrame
        df_results = DataFrame(results_array)
        df_results.num_params = fill(num_params, nrow(df_results))
        df_results.nt_param_removal_codes = nt_param_removal_codes
        df_train_results = vcat(df_train_results, df_results)

        # Optinally consider saving results to csv file for long running calculation of cluster
        # CSV.write(
        #     "$(Dates.format(now(),"mmddyy"))_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
        #     df_results,
        # )
        #store top 10% for next loop as `previous_param_removal_codes`
        filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)
        previous_param_removal_codes = values.(df_results.nt_param_removal_codes)

        #calculate loocv test loss for top subset for each `num_params`
        #TODO: change to pmap
        best_nt_param_removal_code =
            df_results.nt_param_removal_codes[argmin(df_results.train_loss)]
        best_subset_rescaled_params = df_results.params

        test_loss = test_rate_equation(
                general_rate_equation,
                test_data,
                best_subset_rescaled_params,
                metab_names,
                param_names
        )

        df_results =  DataFrame(
            test_loss = test_loss,          
            num_params = num_params,          
            nt_param_removal_code =best_nt_param_removal_code, 
            test_fig =test_fig,
            params = best_subset_rescaled_params           
        )
            
        df_test_results = vcat(df_test_results, df_results)
    end

    return (train_results = df_train_results, test_results = df_test_results)

end


function fit_rate_equation_selection_all_subsets(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    n_repetiotions_opt::Int,
    maxiter_opt::Int,
    )

    figs = unique(data.source)

    # Initialize an empty list for the combined results
    all_subsets_figs_to_fit = []
    lengths = []

    for (n_params, subsets) in param_subsets_per_n_params

        nt_param_subsets = [
            NamedTuple{param_removal_code_names}(x) for
            x in unique(param_removal_codes)
        ]

        # Create the product for this particular number of parameters
        temp_product = collect(Iterators.product(nt_param_subsets, figs))
        # Append the product to the main list
        append!(all_subsets_figs_to_fit, temp_product)
        # Record the length of the product
        push!(lengths, length(temp_product))
    end
    
    # Create the parameter mapping using the recorded lengths
    n_params_mapping = Int[]
    for (n_params, length) in zip(keys(param_subsets_per_n_params), lengths)
        append!(n_params_mapping, fill(n_params, length))
    end

    results_array = pmap(
        subset_fig_to_fit -> loocv_rate_equation(
            subset_fig_to_fit[2],
            general_rate_equation,
            data,
            metab_names,
            param_names;
            n_iter = n_repetiotions_opt,
            maxiter_opt = maxiter_opt,
            nt_param_removal_code = subset_fig_to_fit[1],
        ),
        all_subsets_figs_to_fit, 
    )

    df_results = DataFrame(results_array)
    df_results.num_params = n_params_mapping
    all_subsets  = [item[1] for item in all_subsets_figs_to_fit]
    df_results.nt_param_removal_codes = all_subsets
       
    return (train_test_results = df_results)

end


"function to calculate train loss without a figure and test loss on removed figure"
function loocv_rate_equation(
    fig,
    rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    n_iter = 20,
    maxiter_opt = 50_000,
    nt_param_removal_code = nothing,
)
    # Drop selected figure from data
    train_data = data[data.source.!=fig, :]
    test_data = data[data.source.==fig, :]
    # Calculate fit
    train_res = train_rate_equation(
        rate_equation,
        train_data,
        metab_names,
        param_names;
        n_iter = n_iter,
        maxiter_opt = maxiter_opt,
        nt_param_removal_code = nt_param_removal_code,
    )
    test_loss = test_rate_equation(
        rate_equation,
        test_data,
        train_res.params,
        metab_names,
        param_names,
    )
    return (
        dropped_fig = fig,
        train_loss_wo_fig = train_res.train_loss,
        test_loss_leftout_fig = test_loss,
        params = train_res.params,
    )
end

"""Function to calculate loss for a given `rate_equation` and `nt_fitted_params` on `data` that was not used for training"""
function test_rate_equation(
    rate_equation::Function,
    data,
    nt_fitted_params::NamedTuple,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
)
    filtered_data = data[.!isnan.(data.Rate), [:Rate, metab_names..., :source]]
    # Add a new column to data to assign an integer to each source/figure from publication
    filtered_data.fig_num = vcat(
        [
            i * ones(
                Int64,
                count(==(unique(filtered_data.source)[i]), filtered_data.source),
            ) for i = 1:length(unique(filtered_data.source))
        ]...,
    )
    # Add a column containing indexes of points corresponding to each figure
    fig_point_indexes =
        [findall(filtered_data.fig_num .== i) for i in unique(filtered_data.fig_num)]
    # Convert DF to NamedTuple for better type stability / speed
    rate_data_nt = Tables.columntable(filtered_data)

    fitted_params = values(nt_fitted_params)
    test_loss = loss_rate_equation(
        fitted_params,
        rate_equation::Function,
        rate_data_nt::NamedTuple,
        param_names::Tuple{Symbol,Vararg{Symbol}},
        fig_point_indexes::Vector{Vector{Int64}};
        rescale_params_from_0_10_scale = false,
        nt_param_removal_code = nothing,
    )
    return test_loss
end


function calculate_number_of_parameters(x,n, num_alpha_params)
    return n - num_alpha_params - sum(x[1:end-num_alpha_params] .> 0)
end

"""Generate all possibles codes for ways that mirror params for a and i states of MWC enzyme can be removed from the rate equation"""
function calculate_all_parameter_removal_codes(param_names::Tuple{Symbol,Vararg{Symbol}}, range_number_params::Tuple{Int,Int})
    feasible_param_subset_codes = ()
    for param_name in param_names
        param_name_str = string(param_name)
        if param_name == :L
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(param_name_str, "Vmax_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1,2])
        elseif startswith(param_name_str, "K_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1,2,3])
        elseif startswith(param_name_str, "K_") &&
               !startswith(param_name_str, "K_i") &&
               !startswith(param_name_str, "K_a") &&
               length(split(param_name_str, "_")) == 2
               feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(param_name_str, "K_") &&
               !startswith(param_name_str, "K_i") &&
               !startswith(param_name_str, "K_a") &&
               length(split(param_name_str, "_")) > 2
               feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1,2])
        elseif startswith(string(param_name), "alpha")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        end
    end

    all_param_removal_codes = Iterators.product(feasible_param_subset_codes...)
    n_param_subset = length(first(all_param_removal_codes))
    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    n = length(param_names)

    total_elements = prod(length.(feasible_param_subset_codes))
    # keep for each number of params: all the subsets with this number
    # TODO: TRY FIX THIS
    param_subsets_per_n_params = Dict{Int, Vector{NTuple{n_param_subset, Int}}}()
    println("before param subsets per n params")
    for x in Iterators.take(all_param_removal_codes, 30000)
    # for x in all_param_removal_codes
        n_param = n - num_alpha_params - sum(x[1:end-num_alpha_params] .> 0)
        #param_subset = values(x)
        # Organize into the dictionary
        if !haskey(param_subsets_per_n_params, n_param)
            param_subsets_per_n_params[n_param] = Vector{NTuple{n_param_subset, Int}}()
        end
        push!(param_subsets_per_n_params[n_param], x)
    end
    println("Memory usage of dictionary: ", Base.summarysize(param_subsets_per_n_params) / (1024^3), " GiB")

    println("after param_subsets_per_n_params")
    # param_subsets_tuple = [(
    #     length(param_names) - num_alpha_params - sum(x[1:end-num_alpha_params] .> 0),
    #     values(x) 
    # ) for x in all_param_removal_codes]
   
    # param_subsets_per_n_params = Dict{Int, Vector}()
    # for (key, value) in param_subsets_tuple
    #     if haskey(param_subsets_per_n_params, key)
    #         push!(param_subsets_per_n_params[key], value)
    #     else
    #         param_subsets_per_n_params[key] = [value]
    #     end
    # end

    #check that range_number_params within bounds of minimal and maximal number of parameters
    # TODO: uncomment these lines after debugging
    # @assert range_number_params[1] >=
    # length(param_names) - maximum([sum(x .> 0) for x in all_param_removal_codes]) "starting range_number_params cannot be below $(length(param_names) - maximum([sum(x .> 0) for x in all_param_removal_codes]))"
    # @assert range_number_params[2] <= length(param_names) "ending range_number_params cannot be above $(length(param_names))"

    return param_subsets_per_n_params
end

"""
Function to convert parameter vector to vector where some params are equal to 0, Inf or each other based on nt_param_removal_code
"""
function param_subset_select_denis(params, param_names, nt_param_removal_code)
    @assert length(params) == length(param_names)
    params_dict =
        Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

    for param_choice in keys(nt_param_removal_code)
        if startswith(string(param_choice), "L") && nt_param_removal_code[param_choice] == 1
            params_dict[:L] = 0.0
        elseif startswith(string(param_choice), "Vmax") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[:Vmax_i] = params_dict[:Vmax_a]
        elseif startswith(string(param_choice), "Vmax") &&
               nt_param_removal_code[param_choice] == 2
            global params_dict[:Vmax_i] = 0.0
        elseif startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 1
            K_i = Symbol("K_i_" * string(param_choice)[8:end])
            K_a = Symbol("K_a_" * string(param_choice)[8:end])
            params_dict[K_i] = params_dict[K_a]
        elseif startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 2
            K_a = Symbol("K_a_" * string(param_choice)[8:end])
            params_dict[K_a] = Inf
        elseif startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 3
            K_i = Symbol("K_i_" * string(param_choice)[8:end])
            params_dict[K_i] = Inf
        elseif startswith(string(param_choice), "K_") &&
               !startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[param_choice] = Inf
        elseif startswith(string(param_choice), "K_") &&
               !startswith(string(param_choice), "K_allo") &&
               length(split(string(param_choice), "_")) > 2 &&
               nt_param_removal_code[param_choice] == 2
            params_dict[param_choice] =
                prod([
                    params_dict[Symbol("K_" * string(metab))] for
                    metab in split(string(param_choice), "_")[2:end]
                ])^(1 / (length(split(string(param_choice), "_")[2:end])))
        elseif startswith(string(param_choice), "alpha") &&
               nt_param_removal_code[param_choice] == 0
            params_dict[param_choice] = 0.0
        elseif startswith(string(param_choice), "alpha") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[param_choice] = 1.0
        end
    end

    new_params_sorted = [params_dict[param_name] for param_name in param_names]
    return new_params_sorted
end

function param_subset_select(params, param_names, nt_param_removal_code)
    @assert length(params) == length(param_names)
    params_dict =
        Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

    # for param_choice in keys(nt_param_removal_code)
    for (param, choice) in pairs(nt_param_removal_code)
        param_str = string(param)

        # handle K params
        if startswith(param_str, "K_allo")
            param_name = split(param_str, "K_allo_")[2]
            K_i = Symbol("K_i_" * param_name)
            K_a = Symbol("K_a_" * param_name)
      
            if choice > 0
                if choice == 1
                    params_dict[K_i] = params_dict[K_a] 

                elseif choice == 2
                    params_dict[K_a] = Inf

                elseif choice == 3
                    params_dict[K_i] = Inf
                end
            end

        elseif startswith(param_str, "K_") && !startswith(param_str, "K_allo")
            if choice == 1
                params_dict[Symbol(param_str)] = Inf
            elseif length(split(param_str, "_")) > 2 && choice == 2
                metabs = split(param_str, "_")[2:end]
                params_dict[Symbol(param_str)] = prod(params_dict[Symbol("K_" * metab)] for metab in metabs) ^ (1 / length(metabs))
            end
        
        elseif startswith(param_str, "alpha")
            if choice == 0
                params_dict[Symbol(param_str)] = 0.0
            elseif choice == 1
                params_dict[Symbol(param_str)] = 1.0
            end
            
        elseif param_str == "Vmax"
            if choice == 1
                params_dict[:Vmax_i] = params_dict[:Vmax_a]
            elseif choice == 2
                #TODO: check why it's appear with global in denis's code
                params_dict[:Vmax_i] = 0.0
            end

        elseif startswith(param_str, "L")
            if choice == 1
                params_dict[:L] = 0.0
            end

        end
    end
    new_params_sorted = [params_dict[param_name] for param_name in param_names]
    return new_params_sorted
end
"""
Calculate `nt_param_removal_codes` with `num_params` including non-zero term combinations for codes (excluding alpha terms) in each `previous_param_removal_codes` that has `num_params-1`
"""
function forward_selection_next_param_removal_codes(
    param_subsets_per_n_params,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
    )
    
    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    @assert all([
        (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params + 1
        ) || (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params
        ) for param_removal_code in previous_param_removal_codes
    ])
    
    previous_param_subset_masks = unique([
        (
            mask = (
                (previous_param_removal_code[1:(end-num_alpha_params)] .== 0)...,
                zeros(Int64, num_alpha_params)...,
            ),
            non_zero_params = previous_param_removal_code .*
                              (previous_param_removal_code .!= 0),
        ) for previous_param_removal_code in previous_param_removal_codes
    ])

    #select all param_removal_codes that yield equations with `num_params` number of parameters
    all_param_codes_w_num_params = param_subsets_per_n_params[num_params]

    #choose param_removal_codes with n_removed_params number of parameters removed that also contain non-zero elements from previous_param_removal_codes
    param_removal_codes = []
    for previous_param_subset_mask in previous_param_subset_masks
        push!(
            param_removal_codes,
            unique([
                param_code_w_num_params .* previous_param_subset_mask.mask .+
                previous_param_subset_mask.non_zero_params for
                param_code_w_num_params in all_param_codes_w_num_params
            ])...,
        )
    end
    nt_param_removal_codes = [
        NamedTuple{param_removal_code_names}(x) for
        x in unique(param_removal_codes) if (
            length(param_names) - num_alpha_params - sum(x[1:(end-num_alpha_params)] .> 0)
        ) == num_params
    ]
    return nt_param_removal_codes
end

"""
Calculate `param_removal_codes` with `num_params` including zero term combinations for codes (excluding alpha terms) in each `previous_param_removal_codes` that has `num_params+1`
"""
function reverse_selection_next_param_removal_codes(
    param_subsets_per_n_params,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
)

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    @assert all([
        (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params - 1
        ) || (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params
        ) for param_removal_code in previous_param_removal_codes
    ])
    previous_param_subset_masks = unique([
        (
            mask = [
                (previous_param_removal_code[1:(end-num_alpha_params)] .== 0)...,
                zeros(Int64, num_alpha_params)...,
            ],
            non_zero_params = previous_param_removal_code .*
                              (previous_param_removal_code .!= 0),
        ) for previous_param_removal_code in previous_param_removal_codes
    ])

    #select all codes that yield equations with `num_params` number of parameters
    all_param_codes_w_num_params = param_subsets_per_n_params[num_params]

    #choose param_removal_codes with n_removed_params number of parameters removed that also contain non-zero elements from previous_param_removal_codes
    param_removal_codes = []
    for previous_param_subset_mask in previous_param_subset_masks
        push!(
            param_removal_codes,
            unique([
                previous_param_subset_mask.non_zero_params .*
                (param_code_w_num_params .!= 0) for
                param_code_w_num_params in all_param_codes_w_num_params
            ])...,
        )
    end
    nt_param_removal_codes = [
        NamedTuple{param_removal_code_names}(x) for
        x in unique(param_removal_codes) if (
            length(param_names) - num_alpha_params - sum(x[1:(end-num_alpha_params)] .> 0)
        ) == num_params
    ]
    return nt_param_removal_codes
end

# Compare model performances of different number of parameters based on test losses using the Wilcoxon signed-rank test.
function compare_models_wilcoxon(df::DataFrame, method::Symbol)
    # Sort the DataFrame by the number of parameters
    sort!(df, :num_params)
    
    # Group data by number of parameters and collect test losses
    grouped = groupby(df, :num_params)
    losses = [group[!, :test_loss] for group in grouped]
    
    n = length(losses)
    results = []

    if method == :all_pairs
        # Comparing all pairs of models
        for i in 1:n
            for j in i+1:n
                test_result = SignedRankTest(losses[i], losses[j])
                push!(results, (model_a_num_params = grouped[i][1, :num_params], 
                                model_b_num_params = grouped[j][1, :num_params], 
                                p_value = pvalue(test_result)))
            end
        end
    elseif method == :forward_stepwise
        # Comparing each model with the next one (increasing number of parameters)
        for i in 1:n-1
            test_result = SignedRankTest(losses[i], losses[i+1])
            push!(results, (model_a_num_params = grouped[i][1, :num_params], 
                            model_b_num_params = grouped[i+1][1, :num_params], 
                            p_value = pvalue(test_result)))
        end
    elseif method == :backward_stepwise
        # Comparing each model with the previous one (decreasing number of parameters)
        for i in n:-1:2
            test_result = SignedRankTest(losses[i], losses[i-1])
            push!(results, (model_a_num_params = grouped[i][1, :num_params], 
                            model_b_num_params = grouped[i-1][1, :num_params], 
                            p_value = pvalue(test_result)))
        end
    else
        error("Invalid method specified. Choose :all_pairs, :forward_stepwise, or :backward_stepwise")
    end
    
    return DataFrame(results)
end
