using Dates, CSV, DataFrames, Distributed, HypothesisTests

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


"""
    data_driven_rate_equation_selection(
        general_rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}};
        range_number_params::Union{Nothing, Tuple{Int,Int}} = nothing,
        forward_model_selection::Bool = true,
        max_zero_alpha::Int = 1 + ceil(Int, length(metab_names) / 2),
        n_reps_opt::Int = 20, 
        maxiter_opt::Int = 50_000,
        model_selection_method::String = "current_subsets_filtering",
        p_val_threshold::Float64 = 0.4,
        save_train_results::Bool = false,
        enzyme_name::String = "Enzyme",
        subsets_min_limit::Int = 1, 
        subsets_max_limit::Union{Int, Nothing}=nothing,
        subsets_filter_threshold::Float64=0.1,
    )

This function is used to perform data-driven rate equation selection using a general rate equation and data. 

There are three model_selection methods:

1. current_subsets_filtering: 
This method iteratively fits models that are subsets of the top 10% from the previous iteration,
saving the best model for each n params based on training loss. Optimal number of parameters are selected using 
the Wilcoxon test on test scores from LOOCV, and the best equation is the best model with this optimal number.

2. cv_subsets_filtering:
This method implements current_subsets_filtering separately for each figure,
leaving one figure out as a test set while training on the remaining data.
For each number of parameters, it saves the test loss of the best subset for that figure.
It uses the Wilcoxon test across all figures' results to select the optimal number of parameters. 
Then, for the chosen number, it trains all subset with this n params on the entire dataset and selects the best
rate equation based on minimal training loss.

3. cv_all_subsets: 
This method fits all subsets for each figure, using the others as training data and the left-out figure as the test set.
It selects the best model for each number of parameters and figure based on training error and computes LOOCV test scores. 
The optimal n params is determined by the Wilcoxon test across all figures' test scores.
The best equation is the subset with minimal training loss for this optimal n params when trained on the entire dataset.    

# Arguments
- `general_rate_equation::Function`: Function that takes a NamedTuple of metabolite concentrations (with `metab_names` keys) and parameters (with `param_names` keys) and returns an enzyme rate.
- `data::DataFrame`: DataFrame containing the data with column `Rate` and columns for each `metab_names` where each row is one measurement. It also needs to have a column `source` that contains a string that identifies the source of the data. This is used to calculate the weights for each figure in the publication.
- `metab_names::Tuple`: Tuple of metabolite names that correspond to the metabolites of `rate_equation` and column names in `data`.
- `param_names::Tuple`: Tuple of parameter names that correspond to the parameters of `rate_equation`.

# Keyword Arguments
- `save_train_results::Bool`: A boolean indicating whether to save the results of the training for each number of parameters as a csv file.
- `enzyme_name::String`: A string for enzyme name that is used to name the csv files that are saved.
- `range_number_params::Tuple{Int,Int}`: A tuple of integers representing the range of the number of parameters of general_rate_equation to search over.
- `forward_model_selection::Bool`: A boolean indicating whether to use forward model selection (true) or reverse model selection (false).
- `max_zero_alpha::Int`: An integer representing the maximum number of alpha parameters that can be set to 0.
- `n_reps_opt::Int` n repetitions of optimization  
- `maxiter_opt::Int` max iterations of optimization algorithm
-  `model_selection_method::String` - which model selection to find best rate equation (default is current_subsets_filtering)
-  `p_val_threshold::Float64` - pval threshold for Wilcoxon test
- `save_train_results::Bool`: A boolean indicating whether to save the results of the training for each number of parameters as a csv file.
- `enzyme_name::String`: A string for enzyme name that is used to name the csv files that are saved.
- `subsets_min_limit::Int` - The minimum number of filtered subsets (those with training loss within 10% of the minimum) 
that must be kept for each number of parameters. These subsets are used to generate the subsets for the next iteration (only subsets of these are considered).
Relevant to model selection methods current_subsets_filtering or cv_subsets_filtering.
- `subsets_max_limit::Union{Int, Nothing}` - The maximum number of filtered subsets (those with training loss within 10% of the minimum) 
that must be kept for each number of parameters. These subsets are used to generate the subsets for the next iteration (only subsets of these are considered).
Relevant to model selection methods current_subsets_filtering or cv_subsets_filtering.
- `subsets_filter_threshold::Float64` - This sets the percentage limit for filtering subsets in each iteration.
 Only the subsets with a training loss close to the best (within this percentage) are kept.
 Relevant to model selection methods current_subsets_filtering or cv_subsets_filtering.

# Returns
- `NamedTuple`: A named tuple with the following fields:
  - `results`: df with train and test results
  - `best_n_params`: optimal number of parameters
  - `best_subset_row`: row of the best rate equation selected - includes fitted params 
"""
function data_driven_rate_equation_selection(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    range_number_params::Union{Nothing, Tuple{Int,Int}} = nothing,
    forward_model_selection::Bool = true,
    max_zero_alpha::Int = 1 + ceil(Int, length(metab_names) / 2),
    n_reps_opt::Int = 20, 
    maxiter_opt::Int = 50_000,
    model_selection_method::String = "current_subsets_filtering",
    p_val_threshold::Float64 = 0.05,
    save_train_results::Bool = false,
    enzyme_name::String = "Enzyme",
    subsets_min_limit::Int = 1, 
    subsets_max_limit::Union{Int, Nothing}=nothing,
    subsets_filter_threshold::Float64=0.1,
    train_loss = "sse",
    test_loss = "sse"
 )
    
    data = prepare_data(data, metab_names)
    
    #generate param_removal_code_names by converting each mirror parameter for a and i into one name
    #(e.g. K_a_Metabolite1 and K_i_Metabolite1 into K_allo_Metabolite1)
    param_removal_code_names = (
        [
            Symbol(replace(string(param_name), "_a_" => "_allo_", "Vmax_a" => "Vmax_allo")) for param_name in param_names if !contains(string(param_name), "_i") &&
            param_name != :Vmax &&
            param_name != :L
        ]...,
    )

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))

    if isnothing(range_number_params)
        if :L in param_names
            range_number_params =
                (length(metab_names) + 2, length(param_names) - num_alpha_params)
        else
            range_number_params =
                (length(metab_names) + 1, length(param_names) - num_alpha_params)
        end
    end

    #calculate starting_param_removal_codes parameters
    practically_unidentifiable_params = find_practically_unidentifiable_params(data, param_names)
    all_param_removal_codes = calculate_all_parameter_removal_codes(
        param_names,
        practically_unidentifiable_params,
    )

    if model_selection_method == "current_subsets_filtering"
        results = fit_rate_equation_selection_current(
            general_rate_equation,
            data,
            metab_names,
            param_names,
            param_removal_code_names, 
            range_number_params,
            forward_model_selection,
            max_zero_alpha,
            n_reps_opt,
            maxiter_opt,
            all_param_removal_codes, 
            practically_unidentifiable_params, 
            save_train_results, 
            enzyme_name,
            subsets_min_limit,
            subsets_max_limit,
            subsets_filter_threshold,
            train_loss,
            test_loss
            )

        best_n_params = find_optimal_n_params(results.test_results, p_val_threshold)
        best_subset = get_nt_subset(results.test_results, best_n_params)
        println("Best subset")
        println(best_subset)

        # find best_subset row in train_results
        best_subset_row = filter(row -> row.nt_param_removal_codes == best_subset, results.train_results)
        println("best subset row")
        println(best_subset_row)


    elseif model_selection_method == "cv_subsets_filtering"
        figs = unique(data.source) 
        results_figs_df = map(
            dropped_fig -> fit_rate_equation_selection_per_fig(
                general_rate_equation,
                data,
                metab_names,
                param_names,
                param_removal_code_names,
                range_number_params,
                forward_model_selection,
                max_zero_alpha,
                n_reps_opt,
                maxiter_opt,
                all_param_removal_codes,
                practically_unidentifiable_params, 
                dropped_fig, 
                subsets_min_limit,
                subsets_max_limit, 
                subsets_filter_threshold
                ), 
            figs
        )
        train_results = [res.train_results for res in results_figs_df]
        test_results = [res.test_results for res in results_figs_df]
        combined_train_results = vcat(train_results...)
        combined_test_results = vcat(test_results...)
        results = (train_results =combined_train_results, test_results =combined_test_results )

        best_n_params = find_optimal_n_params(results.test_results, p_val_threshold)

        best_subset_row = train_and_choose_best_subset(
            general_rate_equation, 
            data,
            all_param_removal_codes, 
            best_n_params,
            metab_names, 
            param_names, 
            param_removal_code_names, 
            practically_unidentifiable_params, 
            max_zero_alpha,
            n_reps_opt, 
            maxiter_opt, 
            save_train_results, 
            enzyme_name
        )
        println("best subset row")
        println(best_subset_row)

    elseif model_selection_method == "cv_all_subsets"

        results = fit_rate_equation_selection_all_subsets(
            general_rate_equation,
            data,
            all_param_removal_codes, 
            practically_unidentifiable_params, 
            max_zero_alpha,
            metab_names,
            param_names,
            param_removal_code_names,
            n_reps_opt, 
            maxiter_opt
        )

        # This code groups results by dropped_fig and num_params, finds the row with the minimum train_loss in each group,
        # and creates a new DataFrame with dropped_fig, test_loss, and num_params.
        grouped = groupby(results.train_test_results, [:dropped_fig, :num_params])
        agg_results = combine(grouped) do subdf
            idx = argmin(subdf.train_loss)
            subdf[idx, [:dropped_fig, :test_loss, :num_params]]
        end

        best_n_params = find_optimal_n_params(agg_results, p_val_threshold)

        best_subset_row = train_and_choose_best_subset(
            general_rate_equation, 
            data,
            all_param_removal_codes, 
            best_n_params,
            metab_names, 
            param_names, 
            param_removal_code_names, 
            practically_unidentifiable_params, 
            max_zero_alpha,
            n_reps_opt, 
            maxiter_opt, 
            save_train_results, 
            enzyme_name
        )
        println("best subset row")
        println(best_subset_row)
    else
       throw(ArgumentError("Invalid model selection method $(model_selection_method)"))
    end
    return (results = results, best_n_params = best_n_params, best_subset_row = best_subset_row)
end

function get_nt_subset(df, num)
    # Filter the DataFrame where n_params equals num
    filtered_df = filter(row -> row.num_params == num, df)

    return filtered_df.nt_param_removal_codes[1]

end

"""
    find_best_n_params_wilcoxon(
        df_results::DataFrame, 
        avg_losses::DataFrame, 
        p_value_threshold::Float64, 
        n_param_minimal_loss::Int, 
        losses_minimal_loss::Vector{Float64}
    ) :: Int

This function identifies the best number of parameters for a model based on Wilcoxon signed-rank tests.

## Parameters:
- `df_results::DataFrame`: A DataFrame containing test results, including the number of parameters and corresponding log test losses.
- `avg_losses::DataFrame`: A DataFrame with the average test log losses for different numbers of parameters.
- `p_value_threshold::Float64`: The threshold p-value for determining if the difference in losses is statistically significant.
- `n_param_minimal_loss::Int`: The number of parameters for the model with the minimal average test log loss.
- `losses_minimal_loss::Vector{Float64}`: The log test losses for the model with the minimal average test log loss.

## Returns:
- `best_n_params::Int`: The number of parameters that provides the best model performance based on the Wilcoxon test results.

## Description:
The function starts by finding the model with the minimal average test log loss and then iterates through models with fewer parameters. 
    For each model, it performs a Wilcoxon signed-rank test comparing its test losses to those of the minimal loss model. 
    The p-values from these tests are stored, and the function determines the best number of parameters based on the smallest model that has a p-value above the given threshold.
"""
function find_best_n_params_wilcoxon(
    df_results::DataFrame, 
    avg_losses::DataFrame, 
    p_value_threshold::Float64, 
    n_param_minimal_loss::Int, 
    losses_minimal_loss::Vector{Float64}
) :: Int
    idx_min_loss = argmin(avg_losses.avg_test_log_loss)
    
    wilcoxon_results = DataFrame(num_params = Int[], pval= Float64[])
    current_n_params = n_param_minimal_loss
    # Start checking from the model just below the minimal average loss model downwards
    for i in idx_min_loss-1:-1:1        
        current_n_params = avg_losses[i, :num_params]
        # Perform Wilcoxon signed-rank test on test losses
        losses_current = filter(row -> row.num_params == current_n_params, df_results).log_test_loss
        # compare with best n params: 
        test_result = ExactSignedRankTest(losses_current, losses_minimal_loss)
        pval = pvalue(test_result)
        push!(wilcoxon_results, (current_n_params, pval))
    end
    println(wilcoxon_results)
    best_n_params = n_param_minimal_loss
    
    if !isempty(wilcoxon_results)
        above_threshold = wilcoxon_results[wilcoxon_results.pval .> p_value_threshold, :]
        if !isempty(above_threshold)
        best_n_params = minimum(above_threshold.num_params)
        end
    end
    println("Best n params Wilcoxon: $(best_n_params)")
    return best_n_params
end


"""
    find_best_n_params_within_one_se(
        losses_minimal_loss::Vector{Float64}, 
        avg_losses::DataFrame, 
        n_param_minimal_loss::Int
    ) :: Int

This function identifies the best number of parameters based on the "one standard error" (1-SE) rule.

## Parameters:
- `losses_minimal_loss::Vector{Float64}`: The log test losses for the model with the minimal average test log loss.
- `avg_losses::DataFrame`: A DataFrame containing the average test log losses for different numbers of parameters.
- `n_param_minimal_loss::Int`: The number of parameters for the model with the minimal average test log loss.

## Returns:
- `best_n_prams_se::Int`: The number of parameters that provides the best model performance within one standard error of the minimal average loss.

## Description:
This function applies the 1-SE rule to select the best number of parameters. 
    It calculates the average test log loss and the standard error (SE) for the model with the minimal loss. 
        It then filters models that have fewer or equal parameters and whose average loss is within one standard error of the minimal average loss. 
        Finally, it selects and returns the model with the fewest parameters that satisfies this condition.
"""
function find_best_n_params_within_one_se(
    losses_minimal_loss::Vector{Float64}, 
    avg_losses::DataFrame, 
    n_param_minimal_loss::Int
) :: Int
    best_log_avg_loss = mean(losses_minimal_loss)
    log_best_se = std(losses_minimal_loss) / sqrt(length(losses_minimal_loss))
    println("Best log avg loss: $(best_log_avg_loss), avg+se: $(best_log_avg_loss+log_best_se)")
    avg_losses_filter = filter(row -> row.num_params <= n_param_minimal_loss, avg_losses)
    avg_losses_filter[:, :within_one_se] = avg_losses_filter[:, :avg_test_log_loss] .<= best_log_avg_loss + log_best_se
    println(avg_losses_filter)
    best_n_prams_se = minimum(avg_losses_filter[avg_losses_filter.within_one_se .== true, :num_params])
    println("Best n params SE: $(best_n_prams_se)")
    
    return best_n_prams_se
end

"""
    find_optimal_n_params(
        df_results::DataFrame, 
        p_value_threshold::Float64
    ) :: Int

This function determines the optimal number of parameters for a model by considering both the Wilcoxon signed-rank test and the "one standard error" (1-SE) rule.

## Parameters:
- `df_results::DataFrame`: A DataFrame containing test results, including the number of parameters and corresponding test losses.
- `p_value_threshold::Float64`: The threshold p-value used in the Wilcoxon signed-rank test to assess statistical significance.

## Returns:
- `best_n_params::Int`: The optimal number of parameters based on the combined results of the Wilcoxon test and the 1-SE rule.

## Description:
The function first groups the results by the number of parameters and calculates the average log test loss for each group. It then identifies the model with the minimal average test loss and calculates its corresponding test losses.

Two methods are applied to determine the best number of parameters:
1. **Wilcoxon signed-rank test**: This method compares models with fewer parameters to the model with the minimal average loss, selecting the best model based on statistical significance.
2. **One standard error (1-SE) rule**: This method selects the model with the fewest parameters whose average loss is within one standard error of the minimal average loss.

The function returns the smaller of the two best numbers of parameters from these methods.
"""
function find_optimal_n_params(df_results::DataFrame, p_value_threshold::Float64) :: Int
    # Group by number of parameters and calculate average test loss
    df_results[!, :log_test_loss] = log.(df_results.test_loss)
    grouped = groupby(df_results, :num_params)
    avg_losses = combine(grouped, :log_test_loss => mean => :avg_test_log_loss)
    # Sort by number of parameters
    sort!(avg_losses, :num_params)
    println("Avg LOG CV error for each n params:")
    println(avg_losses)
    # Find the row with the minimum average test loss
    
    idx_min_loss = argmin(avg_losses.avg_test_log_loss)
    n_param_minimal_loss = avg_losses[idx_min_loss, :num_params]
    losses_minimal_loss = filter(row -> row.num_params == n_param_minimal_loss, df_results).log_test_loss

    best_n_params_wilcoxon = find_best_n_params_wilcoxon(
        df_results, 
        avg_losses,
        p_value_threshold,
        n_param_minimal_loss, 
        losses_minimal_loss
    )

    best_n_params_se = find_best_n_params_within_one_se(
        losses_minimal_loss,
        avg_losses, 
        n_param_minimal_loss
    )

    return min(best_n_params_wilcoxon,best_n_params_se)
end


"""
This function iteratively fits models that are subsets of the top 10% from the previous iteration (loop over range num params), saving the best model for each
n params based on training loss, and compute LOOCV test scores for best models.
"""
function fit_rate_equation_selection_current(
        general_rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}},
        param_removal_code_names, 
        range_number_params::Tuple{Int,Int},
        forward_model_selection::Bool,
        max_zero_alpha::Int,
        n_repetiotions_opt::Int,
        maxiter_opt::Int,
        all_param_removal_codes, 
        practically_unidentifiable_params, 
        save_train_results::Bool, 
        enzyme_name::String,
        subsets_min_limit::Int, 
        subsets_max_limit::Union{Int, Nothing},
        subsets_filter_threshold::Float64,
        train_loss,
        test_loss
        )

        num_alpha_params = count(occursin.("alpha", string.([param_names...])))

        if forward_model_selection
            num_param_range = (range_number_params[2]):-1:range_number_params[1]
        elseif !forward_model_selection
            num_param_range = (range_number_params[1]):1:range_number_params[2]
        end

        starting_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
            num_param_range[1],
            all_param_removal_codes,
            param_names,
            param_removal_code_names,
            metab_names,
            practically_unidentifiable_params,
            num_alpha_params,
            max_zero_alpha,
        )

        while isempty(starting_param_removal_codes)
            num_param_range = ifelse(
                forward_model_selection,
                (num_param_range[1]-1:-1:num_param_range[end]),
                (num_param_range[1]+1:+1:num_param_range[end]),
            )
            if ifelse(
                forward_model_selection,
                (num_param_range[1] < num_param_range[end]),
                (num_param_range[1] > num_param_range[end]),
            )
                @error "Could not find any fesible equations for this enzyme within range_number_params"
            end
            println("Trying new range_number_params: $num_param_range")
            starting_param_removal_codes =
                DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes_w_num_params(
                    num_param_range[1],
                    all_param_removal_codes,
                    param_names,
                    param_removal_code_names,
                    metab_names,
                    practically_unidentifiable_params,
                    num_alpha_params,
                    max_zero_alpha,
                )
        end

        nt_param_removal_codes = starting_param_removal_codes
        nt_previous_param_removal_codes = similar(nt_param_removal_codes)
        println("About to start loop with num_params: $num_param_range")
        df_train_results = DataFrame()
        df_test_results = DataFrame()
        for num_params in num_param_range
            println("Running loop with num_params: $num_params")
    
            #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
            if num_params != num_param_range[1]
                if forward_model_selection
                    nt_param_removal_codes = forward_selection_next_param_removal_codes(
                        nt_previous_param_removal_codes,
                        metab_names,
                        practically_unidentifiable_params,
                        num_alpha_params,
                        max_zero_alpha,
                    )
                elseif !forward_model_selection
                    nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                        nt_previous_param_removal_codes,
                        metab_names,
                        practically_unidentifiable_params,
                        num_alpha_params,
                        max_zero_alpha,
                    )
                end
            end

            if isempty(nt_param_removal_codes)
                println(
                    "Stoping the search early as no feasible equations for this enzyme with $num_params parameters could be found.",
                )
                break
            end

            #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
            results_array = map(
                nt_param_removal_code -> train_rate_equation(
                    general_rate_equation,
                    data,
                    metab_names,
                    param_names;
                    n_iter = n_repetiotions_opt,
                    maxiter_opt = maxiter_opt,
                    nt_param_removal_code = nt_param_removal_code,
                    train_loss = train_loss,
                    test_loss = test_loss
                ),
                nt_param_removal_codes,
            )
            #convert results_array to DataFrame
            df_results = DataFrame(results_array)
            df_results.num_params = fill(num_params, nrow(df_results))
            df_results.nt_param_removal_codes = nt_param_removal_codes
            df_train_results = vcat(df_train_results, df_results)
    
             # Optinally consider saving results to csv file for long running calculation of cluster
            if save_train_results
                CSV.write(
                    "$(Dates.format(now(),"mmddyy"))_$(enzyme_name)_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
                    df_results,
                )
            end   

            #if all train_loss are Inf, then skip to next loop
            if all(df_results.train_loss .== Inf)
                nt_previous_param_removal_codes = [
                    NamedTuple{param_removal_code_names}(x) for
                    x in values.(df_results.nt_param_removal_codes)
                ]
                continue
            end

            #store top 10% for next loop as `previous_param_removal_codes`
            df_results = filter_and_limit_rows(df_results, :train_loss, subsets_min_limit, subsets_max_limit, subsets_filter_threshold)

            # previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
            nt_previous_param_removal_codes = [
                NamedTuple{param_removal_code_names}(x) for
                x in values.(df_results.nt_param_removal_codes)
            ]

            # save best subset for each `num_params` (afterwards loocv test loss will be calculated)
            best_nt_param_removal_code =
                df_results.nt_param_removal_codes[argmin(df_results.train_loss)]
            
            df_results = DataFrame(:num_params => [num_params], :nt_param_removal_codes => [best_nt_param_removal_code])
            df_test_results = vcat(df_test_results, df_results)
        end
       
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
                train_loss = train_loss,
                test_loss = test_loss
            ), 
        subsets_to_fit
        )
        # arrange test result ds
        result_dfs = DataFrame[]
        for (res, subset) in zip(results, subsets_to_fit)
            res_df = DataFrame([res])
            res_df[!, :nt_param_removal_codes] = [subset[1]]
            res_df[!, :num_params] = [subset[3]]
            push!(result_dfs, res_df)
        end

        df_test_results = vcat(result_dfs...)

    return (train_results = df_train_results, test_results = df_test_results, practically_unidentifiable_params = practically_unidentifiable_params)
end

"""
This function takes a given figure, splits it into training data (all other figures) and a test set (the figure itself). 
It then iteratively fits models that are subsets of the top 10% from the previous iteration (loop over range num params), 
saving the best model for each number of parameters based on training loss. 
Finally, it computes LOOCV test scores for the best models.
"""
function fit_rate_equation_selection_per_fig(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    range_number_params::Tuple{Int,Int},
    forward_model_selection::Bool,
    max_zero_alpha::Int, 
    n_repetiotions_opt::Int,
    maxiter_opt::Int,
    all_param_removal_codes,
    practically_unidentifiable_params, 
    test_fig, 
    subsets_min_limit, 
    subsets_max_limit, 
    subsets_filter_threshold
    )

    train_data = data[data.source.!=test_fig, :]
    test_data = data[data.source.==test_fig, :]

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))

    if forward_model_selection
        num_param_range = (range_number_params[2]):-1:range_number_params[1]
    elseif !forward_model_selection
        num_param_range = (range_number_params[1]):1:range_number_params[2]
    end

    starting_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
        num_param_range[1],
        all_param_removal_codes,
        param_names,
        param_removal_code_names,
        metab_names,
        practically_unidentifiable_params,
        num_alpha_params,
        max_zero_alpha,
    )

    while isempty(starting_param_removal_codes)
        num_param_range = ifelse(
            forward_model_selection,
            (num_param_range[1]-1:-1:num_param_range[end]),
            (num_param_range[1]+1:+1:num_param_range[end]),
        )
        if ifelse(
            forward_model_selection,
            (num_param_range[1] < num_param_range[end]),
            (num_param_range[1] > num_param_range[end]),
        )
            @error "Could not find any fesible equations for this enzyme within range_number_params"
        end
        println("Trying new range_number_params: $num_param_range")
        starting_param_removal_codes =
            DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes_w_num_params(
                num_param_range[1],
                all_param_removal_codes,
                param_names,
                param_removal_code_names,
                metab_names,
                practically_unidentifiable_params,
                num_alpha_params,
                max_zero_alpha,
            )
    end

    nt_param_removal_codes = starting_param_removal_codes
    nt_previous_param_removal_codes = similar(nt_param_removal_codes)
    println("Leftout figure: $(test_fig), About to start loop with num_params: $num_param_range")
    df_train_results = DataFrame()
    df_test_results = DataFrame()

    for num_params in num_param_range
        println("Running loop with num_params: $num_params")
        
        #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
        if num_params != num_param_range[1]
            if forward_model_selection
                nt_param_removal_codes = forward_selection_next_param_removal_codes(
                    nt_previous_param_removal_codes,
                    metab_names,
                    practically_unidentifiable_params,
                    num_alpha_params,
                    max_zero_alpha,
                )
            elseif !forward_model_selection
                nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                    nt_previous_param_removal_codes,
                    metab_names,
                    practically_unidentifiable_params,
                    num_alpha_params,
                    max_zero_alpha,
                )
            end
        end

        if isempty(nt_param_removal_codes)
            println(
                "Stoping the search early as no feasible equations for this enzyme with $num_params parameters could be found.",
            )
            break
        end

        #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
        results_array = pmap(
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
        df_results.dropped_fig = fill(test_fig, nrow(df_results))
        df_results.nt_param_removal_codes = nt_param_removal_codes
        df_train_results = vcat(df_train_results, df_results)
        
        #if all train_loss are Inf, then skip to next loop
        if all(df_results.train_loss .== Inf)
            nt_previous_param_removal_codes = [
                NamedTuple{param_removal_code_names}(x) for
                x in values.(df_results.nt_param_removal_codes)
            ]
            continue
        end

        #store top 10% for next loop as `previous_param_removal_codes`
        df_results = filter_and_limit_rows(df_results, :train_loss, subsets_min_limit, subsets_max_limit, subsets_filter_threshold)
        # previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
        nt_previous_param_removal_codes = [
            NamedTuple{param_removal_code_names}(x) for
            x in values.(df_results.nt_param_removal_codes)
        ]

        # Save the best subset for each num_params. afterwards, test loss will be calculated using test_fig
        idx_min_row = argmin(df_results.train_loss)
        best_nt_param_removal_code = df_results[idx_min_row, :nt_param_removal_codes]
        best_subset_rescaled_params =  df_results[idx_min_row, :params]

        df_results = DataFrame(:num_params => num_params,
        :nt_param_removal_codes => best_nt_param_removal_code,
        :params => best_subset_rescaled_params)
        
        df_test_results = vcat(df_test_results, df_results)
    end

    # calculate test loss for top subsets:
    # Prepare the data for pmap
    subsets_to_test = [(row.params, row.nt_param_removal_codes,row.num_params) for row in eachrow(df_test_results)]

    test_results = pmap(
        best_subset_params -> test_rate_equation(
            general_rate_equation,
            test_data,
            best_subset_params[1], #rescaled params 
            metab_names, 
            param_names
        ), 
        subsets_to_test
    )

    result_dfs = DataFrame[]
    for (res, subset) in zip(test_results, subsets_to_test)
        res_df = DataFrame(
            test_loss = res,          
            num_params = subset[3],          
            nt_param_removal_code =subset[2], 
            test_fig =test_fig,
            params = subset[1]           
        )
        push!(result_dfs, res_df)
    end

    df_test_results = vcat(result_dfs...)
    return (
        train_results = df_train_results, 
        test_results = df_test_results, 
        practically_unidentifiable_params = practically_unidentifiable_params
        )
end

"""
This function fits all subsets for each figure, and computes LOOCV test scores for each.
"""
function fit_rate_equation_selection_all_subsets(
    general_rate_equation::Function,
    data::DataFrame,
    all_param_removal_codes, 
    practically_unidentifiable_params, 
    max_zero_alpha,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    n_reps_opt::Int,
    maxiter_opt::Int,
    )

    # create param_subsets_per_n_params
    len_param_subset = length(first(all_param_removal_codes))
    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    n = length(param_names)

    # keep for each number of params: all the subsets with this number
    param_subsets_per_n_params = Dict{Int, Vector{NTuple{len_param_subset, Int}}}()
    # for x in Iterators.take(all_param_removal_codes, 500)
    for x in all_param_removal_codes
        n_param = n - num_alpha_params - sum(x[1:end-num_alpha_params] .> 0)
        if !haskey(param_subsets_per_n_params, n_param)
            param_subsets_per_n_params[n_param] = Vector{NTuple{len_param_subset, Int}}()
        end
        push!(param_subsets_per_n_params[n_param], x)
    end

    figs = unique(data.source)

    # Initialize an empty list for the combined results
    all_subsets_figs_to_fit = []
    lengths = []

    for (n_params, subsets) in param_subsets_per_n_params
        nt_param_removal_codes = [
            NamedTuple{param_removal_code_names}(x) for
            x in unique(subsets)
        ]
        if isempty(nt_param_removal_codes)
            filtered_nt_param_removal_codes = NamedTuple[]
        else
            filtered_nt_param_removal_codes =
                filter_param_removal_codes_to_prevent_wrong_param_combos(
                    nt_param_removal_codes,
                    metab_names,
                )
        end
        if isempty(filtered_nt_param_removal_codes)
            filtered_nt_param_removal_codes_max_alpha = NamedTuple[]
        else
            filtered_nt_param_removal_codes_max_alpha =
                filter_param_removal_codes_for_max_zero_alpha(
                    filtered_nt_param_removal_codes,
                    practically_unidentifiable_params,
                    max_zero_alpha,
                )
        end
        nt_param_subsets = unique(filtered_nt_param_removal_codes_max_alpha)
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
            n_iter = n_reps_opt,
            maxiter_opt = maxiter_opt,
            nt_param_removal_code = subset_fig_to_fit[1],
        ),
        all_subsets_figs_to_fit, 
    )

    df_results = DataFrame(results_array)
    df_results.num_params = n_params_mapping
    all_subsets  = [item[1] for item in all_subsets_figs_to_fit]
    df_results.nt_param_removal_codes = all_subsets
       
    return (
        train_test_results = df_results,
        practically_unidentifiable_params = practically_unidentifiable_params
        )
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
    train_loss = "sse",
    test_loss = "sse"
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
        train_loss = train_loss,
        test_loss = test_loss
    )
    if test_loss == "sse"
        test_loss = test_rate_equation(
            rate_equation,
            test_data,
            train_res.params,
            metab_names,
            param_names;
            test_loss = test_loss 
        )
    elseif test_loss == "likelihhod"
        test_loss = test_rate_equation(
            rate_equation,
            test_data,
            train_res.params,
            metab_names,
            param_names;
            test_loss = test_loss,
            sigma = train_res.sigma_nt.sigma,
            sigma_re = train_res.sigme_nt.sigma_re
        )
    end 
    # TODO: check if I need to return sigma here 
    return (
        dropped_fig = fig,
        train_loss = train_res.train_loss,
        test_loss = test_loss,
        params = train_res.params,
    )
end

"""Function to calculate loss for a given `rate_equation` and `nt_fitted_params` on `data` that was not used for training"""
function test_rate_equation(
    rate_equation::Function,
    data::DataFrame,
    nt_fitted_params::NamedTuple,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    test_loss = "sse",
    sigma = nothing,
    sigma_re = nothing
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
    loss_rate_equation = get_loss_function(test_loss)

    if test_loss == "sse"
        test_loss = loss_rate_equation(
            fitted_params,
            rate_equation::Function,
            rate_data_nt::NamedTuple,
            param_names::Tuple{Symbol,Vararg{Symbol}},
            fig_point_indexes::Vector{Vector{Int64}};
            rescale_params_from_0_10_scale = false,
            nt_param_removal_code = nothing,
        )
    elseif test_loss == "likelihood"
        test_loss = loss_rate_equation(
            fitted_params,
            rate_equation::Function,
            rate_data_nt::NamedTuple,
            param_names::Tuple{Symbol,Vararg{Symbol}},
            fig_point_indexes::Vector{Vector{Int64}};
            rescale_params_from_0_10_scale = false,
            nt_param_removal_code = nothing,
            sigma = sigma,
            sigma_re = sigma_re,
            is_test = true
        )
    end

    return test_loss
end

"""Generate all possibles codes for ways that params can be removed from the rate equation"""
function calculate_all_parameter_removal_codes(
    param_names::Tuple{Symbol,Vararg{Symbol}},
    practically_unidentifiable_params::Tuple{Vararg{Symbol}},
)
    feasible_param_subset_codes = ()
    for param_name in param_names
        if startswith(string(param_name), "Vmax_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2])
        elseif startswith(string(param_name), "K_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2, 3])
        elseif startswith(string(param_name), "K_") &&
               !startswith(string(param_name), "K_i") &&
               !startswith(string(param_name), "K_a") &&
               length(split(string(param_name), "_")) == 2
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(string(param_name), "K_") &&
               !startswith(string(param_name), "K_i") &&
               !startswith(string(param_name), "K_a") &&
               length(split(string(param_name), "_")) > 2
            if param_name in practically_unidentifiable_params
                feasible_param_subset_codes = (feasible_param_subset_codes..., [1])
            else
                feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2])
            end
        elseif startswith(string(param_name), "alpha")
            if param_name in practically_unidentifiable_params
                feasible_param_subset_codes = (feasible_param_subset_codes..., [1])
            else
                feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
            end
        end
    end
    # return collect(Iterators.product(feasible_param_subset_codes...))
    return Iterators.product(feasible_param_subset_codes...)
end

"""Find parameters that cannot be identified based on data and they are in front of products of metabolites concentrations that are always zero as these combinations of metabolites are absent in the data."""
function find_practically_unidentifiable_params(
    data::DataFrame,
    param_names::Tuple{Symbol,Vararg{Symbol}},
)
    practically_unidentifiable_params = []
    for param_name in param_names
        if startswith(string(param_name), "K_") &&
           !startswith(string(param_name), "K_i") &&
           !startswith(string(param_name), "K_a") &&
           length(split(string(param_name), "_")) >= 3
            if all([
                prod(row) == 0 for
                row in eachrow(data[:, Symbol.(split(string(param_name), "_")[2:end])])
            ])
                push!(practically_unidentifiable_params, param_name)
            end
        elseif startswith(string(param_name), "alpha_")
            metabs_in_param_name = Symbol.(split(string(param_name), "_")[2:3])
            if all([prod(row) == 0 for row in eachrow(data[:, metabs_in_param_name])])
                push!(practically_unidentifiable_params, param_name)
            end
        end
    end
    return Tuple(practically_unidentifiable_params)
end

"""Generate NamedTuple of codes for ways that params can be removed from the rate equation but still leave `num_params`"""
function calculate_all_parameter_removal_codes_w_num_params(
    num_params::Int,
    all_param_removal_codes,
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names::Tuple{Symbol,Vararg{Symbol}},
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    practically_unidentifiable_params::Tuple{Vararg{Symbol}},
    num_alpha_params::Int,
    max_zero_alpha::Int,
)
    codes_with_num_params = Tuple[]
    num_non_zero_in_each_code = Int[]
    for code in all_param_removal_codes
        sum_non_zero = 0
        for i = 1:(length(code)-num_alpha_params)
            if code[i] > 0
                sum_non_zero += 1
            end
        end
        push!(num_non_zero_in_each_code, sum_non_zero)
    end
    num_params_in_each_code =
        length(param_names) .- num_alpha_params .- num_non_zero_in_each_code
    for (i, code) in enumerate(all_param_removal_codes)
        if num_params_in_each_code[i] == num_params
            push!(codes_with_num_params, code)
        end
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(codes_with_num_params)]
    if isempty(nt_param_removal_codes)
        filtered_nt_param_removal_codes = NamedTuple[]
    else
        filtered_nt_param_removal_codes =
            filter_param_removal_codes_to_prevent_wrong_param_combos(
                nt_param_removal_codes,
                metab_names,
            )
    end
    if isempty(filtered_nt_param_removal_codes)
        filtered_nt_param_removal_codes_max_alpha = NamedTuple[]
    else
        filtered_nt_param_removal_codes_max_alpha =
            filter_param_removal_codes_for_max_zero_alpha(
                filtered_nt_param_removal_codes,
                practically_unidentifiable_params,
                max_zero_alpha,
            )
    end
    return unique(filtered_nt_param_removal_codes_max_alpha)
end

# """
# Function to convert parameter vector to vector where some params are equal to 0, Inf or each other based on nt_param_removal_code
# """
# function param_subset_select_denis(params, param_names, nt_param_removal_code)
#     @assert length(params) == length(param_names)
#     params_dict =
#         Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

#     for param_choice in keys(nt_param_removal_code)
#         if startswith(string(param_choice), "L") && nt_param_removal_code[param_choice] == 1
#             params_dict[:L] = 0.0
#         elseif startswith(string(param_choice), "Vmax") &&
#                nt_param_removal_code[param_choice] == 1
#             params_dict[:Vmax_i] = params_dict[:Vmax_a]
#         elseif startswith(string(param_choice), "Vmax") &&
#                nt_param_removal_code[param_choice] == 2
#             global params_dict[:Vmax_i] = 0.0
#         elseif startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 1
#             K_i = Symbol("K_i_" * string(param_choice)[8:end])
#             K_a = Symbol("K_a_" * string(param_choice)[8:end])
#             params_dict[K_i] = params_dict[K_a]
#         elseif startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 2
#             K_a = Symbol("K_a_" * string(param_choice)[8:end])
#             params_dict[K_a] = Inf
#         elseif startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 3
#             K_i = Symbol("K_i_" * string(param_choice)[8:end])
#             params_dict[K_i] = Inf
#         elseif startswith(string(param_choice), "K_") &&
#                !startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 1
#             params_dict[param_choice] = Inf
#         elseif startswith(string(param_choice), "K_") &&
#                !startswith(string(param_choice), "K_allo") &&
#                length(split(string(param_choice), "_")) > 2 &&
#                nt_param_removal_code[param_choice] == 2
#             params_dict[param_choice] =
#                 prod([
#                     params_dict[Symbol("K_" * string(metab))] for
#                     metab in split(string(param_choice), "_")[2:end]
#                 ])^(1 / (length(split(string(param_choice), "_")[2:end])))
#         elseif startswith(string(param_choice), "alpha") &&
#                nt_param_removal_code[param_choice] == 0
#             params_dict[param_choice] = 0.0
#         elseif startswith(string(param_choice), "alpha") &&
#                nt_param_removal_code[param_choice] == 1
#             params_dict[param_choice] = 1.0
#         end
#     end

#     new_params_sorted = [params_dict[param_name] for param_name in param_names]
#     return new_params_sorted
# end

"""
Function to convert parameter vector to vector where some params are equal to 0, Inf or each other based on nt_param_removal_code
"""
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
            
        elseif startswith(param_str, "Vmax")
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
Calculate `nt_param_removal_codes` with `num_params` including non-zero term combinations for codes (excluding alpha terms) in each `nt_previous_param_removal_codes` that has `num_params-1`
"""
function forward_selection_next_param_removal_codes(
    nt_previous_param_removal_codes::Vector{T} where {T<:NamedTuple},
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    practically_unidentifiable_params::Tuple{Vararg{Symbol}},
    num_alpha_params::Int,
    max_zero_alpha::Int,
)
    param_removal_code_names = keys(nt_previous_param_removal_codes[1])
    next_param_removal_codes = Vector{Vector{Int}}()
    for previous_param_removal_code in nt_previous_param_removal_codes
        i_cut_off = length(previous_param_removal_code) - num_alpha_params
        for (i, code_element) in enumerate(previous_param_removal_code)
            if i <= i_cut_off && code_element == 0
                if startswith(string(param_removal_code_names[i]), "Vmax_allo")
                    feasible_param_subset_codes = [1, 2]
                elseif startswith(string(param_removal_code_names[i]), "K_allo")
                    feasible_param_subset_codes = [1, 2, 3]
                elseif startswith(string(param_removal_code_names[i]), "K_") &&
                       !startswith(string(param_removal_code_names[i]), "K_allo") &&
                       length(split(string(param_removal_code_names[i]), "_")) == 2
                    feasible_param_subset_codes = [1]
                elseif startswith(string(param_removal_code_names[i]), "K_") &&
                       !startswith(string(param_removal_code_names[i]), "K_allo") &&
                       length(split(string(param_removal_code_names[i]), "_")) > 2
                    if param_removal_code_names[i] in practically_unidentifiable_params
                        feasible_param_subset_codes = [1]
                    else
                        feasible_param_subset_codes = [1, 2]
                    end
                end
                for code_element in feasible_param_subset_codes
                    next_param_removal_code = collect(Int, previous_param_removal_code)
                    next_param_removal_code[i] = code_element
                    push!(next_param_removal_codes, next_param_removal_code)
                end
            end
        end
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(next_param_removal_codes)]
    if isempty(nt_param_removal_codes)
        filtered_nt_param_removal_codes = NamedTuple[]
    else
        filtered_nt_param_removal_codes =
            filter_param_removal_codes_to_prevent_wrong_param_combos(
                nt_param_removal_codes,
                metab_names,
            )
    end
    if isempty(filtered_nt_param_removal_codes)
        filtered_nt_param_removal_codes_max_alpha = NamedTuple[]
    else
        filtered_nt_param_removal_codes_max_alpha =
            filter_param_removal_codes_for_max_zero_alpha(
                filtered_nt_param_removal_codes,
                practically_unidentifiable_params,
                max_zero_alpha,
            )
    end
    return unique(filtered_nt_param_removal_codes_max_alpha)
end

"""
Use `nt_previous_param_removal_codes` to calculate `nt_next_param_removal_codes` that have one additional zero elements except for for elements <= `num_alpha_params` from the end
"""
function reverse_selection_next_param_removal_codes(
    nt_previous_param_removal_codes::Vector{T} where {T<:NamedTuple},
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    practically_unidentifiable_params::Tuple{Vararg{Symbol}},
    num_alpha_params::Int,
    max_zero_alpha::Int,
)
    param_removal_code_names = keys(nt_previous_param_removal_codes[1])
    next_param_removal_codes = Vector{Vector{Int}}()
    for previous_param_removal_code in nt_previous_param_removal_codes
        i_cut_off = length(previous_param_removal_code) - num_alpha_params
        for (i, code_element) in enumerate(previous_param_removal_code)
            if i <= i_cut_off &&
               code_element != 0 &&
               param_removal_code_names[i] ∉ practically_unidentifiable_params
                next_param_removal_code = collect(Int, previous_param_removal_code)
                next_param_removal_code[i] = 0
                push!(next_param_removal_codes, next_param_removal_code)
            end
        end
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(next_param_removal_codes)]
    if isempty(nt_param_removal_codes)
        filtered_nt_param_removal_codes = NamedTuple[]
    else
        filtered_nt_param_removal_codes =
            filter_param_removal_codes_to_prevent_wrong_param_combos(
                nt_param_removal_codes,
                metab_names,
            )
    end
    if isempty(filtered_nt_param_removal_codes)
        filtered_nt_param_removal_codes_max_alpha = NamedTuple[]
    else
        filtered_nt_param_removal_codes_max_alpha =
            filter_param_removal_codes_for_max_zero_alpha(
                filtered_nt_param_removal_codes,
                practically_unidentifiable_params,
                max_zero_alpha,
            )
    end
    return unique(filtered_nt_param_removal_codes_max_alpha)
end

"""Filter removal codes to ensure that if K_S1 = Inf then all K_S1_S2 and all other K containing S1 in qssa cannot be 2, which stands for (K_S1_S2)^2 = K_S1 * K_S2"""
function filter_param_removal_codes_to_prevent_wrong_param_combos(
    nt_param_removal_codes,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
)
    # ensure that if K_S1 = Inf then all K_S1_S2 and all other K containing S1 in qssa cannot be 2, which stands for (K_S1_S2)^2 = K_S1 * K_S2
    if any([occursin("allo", string(key)) for key in keys(nt_param_removal_codes[1])])
        filtered_nt_param_removal_codes = nt_param_removal_codes
    else
        filtered_nt_param_removal_codes = NamedTuple[]
        for nt_param_removal_code in nt_param_removal_codes
            if all(
                nt_param_removal_code[Symbol("K_" * string(metab))] != 1 for
                metab in metab_names
            )
                push!(filtered_nt_param_removal_codes, nt_param_removal_code)
            else
                one_metab_codes = metab_names[findall(
                    nt_param_removal_code[Symbol("K_" * string(metab))] == 1 for
                    metab in metab_names
                )]
                if all(
                    nt_param_removal_code[param_name] != 2 for
                    param_name in keys(nt_param_removal_code) if
                    any(occursin.(string.(one_metab_codes), string(param_name)))
                )
                    push!(filtered_nt_param_removal_codes, nt_param_removal_code)
                end
            end
        end
    end
    return filtered_nt_param_removal_codes
end

"""Filter removal codes to ensure that number of alpha that are 0 is max_zero_alpha"""
function filter_param_removal_codes_for_max_zero_alpha(
    nt_param_removal_codes,
    practically_unidentifiable_params::Tuple{Vararg{Symbol}},
    max_zero_alpha::Int,
)
    practically_unidentifiable_alphas = [
        param for
        param in practically_unidentifiable_params if occursin("alpha", string(param))
    ]
    alpha_keys = [
        key for key in keys(nt_param_removal_codes[1]) if
        occursin("alpha", string(key)) && key ∉ practically_unidentifiable_alphas
    ]

    if isempty(alpha_keys)
        filtered_nt_param_removal_codes = nt_param_removal_codes
    else
        filtered_nt_param_removal_codes = NamedTuple[]
        for nt_param_removal_code in nt_param_removal_codes
            if sum([nt_param_removal_code[key] == 0 for key in alpha_keys]) <= max_zero_alpha
                push!(filtered_nt_param_removal_codes, nt_param_removal_code)
            end
        end
    end
    return filtered_nt_param_removal_codes
end

"""
This function taked the best number of parameters, trains all possible subsets of these num parameters on the entire dataset, 
and then chooses the best subset as the one with the minimal training loss.
"""
function train_and_choose_best_subset(
    general_rate_equation::Function,
    data::DataFrame,
    all_param_removal_codes, 
    best_n_params::Int, 
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    practically_unidentifiable_params, 
    max_zero_alpha,
    n_reps_opt::Int, 
    maxiter_opt::Int, 
    save_train_results::Bool, 
    enzyme_name::String;
    train_loss = "sse"
)
    num_alpha_params = count(occursin.("alpha", string.([param_names...])))

    nt_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
        best_n_params,
        all_param_removal_codes,
        param_names,
        param_removal_code_names,
        metab_names,
        practically_unidentifiable_params, 
        num_alpha_params,
        max_zero_alpha, 
    )

    results_array = pmap(
        nt_param_removal_code -> train_rate_equation(
            general_rate_equation,
            data,
            metab_names,
            param_names;
            n_iter = n_reps_opt,
            maxiter_opt = maxiter_opt,
            nt_param_removal_code = nt_param_removal_code,
            train_loss = train_loss
        ),
        nt_param_removal_codes,
    )

    #convert results_array to DataFrame
    df_results = DataFrame(results_array)
    df_results.num_params = fill(best_n_params, nrow(df_results))
    df_results.nt_param_removal_codes = nt_param_removal_codes

    # Optinally consider saving results to csv file for long running calculation of cluster
    if save_train_results
        CSV.write(
            "$(Dates.format(now(),"mmddyy"))_$(enzyme_name)_training_results_for_all_subsets_with_best_num_params_$(best_n_params).csv",
            df_results,
        )
    end  

    best_param_subset = DataFrame(df_results[argmin(df_results.train_loss),:])
    println("Best subset: $(best_param_subset.nt_param_removal_codes)")

    return best_param_subset
end

"""
    filter_and_limit_rows(
        df::DataFrame, 
        train_loss_col::Symbol, 
        min_limit::Int, 
        max_limit::Union{Int, Nothing}=nothing, 
        filter_threshold::Float64=0.1
    ) :: DataFrame

This function filters and limits the rows of a DataFrame based on a specified training loss column, ensuring the number of rows falls within a given range.

## Parameters:
- `df::DataFrame`: The input DataFrame containing the data to be filtered and limited.
- `train_loss_col::Symbol`: The column symbol representing the training loss in the DataFrame.
- `min_limit::Int`: The minimum number of rows to keep in the filtered DataFrame.
- `max_limit::Union{Int, Nothing}`: The maximum number of rows to keep in the filtered DataFrame. If `Nothing`, no maximum limit is applied.
- `filter_threshold::Float64`: The percentage threshold used to filter rows based on the minimum training loss. Rows with a training loss within this threshold of the minimum loss are kept.

## Returns:
- `filtered_df::DataFrame`: The filtered and limited DataFrame based on the specified criteria.

## Description:
This function performs the following steps:
1. **Sorting**: The DataFrame is sorted by the specified training loss column.
2. **Filtering**: Rows are filtered based on a threshold percentage of the minimum training loss. Only rows with a training loss within this percentage of the minimum are retained.
3. **Ensuring Minimum Rows**: If the filtered DataFrame has fewer rows than `min_limit`, additional rows are taken from the sorted DataFrame to meet the minimum requirement.
4. **Applying Maximum Limit**: If `max_limit` is specified and the filtered DataFrame exceeds this limit, the DataFrame is truncated to `max_limit` rows.

The function ensures that the resulting DataFrame has at least `min_limit` rows and, if applicable, no more than `max_limit` rows.
"""
function filter_and_limit_rows(
    df::DataFrame, 
    train_loss_col::Symbol, 
    min_limit::Int, 
    max_limit::Union{Int, Nothing}=nothing, 
    filter_threshold::Float64=0.1
) :: DataFrame
# function filter_and_limit_rows(df::DataFrame, train_loss_col::Symbol, min_limit::Int, max_limit::Union{Int, Nothing}=nothing, filter_threshold::Float64=0.1)
    # Check if min_limit is greater than the size of the original df
    if min_limit > nrow(df)
        println("min_limit ($(min_limit)) is greater than the number of rows in the dataframe ($(nrow(df))). Using all available rows.")
        return df
    end

    # Sort the dataframe by train loss
    sorted_df = sort(df, train_loss_col)
    
    # Get the minimum train loss
    min_loss = minimum(df[!, train_loss_col])
    
    # Calculate the threshold value based on the filter_threshold percentage
    threshold_value = min_loss * (1 + filter_threshold)
    
    # Filter rows where train loss is less than or equal to the threshold value
    filtered_df = filter(row -> row[train_loss_col] <= threshold_value, sorted_df)
    
    # Get the number of rows in the filtered dataframe
    num_rows = nrow(filtered_df)
    
    # Handle the min limit
    if num_rows < min_limit
        # If we have fewer rows than the min limit, take more rows from the sorted dataframe
        additional_rows = min(min_limit - num_rows, nrow(sorted_df) - num_rows)
        filtered_df = sorted_df[1:(num_rows + additional_rows), :]
    end
    
    # Handle the max limit only if it's not Nothing
    if !isnothing(max_limit) && nrow(filtered_df) > max_limit
        filtered_df = filtered_df[1:max_limit, :]
    end
    
    return filtered_df
end
