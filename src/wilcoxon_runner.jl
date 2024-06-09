using HypothesisTests, Random, DataFrames, Statistics


function compare_models(df::DataFrame, method::Symbol)
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

function find_best_n_params(results_df::DataFrame, p_value_threshold::Float64, comparison_direction::Symbol) :: Int
    # Determine the key column based on the direction of comparison
    key_column = comparison_direction == :forward ? :model_b_num_params : :model_a_num_params
    
    # Filter results where the p-value indicates no significant difference
    no_significant_difference = filter(row -> row.p_value > p_value_threshold, results_df)
    
    # Find the optimal model depending on the comparison direction
    if nrow(no_significant_difference) > 0
        best_model = minimum(no_significant_difference[!, key_column])
    else
        # If all comparisons are significant, choose based on the safest approach to avoid overfitting
        best_model = comparison_direction == :forward ? minimum(results_df[!, :model_a_num_params]) :
                                                        maximum(results_df[!, :model_b_num_params])
    end

    return best_model
end

function find_optimal_n_params(df_results::DataFrame, p_value_threshold::Float64) :: Int
    # Group by number of parameters and calculate average test loss
    grouped = groupby(df_results, :num_params)
    avg_losses = combine(grouped, :test_loss => mean => :avg_test_loss)
    # Sort by number of parameters
    sort!(avg_losses, :num_params)
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
    
    return current_n_params
end

Random.seed!(1353)
test_results  = DataFrame(
    num_params = repeat([1, 2, 3, 4, 5], inner = 6),
    removed_fig = repeat(1:6, outer = 5),
    test_loss = rand(30)  # Random test losses
)
println(test_results)
println(compare_models(test_results, :all_pairs))
println(find_optimal_n_params(test_results, 0.05)) 



# Run the comparison
# results_df = compare_models(test_results, :forward_stepwise)
# println(results_df)

# best_n_params = find_best_n_params(results_df, 0.05, :forward)
# print(best_n_params)
# # Example data
# losses_modelA = [0.1, 0.2, 0.15, 0.18, 0.16]
# losses_modelB = [0.12, 0.19, 0.17, 0.16, 0.15]

# # Calculate differences
# differences = losses_modelA .- losses_modelB

# # Apply Wilcoxon signed-rank test
# test_result = SignedRankTest(differences)

# # Output the result
# println(test_result)
# println(pvalue(test_result))