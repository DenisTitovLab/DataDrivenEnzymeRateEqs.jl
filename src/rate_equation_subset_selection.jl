
mytupletype = Tuple{Symbol, Vararg{Symbol}}
(:K,:L) isa mytupletype

function rate_equations_subset_selection(
    rate_equation::Function,
    metab_names::Tuple{Symbol, Vararg{Symbol}},
    param_names::Tuple{Symbol, Vararg{Symbol}},
    range_number_params::Tuple{Int, Int},
    forward_subset_selection::Bool,

)
    #generate `all_param_subsets` at `range_number_params`

    #loop over `all_param_subsets`
    for n in range_number_params[1]:range_number_params[2]
        #calculate param_subsets for `n` given `all_param_subsets` and fixed params from previous `n`

        #pmap over param_subsets for a given `n` return rescaled and nt_param_subset added

        #calculate test loss for top 10% subsets for each `n`

        #store top 10% for next loop

        #store rescaled results

    end

    #return train loss and params for all tested subsets, test loss for all tested subsets
end



function test_rate_equation(
    data,
    fitted_params,
    metab_names,
    param_names;
    nt_param_choice = zero_default_nt_param_choice,
    # rate_equation::Function = rate_equation,
)
    # Add a new column to data to assign an integer to each source/figure from publication
    data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for
            i = 1:length(unique(data.source))
        ]...,
    )

    # Convert DF to NamedTuple for better type stability / speed
    rate_data_nt = Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]])

    # Make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]

    # Check if nt_param_choice makes loss returns NaN and abort early if it does. The latter could happens due to nt_param_choice making fitted_params=Inf
    if isnan(
        loss_rate_equation(
            rate_equation::Function,
            5 .* ones(length(param_names)),
            rate_data_nt,
            param_names,
            fig_point_indexes;
            rescale_params_from_0_10_scale = true,
            nt_param_choice = nt_param_choice,
        ),
    )
        println("Loss returns NaN for this param combo")
        return Inf
    end
    # Check if one of the parameters is NaN
    if any(isnan.(fitted_params))
        println("One of the kinetic parameters used by test_loss() is NaN")
        return Inf
    end
    test_loss = loss_rate_equation(
        rate_equation::Function,
        fitted_params,
        rate_data_nt,
        param_names,
        fig_point_indexes;
        rescale_params_from_0_10_scale = false,
        nt_param_choice = nt_param_choice,
    )
    return test_loss
end
