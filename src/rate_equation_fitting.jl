#=
CODE FOR RATE EQUATION FITTING
=#
using CMAEvolutionStrategy, DataFrames
"""
Fit `rate_equation` to `data` and report the loss and best fit parameters.
"""
function train_rate_equation(
        rate_equation::Function,
        data::DataFrame,
        metab_names::Vector{Symbol},
        param_names::Vector{Symbol};
        n_iter = 20,
        nt_param_choice = nothing,
        # optimization_kwargs = optimization_kwargs
)
    # Add a new column to data to assign an integer to each source/figure from publication
    data.fig_num = vcat(
        [i * ones(Int64, count(==(unique(data.source)[i]), data.source))
         for
         i in 1:length(unique(data.source))]...,
    )

    # Convert DF to NamedTuple for better type stability / speed
    #TODO: add fig_point_indexes to rate_data_nt to avoid passing it as an argument to loss_rate_equation
    rate_data_nt = Tables.columntable(data[
        .!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]])

    # Make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]

    # Check if nt_param_choice makes loss returns NaN and abort early if it does. The latter
    # could happens due to nt_param_choice making params=Inf
    if isnan(
        loss_rate_equation(
        5 .* ones(length(param_names)),
        rate_equation,
        rate_data_nt,
        param_names,
        fig_point_indexes;
        nt_param_choice = nt_param_choice
    ),
    )
        println("Loss returns NaN for this param combo")
        return Inf, fill(NaN, length(param_names))
    end

    solns = []
    for i in 1:n_iter
        x0 = 10 .* rand(length(param_names))
        sol = try
            minimize(
                x -> loss_rate_equation(
                    x,
                    rate_equation,
                    rate_data_nt,
                    param_names,
                    fig_point_indexes;
                    nt_param_choice = nt_param_choice
                ),
                x0,
                0.01,
                lower = repeat([0.0], length(x0)),
                upper = repeat([10.0], length(x0)),
                popsize = 4 * (4 + floor(Int, 3 * log(length(x0)))),
                maxiter = 50_000,
                verbosity = 0,
                ftol = 1e-10
            )
        catch error
            # bypass rare errors (~1 in 10,000 runs) where the minimize() fails to converge with "ArgumentError: matrix contains Infs or NaNs"
            if isa(error, ArgumentError)
                println(error)
                sol = Inf
            end
        end
        push!(solns, sol)
    end
    filter!(sol -> sol != Inf ? !isinf(fbest(sol)) : !isnan(fbest(sol)), solns)
    filter!(sol -> sol != NaN ? !isinf(fbest(sol)) : !isnan(fbest(sol)), solns)
    # filter!(!isinf, solns)
    # filter!(!isnan, solns)
    if isempty(solns)
        println("All of the iterations of fits for this param combo return NaN or Inf")
        return Inf, fill(NaN, length(param_names))
    end
    index_best_sol = argmin([fbest(sol) for sol in solns])
    best_sol = try
        minimize(
            x -> loss_rate_equation(
                x,
                rate_equation::Function,
                rate_data_nt,
                param_names,
                fig_point_indexes;
                nt_param_choice = nt_param_choice
            ),
            xbest(solns[index_best_sol]),
            0.001,
            lower = repeat([0.0], length(xbest(solns[index_best_sol]))),
            upper = repeat([10.0], length(xbest(solns[index_best_sol]))),
            popsize = 4 * (4 + floor(Int, 3 * log(length(xbest(solns[index_best_sol]))))),
            maxiter = 50_000,
            verbosity = 0,
            ftol = 1e-14
        )
    catch error
        # bypass rare errors where the minimize() fails to converge with "ArgumentError: matrix contains Infs or NaNs"
        if isa(error, ArgumentError)
            println(error)
            best_sol = solns[index_best_sol]
        end
    end
    return fbest(best_sol), xbest(best_sol)
end

# Base.isnan(sol::CMAEvolutionStrategy.Optimizer) = isnan(fbest(sol))
# Base.isinf(sol::CMAEvolutionStrategy.Optimizer) = isinf(fbest(sol))

"Loss function used for fitting that calculate log of ratio of rate equation predicting of rate and rate data"
function loss_rate_equation(
        params,
        rate_equation::Function,
        rate_data_nt::NamedTuple,
        param_names,
        fig_point_indexes::Vector{Vector{Int}};
        nt_param_choice = nothing
)
    kinetic_params = param_rescaling(params, param_names)
    if !isnothing(nt_param_choice)
        kinetic_params .= param_subset_select(kinetic_params, nt_param_choice)
    end

    #precalculate log_pred_vs_data_ratios for all points as it is expensive and reuse it for weights and loss
    #convert kinetic_params to NamedTuple with field names from param_names for better type stability
    kinetic_params_nt = NamedTuple{param_names}(kinetic_params)
    log_pred_vs_data_ratios = log_ratio_predict_vs_data(
        rate_equation, rate_data_nt, kinetic_params_nt)

    #calculate figures weights and loss on per figure basis
    loss = zero(eltype(kinetic_params))
    for i in 1:maximum(rate_data_nt.fig_num)
        # calculate Vmax weights for each figure which have analytical solution as ratio of gemetric means of data vs prediction
        log_weight = mean(-log_pred_vs_data_ratios[fig_point_indexes[i]])
        loss += sum(abs2.(log_weight .+ log_pred_vs_data_ratios[fig_point_indexes[i]]))
    end
    return loss / length(rate_data_nt.Rate)
end

function log_ratio_predict_vs_data(
        rate_equation::Function,
        rate_data_nt::NamedTuple,
        kinetic_params_nt::NamedTuple;
)
    log_pred_vs_data_ratios = 10 .* ones(Float64, length(rate_data_nt.Rate))
    #TODO: maybe convert this to broacasting calculation of rate using rate_equation.(row, Ref(kinetic_params_nt))
    # and then log.(pred ./ data) instead of a loop.
    for (i, row) in enumerate(Tables.namedtupleiterator(rate_data_nt))
        log_pred_vs_data_ratios[i] = log(rate_equation(row, kinetic_params_nt) / row.Rate)
    end
    return log_pred_vs_data_ratios
end

# TODO: need to add an option to set different ranges for L, Vmax, K and alpha
"Rescaling of fitting parameters from [0, 10] scale that optimizer uses to actual values"
function param_rescaling(p, param_names)
    @assert length(p) == length(param_names)
    new_p = similar(p)
    for i in eachindex(p)
        if param_names[i] == :L
            new_p[i] = 10^(-5) * 10^(10 * p[i] / 10)
        elseif startswith(string(param_names[i]), "Vmax_")
            new_p[i] = 10^(-3) * 10^(3 * p[i] / 10)
        elseif startswith(string(param_names[i]), "K_")
            new_p[i] = 10^(-10) * 10^(13 * p[i] / 10)
        elseif startswith(string(param_names[i]), "alpha_")
            p[i] >= 5.0 ? new_p[i] = 1.0 : new_p[i] = 0.0
        else
            error("Cannot rescale unknown parameter name $(string(param_names[i])) using `param_rescaling()`")
        end
    end
    return new_p
end
