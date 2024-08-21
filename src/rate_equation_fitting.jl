#=
CODE FOR RATE EQUATION FITTING
=#
using CMAEvolutionStrategy, DataFrames, Statistics, MixedModels, Distributions, LinearAlgebra

#TODO; add optimization_kwargs and use Optimization.jl
#TODO: add an option to set different ranges for L, Vmax, K and alpha
#TODO: add an option to fit real Vmax values instead of fixing Vmax=1.0
"""
    fit_rate_equation(
        rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol, Vararg{Symbol}},
        param_names::Tuple{Symbol, Vararg{Symbol}};
        n_iter = 20
)

Fit `rate_equation` to `data` and return loss and best fit parameters.

# Arguments
- `rate_equation::Function`: Function that takes a NamedTuple of metabolite concentrations (with `metab_names` keys) and parameters (with `param_names` keys) and returns an enzyme rate.
- `data::DataFrame`: DataFrame containing the data with column `Rate` and columns for each `metab_names` where each row is one measurement. It also needs to have a column `source` that contains a string that identifies the source of the data. This is used to calculate the weights for each figure in the publication.
- `metab_names::Tuple{Symbol, Vararg{Symbol}}`: Tuple of metabolite names that correspond to the metabolites of `rate_equation` and column names in `data`.
- `param_names::Tuple{Symbol, Vararg{Symbol}}`: Tuple of parameter names that correspond to the parameters of `rate_equation`.
- `n_iter::Int`: Number of iterations to run the fitting process.

# Returns
- `loss::Float64`: Loss of the best fit.
- `params::NamedTuple`: Best fit parameters with `param_names` keys

# Example
```julia
using DataFrames
data = DataFrame(
    Rate = [1.0, 2.0, 3.0],
    A = [1.0, 2.0, 3.0],
    source = ["Figure 1", "Figure 1", "Figure 2"]
)
rate_equation(metabs, params) = params.Vmax * metabs.S / (1 + metabs.S / params.K_S)
fit_rate_equation(rate_equation, data, (:A,), (:Vmax, :K_S))
```
"""
function fit_rate_equation(
    rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    n_iter = 20,
    maxiter_opt = 50_000,
 )
    train_results = train_rate_equation(
        rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}};
        n_iter = n_iter,
        maxiter_opt = maxiter_opt,
        nt_param_removal_code = nothing,
    )
    # rescaled_params = param_rescaling(train_results[2], param_names)
    # return (loss = train_results[1], params = NamedTuple{param_names}(rescaled_params))
    return (train_loss = train_results.train_loss, params = train_results.params)
end

function train_rate_equation(
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
    sigma_nt = nothing
    # Add a new column to data to assign an integer to each source/figure from publication
    data.fig_num = vcat(
        [
            i * ones(
                Int64,
                count(==(unique(data.source)[i]), data.source),
            ) for i = 1:length(unique(data.source))
        ]...,
    )
    # Add a column containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]
    # Convert DF to NamedTuple for better type stability / speed
    rate_data_nt = Tables.columntable(data)

    loss_rate_equation = get_loss_function(train_loss)
    # Check if nt_param_removal_code makes loss returns NaN and abort early if it does. The latter
    # could happens due to nt_param_removal_code making params=Inf
    if isnan(
        loss_rate_equation(
            5 .* ones(length(param_names)),
            rate_equation,
            rate_data_nt,
            param_names,
            fig_point_indexes;
            rescale_params_from_0_10_scale = true,
            nt_param_removal_code = nt_param_removal_code,
        ),
    )
        # @warn "Loss returns NaN for this param combo in train_rate_equation() before minimization"
        return (
            train_loss = Inf,
            params = NamedTuple{param_names}(Tuple(fill(NaN, length(param_names)))),
            sigma_nt = sigma_nt
        )
    end

    solns = []
    for i = 1:n_iter
        x0 = 10 .* rand(length(param_names))
        sol = try
            minimize(
                x -> loss_rate_equation(
                    x,
                    rate_equation,
                    rate_data_nt,
                    param_names,
                    fig_point_indexes;
                    rescale_params_from_0_10_scale = true,
                    nt_param_removal_code = nt_param_removal_code,
                ),
                x0,
                0.01,
                lower = repeat([0.0], length(x0)),
                upper = repeat([10.0], length(x0)),
                popsize = 4 * (4 + floor(Int, 3 * log(length(x0)))),
                maxiter = maxiter_opt,
                verbosity = 0,
                ftol = 1e-10,
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
    filter!(sol -> sol != Inf, solns)
    filter!(sol -> fbest(sol) != NaN, solns)

    if isempty(solns)
        @warn "All of the iterations of fits for this param combo return NaN or Inf in train_rate_equation() before minimization"
        return (
            train_loss = Inf,
            params = NamedTuple{param_names}(Tuple(fill(NaN, length(param_names)))),
            sigma_nt = sigma_nt
        )
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
                rescale_params_from_0_10_scale = true,
                nt_param_removal_code = nt_param_removal_code,
            ),
            xbest(solns[index_best_sol]),
            0.001,
            lower = repeat([0.0], length(xbest(solns[index_best_sol]))),
            upper = repeat([10.0], length(xbest(solns[index_best_sol]))),
            popsize = 4 * (4 + floor(Int, 3 * log(length(xbest(solns[index_best_sol]))))),
            maxiter = maxiter_opt,
            verbosity = 0,
            ftol = 1e-14,
        )
    catch error
        # bypass rare errors where the minimize() fails to converge with "ArgumentError: matrix contains Infs or NaNs"
        if isa(error, ArgumentError)
            println(error)
            best_sol = solns[index_best_sol]
        end
    end
    rescaled_params = param_rescaling(xbest(best_sol), param_names)
    if !isnothing(nt_param_removal_code)
        rescaled_params =
            param_subset_select(rescaled_params, param_names, nt_param_removal_code)
    end

    if test_loss == "likelihood"
        sigma_nt = get_sigmas_of_bes_sol(
            rescaled_params,
            rate_equation, 
            rate_data_nt,
            param_names
        )
    end
    return (train_loss = fbest(best_sol), params = NamedTuple{param_names}(rescaled_params), sigma_nt = sigma_nt)
end

"Loss function used for fitting that calculate log of ratio of rate equation predicting of rate and rate data"
# function loss_rate_equation(
#     params,
#     rate_equation::Function,
#     rate_data_nt::NamedTuple,
#     param_names::Tuple{Symbol,Vararg{Symbol}},
#     fig_point_indexes::Vector{Vector{Int}};
#     rescale_params_from_0_10_scale = true,
#     nt_param_removal_code = nothing,
#  )
#     if rescale_params_from_0_10_scale
#         kinetic_params = param_rescaling(params, param_names)
#     else
#         !rescale_params_from_0_10_scale
#         kinetic_params = params
#     end
#     if !isnothing(nt_param_removal_code)
#         kinetic_params .=
#             param_subset_select(kinetic_params, param_names, nt_param_removal_code)
#     end

#     #precalculate log_pred_vs_data_ratios for all points as it is expensive and reuse it for weights and loss
#     kinetic_params_nt =
#         NamedTuple{param_names}(ntuple(i -> kinetic_params[i], Val(length(kinetic_params))))
#     log_pred_vs_data_ratios =
#         log_ratio_predict_vs_data(rate_equation, rate_data_nt, kinetic_params_nt)

#     #calculate figures weights and loss on per figure basis
#     loss = zero(eltype(log_pred_vs_data_ratios))
#     for i = 1:maximum(rate_data_nt.fig_num)
#         log_fig_weight = zero(eltype(log_pred_vs_data_ratios))
#         counter = 0
#         for j in fig_point_indexes[i]
#             log_fig_weight += log_pred_vs_data_ratios[j]
#             counter += 1
#         end
#         log_fig_weight /= counter
#         for j in fig_point_indexes[i]
#             loss += abs2(log_fig_weight - log_pred_vs_data_ratios[j])
#         end
#     end
#     return loss / length(rate_data_nt.Rate)
# end


"Loss function used for fitting that calculate log of ratio of rate equation predicting of rate and rate data"
function loss_sse_rate_equation(
    params,
    rate_equation::Function,
    rate_data_nt::NamedTuple,
    param_names::Tuple{Symbol,Vararg{Symbol}},
    fig_point_indexes::Vector{Vector{Int}};
    rescale_params_from_0_10_scale = true,
    nt_param_removal_code = nothing,
 )
    if rescale_params_from_0_10_scale
        kinetic_params = param_rescaling(params, param_names)
    else
        !rescale_params_from_0_10_scale
        kinetic_params = params
    end
    if !isnothing(nt_param_removal_code)
        kinetic_params .=
            param_subset_select(kinetic_params, param_names, nt_param_removal_code)
    end

    #precalculate log_pred_vs_data_ratios for all points as it is expensive and reuse it for weights and loss
    kinetic_params_nt =
        NamedTuple{param_names}(ntuple(i -> kinetic_params[i], Val(length(kinetic_params))))
    log_pred_vs_data_ratios =
        log_ratio_predict_vs_data(rate_equation, rate_data_nt, kinetic_params_nt)

    #calculate figures weights and loss on per figure basis
    loss = zero(eltype(log_pred_vs_data_ratios))
    for i = 1:maximum(rate_data_nt.fig_num)
        log_fig_weight = zero(eltype(log_pred_vs_data_ratios))
        counter = 0
        for j in fig_point_indexes[i]
            log_fig_weight += log_pred_vs_data_ratios[j]
            counter += 1
        end
        log_fig_weight /= counter
        for j in fig_point_indexes[i]
            loss += abs2(log_fig_weight - log_pred_vs_data_ratios[j])
        end
    end
    return loss / length(rate_data_nt.Rate)
end

function log_likelihood_calc(sigma,sigma_re, data)
    
    log_likelihood_total = 0.0

    for j in 1:maximum(data.figure)
        figure = data.figure[j]
        rates_for_figure = data.log_ratios[data.figure .== j]

        n_j = length(rates_for_figure)

        # Construct the covariance matrix for the figure
        Sigma_j = (sigma^2) * I(n_j) + (sigma_re^2) *  ones(n_j, n_j)

        # Calculate the log-likelihood for the figure using the multivariate normal distribution
        dist_j = MvNormal(zeros(n_j), Sigma_j) 
        log_likelihood_total += logpdf(dist_j, rates_for_figure)
    end

    return log_likelihood_total  
end

function get_sigmas_of_bes_sol(
    kinetic_params,
    rate_equation::Function,
    rate_data_nt::NamedTuple,
    param_names::Tuple{Symbol,Vararg{Symbol}},
 )
    #precalculate log_pred_vs_data_ratios for all points as it is expensive and reuse it for weights and loss
    kinetic_params_nt =
    NamedTuple{param_names}(ntuple(i -> kinetic_params[i], Val(length(kinetic_params))))
    log_pred_vs_data_ratios =
        log_ratio_predict_vs_data(rate_equation, rate_data_nt, kinetic_params_nt)

    df = DataFrame(figure=rate_data_nt.fig_num, log_ratios=log_pred_vs_data_ratios)
    if length(unique(df.log_ratios))==1
        return (sigma = nothing, sigma_re = nothing)
    end
    
    model_formula = @formula(-log_ratios ~ 0+ (1|figure)) # - since it suppose to be log(actual/pred)
    model = fit!(LinearMixedModel(model_formula, df))
    fitted_sigma = model.σ
    fitted_sigma_rand = model.σs.figure[1]
    return (sigma = fitted_sigma, sigma_re = fitted_sigma_rand)
end

function loss_likelihood_rate_equation(
    params,
    rate_equation::Function,
    rate_data_nt::NamedTuple,
    param_names::Tuple{Symbol,Vararg{Symbol}},
    fig_point_indexes::Vector{Vector{Int}};
    rescale_params_from_0_10_scale = true,
    nt_param_removal_code = nothing,
    sigma = nothing,
    sigma_re = nothing,
    is_test = false
 )
    if rescale_params_from_0_10_scale
        kinetic_params = param_rescaling(params, param_names)
    else
        !rescale_params_from_0_10_scale
        kinetic_params = params
    end
    if !isnothing(nt_param_removal_code)
        kinetic_params .=
            param_subset_select(kinetic_params, param_names, nt_param_removal_code)
    end

    #precalculate log_pred_vs_data_ratios for all points as it is expensive and reuse it for weights and loss
    kinetic_params_nt =
        NamedTuple{param_names}(ntuple(i -> kinetic_params[i], Val(length(kinetic_params))))
    log_pred_vs_data_ratios =
        log_ratio_predict_vs_data(rate_equation, rate_data_nt, kinetic_params_nt)

    df = DataFrame(figure=rate_data_nt.fig_num, log_ratios=log_pred_vs_data_ratios)
    
    if is_test == false
        model_formula = @formula(-log_ratios ~ 0+ (1|figure)) # - since it suppose to be log(actual/pred)

        if length(unique(df.log_ratios))==1
            return Inf
        else
            model = fit!(LinearMixedModel(model_formula, df))
            return -loglikelihood(model)
        end
    elseif is_test == true
        return log_likelihood_calc(sigma,sigma_re, df)
    end
end

function get_loss_function(loss_choice = "sse")
    if loss_choice == "sse"
        return loss_sse_rate_equation
    elseif loss_choice == "likelihood"
        return loss_likelihood_rate_equation
    end
end

function log_ratio_predict_vs_data(
    rate_equation::Function,
    rate_data_nt::NamedTuple,
    kinetic_params_nt::NamedTuple;
)
    log_pred_vs_data_ratios = ones(Float64, length(rate_data_nt.Rate))
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
        elseif startswith(string(param_names[i]), "Vmax")
            new_p[i] = 10^(-3) * 10^(3 * p[i] / 10)
        elseif startswith(string(param_names[i]), "K_")
            new_p[i] = 10^(-10) * 10^(13 * p[i] / 10)
        elseif startswith(string(param_names[i]), "alpha_")
            p[i] >= 5.0 ? new_p[i] = 1.0 : new_p[i] = 0.0
        else
            error(
                "Cannot rescale unknown parameter name $(string(param_names[i])) using `param_rescaling()`",
            )
        end
    end
    return new_p
end
