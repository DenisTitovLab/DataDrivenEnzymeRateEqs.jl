# DataDrivenEnzymeRateEqs

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://denistitovlab.github.io/DataDrivenEnzymeRateEqs.jl/dev/)
[![Build Status](https://github.com/denistitovlab/DataDrivenEnzymeRateEqs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/denistitovlab/DataDrivenEnzymeRateEqs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/DenisTitovLab/DataDrivenEnzymeRateEqs.jl/graph/badge.svg?token=3RBUX547CQ)](https://codecov.io/gh/DenisTitovLab/DataDrivenEnzymeRateEqs.jl)

## Package status

This package is under rapid development. Expect frequent breaking changes.

## Description

DataDrivenEnzymeRateEqs.jl is a Julia package for identifying metabolic enzyme rate equations based on *in vitro* kinetic data.  

The user provides the following info about enzyme:
- names of substrates, products, and regulators of the enzyme
- number of metabolite binding sites and list of metabolites bound to each site
- oligomeric state of the enzyme
- *in vitro* kinetic data

DataDrivenEnzymeRateEqs.jl outputs the following info about enzyme:
- enzyme rate equation based on the Monod-Wyman-Changeux model that is most consistent with the user provided data (i.e., rate equation that is best at predicting *in vitro* kinetic data not used for fitting measured with lowest leave one out crossvalidation score)
- kinetic parameters and their uncertainty (*work in progress*)

## Example

#### Make up an enzyme rate equation and generate noisy data
```julia
using DataDrivenEnzymeRateEqs, DataFrames

data_gen_rate_equation_Keq = 1.0
data_gen_rate_equation(metabs, params) = params.Vmax * (metabs.S / params.K_S - (1 / data_gen_rate_equation_Keq) * metabs.P / params.K_P) / (1 + metabs.S / params.K_S + metabs.P / params.K_P)
param_names = (:Vmax, :K_S, :K_P)
metab_names = (:S, :P)
params = (Vmax=10.0, K_S=1e-3, K_P=5e-3)
```
#### Generate noisy data using `data_gen_rate_equation`
```julia
num_datapoints = 10
num_figures = 4
S_concs = Float64[]
P_concs = Float64[]
sources = String[]

for i in 1:num_figures
    if i < num_figures ÷ 2
        for S in range(0, rand(1:10) * params.K_S, rand(num_datapoints ÷ 2 : num_datapoints * 2))
            push!(S_concs, S)
            push!(P_concs, 0.0)
            push!(sources, "Figure$i")
        end
    else
        for P in range(0, rand(1:10) * params.K_P, rand(num_datapoints ÷ 2 : num_datapoints * 2))
            push!(S_concs, 0.0)
            push!(P_concs, P)
            push!(sources, "Figure$i")
        end
    end
end
data = DataFrame(S=S_concs, P=P_concs, source=sources)
noise_sd = 0.2
data.Rate = [data_gen_rate_equation(row, params) * (1 + noise_sd * randn()) for row in eachrow(data)]
```
#### Use the `data` above to identify the rate equation and check that it is same as `data_gen_rate_equation`  
```julia
#Define enzyme
enzyme_parameters = (; substrates=[:S,], products=[:P], regulators=[], Keq=1.0, oligomeric_state=1, rate_equation_name=:mwc_derived_rate_equation)
metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)

#Find the best rate equation
mwc_derived_rate_equation_no_Keq(nt_metabs, nt_params) = mwc_derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(mwc_derived_rate_equation_no_Keq, data, metab_names, param_names)

#Display best equation with 3 parameters. Compare sym_rate_equation with data_gen_rate_equation with Vmax=1
nt_param_removal_code = filter(x -> x.num_params .== 3, selection_result.test_results).nt_param_removal_codes[1]
sym_rate_equation = display_rate_equation(mwc_derived_rate_equation, metab_names, param_names; nt_param_removal_code=nt_param_removal_code)
```

## Documentation

Check the [Documentation](https://denistitovlab.github.io/DataDrivenEnzymeRateEqs.jl/dev/) for more info

## Plans for the future
- add plotting functions
- add function for bootstrapping for calculation of kinetic parameter uncertainty based on log-normal distribution
- update rate equations based on Quasy-Steady-State approximation to choose actual kinetic constants instead of aggregate constants like K_S1_P1_P2


