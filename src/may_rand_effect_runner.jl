
using Pkg
package_path = "/home/ec2-user/code/DataDrivenEnzymeRateEqs.jl"
Pkg.activate(package_path)
using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools
# include("rate_equation_selection.jl")

file_path = joinpath(package_path, "test/Data_for_tests/PKM2_data.csv")
data = CSV.read(file_path, DataFrame)


model_selection_method = "current_subsets_filtering" # cv_subsets_filtering / cv_all_subsets /current_subsets_filtering
n_reps_opt=1  # n repeats optimization
maxiter_opt=20 # n iteration opt algorithm
p_val_threshold = 0.05
train_loss = "likelihood"
test_loss = "likelihood"
enzyme_name = "PKM2"
forward = true


# enzyme_parameters = (; 
# substrates=[:PEP,:ADP], 
# products=[:Pyruvate, :ATP], 
# cat1=[:PEP, :Pyruvate],
# cat2 = [:ADP, :ATP], 
# reg1=[:F16BP], reg2=[:Phenylalanine], 
# Keq=20_000, oligomeric_state=4,
# rate_equation_name=:derived_rate_equation)

PKM2_enzyme = (;
    substrates=[:PEP, :ADP],
    products=[:Pyruvate, :ATP],
    regulators=[:F16BP, :Phenylalanine],
    Keq=20_000.0,
    oligomeric_state=4,
    rate_equation_name=:pkm2_rate_equation,
)
metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)
pkm2_rate_equation_no_Keq(metabs, p) = pkm2_rate_equation(metabs, p, 20000.0)

# metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
# derived_rate_equation_no_Keq(nt_metabs, nt_params) = derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(pkm2_rate_equation_no_Keq,
 data,
  metab_names, 
  param_names;
  forward_model_selection = forward,
  n_reps_opt=n_reps_opt, # n repeats optimization
  maxiter_opt=maxiter_opt,# n iteration opt algorithm
  model_selection_method = model_selection_method,
  p_val_threshold =p_val_threshold,
  save_train_results = false,
  enzyme_name = enzyme_name, 
  train_loss = train_loss,
  test_loss = test_loss
  )

