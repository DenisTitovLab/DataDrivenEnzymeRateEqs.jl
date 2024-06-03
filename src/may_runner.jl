using Pkg
package_path = "/home/ec2-user/code/DataDrivenEnzymeRateEqs.jl"
Pkg.activate(package_path)

using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools
include("rate_equation_selection.jl")

file_path = joinpath(package_path, "test/Data_for_tests/PKM2_data.csv")
data = CSV.read(file_path, DataFrame)

enzyme_parameters = (; 
substrates=[:PEP,:ADP], 
products=[:Pyruvate, :ATP], 
cat1=[:PEP, :Pyruvate],
cat2 = [:ADP, :ATP], 
reg1=[:F16BP], reg2=[:Phenylalanine], 
Keq=20_000, oligomeric_state=4,
rate_equation_name=:derived_rate_equation)
metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
derived_rate_equation_no_Keq(nt_metabs, nt_params) = derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(derived_rate_equation_no_Keq,
 data,
  metab_names, 
  param_names, 
  (7, 15), 
  true,
  1, # n repeats optimization
  100 # n iteration opt algorithm
  )