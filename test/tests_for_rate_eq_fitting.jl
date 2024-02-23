using EnzymeFitting, Test, BenchmarkTools, CMAEvolutionStrategy, DataFrames, CSV

#Load and process data
PKM2_data_for_fit = CSV.read("Data_for_tests/PKM2_data.csv", DataFrame)

#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source .= PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig
