using CairoMakie, CSV, DataFrames, Printf, Statistics

"""
    plot_fit_on_data(
        rate_equation,
        nt_kinetic_params,
        metab_names,
        data;
        enzyme_name::String = "",
        num_col::Int = 5,
        scaler = 4.0
    )

Plot the fit of `rate_equation` with `nt_kinetic_params` on `data`.

# Arguments
- `rate_equation::Function`: Function that calculates the rate of the enzyme.
- `nt_kinetic_params::NamedTuple`: NamedTuple containing the kinetic parameters of the enzyme.
- `metab_names::Vector{Symbol}`: Vector of metabolite names.
- `data::DataFrame`: DataFrame containing the data to plot.

# Optional Arguments
- `enzyme_name::String = ""`: Name of the enzyme to be displayed on the plot.
- `num_col::Int = 5`: Number of columns in the plot.
- `scaler::Float = 4.0`: Scaling factor for the plot.

# Returns
- `fig::Figure`: Figure containing the plot.
"""
function plot_fit_on_data(
        rate_equation,
        nt_kinetic_params,
        metab_names,
        data;
        enzyme_name::String = "",
        num_col::Int = 5,
        scaler = 4.0
)
    fontsize = scaler * 5
    markersize = scaler * 3
    linewidth = scaler * 1
    set_theme!(
        Theme(
        fontsize = fontsize,
        Axis = (
            titlefont = :regular,
            xticksize = scaler * 1,
            yticksize = scaler * 1,
            xlabelfont = :bold,
            ylabelfont = :bold,
            yticklabelpad = scaler * 1,
            ylabelpadding = scaler * 3
        ),
        Legend = (titlefont = :regular,)
    ),
    )
    num_rows = ceil(Int, length(unique(data.source)) / num_col)
    size_inches = scaler .* (7, num_rows)
    size_pt = 72 .* size_inches

    fig = Figure(size = size_pt)
    for (i, figure) in enumerate(unique(data.source))
        grid_layout = fig[i % num_col == 0 ? i ÷ num_col : 1 + i ÷ num_col, i % num_col == 0 ? num_col : i % num_col] = GridLayout()

        fig_data_points = data[data.source .== figure, :]
        #assert that all data points in the figure have the same X_axis_label
        @assert all(fig_data_points[:, :X_axis_label] .==
                    fig_data_points[:, :X_axis_label][1])

        x_axis_metabolite = fig_data_points[:, :X_axis_label][1]
        other_metabolites = setdiff(metab_names, [Symbol(x_axis_metabolite)])
        fig_Vmax_for_fit = exp(mean(-DataDrivenEnzymeRateEqs.log_ratio_predict_vs_data(
            rate_equation,
            Tables.columntable(fig_data_points),
            nt_kinetic_params
        )))
        #make Axis for plotting
        ax = Axis(
            grid_layout[1, 1],
            title = replace(figure, "_" => " "),
            xticks = LinearTicks(3),
            yticks = LinearTicks(3),
            limits = begin
                maximum(fig_data_points.Rate) > 0.0 ?
                (nothing, (-0.05 * maximum(fig_data_points.Rate), nothing)) :
                (nothing, (nothing, -0.05 * minimum(fig_data_points.Rate)))
            end,
            ytickformat = ys -> ["$(round(x/maximum(abs.(ys)), sigdigits=1))" for x in ys],
            xtickformat = xs -> ["$(round(x/1e-3, sigdigits=2))" for x in xs],
            yticklabelrotation = π / 2,
            xlabelpadding = scaler * 0,
            ylabelpadding = scaler * 0,
            xticklabelpad = 0,
            yticklabelpad = 0,
            xlabel = crop_to_vowel_end(string(x_axis_metabolite), 5) * ", mM",
            ylabel = begin
                i % num_col == 1 ? "$(enzyme_name) Rate" : ""
            end
        )
        #find unque and changing metabolites concs in fig_data_points except x_axis_metabolite
        const_metab_concs = [metab
                             for metab in other_metabolites
                             if length(unique(fig_data_points[!, metab])) == 1 &&
                                unique(fig_data_points[!, metab])[1] != 0.0]
        changing_metab_concs = [metab
                                for metab in other_metabolites
                                if length(unique(fig_data_points[!, metab])) > 1]
        #loop over unique values of metabolites other than x_axis_metabolite and plot data and fit
        for (i, unique_metab) in enumerate(eachrow(unique(fig_data_points[
            :, other_metabolites])))

            #get data points for the unique_metab
            data_for_scatter = filter(
                row -> all([row[metab] == unique_metab[metab]
                            for metab in other_metabolites]),
                fig_data_points)
            #get the fit for the unique_metab using nt_kinetic_params
            function fit(x)
                fig_Vmax_for_fit * rate_equation(
                    merge((; Symbol(x_axis_metabolite) => x,), NamedTuple(unique_metab)),
                    nt_kinetic_params)
            end
            #make a label of changing metabolites conc
            metab_conc_label = join(
                [unit_conc_convertion(unique_metab[metab])
                 for metab in changing_metab_concs],
                ", "
            )
            metab_conc_label = fallback_if_blank(metab_conc_label; fallback = "No varying metabolites")
            #plot data and fit
            scatter!(ax, data_for_scatter[!, x_axis_metabolite],
                data_for_scatter.Rate, markersize = markersize, label = metab_conc_label)
            lines!(ax, 0 .. maximum(fig_data_points[!, x_axis_metabolite]),
                x -> fit(x), linewidth = linewidth, label = metab_conc_label)
        end
        #make legend title
        legend_title = begin
            str = ""
            for (i, metab) in enumerate(const_metab_concs)
                str = str * crop_to_vowel_end(string(metab), 5) * "=" *
                      unit_conc_convertion(unique(fig_data_points[!, metab])[1])
                str = str * "\n"
            end
            for (i, metab) in enumerate(changing_metab_concs)
                str = str * crop_to_vowel_end(string(metab), 5) * ""
                if i < length(changing_metab_concs)
                    str = str * ", "
                end
            end
            str = fallback_if_blank(str; fallback = "No varying metabolites")
            str
        end
        leg = Legend(grid_layout[1, 2],
            ax,
            legend_title,
            merge = true,
            unique = true,
            labelsize = fontsize,
            titlesize = fontsize,
            patchsize = scaler .* (6.0f0, 5.0f0),
            patchlabelgap = scaler * 2,
            padding = scaler .* (0.5f0, 0.0f0, 0.0f0, 0.0f0),
            framevisible = false,
            rowgap = 0,
            titlegap = 0,
            titlevalign = :top,
            titlehalign = :left,
            valign = :top,
            halign = :left
        )
        colgap!(grid_layout, 1)
    end
    colgap!(fig.layout, scaler * 5)
    rowgap!(fig.layout, scaler * 5)
    fig
end

function crop_to_vowel_end(str, max_length::Int = 5)
    vowels = Set(['a', 'e', 'i', 'o', 'u', 'y', 'A', 'E', 'I', 'O', 'U', 'Y'])
    if length(str) > max_length
        str = str[1:max_length]
    end
    while length(str) > 0 && (last(str) in vowels)
        str = str[1:(end - 1)]
    end
    return str
end

function unit_conc_convertion(conc::Number)
    if conc == 0.0
        str = "––––"
    elseif conc >= 0.1
        str = @sprintf("%.2f", round(conc, sigdigits = 2)) * "M"
    elseif conc >= 10e-3
        str = @sprintf("%.0f", round(conc, sigdigits = 2)/1e-3) * "mM"
    elseif conc >= 1e-3
        str = @sprintf("%.1f", round(conc, sigdigits = 2)/1e-3) * "mM"
    elseif conc >= 100e-6
        str = @sprintf("%.0f", round(conc, sigdigits = 2)/1e-6) * "µM"
    elseif conc >= 10e-6
        str = @sprintf("%.0f", round(conc, sigdigits = 2)/1e-6) * "µM"
    elseif conc >= 1e-6
        str = @sprintf("%.1f", round(conc, sigdigits = 2)/1e-6) * "µM"
    elseif conc >= 100e-9
        str = @sprintf("%.0f", round(conc, sigdigits = 2)/1e-9) * "nM"
    elseif conc >= 10e-9
        str = @sprintf("%.0f", round(conc, sigdigits = 2)/1e-9) * "nM"
    elseif conc >= 1e-9
        str = @sprintf("%.1f", round(conc, sigdigits = 2)/1e-9) * "nM"
    elseif conc >= 0.1e-9
        str = @sprintf("%.0f", round(conc, sigdigits = 2)/1e-12) * "pM"
    elseif conc == 0.0
        str = "0.0     "
    else
        str = "Number Out of Scale"
    end
    return str
end

@inline function fallback_if_blank(text::AbstractString; fallback::AbstractString = "No varying metabolites")
    isempty(strip(text)) && return fallback
    return text
end
