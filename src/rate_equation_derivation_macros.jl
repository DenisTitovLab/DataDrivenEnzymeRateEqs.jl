#=
CODE FOR RATE EQUATION DERIVATION
=#

function general_mwc_rate_equation(
        S1::T,
        S2::T,
        P1::T,
        P2::T,
        R1_reg1::T,
        R2_reg1::T,
        R3_reg1::T,
        R1_reg2::T,
        R2_reg2::T,
        R3_reg2::T,
        R1_reg3::T,
        R2_reg3::T,
        R3_reg3::T,
        L::T,
        Vmax_a::T,
        Vmax_i::T,
        K_a_S1_cat::T,
        K_i_S1_cat::T,
        K_a_S2_cat::T,
        K_i_S2_cat::T,
        K_a_P1_cat::T,
        K_i_P1_cat::T,
        K_a_P2_cat::T,
        K_i_P2_cat::T,
        K_a_R1_reg1::T,
        K_i_R1_reg1::T,
        K_a_R2_reg1::T,
        K_i_R2_reg1::T,
        K_a_R3_reg1::T,
        K_i_R3_reg1::T,
        K_a_R1_reg2::T,
        K_i_R1_reg2::T,
        K_a_R2_reg2::T,
        K_i_R2_reg2::T,
        K_a_R3_reg2::T,
        K_i_R3_reg2::T,
        K_a_R1_reg3::T,
        K_i_R1_reg3::T,
        K_a_R2_reg3::T,
        K_i_R2_reg3::T,
        K_a_R3_reg3::T,
        K_i_R3_reg3::T,
        alpha_S1_P2::T,
        alpha_S2_P1::T,
        Keq::T
) where {T <: Float64}
    Vmax_a = 1.0

    Vmax_a_rev = ifelse(
        (K_a_P2_cat * K_a_P1_cat) != Inf,
        Vmax_a * K_a_P2_cat * K_a_P1_cat / (Keq * K_a_S1_cat * K_a_S2_cat),
        0.0
    )
    Vmax_i_rev = ifelse(
        (K_i_P2_cat * K_i_P1_cat) != Inf,
        Vmax_i * K_i_P2_cat * K_i_P1_cat / (Keq * K_i_S1_cat * K_i_S2_cat),
        0.0
    )

    Z_a_cat = (
        1 +
        (S1 / K_a_S1_cat) +
        (S2 / K_a_S2_cat) +
        (P1 / K_a_P1_cat) +
        (P2 / K_a_P2_cat) +
        (S1 / K_a_S1_cat) * (S2 / K_a_S2_cat) +
        (P1 / K_a_P1_cat) * (P2 / K_a_P2_cat) +
        alpha_S1_P2 * (S1 / K_a_S1_cat) * (P2 / K_a_P2_cat) +
        alpha_S2_P1 * (P1 / K_a_P1_cat) * (S2 / K_a_S2_cat)
    )
    Z_i_cat = (
        1 +
        (S1 / K_i_S1_cat) +
        (S2 / K_i_S2_cat) +
        (P1 / K_i_P1_cat) +
        (P2 / K_i_P2_cat) +
        (S1 / K_i_S1_cat) * (S2 / K_i_S2_cat) +
        (P1 / K_i_P1_cat) * (P2 / K_i_P2_cat) +
        alpha_S1_P2 * (S1 / K_i_S1_cat) * (P2 / K_i_P2_cat) +
        alpha_S2_P1 * (P1 / K_i_P1_cat) * (S2 / K_i_S2_cat)
    )
    Z_a_reg = (
        (1 + R1_reg1 / K_a_R1_reg1 + R2_reg1 / K_a_R2_reg1 + R3_reg1 / K_a_R3_reg1) *
        (1 + R1_reg2 / K_a_R1_reg2 + R2_reg2 / K_a_R2_reg2 + R3_reg2 / K_a_R3_reg2) *
        (1 + R1_reg3 / K_a_R1_reg3 + R2_reg3 / K_a_R2_reg3 + R3_reg3 / K_a_R3_reg3)
    )
    Z_i_reg = (
        (1 + R1_reg1 / K_i_R1_reg1 + R2_reg1 / K_i_R2_reg1 + R3_reg1 / K_i_R3_reg1) *
        (1 + R1_reg2 / K_i_R1_reg2 + R2_reg2 / K_i_R2_reg2 + R3_reg2 / K_i_R3_reg2) *
        (1 + R1_reg3 / K_i_R1_reg3 + R2_reg3 / K_i_R2_reg3 + R3_reg3 / K_i_R3_reg3)
    )

    Rate = (
        (
            Vmax_a * (S1 / K_a_S1_cat) * (S2 / K_a_S2_cat) -
            Vmax_a_rev * (P1 / K_a_P1_cat) * (P2 / K_a_P2_cat)
        ) *
        (Z_a_cat^3) *
        (Z_a_reg^4) +
        L *
        (
            Vmax_i * (S1 / K_i_S1_cat) * (S2 / K_i_S2_cat) -
            Vmax_i_rev * (P1 / K_i_P1_cat) * (P2 / K_i_P2_cat)
        ) *
        (Z_i_cat^3) *
        (Z_i_reg^4)
    ) / ((Z_a_cat^4) * (Z_a_reg^4) + L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

Base.@kwdef struct MonodWymanChangeuxEnzyme
    Keq::Float64
    substrates::Union{NTuple{2, Symbol}, Nothing} = nothing
    products::Union{NTuple{2, Symbol}, Nothing} = nothing
    reg1::Union{NTuple{1, Symbol}, NTuple{2, Symbol}, NTuple{3, Symbol}, Nothing} = nothing
    reg2::Union{NTuple{1, Symbol}, NTuple{2, Symbol}, NTuple{3, Symbol}, Nothing} = nothing
    reg3::Union{NTuple{1, Symbol}, NTuple{2, Symbol}, NTuple{3, Symbol}, Nothing} = nothing
end

Base.@kwdef struct MWCEnzymeTemplate
    Keq::Float64
    # oligomeric_state::Int
    S1::Union{Symbol, Nothing} = nothing
    S2::Union{Symbol, Nothing} = nothing
    P1::Union{Symbol, Nothing} = nothing
    P2::Union{Symbol, Nothing} = nothing
    #rename R1_reg1 to R1_site1 and add R1_site2 to R4_site1 etc
    R1_reg1::Union{Symbol, Nothing} = nothing
    R2_reg1::Union{Symbol, Nothing} = nothing
    R3_reg1::Union{Symbol, Nothing} = nothing
    R1_reg2::Union{Symbol, Nothing} = nothing
    R2_reg2::Union{Symbol, Nothing} = nothing
    R3_reg2::Union{Symbol, Nothing} = nothing
    R1_reg3::Union{Symbol, Nothing} = nothing
    R2_reg3::Union{Symbol, Nothing} = nothing
    R3_reg3::Union{Symbol, Nothing} = nothing
end

function make_mwc_enzyme_template(Enzyme::MonodWymanChangeuxEnzyme)
    return MWCEnzymeTemplate(
        Keq = Enzyme.Keq,
        #rewrite substrate and products using get(substrate, 1, nothing) and get(substrate, 2, nothing) etc
        S1 = !isnothing(Enzyme.substrates) ? Enzyme.substrates[1] : nothing,
        S2 = !isnothing(Enzyme.substrates) ? Enzyme.substrates[2] : nothing,
        P1 = !isnothing(Enzyme.products) ? Enzyme.products[1] : nothing,
        P2 = !isnothing(Enzyme.products) ? Enzyme.products[2] : nothing,
        R1_reg1 = !isnothing(Enzyme.reg1) ? get(Enzyme.reg1, 1, nothing) : nothing,
        R2_reg1 = !isnothing(Enzyme.reg1) ? get(Enzyme.reg1, 2, nothing) : nothing,
        R3_reg1 = !isnothing(Enzyme.reg1) ? get(Enzyme.reg1, 3, nothing) : nothing,
        R1_reg2 = !isnothing(Enzyme.reg2) ? get(Enzyme.reg2, 1, nothing) : nothing,
        R2_reg2 = !isnothing(Enzyme.reg2) ? get(Enzyme.reg2, 2, nothing) : nothing,
        R3_reg2 = !isnothing(Enzyme.reg2) ? get(Enzyme.reg2, 3, nothing) : nothing,
        R1_reg3 = !isnothing(Enzyme.reg3) ? get(Enzyme.reg3, 1, nothing) : nothing,
        R2_reg3 = !isnothing(Enzyme.reg3) ? get(Enzyme.reg3, 2, nothing) : nothing,
        R3_reg3 = !isnothing(Enzyme.reg3) ? get(Enzyme.reg3, 3, nothing) : nothing
    )
end

macro derive_rate_eq(MWC_enzyme_name)
    MWCenzyme = eval(MWC_enzyme_name)
    @assert enzyme isa MonodWymanChangeuxEnzyme "@derive_rate_eq requires an argument of type ::MWCEnzyme"
    enzyme = make_mwc_enzyme_template(MWCenzyme)
    return :(
        function rate_equation(metabs, params)
        general_mwc_rate_equation(
            $(enzyme.S1 isa Symbol) ? metabs.$(enzyme.S1) : 1.0,
            $(enzyme.S2 isa Symbol) ? metabs.$(enzyme.S2) : 1.0,
            $(enzyme.P1 isa Symbol) ? metabs.$(enzyme.P1) : 1.0,
            $(enzyme.P2 isa Symbol) ? metabs.$(enzyme.P2) : 1.0,
            $(enzyme.R1_reg1 isa Symbol) ? metabs.$(enzyme.R1_reg1) : 1.0,
            $(enzyme.R2_reg1 isa Symbol) ? metabs.$(enzyme.R2_reg1) : 1.0,
            $(enzyme.R3_reg1 isa Symbol) ? metabs.$(enzyme.R3_reg1) : 1.0,
            $(enzyme.R1_reg2 isa Symbol) ? metabs.$(enzyme.R1_reg2) : 1.0,
            $(enzyme.R2_reg2 isa Symbol) ? metabs.$(enzyme.R2_reg2) : 1.0,
            $(enzyme.R3_reg2 isa Symbol) ? metabs.$(enzyme.R3_reg2) : 1.0,
            $(enzyme.R1_reg3 isa Symbol) ? metabs.$(enzyme.R1_reg3) : 1.0,
            $(enzyme.R2_reg3 isa Symbol) ? metabs.$(enzyme.R2_reg3) : 1.0,
            $(enzyme.R3_reg3 isa Symbol) ? metabs.$(enzyme.R3_reg3) : 1.0,
            params.L,
            params.Vmax_a,
            params.Vmax_i,
            $(enzyme.S1 isa Symbol) ? params.$(Symbol("K_a_$(enzyme.S1)_cat")) : 1.0,
            $(enzyme.S1 isa Symbol) ? params.$(Symbol("K_i_$(enzyme.S1)_cat")) : 1.0,
            $(enzyme.S2 isa Symbol) ? params.$(Symbol("K_a_$(enzyme.S2)_cat")) : 1.0,
            $(enzyme.S2 isa Symbol) ? params.$(Symbol("K_i_$(enzyme.S2)_cat")) : 1.0,
            $(enzyme.P1 isa Symbol) ? params.$(Symbol("K_a_$(enzyme.P1)_cat")) : 1.0,
            $(enzyme.P1 isa Symbol) ? params.$(Symbol("K_i_$(enzyme.P1)_cat")) : 1.0,
            $(enzyme.P2 isa Symbol) ? params.$(Symbol("K_a_$(enzyme.P2)_cat")) : 1.0,
            $(enzyme.P2 isa Symbol) ? params.$(Symbol("K_i_$(enzyme.P2)_cat")) : 1.0,
            $(enzyme.R1_reg1 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R1_reg1)_reg1")) : Inf,
            $(enzyme.R1_reg1 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R1_reg1)_reg1")) : Inf,
            $(enzyme.R2_reg1 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R2_reg1)_reg1")) : Inf,
            $(enzyme.R2_reg1 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R2_reg1)_reg1")) : Inf,
            $(enzyme.R3_reg1 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R3_reg1)_reg1")) : Inf,
            $(enzyme.R3_reg1 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R3_reg1)_reg1")) : Inf,
            $(enzyme.R1_reg2 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R1_reg2)_reg2")) : Inf,
            $(enzyme.R1_reg2 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R1_reg2)_reg2")) : Inf,
            $(enzyme.R2_reg2 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R2_reg2)_reg2")) : Inf,
            $(enzyme.R2_reg2 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R2_reg2)_reg2")) : Inf,
            $(enzyme.R3_reg2 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R3_reg2)_reg2")) : Inf,
            $(enzyme.R3_reg2 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R3_reg2)_reg2")) : Inf,
            $(enzyme.R1_reg3 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R1_reg3)_reg3")) : Inf,
            $(enzyme.R1_reg3 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R1_reg3)_reg3")) : Inf,
            $(enzyme.R2_reg3 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R2_reg3)_reg3")) : Inf,
            $(enzyme.R2_reg3 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R2_reg3)_reg3")) : Inf,
            $(enzyme.R3_reg3 isa Symbol) ?
            params.$(Symbol("K_a_$(enzyme.R3_reg3)_reg3")) : Inf,
            $(enzyme.R3_reg3 isa Symbol) ?
            params.$(Symbol("K_i_$(enzyme.R3_reg3)_reg3")) : Inf,
            params.alpha_PEP_ATP,
            params.alpha_ADP_Pyruvate,
            enzyme.Keq
        )
    end
    )
end
