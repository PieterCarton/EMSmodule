using JuMP, InfiniteOpt
using makeEMSobjs

abstract type complement_formulation end
struct fully_relaxed      <: complement_formulation end
struct native             <: complement_formulation end
struct indicator          <: complement_formulation end
struct scholtes           <: complement_formulation end
struct lin_fukushima      <: complement_formulation end
struct fischer_burmeister <: complement_formulation end

abstract type battery_model_relaxation end
struct exact             <: battery_model_relaxation end
struct nazir_almassalkhi <: battery_model_relaxation end
struct extn_lp           <: battery_model_relaxation end
struct to_lp             <: battery_model_relaxation end

struct FormulationSettings
    battery_model_relaxatin
    complements::complement_formulation
    relaxation::Float32
end

const default_formulation_settings = FormulationSettings(scholtes(), 1e-5)

function complement!(model::InfiniteModel, xpos, xneg)
    return complement!(model, xpos, xneg,
        default_formulation_settings.relaxation,
        default_formulation_settings.complements
    )
end

function complement!(model::InfiniteModel, xpos, xneg, formulation_settings::FormulationSettings)
    return complement!(model, xpos, xneg, formulation_settings.relaxation, formulation_settings.complements)
end

function complement!(model::InfiniteModel, xpos, xneg, relaxation, complement_formulation::fully_relaxed)
end

function complement!(model::InfiniteModel, xpos, xneg, relaxation, complement_formulation::native)
    return @constraint(model, xpos ⟂ xneg)
end

function complement!(model::InfiniteModel, xpos, xneg, relaxation, complement_formulation::scholtes)
    @constraint(model, xpos * xneg <= relaxation)
end

function complement!(model::InfiniteModel, xpos, xneg, relaxation, complement_formulation::lin_fukushima)
    @constraint(model, xpos * xneg <= relaxation^2)
    @constraint(model, (xpos + relaxation) * (xneg + relaxation) >= relaxation^2)
end

function complement!(model::InfiniteModel, xpos, xneg, relaxation, complement_formulation::fischer_burmeister)
    @constraint(model, xpos + xneg - (xpos^2 + xneg^2 + relaxation)^(1/2) == 0)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::nazir_almassalkhi)
    PbessMax = upper_bound(variable_by_name(model, "PbessPos"))
    PbessMin = upper_bound(variable_by_name(model, "PbessNeg"))
    Pmax = max(PbessMax, PbessMin)

    @constraint(model, PbessPos + PbessNeg <= Pmax)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::nazir_almassalkhi, genInfo::Generic)
    # Extract data
    @unpack PowerLim, P0, SoCLim, SoC0, ηC, termCond = data["BESS"].GenInfo
    PbessMax = PowerLim[2]; # Max power [kW]
    PbessMin = PowerLim[1]; # Min power [kW]
    SoCbessMin = SoCLim[1]; # Min State of Charge [p.u.]
    SoCbessMax = SoCLim[2]; # Max State of Charge [p.u.]
    SoCbess0 = SoC0; # Initial SoC [p.u.]

    # PbessMax = upper_bound(variable_by_name(model, "PbessPos"))
    # PbessMin = upper_bound(variable_by_name(model, "PbessNeg"))
    Pmax = max(PbessMax, PbessMin)

    @constraint(model, PbessPos + PbessNeg <= Pmax)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::extn_lp, genInfo::Generic)
    Δt=supports(t)[2]-supports(t)[1]
    
    # Extract data
    @unpack PowerLim, P0, SoCLim, SoC0, ηC, termCond = data["BESS"].GenInfo
    PbessMax = PowerLim[2]; # Max power [kW]
    PbessMin = PowerLim[1]; # Min power [kW]
    SoCbessMin = SoCLim[1]; # Min State of Charge [p.u.]
    SoCbessMax = SoCLim[2]; # Max State of Charge [p.u.]
    SoCbess0 = SoC0; # Initial SoC [p.u.]

    @unpack initQ, SoHQ = GenInfo
    Qbess0 = initQ*SoHQ;

    # 1g
    @constraint(model, PbessPos <= (PbessMax - (PbessMax/PbessMin)) * PbessNeg)
    # 2e
    @constraint(model, ibess⁻(t) <= Qbess(t)*SoCbessMax - Qbess0*SoCbess(t))
    @constraint(model, ibess⁺ <= (Qbess*SoCbess - (PbessMax/PbessMin)) * PbessNeg)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::extn_lp, genInfo::Generic)
    Δt=supports(t)[2]-supports(t)[1]
    
    # Extract data
    @unpack PowerLim, P0, SoCLim, SoC0, ηC, termCond = data["BESS"].GenInfo
    PbessMax = PowerLim[2]; # Max power [kW]
    PbessMin = PowerLim[1]; # Min power [kW]
    SoCbessMin = SoCLim[1]; # Min State of Charge [p.u.]
    SoCbessMax = SoCLim[2]; # Max State of Charge [p.u.]
    SoCbess0 = SoC0; # Initial SoC [p.u.]

    @unpack initQ, SoHQ = GenInfo
    Qbess0 = initQ*SoHQ;

    @variable(model, 0 <= δₚ <= 1, Infinite(t))

    # 1c
    @constraint(model, PbessNeg <= -PbessMin * δₚ)
    # 1d
    @constraint(model, PbessPos <= PbessMax *(1 - δₚ))

    # 4a
    @constraint(model, SoCbess(t-1) * Qbess0 >= SoCbessMin * Qbess0 + ibess⁺(t) * Δt)
    # 4b
    @constraint(model, SoCbess(t-1) * Qbess0 <= SoCbessMax * Qbess0 - ibess⁻(t) * Δt)
end

# function complement!(model::InfiniteModel, xpos, xneg, t, complement_formulation::indicator)
#     # create
#     xpos_indicator = @variable(model, Infinite(t), base_name=name(xpos)*"_indicator", binary=true)
#     xneg_indicator = @variable(model, Infinite(t), base_name=name(xneg)*"_indicator", binary=true)

#     @constraint(model, xpos_indicator + xneg_indicator = 0)
#     @constraint(model, xneg_indicator * xpos = 0)
#     @constraint(model, xpos_indicator * xneg = 0)

#     @constraint(model, xpos * xneg <= epsilon)
# end