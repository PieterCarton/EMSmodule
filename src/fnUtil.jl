using JuMP, InfiniteOpt

abstract type complement_formulation end
struct fully_relaxed      <: complement_formulation end
struct native             <: complement_formulation end
struct indicator          <: complement_formulation end
struct scholtes           <: complement_formulation end
struct lin_fukushima      <: complement_formulation end
struct fischer_burmeister <: complement_formulation end

abstract type battery_model_relaxation end
struct exact             <: battery_model_relaxation end
struct simple_lp         <: battery_model_relaxation end
struct nazir_almassalkhi <: battery_model_relaxation end
struct extn_lp           <: battery_model_relaxation end
struct to_lp             <: battery_model_relaxation end

struct FormulationSettings
    battery_model_relaxation
    complements::complement_formulation
    relaxation::Float32
end

const default_formulation_settings = FormulationSettings(exact(), scholtes(), 1e-5)

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
    # lower and upper bounds not needed!
    delete_lower_bound(xpos)
    delete_lower_bound(xneg)
    @constraint(model, xpos * xneg <= relaxation^2)
    @constraint(model, (xpos + relaxation) * (xneg + relaxation) >= relaxation^2)
end

function complement!(model::InfiniteModel, xpos, xneg, relaxation, complement_formulation::fischer_burmeister)
    @constraint(model, xpos + xneg - (xpos^2 + xneg^2 + relaxation)^(1/2) == 0)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::exact, genInfo::Generic)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::simple_lp, genInfo::Generic)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::nazir_almassalkhi, genInfo::Generic)
    # parameters
    # @unpack PowerLim, vLim, Np, Ns, ηC = genInfo
    # imax     =1e3*PbessMax/Np/Ns/vLim[1]

    # # variables
    # ibess⁻   = model[:ibess⁻]
    # ibess⁺   = model[:ibess⁺]

    # # parameters
    # PbessMax = upper_bound(PbessPos)
    # PbessMin = upper_bound(PbessNeg)
    # Pmax = max(PbessMax, PbessMin)

    # @constraint(model, ibess, ibess⁺ - ibess⁻)
    # @constraint(model, ibess⁺ - ibess⁻ <= Pmax)
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::extn_lp, genInfo::Generic)
    t=model[:t];
    Δt=Float64(supports(t)[2]-supports(t)[1])

    # variables
    PbessPos = model[:PbessPos]
    PbessNeg = model[:PbessNeg]
    ibess⁻   = model[:ibess⁻]
    ibess⁺   = model[:ibess⁺]
    Qbess    = model[:Qbess]
    SoCbess  = model[:SoCbess]
    
    # Extract data
    @unpack PowerLim, SoCLim, initQ, SoHQ, η = genInfo
    PbessMax = PowerLim[2]; # Max power [kW]
    PbessMin = PowerLim[1]; # Min power [kW]
    SoCbessMin = SoCLim[1]; # Min State of Charge [p.u.]
    SoCbessMax = SoCLim[2]; # Max State of Charge [p.u.]
    Qbess0 = initQ*SoHQ;    # Intial capacity [Ah]

    # 2e
    for (t_minus_one, t) in zip(supports(t)[1:end-1], supports(t)[2:end])
        @constraint(model, ibess⁻(t) * Δt * (1/3600) <= (Qbess(t) * SoCbessMax - Qbess0 * SoCbess(t_minus_one))/η)
    end
    
    # 2f
    for (t_minus_one, t) in zip(supports(t)[1:end-1], supports(t)[2:end])
        @constraint(model, ibess⁺(t) * Δt * (1/3600) <= Qbess0 * SoCbess(t_minus_one) - Qbess(t) * SoCbessMin)
    end

    # 1g
    @constraint(model, PbessPos <= (PbessMax - (PbessMax/PbessMin)) * PbessNeg)

    # println(model)
    # throw("STOPPED")
end

function battery_model!(model::InfiniteModel, battery_model_relaxation::to_lp, genInfo::Generic)
    t=model[:t];
    Δt=Float64(supports(t)[2]-supports(t)[1])

    # variables
    PbessPos = model[:PbessPos]
    PbessNeg = model[:PbessNeg]
    ibess⁻   = model[:ibess⁻]
    ibess⁺   = model[:ibess⁺]
    Qbess    = model[:Qbess]
    SoCbess  = model[:SoCbess]

    @variable(model, 0 <= δₚ <= 1, Infinite(t))
    
    # Extract data
    @unpack PowerLim, SoCLim, initQ, SoHQ = genInfo
    PbessMax = PowerLim[2]; # Max power [kW]
    PbessMin = PowerLim[1]; # Min power [kW]
    SoCbessMin = SoCLim[1]; # Min State of Charge [p.u.]
    SoCbessMax = SoCLim[2]; # Max State of Charge [p.u.]
    Qbess0 = initQ*SoHQ*3600;    # Intial capacity in As

    # # 4a
    for (t_minus_one, t) in zip(supports(t)[1:end-1], supports(t)[2:end])
        @constraint(model, SoCbess(t_minus_one) * Qbess0 >= SoCbessMin * Qbess0 + ibess⁺(t) * Δt)
    end
    
    # 4b
    for (t_minus_one, t) in zip(supports(t)[1:end-1], supports(t)[2:end])
        @constraint(model, SoCbess(t_minus_one) * Qbess0 <= SoCbessMax * Qbess0 - ibess⁻(t) * Δt)
    end

    # 1c
    @constraint(model, PbessNeg <= -PbessMin * δₚ)
    # 1d
    @constraint(model, PbessPos <= PbessMax *(1 - δₚ))
    # println(model)
    # throw("STOPPED")
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