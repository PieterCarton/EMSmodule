using JuMP, InfiniteOpt

abstract type complement_formulation end
struct no_complement <: complement_formulation end
struct native <: complement_formulation end
struct indicator <: complement_formulation end
struct multiplicative <: complement_formulation end

struct FormulationSettings
    complements::complement_formulation
end

const default_formulation_settings = FormulationSettings(multiplicative())

function complement!(model::InfiniteModel, xpos, xneg)
    return complement!(model, xpos, xneg, default_formulation_settings.complements)
end

function complement!(model::InfiniteModel, xpos, xneg, formulation_settings::FormulationSettings)
    return complement!(model, xpos, xneg, formulation_settings.complements)
end

function complement!(model::InfiniteModel, xpos, xneg, complement_formulation::native)
    return @constraint(model, xpos ⟂ xneg)
end

function complement!(model::InfiniteModel, xpos, xneg, complement_formulation::multiplicative)
    epsilon = 0.001
    @constraint(model, xpos * xneg <= epsilon)
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