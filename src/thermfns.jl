## Thermal carrier functions
# This file contains the functions used for the thermal devices (HP, TESS and loads) of the residential Energy Hub.
# The functions are made for modelling the thermal device performance based on its thermodynamics.
#
# By: Darío Slaifstein, PhD-student @TU Delft, DCES.
# Branch: thermal_models
# Version: 1.0
# Date: 07/10/2024
# import EMSmodule.simulate_storage_asset!

function occupancy_profiles(γ::Array{Float64,1},
    fs::Int64 = 4, # samples per hour
    )
    # This function generates the occupancy profiles for the residential buildings.
    # The profiles are generated based on the availability of EV, γ.
    # The occupancy 𝒪 is the timeseries of a person being at home at a given time t.
    ns = length(γ) # total number of samples of the timeseries
    ndays = Int(floor(ns/fs/24)); # number of days
    𝒪 = ones(ns); # initialize the occupancy profile
    # generate the time array t
    t = collect(1/fs:1/fs:ndays*24*fs);
    # Identify arrival and departure times
    arrIdx_γ = findall(diff(γ) .== 1)
    depIdx_γ = findall(diff(γ) .== -1)
    length(arrIdx_γ) == ndays ? nothing : arrIdx_γ = [arrIdx_γ; ns]

    # Create the occupancy signal
    for day in 1:ndays
        # Determine the time indices corresponding to arrival and departure for this day
        t = collect(0:1/fs:(24-1/fs))
        # create noisy arrival and departure times
        εArr = round.(Int, rand(Normal(0, 1)))
        εDep = round.(Int, rand(Normal(0, 1)))
        # daily arrival and departure times
        depIdx = depIdx_γ[day]
        arrIdx = arrIdx_γ[day]
        # check if the arrival and departure times are within the timeseries
        if (any((depIdx .+ εDep) .> ns))
            εDep = 0;
        elseif (any((depIdx .+ εDep) .< 1))
            εDep = 0;
        elseif (any((arrIdx .+ εArr) .> ns))
            εArr = 0;
        elseif (any((arrIdx .+ εArr) .< 1))
            εArr = 0;
        end
        depIdx += εDep
        arrIdx += εArr
        # # Correct for the day
        # depIdx = depIdx + (day-1)*24*fs
        # arrIdx = arrIdx + (day-1)*24*fs
        # Mark the time series as parked during the parked interval for this day
        if depIdx <= arrIdx
            𝒪[depIdx:arrIdx] .= 0
        end
    end
    return 𝒪;
end

# Thermal devices
"""
    st(model, data)
    
The `st` function calculates the solar thermal power based on the given model and data.
# Arguments
- `model::InfiniteModel`: The model object representing the system.
- `data::Dict`: A dictionary containing the necessary data for the calculation.
# Returns
- `model::InfiniteModel`: The updated model object.
Note: The curtailment part of the code is currently commented out and left for future implementation.
"""
function st!(model::InfiniteModel, data::Dict; add_noise::Bool = true) # solar thermal 
    # The thermal pv/heatpipes are only the converted measurement of the irradiance.
    t = model[:t]; 
    Dt = supports(t);
    t0 = supports(t)[1]; Δt = supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(supports(t))-1;
    # Extract params
    # PstRated = data["ST"].RatedPower; # rated power of the panels
    ηST = data["ST"].η; # conversion factor from Electric PV to thermal
    # Extract data
    MPPTmeas = data["SPV"].MPPTData[it0:itend];
    # simulated forecast
    MPPTnoisy = copy(MPPTmeas);
    if add_noise
        MPPTnoisy[MPPTmeas .> 0] = [rand(Uniform(0.5*MPPTmeas[tt],1.1*MPPTmeas[tt])) for tt ∈ findall(MPPTmeas .> 0)]
    end
    @parameter_function(model, Pst == (t) -> ηST*MPPTnoisy[zero_order_time_index(Dt, t)])

    #= Curtailment is left for the future
        @variable(model, 0 ≤ Pst ≤ PstRated, Infinite(t));
        @constraint(model, Pst == ηST*(model[:Ppve]+model[:PpvDwn]));
    =#
    return model;
end;

"""
    st(model, data)
    
The `st_T` function calculates the solar thermal power based on the given model and data.
# Arguments
- `model::InfiniteModel`: The model object representing the system.
- `data::Dict`: A dictionary containing the necessary data for the calculation.
# Returns
- `model::InfiniteModel`: The updated model object.
Note: The curtailment part of the code is currently commented out and left for future implementation.
"""
function st_T!(model::InfiniteModel, data::Dict) # solar thermal 
    # The thermal pv/heatpipes are only the converted measurement of the irradiance.
    t = model[:t]; 
    Dt = supports(t);
    t0 = supports(t)[1]; Δt = supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(supports(t))-1;
    # Extract params
    # PstRated = data["ST"].RatedPower; # rated power of the panels
    ηST = data["ST"].η; # conversion factor from Electric PV to thermal

    # Extract data
    MPPTmeas = data["SPV"].MPPTData[it0:itend];
    radiation = data["grid"].Gir[it0:itend]; # CHECK add to object
    # simulated forecast
    MPPTnoisy = zeros(size(MPPTmeas));
    MPPTnoisy[MPPTmeas .> 0] = [rand(Uniform(0.5*MPPTmeas[tt],1.1*MPPTmeas[tt])) for tt ∈ findall(MPPTmeas .> 0)]
    @parameter_function(model, Gir == (t) -> radiation[zero_order_time_index(Dt, t)])
    @parameter_function(model, Pst == (t) -> ηST*MPPTnoisy[zero_order_time_index(Dt, t)])
    @variables(model, begin
        0 ≤ Q̇ₛₜ, Infinite(t)
        0 ≤ Q̇ₛₜᴰ, Infinite(t)
        0 ≤ Q̇ₛₜᵗᵉˢˢ, Infinite(t)
        bST, Infinite(t), Bin # ON-OFF variable for the solar thermal
        # 273 ≤ TstIn ≤ 273+100, Infinite(t)
        # 273 ≤ TstOut ≤ 273+100, Infinite(t)
        0 ≤ ΔTₛₜ ≤ 100, Infinite(t)
    end)
    Aₛₜ = 10; # area of the PVT solar collector [m²]
    h_air_window = 25 # Indoor air -> windows [W/m².K]
    h_window_atm = 32 # Windows -> atmosphere [W/m².K]
    LWindow = 0.004 # Thickness of a single window pane [m]
    LCavity = 0.014 # Thickness of the cavity between the double glass window [m]  
    kWindow = 0.8 # Thermal conductivity of glass [W/mK]
    kAir = 0.0257 # Thermal conductivity of air at 293 K [W/mK]
    Uₛₜ = 1/((1/h_air_window) + (LWindow/kWindow) + (LCavity/kAir) + (LWindow/kWindow) + (1/h_window_atm)); # heat loss coefficient [W/m².K]
    ΔTₛₜ = Gir/Uₛₜ/Aₛₜ
    @constraints(model, begin
        Q̇ₛₜ .== Q̇ₛₜᴰ + Q̇ₛₜᵗᵉˢˢ
        Q̇ₛₜ .* bST .== ηₛₜ * Uₛₜ * Aₛₜ  * ΔTₛₜ
    end);

    return model;
end;
"""
    tess(model::InfiniteModel, data::Dict)

The `tess` function models a Thermal Energy Storage System (TESS) as bucket model.
It adds variables and constraints to the model to represent the TESS's state of charge,
thermal power, and binary variable for TESS power. The function also includes initial
conditions and a bucket model constraint.

# Arguments
- `model::InfiniteModel`: The InfiniteModel to which the TESS variables and constraints will be added.
- `data::Dict`: A dictionary containing the TESS data, including capacity, power limits, initial state of charge, and thermal efficiency.

# Returns
- `model::InfiniteModel`: The updated InfiniteModel with the TESS variables and constraints added.
"""
function tess!(model::InfiniteModel, data::Dict) # thermal energy storage buffer
    # Thermal Energy Storage System
    t = model[:t];
    t0 = supports(t)[1];
    # Bucket model
    # Extract data
    Qtess = data["TESS"].Q; # Capacity [kWh]
    PtessMin = data["TESS"].PowerLim[1]; # Min power [kW]
    PtessMax = data["TESS"].PowerLim[2]; # Max power [kW]
    SoCtessMin = data["TESS"].SoCLim[1]; # Min State of Charge [p.u.]
    SoCtess0 = data["TESS"].SoC0; # Initial State of Charge [p.u.]
    ηtess = data["TESS"].η; # thermal efficiency
       
    # Add variables
    @variables(model, begin
        # SoCtessMin ≤ SoCtess ≤ SoCtessMax, Infinite(t) # State of Charge
        SoCtess ≥ SoCtessMin, Infinite(t) # State of Charge
        0 ≤ PtessPos ≤ PtessMax, Infinite(t) # Ptess^+ out power
        0 ≤ PtessNeg ≤ -PtessMin, Infinite(t)
    end);
    # Dummy variables for bidirectional power flow, ensuring only export or import
    # PtessNeg + PtessPos .== Ptess
    @expression(model, Ptess, PtessPos * ηtess .- PtessNeg * (1/ηtess))
    @constraints(model, begin
        # MPEC with ⟂
        PtessNeg ⟂ PtessPos
    end);
    # Initial conditions
    @constraint(model, SoCtess(t0) .== SoCtess0)
    # Bucket model
    @constraint(model, ∂.(SoCtess, t) .== -Ptess/Qtess/3600);
    return model;
end;

"""
    tess_T(model::InfiniteModel, data::Dict)

The `tess_T` function models a Thermal Energy Storage System (TESS) with the temperature-based approach.
It adds variables and constraints to the model to represent the TESS's state of charge,
thermal power, and binary variable for TESS power. The function also includes initial
conditions and a bucket model constraint.

# Arguments
- `model::InfiniteModel`: The InfiniteModel to which the TESS variables and constraints will be added.
- `data::Dict`: A dictionary containing the TESS data, including capacity, power limits, initial state of charge, and thermal efficiency.

# Returns
- `model::InfiniteModel`: The updated InfiniteModel with the TESS variables and constraints added.
"""
function tess_T!(model::InfiniteModel, data::Dict) # thermal energy storage buffer
    # Thermal Energy Storage System
    t = model[:t]; Q̇ₕₚᵗᵉˢˢ = model[:Q̇ₕₚᵗᵉˢˢ]; 
    # Q̇ₚᵥₜᵗᵉˢˢ = model[:Q̇ₚᵥₜᵗᵉˢˢ];
    t0 = supports(t)[1]; Δt = supports(t)[2]-supports(t)[1];
    # Bucket model
    # Extract data
    # Ptess0 = data["TESS"].P0; # Initial Power [p.u.]
    # SoCtessMin = data["TESS"].SoCLim[1]; # Min State of Charge [p.u.]
    # SoCtessMax = data["TESS"].SoCLim[2]; # Max State of Charge [p.u.]
    # SoCtess0 = data["TESS"].SoC0; # Initial State of Charge [p.u.]
    @unpack_TESSData_T data["TESS"]
    Qtess = Q; # Capacity [kWh]
    Q̇tessMin = QLim[1]; # Min power [kW]
    Q̇tessMax = QLim[2]; # Max power [kW]
    TtessMin = TLim[1]; # Min T [°C]
    TtessMax = TLim[2]; # Max T [°C]
    Ttess0 = T0; # Initial T [°C]
    ηtess = η; # thermal efficiency
    mtess = m; # tank mass [kg]
    ctess = c; # tank specific heat capacity [J/kg.K]
    ṁf = data["HP"].mdot; # mass flow rate [kg/s]
    cf = data["HP"].c; # specific heat capacity [J/kg.K]
    ṁf *= 3600; # to [kg/h]
    cf /= (3600 * 1e3); # to [kWh/kg.K]
    # ctess /= (3600 * 1e3); # to [kWh/kg.K]
    # cf /= 1e3; # to [kWs/kg.K]
    ctess /= 1e3; # to [kWs/kg.K]
    
    # Add variables
    @variables(model, begin
        # SoCtessMin ≤ SoCtess ≤ SoCtessMax, Infinite(t) # State of Charge
        273 + TtessMin ≤ Ttess ≤ 273 + TtessMax, Infinite(t) # TESS internal Temperature
        273 + 30 ≤ TtessIn ≤ 273 + 90, Infinite(t) # Inlet T of the heat exchanger
        273 + 30 ≤ TtessOut ≤ 273 + 90, Infinite(t) # Outlet T of the heat exchanger
        # Ptess, Infinite(t) # Thermal power
        # bPtess, Infinite(t), Bin # Binary variable for TESS power
        0. ≤ Q̇ₜₑₛₛᴰ ≤ Q̇tessMax, Infinite(t) # TESS Heat to the demand
        # bQtess, Infinite(t), Bin # Binary variable for TESS heat
        # 0 ≤ Q̇ₜₑₛₛᴰ⁺, Infinite(t) # Qtess⁺ out heat
        # Q̇ₜₑₛₛᴰ⁻ ≤ 0, Infinite(t) # Qtess⁻ in heat
    end);
    # Dummy variables for bidirectional power flow, ensuring only export or import
    # @constraints(model, begin
    #     bQtess*Q̇tessMin ≤ Q̇ₜₑₛₛᴰ⁻
    #     Q̇ₜₑₛₛᴰ⁻ + Q̇ₜₑₛₛᴰ⁺ .== Q̇ₜₑₛₛᴰ
    #     Q̇ₜₑₛₛᴰ⁺ ≤ (1-bQtess)*Q̇tessMax
    # end);
    fix(TtessOut, 55. + 273, force = true) # fix the outlet temperature of the TESS
    # Initial conditions
    @constraint(model, Ttess(t0) .== Ttess0 + 273)

    # self-discharge loss CHECK
    # @parameter_function(model, Q̇sd == (t) -> -0.01*Qtess); # 1% of the capacity per hour
    Q̇sd = 0.001*Qtess # 0.1% of the capacity per hour
    mf = 50
    # Bucket model
    @constraints(model, begin
                # ∂.(Ttess, t) * mtess * ctess .== Q̇ₕₚᵗᵉˢˢ + Q̇ₚᵥₜᵗᵉˢˢ - Q̇ₜₑₛₛᴰ - Q̇sd
                ∂.(Ttess, t) * mtess * ctess .== Q̇ₕₚᵗᵉˢˢ - Q̇ₜₑₛₛᴰ - Q̇sd
                # ∂.(TtessIn, t) * mf * cf .== Q̇ₜₑₛₛᴰ - ηtess * ṁf * cf * (TtessIn - TtessOut) # CHECK
                Q̇ₜₑₛₛᴰ - ηtess * ṁf * cf * (TtessOut - TtessIn) .== 0 
                # TtessIn(t0) .== 40. + 273.
                # Q̇ₜₑₛₛᴰ = ηTESS * ṁf * cf * (Tsup - T2)
                # goal is to T3->T4->Tsup
                end);
    return model;
end;

"""
    heatpump!(model::InfiniteModel, data::Dict)

Adds a basic heat pump model to the EMS optimization model with electric power variable only.

This is a simplified heat pump model that only adds the electric power consumption variable
without thermal dynamics or COP calculations. For more detailed heat pump models with 
thermal behavior, use `heatpump_nl!` or `heatpump_milp!`.

# Arguments
- `model::InfiniteModel`: The InfiniteOpt model to which the heat pump variables will be added.
- `data::Dict`: A dictionary containing the heat pump data with the following required key:
  - `"HP"`: Heat pump data structure containing:
    - `RatedPower`: Maximum electric power rating of the heat pump [kW]

# Returns
- `model::InfiniteModel`: The updated InfiniteModel with the heat pump electric power variable added.

# Variables Added
- `Phpe`: Electric power consumption of the heat pump [kW], bounded between 0 and rated power

# Example
```julia
# Add basic heat pump to model
model = heatpump!(model, data)
```

# See Also
- [`heatpump_nl!`](@ref): Nonlinear heat pump model with thermal dynamics
- [`heatpump_milp!`](@ref): MILP approximation of heat pump with COP modeling
"""
function heatpump!(model::InfiniteModel, data::Dict) # heat pump
    # The heat pump has a variable (electrical) and a subordinate finite_param (thermal)
    # Extract params
    t = model[:t];
    PhpRated = data["HP"].RatedPower; # rated power of the heat pump
    @variable(model, 0 ≤ Phpe ≤ PhpRated, Infinite(t));  # Electric power
    return model;
end;

"""
    heatpump_nl!(model::InfiniteModel, data::Dict; add_noise::Bool=false)

Adds a nonlinear heat pump model with detailed thermal dynamics to the EMS optimization model.

This function implements a comprehensive heat pump model that includes:
- Electric power consumption with thermal heat generation
- Dual heat output streams (to house demand and TESS)
- Temperature-based heat exchanger modeling with inlet/outlet temperatures
- Ambient temperature consideration with optional noise modeling
- Fixed COP (Coefficient of Performance) approach
- Complementarity constraints ensuring exclusive heat distribution

The model uses two parallel heat exchangers for supplying heat to the building demand
and the Thermal Energy Storage System (TESS), with fixed outlet temperatures and
thermal efficiency considerations.

# Arguments
- `model::InfiniteModel`: The InfiniteOpt model to which the heat pump variables and constraints will be added.
- `data::Dict`: A dictionary containing the required heat pump and building data:
  - `"HP"`: Heat pump data structure containing:
    - `mdot`: Mass flow rate [kg/s]
    - `c`: Specific heat capacity [J/kg·K]
    - `η`: Thermal efficiency of the heat exchanger [dimensionless]
    - `RatedPower`: Maximum electric power rating [kW]
  - `"Building"`: Building data structure containing:
    - `ambT`: Ambient temperature time series [°C]
- `add_noise::Bool=false`: Whether to add Gaussian noise (2°C std dev) to ambient temperature for forecast simulation

# Returns
- `model::InfiniteModel`: The updated InfiniteModel with heat pump variables and constraints added.

# Variables Added
- `Phpe`: Electric power consumption [kW], bounded by rated power
- `Q̇ₕₚᴰ`: Heat flow to house demand [kW], max 12 kW
- `Q̇ₕₚᵗᵉˢˢ`: Heat flow to TESS [kW], max 12 kW  
- `ThpDⁱⁿ`, `ThpTESSⁱⁿ`: Inlet temperatures for demand and TESS heat exchangers [K]
- `ThpDᵒᵘᵗ`, `ThpTESSᵒᵘᵗ`: Outlet temperatures for demand and TESS heat exchangers [K]
- `COPhat`: Heat pump coefficient of performance (fixed at 3.0)

# Parameters Added
- `Tamb`: Ambient temperature parameter function [K]

# Constraints Added
- Complementarity constraint: `Q̇ₕₚᴰ ⟂ Q̇ₕₚᵗᵉˢˢ` (exclusive heat distribution)
- Energy balance: `Q̇ₕₚᴰ + Q̇ₕₚᵗᵉˢˢ = COPhat × Phpe`
- Heat exchanger equations for both demand and TESS circuits
- Fixed outlet temperatures (55°C for demand, 90°C for TESS)

# Example
```julia
# Add nonlinear heat pump model with noise simulation
model = heatpump_nl!(model, data; add_noise=true)
```

# See Also
- [`heatpump!`](@ref): Basic heat pump model with electric power only
- [`heatpump_milp!`](@ref): MILP approximation with temperature-dependent COP
"""
function heatpump_nl!(model::InfiniteModel, data::Dict; add_noise::Bool=false) # heat pump
    # The heat pump has a variable (electrical) and a subordinate finite_param (thermal)
    t = model[:t];
    Dt = supports(t);
    t0 = supports(t)[1]; Δt = supports(t)[2]-supports(t)[1];
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # Extract params
    ṁf = data["HP"].mdot; # mass flow rate
    ṁf *= 3600; # mass flow rate [kg/h]
    cf = data["HP"].c; # specific heat capacity
    cf /= (3600 * 1e3); # specific heat capacity [kWh/kg.K]
    # cf /= 1e3; # specific heat capacity [kWs/kg.K]
    ηₕₚ = data["HP"].η; # thermal efficiency of the heat exchanger
    PhpRated = data["HP"].RatedPower; # rated power of the heat pump
    ambT = data["Building"].ambT[it0:itend] .+ 273; # CHECK add to object
    if add_noise
        εₗᵗʰ= randn(Int(length(ambT)*Δt/3600)) .* 2 # 2 °C noise
        εₗᵗʰ = repeat(εₗᵗʰ, inner = Int(3600/Δt))
        ambT = ambT .+ εₗᵗʰ; # add noise
    end
    @parameter_function(model, Tamb == (t) -> ambT[zero_order_time_index(Dt, t)])
    
    @variables(model, begin
        0 ≤ Phpe ≤ PhpRated, Infinite(t) # Electric power
        0 ≤ Q̇ₕₚᴰ ≤ 12, Infinite(t), (start = 0) # Heat to the house
        0 ≤ Q̇ₕₚᵗᵉˢˢ ≤ 12, Infinite(t), (start = 0) # Heat to the TESS
        273 + 30 ≤ ThpDⁱⁿ ≤ 60 + 273, Infinite(t) # inlet temperature
        273 + 30 ≤ ThpTESSⁱⁿ ≤ 90 + 273, Infinite(t) # inlet temperature
        10 + 273 ≤ ThpDᵒᵘᵗ ≤ 60 + 273, Infinite(t) # outlet temperature of the heat exchanger to the demand
        10 + 273 ≤ ThpTESSᵒᵘᵗ ≤ 90 + 273, Infinite(t) # outlet temperature of the heat exchanger to the TESS
        COPhat, Infinite(t)
    end);
    # fix the outlet temperatures
    fix(COPhat, 3., force = true) # fixed COP
    fix(ThpTESSᵒᵘᵗ, 90. + 273, force = true) 
    fix(ThpDᵒᵘᵗ, 55. + 273, force = true)
    # fix(ThpOut, 50., force = true) # fix the outlet temperature DHW
    # fix(ThpOut, 18., force = true) # fix the outlet temperature cooling
    # CHECK WITH JOEL
    # the 2 heat exchangers are in parallel and at least half as the main pipeline
    # we can try with other values
    mf1 = 500; # total fluid mass in the pipeline [kg]
    mf2 = mf1;
    @constraints(model, begin
        Q̇ₕₚᴰ ⟂ Q̇ₕₚᵗᵉˢˢ
        Q̇ₕₚᴰ + Q̇ₕₚᵗᵉˢˢ  .== COPhat*Phpe # produced heat
        Q̇ₕₚᴰ - ηₕₚ * ṁf * cf * (ThpDᵒᵘᵗ - ThpDⁱⁿ) .== 0 # CHECK
        Q̇ₕₚᵗᵉˢˢ - ηₕₚ * ṁf * cf * (ThpTESSᵒᵘᵗ - ThpTESSⁱⁿ) .== 0 # CHECK
    end);
    return model;
end;

"""
    heatpump_milp!(model::InfiniteModel, data::Dict)

Adds a Mixed-Integer Linear Programming (MILP) heat pump model with temperature-dependent COP to the EMS optimization model.

This function implements an advanced heat pump model using piecewise affine approximation to capture
the nonlinear relationship between the Coefficient of Performance (COP) and ambient temperature.
The model includes:
- Temperature-dependent COP using exponential relationship: COP = 7.90471 × exp(-0.024 × ΔT)
- Piecewise linear approximation with 3 planes for MILP compatibility
- Single heat output stream (to house demand only)
- Heat exchanger modeling with inlet/outlet temperature relationships
- Fixed outlet temperature for floor heating applications

The COP varies based on the temperature difference between inlet and ambient temperatures,
providing more realistic heat pump performance modeling compared to fixed COP approaches.

# Arguments
- `model::InfiniteModel`: The InfiniteOpt model to which the heat pump variables and constraints will be added.
- `data::Dict`: A dictionary containing the required heat pump and grid data:
- `"HP"`: Heat pump data structure containing:
- `mdot`: Mass flow rate [kg/s]
- `c`: Specific heat capacity [J/kg·K]  
- `η`: Thermal efficiency/conversion factor [dimensionless]
- `RatedPower`: Maximum electric power rating [kW]
- `"grid"`: Grid data structure containing:
- `ambT`: Ambient temperature time series [°C]

# Returns
- `model::InfiniteModel`: The updated InfiniteModel with heat pump variables and constraints added.

# Variables Added
- `Phpe`: Electric power consumption [kW], bounded by rated power
- `Q̇ₕₚᴰ`: Heat flow to house demand [kW], unbounded (≥ 0)
- `ThpIn`: Heat pump inlet temperature [K], range 273-372 K
- `ThpOut`: Heat pump outlet temperature [K], range 273-372 K (fixed at 35°C + 273)
- `COPhat`: Temperature-dependent coefficient of performance, approximated by PWA

# Parameters Added
- `Tamb`: Ambient temperature parameter function [K]

# Constraints Added
- Piecewise affine COP constraints using 3 linear planes over temperature range 20-60°C
- Energy balance: `Q̇ₕₚᴰ = COPhat × Phpe`
- Heat exchanger relationship: `ThpIn = ThpOut - Q̇ₕₚᴰ/(η × ṁf × cf)`
- Fixed outlet temperature at 35°C (308 K) for floor heating

# COP Model
The nonlinear COP function `COP(ΔT) = 7.90471 × exp(-0.024 × ΔT)` where ΔT = ThpIn - Tamb
is approximated using a 3-plane piecewise affine function over the range [20°C, 60°C] using
the HiGHS optimizer for the approximation generation.

# Example
```julia
# Add MILP heat pump model with temperature-dependent COP
model = heatpump_milp!(model, data)
```

# See Also
- [`heatpump!`](@ref): Basic heat pump model with electric power only
- [`heatpump_nl!`](@ref): Nonlinear heat pump model with fixed COP and dual outputs

# Notes
- The piecewise affine approximation makes this model suitable for MILP solvers
- Currently optimized for floor heating applications (35°C outlet temperature)
- The model assumes convex COP relationship for the approximation
"""
function heatpump_milp!(model::InfiniteModel, data::Dict) # heat pump
    # The heat pump has a variable (electrical) and a subordinate finite_param (thermal)
    t = model[:t];
    Dt = supports(t);
    t0 = supports(t)[1]; Δt = supports(t)[2]-supports(t)[1];
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;

    # Extract params
    ṁf = data["HP"].mdot; # mass flow rate
    ṁf *= 3600; # to [kg/h]
    cf = data["HP"].c; # specific heat capacity
    cf *= 1e-3; # to [kWh/kg.K]
    ηₕₚ = data["HP"].η; # conversion factor from Electric to thermal heat pump
    t = model[:t]; # t0=supports(t)[1];
    PhpRated = data["HP"].RatedPower; # rated power of the heat pump
    
    ambT = data["grid"].ambT[it0:itend] .+ 273; # CHECK add to object
    @parameter_function(model, Tamb == (t) -> ambT[zero_order_time_index(Dt, t)])

    @variables(model, begin
        0 ≤ Phpe ≤ PhpRated, Infinite(t) # Electric power
        0 ≤ Q̇ₕₚᴰ, Infinite(t) # Heat to the house
        # 0 ≤ Q̇ₕₚᵗᵉˢˢ, Infinite(t) # Heat to the TESS
        273. ≤ ThpIn ≤ 99 + 273, Infinite(t) # inlet temperature
        273. ≤ ThpOut ≤ 99 + 273, Infinite(t) # inlet temperature
        0 ≤ COPhat, Infinite(t) # inlet temperature
        end);  
    fix(ThpOut, 35. + 273, force = true) # fix the outlet temperature floor-heat
    # COP nonlinear relationship
    COPfn(x) = 7.90471 .* exp.(-0.024 .* first(x))
    # COP = COPfn.(20:1:60) # points
    # PWA approximation
    COPpwa=approx(
        COPfn, [(20, 60)],
        Convex(),
        MILP(optimizer = HiGHS.Optimizer,
         planes = 3),
    )
    # add to the model
    # THIS doesn´t work for InfiniteModels
    # pwaffine(model,
    #         sThpIn,
    #         COPpwa;
    #         z = COPhat,
    #         )

    @constraints(model, begin
        COPhat ≥ COPpwa.planes[1].α[1] * (ThpIn - Tamb) + COPpwa.planes[1].β[1]
        COPhat ≥ COPpwa.planes[2].α[1] * (ThpIn - Tamb) + COPpwa.planes[1].β[1]
        COPhat ≥ COPpwa.planes[3].α[1] * (ThpIn - Tamb) + COPpwa.planes[1].β[1]
        Q̇ₕₚᴰ .== COPhat*Phpe # produced heat CHECK COULD BE EXPRESSION
        # Q̇ₕₚᴰ + QtessHP == COPhat*Phpe # produced heat CHECK COULD BE EXPRESSION
        ThpIn .==  ThpOut - Q̇ₕₚᴰ/(ηₕₚ * ṁf * cf) # outlet temperature
        # ThpIn ≤ ThpOut
        end);
    return model;
end;

"""
    gridThermal(model::InfiniteModel, sets::modelSettings, data::Dict, add_noise::Bool=true)

This function adds the thermal balance to the EMS model. Using the flow-based approach.

# Arguments
- `model::InfiniteModel`: The InfiniteModel object representing the grid.
- `sets::modelSettings`: The model settings.
- `data::Dict`: A dictionary containing the necessary data for the calculation.
- `add_noise::Bool`: Whether to add noise to the thermal load for simulating a forecast.

# Returns
- `model::InfiniteModel`: The updated InfiniteModel object with the thermal power balance constraint.

"""
function gridThermal!(model::InfiniteModel, sets::modelSettings, data::Dict; add_noise::Bool = true) # thermal grid
    # Load data
    Dt = sets.dTime;
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # make InfiniteOpt compatible Plt(t) work for arbitrary times
    loadTh = data["grid"].loadTh[it0:itend];
    if add_noise
        # white noise to simulate a forecast
        εₗᵗʰ= randn(Int(length(loadTh)*Δt/3600)) .* 0.05 # 50W
        εₗᵗʰ = repeat(εₗᵗʰ, inner = Int(3600/Δt))
        loadTh = loadTh .+ εₗᵗʰ; # add noise
        loadTh[loadTh .< 0] .= 0
    end
    @parameter_function(model, Plt == (t) -> loadTh[zero_order_time_index(Dt, t)])
    
    # When there´s no FCR then there´s no decision variable Ppve hence Pst can just be a @parameter_function.
    # Extract params
    ηHP=data["HP"].η; # conversion factor from Electric to thermal heat pump
    # Thermal Power balance
    model[:thBalance]=@constraint(model, model[:Pst]+model[:Phpe].*ηHP+model[:Ptess] .==  Plt);
    resLoadTh = loadTh .- data["ST"].η * data["SPV"].MPPTData[it0:itend];
    resLoadTh⁺ = copy(resLoadTh); resLoadTh⁻ = copy(resLoadTh);
    resLoadTh⁺[resLoadTh .< 0] .= 0; resLoadTh⁻[resLoadTh .> 0] .= 0;
    set_start_value_function(model[:PtessPos], t -> resLoadTh⁺[zero_order_time_index(Dt, t)])
    set_start_value_function(model[:PtessNeg], t -> -resLoadTh⁻[zero_order_time_index(Dt, t)])
    return model;
end;

"""
    gridThermal_T(model::InfiniteModel, sets::modelSettings, data::Dict, add_noise::Bool=true)

This function adds the building thermal balance to the EMS model, using the temperature-based approach.

# Arguments
- `model::InfiniteModel`: The InfiniteModel object representing the grid.
- `sets::modelSettings`: The model settings.
- `data::Dict`: A dictionary containing the necessary data for the calculation.
- `add_noise::Bool`: Whether to add noise to the ambient temperature for simulating a forecast.

# Returns
- `model::InfiniteModel`: The updated InfiniteModel object with the thermal power balance constraint.

"""
function gridThermal_T!(model::InfiniteModel, sets::modelSettings, data::Dict; add_noise::Bool = false) # thermal grid
    # Load data
    # Dt = sets.dTime;
    Dt = supports(model[:t]);
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # make InfiniteOpt compatible Wₜ₊₁ work for arbitrary times
    # loadTh = data["grid"].loadTh[it0:itend];
    ṁf = data["HP"].mdot; # mass flow rate [kg/s]
    cf = data["HP"].c; # specific heat capacity [J/kg.K]
    @unpack_BuildingData data["Building"]
    ambT = data["Building"].ambT[it0:itend]; # CHECK add to object
    radiation = data["Building"].Gir[it0:itend]; # CHECK add to object
    𝒮 = length(A) # number of surfaces
    Cair *= 3600 # to [kWs/kg.K]
    Cb *= 3600 # to [kWs/K]
    rb /= 3600 # to [1/s]
    @parameter_function(model, Gir == (t) -> radiation[zero_order_time_index(Dt, t)])
    Tamb = model[:Tamb]; Q̇ₕₚᴰ = model[:Q̇ₕₚᴰ]; Q̇ₜₑₛₛᴰ = model[:Q̇ₜₑₛₛᴰ]; 
    # white noise to simulate a forecast
    if add_noise
        εₗᵗʰ= randn(Int(length(Tamb)*Δt/3600)) .* 2 # 2 °C noise
        εₗᵗʰ = repeat(εₗᵗʰ, inner = Int(3600/Δt))
        Tamb = Tamb .+ εₗᵗʰ; # add noise
    end
    # T3 = model[:T3]; 
    @variables(model, begin
        273 ≤ Tin ≤ 273 + 40., Infinite(t) # Building inside temperature
        # 0. ≤ T1 ≤ 50. + 273., Infinite(t) # temperature 1 of the pipelines
        # 273. ≤ T2 ≤ 80. + 273., Infinite(t) # temperature 2 of the pipelines
        # 273. ≤ T3 ≤ 80. + 273., Infinite(t) # temperature 3 of the pipelines
        # 273. ≤ T4 ≤ 70. + 273., Infinite(t) # temperature 4 of the pipelines
        0. ≤ Q̇vent, Infinite(t) # ventilation losses
        0. ≤ Q̇cond, Infinite(t) # conductive losses
        0. ≤ Q̇loss, Infinite(t) # heat losses
    end);

    # Q̇ir heat from solar irradiance
    Q̇ir = Gir * wb * sb * sum(A[2:3]) # heat from solar irradiance
    # Q̇ir = 0. 
    model[:Q̇ir] = Q̇ir; # save in model
    model[:thBalance]=@constraint(model, ∂.(Tin, t) .== (Q̇ir .+ Q̇ₜₑₛₛᴰ + Q̇ₕₚᴰ - Q̇loss)./(Cb + Vb * ρair * Cair));
    mf = 1000; # total fluid mass in the pipeline [kg]
    η = 0.8; # heat exchanger efficiency
    @constraints(model, begin
        # Temp dynamics
        # ∂.(Tⁱⁿ, t) * mf * cf .== Q̇ - ηₕₚ * ṁf * cf * (Tⁱⁿ - Tᵒᵘᵗ)
        # ∂.(T2, t) * mf * cf .==  Q̇ₕₚᴰ - η * ṁf * cf * (T2 - T3)
        # ∂.(T3, t) * mf * cf .==  Q̇ₜₑₛₛᴰ - η * ṁf * cf * (T3 - T4)
        # ∂.(T4, t) * mf * cf .==  -(Q̇ir .+ Q̇ₜₑₛₛᴰ + Q̇ₕₚᴰ - Q̇loss) - η * ṁf * cf * (T4 - T2)
        # T1 ≤ T2
        # Tin ≤ T2
        # T2 ≤ T3
        # T3 ≤ T4
        # Heat losses
        Q̇vent .== Cair * ρair * rb * Vb * (Tin - Tamb) # ventilation losses
        # Q̇cond .== sum(d[s] * U[s] * A[s] for s in 1:𝒮) * (Tin - Tamb) # conductive losses
        Q̇cond .== sum(U[s] * A[s] for s in 1:𝒮) * (Tin - Tamb) # conductive losses
        Q̇loss .== Q̇cond + Q̇vent # heat losses
        # Initial conditions
        Tin(t0) == Tin0
        # T2(t0) == 35. + 273.
        # T3(t0) == 48. + 273.
        # T4(t0) == 50. + 273.
    end);
    # setup initial guess for TESS
    loss_guess =  Cair * ρair * rb * Vb * (20 .- ambT) + sum(d[s] * U[s] * A[s] for s in 1:𝒮) * (20 .- ambT) # heat losses
    # match freq of the loss_guess to the model
    # if length(loss_guess) > length(radiation)
    #     length(loss_guess) == 360 ? loss_guess = loss_guess[1:4:192] : loss_guess = loss_guess[1:96]
    # end
    bld_netLoad = radiation .- loss_guess  # net load in the house
    set_start_value_function(model[:Q̇ₜₑₛₛᴰ], t -> bld_netLoad[zero_order_time_index(Dt, t)]);
    return model;
end;