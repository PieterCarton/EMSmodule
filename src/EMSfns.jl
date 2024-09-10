### EMS functions
# This file contains the functions used for building a EMS model object.
# The functions are mainly the device models (solar pv, bess, tess, etc.) and other elements (grid balances) of the EMS problem.

# By: Darío Slaifstein, PhD-student @TU Delft, DCES.
# Branch: agingModeling
# Version: 1.5
# Date: 01/09/2023

## Modeling functions
# These functions are used to build the EMS model object. They include the device models, the grid balances and cost function.
# We incorporate data in the form of ZOH approximations. We create an infinite t_supp and from now on we can
# use the data as a continuos function over t.
# The best way to do it would be with an indp function an changing the input data but @parameter_function only
# allows for things whose only arg is par_fn(t::Infinite()).
function zero_order_time_index(ts, t)
    inds = searchsorted(ts, t)
    if !isempty(inds) || inds.start == 1
        return inds.start
    elseif inds.start > length(ts)
        return length(ts)
    else # we don't have an exact match so let's choose the earlier time index
        return inds.stop
    end
end

# Electric devices
"""
    spv(model, data)

The PV panel is composed of a measurement (MPPT power), actual P out, and rolling forecast.
The function translates the measurement into a forecast and adds the forecast to the model.
# Arguments
- model::InfiniteModel: the model object containing the EMS problem.
- data::Dict: the data dictionary containing the Exogenous information (MPPT measurement).
# Returns
- model::InfiniteModel: updated with the added deterministic forecast.
"""
function spv!(model::InfiniteModel, data::Dict) # solar pv panel
    # MPPT measurement
    t = model[:t];
    Dt = supports(t);
    t0 = supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(supports(t))-1;
    
    MPPTmeas = data["SPV"].MPPTData[it0:itend];
    # simulated forecast
    MPPTnoisy = zeros(size(MPPTmeas));
    MPPTnoisy[MPPTmeas .> 0] = [rand(Uniform(0.5*MPPTmeas[tt],1.1*MPPTmeas[tt])) for tt ∈ findall(MPPTmeas .> 0)]
    # Add forecast
    @parameter_function(model, PpvMPPT == (t) -> MPPTnoisy[zero_order_time_index(Dt, t)])

    #= Curtailment is left for the future
        # Extract params
        Dt = sets.dTime; # is the timestamp of the PV measurement
        PpvRated = data["SPV"].RatedPower; # rated power of the panels
        
        # MPPT measurement
        MPPTmeas = data["SPV"].MPPTData[1:length(Dt)]; 
        PpvMPPT = JuMP.Containers.DenseAxisArray(MPPTmeas, Dt);

        # Add variables
        @variable(model, 0 ≤ Ppve ≤ PpvRated, Infinite(t)); # output power
        @variable(model, 0 ≤ PpvDwn ≤ PpvRated, Infinite(t)); # down reg. for FCR
        
        # measurement constraints
        @constraint(model, [tw ∈ Dt], Ppve(tw) ≤ PpvMPPT[tw]);
        @constraint(model, [tw ∈ Dt], PpvDwn(tw) ≤ PpvMPPT[tw]);
        
        # Add power balance btwn FCR and production
        # Ppv+PpvDwn=PpvMPPT
        @constraint(model, [tw ∈ Dt], Ppve(tw) + PpvDwn(tw) == PpvMPPT[tw]);
    =#
    return model;
end;

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
function st!(model::InfiniteModel, data::Dict) # solar thermal 
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
    MPPTnoisy = zeros(size(MPPTmeas));
    MPPTnoisy[MPPTmeas .> 0] = [rand(Uniform(0.5*MPPTmeas[tt],1.1*MPPTmeas[tt])) for tt ∈ findall(MPPTmeas .> 0)]
    @parameter_function(model, Pst == (t) -> ηST*MPPTnoisy[zero_order_time_index(Dt, t)])

    #= Curtailment is left for the future
        @variable(model, 0 ≤ Pst ≤ PstRated, Infinite(t));
        @constraint(model, Pst == ηST*(model[:Ppve]+model[:PpvDwn]));
    =#
    return model;
end;

function heatpump!(model::InfiniteModel, data::Dict) # heat pump
    # The heat pump has a variable (electrical) and a subordinate finite_param (thermal)
    # Extract params
    t = model[:t]; # t0=supports(t)[1];
    PhpRated = data["HP"].RatedPower; # rated power of the heat pump
    # ηHP = data["HP"].η; # electrothermal conv efficiency $/kW
    # Phpe0=data["HP"].P0;
    
    # @variable(model, 0 ≤ Phpe ≤ PhpRated, Infinite(t), start=Phpe0);  # Electric power
    @variable(model, 0 ≤ Phpe ≤ PhpRated, Infinite(t));  # Electric power
    # @constraint(model, Phpe(t0).== Phpe0)
    # @variable(model, - PhpRated ≤ Phpe ≤ PhpRated, Infinite(t), start = Phpe0);  # Electric power
    # @variable(model, - PhpRated*ηHP ≤ Phpt ≤ PhpRated*ηHP, Infinite(t));  # Thermal power
    
    # # Dummy variables for bidirectional flow
    # @variable(model, 1e-4 ≤ PhpPos ≤ PhpRated, Infinite(t), start = Phpe0);
    # @variable(model, -PhpRated ≤ PhpNeg ≤ -1e-4, Infinite(t), start = 0);
    # # Bidirectional power flow, ensuring only export or import
    # @constraint(model, PhpPos + PhpNeg == Phpe);
    # @constraint(model, PhpPos * PhpNeg == 0);
    
    # Thermal power
    # @constraint(model, Phpt == ηHP*Phpe);
    return model;
end;

"""
    tess(model::InfiniteModel, data::Dict)

The `tess` function models a Thermal Energy Storage System (TESS).
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
    # Ptess0 = data["TESS"].P0; # Initial Power [p.u.]
    SoCtessMin = data["TESS"].SoCLim[1]; # Min State of Charge [p.u.]
    # SoCtessMax = data["TESS"].SoCLim[2]; # Max State of Charge [p.u.]
    SoCtess0 = data["TESS"].SoC0; # Initial State of Charge [p.u.]
    ηtess = data["TESS"].η; # thermal efficiency
       
    # Add variables
    @variables(model, begin
        # SoCtessMin ≤ SoCtess ≤ SoCtessMax, Infinite(t) # State of Charge
        SoCtess, Infinite(t) # State of Charge
        Ptess, Infinite(t) # Thermal power
        bPtess, Infinite(t), Bin # Binary variable for TESS power
        0 ≤ PtessPos, Infinite(t) # Ptess^+ out power
        PtessNeg ≤ 0, Infinite(t) # Ptess^- in power
    end);
    # Dummy variables for bidirectional power flow, ensuring only export or import
    @constraints(model, begin
        bPtess*PtessMin ≤ PtessNeg
        PtessNeg + PtessPos .== Ptess
        # PtessNeg * (1/ηtess) + PtessPos * ηtess .== Ptess
        PtessPos ≤ (1-bPtess)*PtessMax
        SoCtessMin ≤ SoCtess
    end);
    
    # Initial conditions
    @constraint(model, SoCtess(t0) .== SoCtess0)
    # @constraint(model, Ptess(t0).== Ptess0)

    # Bucket model
    @constraint(model, ∂.(SoCtess, t) .== -ηtess*Ptess/Qtess/3600);
    # @constraints(model, begin
    #     SoCtess(t0+(Tw+Δt)/2) .- SoCtess(t0) .≤ 0.05
    #     SoCtess(t0+(Tw+Δt)/2) + 0.05 .≤ SoCtessMax
    # end);
    return model;
end;

# other
# Grids

"""
    gridConn(model::InfiniteModel, data::Dict)

The `gridConn` function is used to define a grid connection of a power electronic interface.

# Arguments
- `model::InfiniteModel`: The simulation model.
- `data::Dict`: A dictionary containing the data for the grid connection.

# Returns
- `model::InfiniteModel`: The updated simulation model.
"""
function gridConn!(model::InfiniteModel, data::Dict) # power electronic interface
    t = model[:t];
    # Grid limits
    PgMin = data["grid"].PowerLim[1]; # Min power going out
    PgMax = data["grid"].PowerLim[2]; # Max power coming in
    ηg = data["grid"].η; # converter efficiency

    # Add variables
    @variables(model, begin
        Pg, Infinite(t)  # grid power
        bPg, Infinite(t), Bin # Binary variable for Grid power
        0 ≤ PgPos, Infinite(t) # Pg^+ out/buy power
        PgNeg ≤ 0, Infinite(t), (start=0) # Pg^- in/sell power
    end)
    
    # Dummy variables for bidirectional power flow, ensuring only export or import
    @constraints(model, begin
        PgPos ≤ (1-bPg)*PgMax
        bPg*PgMin ≤ PgNeg
        # PgNeg + PgPos .== Pg
        PgNeg * (1/ηg) + PgPos * ηg .== Pg
    end)
    return model;
end;

"""
    gridThermal(model::InfiniteModel, sets::modelSettings, data::Dict)

This function adds the thermal balance to the EMS model.

# Arguments
- `model::InfiniteModel`: The InfiniteModel object representing the grid.
- `sets::modelSettings`: The model settings.
- `data::Dict`: A dictionary containing the necessary data for the calculation.

# Returns
- `model::InfiniteModel`: The updated InfiniteModel object with the thermal power balance constraint.

"""
function gridThermal!(model::InfiniteModel, sets::modelSettings, data::Dict) # thermal grid
    # Load data
    Dt = sets.dTime;
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # make InfiniteOpt compatible Plt(t) work for arbitrary times
    loadTh = data["grid"].loadTh[it0:itend];
    # white noise to simulate a forecast
    εₗᵗʰ= randn(Int(length(loadTh)*Δt/3600)) .* 0.05 # 50W
    εₗᵗʰ = repeat(εₗᵗʰ, inner = Int(3600/Δt))
    loadTh = loadTh .+ εₗᵗʰ; # add noise
    @parameter_function(model, Plt == (t) -> loadTh[zero_order_time_index(Dt, t)])
    
    # When there´s no FCR then there´s no decision variable Ppve hence Pst can just be a @parameter_function.
    # Extract params
    ηHP=data["HP"].η; # conversion factor from Electric to thermal heat pump
    # Thermal Power balance
    model[:thBalance]=@constraint(model, model[:Pst]+model[:Phpe].*ηHP+model[:Ptess] .==  Plt);
    return model;
end;

"""
    pei(model::InfiniteModel, sets::modelSettings, data::Dict)

This function adds the power balance to the EMS model.

# Arguments
- `model::InfiniteModel`: The InfiniteModel object representing the grid.
- `sets::modelSettings`: The model settings.
- `data::Dict`: A dictionary containing the necessary data for the calculation.

# Returns
- `model::InfiniteModel`: The updated InfiniteModel object with the power balance constraint.

"""
function pei!(model::InfiniteModel, sets::modelSettings, data::Dict) # power electronic interface
    # Extract data
    Dt = sets.dTime; nEV=sets.nEV;
    t = model[:t];
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # Load data
    loadElec = data["grid"].loadE[it0:itend];
    # white noise to simulate a forecast
    εₗᵉ= randn(Int(length(loadElec)*Δt/3600)) .* 0.2 # 200W
    εₗᵉ = repeat(εₗᵉ, inner = Int(3600/Δt))
    loadElec = loadElec .+ εₗᵉ; # add noise
    # loadElec can't be negative
    loadElec[loadElec .< 0] .= 0
    @parameter_function(model, Ple == (t) -> loadElec[zero_order_time_index(Dt, t)])
    
    # Power balance DC busbar
    # If we have EVs
    if any(name.(all_variables(model)) .== "Pev[1]")
        Pev_sum = nEV != 1 ? sum(model[:γf][n].*model[:Pev][n] for n in 1:nEV) : sum(model[:γf].*model[:Pev][n] for n in 1:nEV)
    else
        Pev_sum = 0
    end

    # If we have HP
    Phpe = any(name.(all_variables(model)) .== "Phpe") ? model[:Phpe] : 0

    model[:powerBalance]=@constraint(model, model[:PpvMPPT] + model[:Pbess] + Pev_sum + model[:Pg] .== Ple + Phpe)

    return model;
end;

"""
    costFunction(model, sets::modelSettings, data::Dict)

Adds Objective function to the model. The model may have three components:
``C_{\textrm{grid}}``, ``C_{\textrm{loss}}`` and a penalty for not charging the EVs.
The function includes weights and picks how to build the objective depending on the settings.

# Arguments
- `model`: The model object representing the FLEXINet simulation.
- `sets::modelSettings`: The settings object containing the model settings.
- `data`: The data object containing the simulation data.

# Returns
- `model`: The updated model object.

"""
function costFunction!(model, sets::modelSettings, data::Dict) # Objective function
    W = sets.costWeights;
    Dt = sets.dTime;
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,t0/Δt);
    itend = it0+length(Dt)-1;

    # Note: We need to use hcat() to form the axis cause otherwise λ[:,c] broadcasts into a vector.
    # In order to stay consistent with the other DenseAxisArray we need matrices, so here hcat does that for us.
    priceBuy=hcat(data["grid"].λ[it0:itend, 1].*1e-3/3600) # convert from €/MWh to €/kWs
    priceSell=hcat(data["grid"].λ[it0:itend, 2].*1e-3/3600) # convert from €/MWh to €/kWs
    # simulated forecast
    ελ = randn(Int(length(priceBuy)*Δt/3600)) .* 20*1e-3/3600 # 20 €/MWh noise
    ελ = repeat(ελ, inner = Int(3600/Δt))
    priceBuy = priceBuy .+ ελ; # add noise
    priceSell = priceSell .+ ελ; # add noise
    @parameter_function(model, λbuy == (t) -> priceBuy[zero_order_time_index(Dt, t)])
    @parameter_function(model, λsell == (t) -> priceSell[zero_order_time_index(Dt, t)])
    # Def cost
    # totCAPEX = [sum(data[k].capex .* data[k].RatedPower' for k ∈ keys)]
    # totCAPEX = totCAPEX+data["TESS"].capex*data["TESS"].Q;

    # Grid costs
    Wgrid = W[1]; # regularization factor for grid cost. max(λ)*max(P)
    cgrid = Wgrid .* (model[:PgPos]*λbuy + model[:PgNeg]*λsell);
    
    # Aging costs CHECK
    clossbess = 1.2; # cost of lost capacity EUR/Ah
    Wloss=W[3]; # regularization factor for lost capacity
    if Wloss != 0
        # cQloss = Wloss != 0 ? ( Wloss*(model[:ilossbess]+sum(model[:ilossev][n] for n ∈ 1:sets.nEV))/3600) : 0.0;
        lossBESS = data["BESS"].GenInfo.Ns * data["BESS"].GenInfo.Np * model[:ilossbess];
        lossEV = [data["EV"][n].carBatteryPack.GenInfo.Ns * 
                data["EV"][n].carBatteryPack.GenInfo.Np *
                model[:ilossev][n] for n ∈ 1:sets.nEV]
        cQloss = ∫(( Wloss*(lossBESS+sum(lossEV[n] for n ∈ 1:sets.nEV))/3600),t);
    else
        cQloss = 0.0;
    end

    # # Penalty for excessive actions
    # xπ = Vector{GeneralVariableRef}();
    # for (~,var) ∈ enumerate([:Pbess, :Pev, :Phpe])
    #     append!(xπ, model[var]);
    # end
    # Wπ = W[5]; # penalty for excessive actions
    # cπ = Wπ*∫(sum(xπ[i].^2 for i ∈ eachindex(xπ)),t);

    # Penalty/soft constraint for TESS overcharging
    Wlims = W[4]; # penalty for TESS overcharging
    SoCtess = model[:SoCtess]; SoCtessMax = data["TESS"].SoCLim[2];
    # SoCbess = model[:SoCbess]; SoCbessMin = data["BESS"].GenInfo.SoCLim[1];
    @variable(model, auxTess ≥ 0, Infinite(t));
    @constraint(model, auxTess ≥ SoCtess - SoCtessMax);
    # @variable(model, auxBess ≥ 0, Infinite(t));
    # @constraint(model, auxBess ≥ SoCbessMin - SoCbess);

    # Define penalty for not charging
    WSoCDep = W[2]
    pDep = (isempty(model[:ϵSoC]) ? 0 : WSoCDep*sum(model[:ϵSoC][n]^2 for n ∈ eachindex(model[:ϵSoC])))

    # Define objective function
    if any(name.(all_variables(model)) .== "ilossbess")
        @objective(model, Min, ∫(cgrid,t)/sum(Dt) + pDep + Wlims*∫(auxTess,t)/sum(Dt) +
            # clossbess*∫(iloss,t)/sum(Dt))
            clossbess*cQloss/sum(Dt))
        # @objective(model, Max, -∫(cgrid,t)/sum(Dt) - pDep - Wlims*∫(auxTess + auxBess,t)/sum(Dt) -
        #     clossbess*cQloss/sum(Dt))
    else
        @objective(model, Min, ∫(cgrid,t)/sum(Dt) + pDep + Wlims*∫(auxTess,t)/sum(Dt))
        # @objective(model, Min, ∫(cgrid,t)/sum(Dt) + pDep + Wlims*∫(auxTess + auxBess,t)/sum(Dt) + cπ/sum(Dt))
    end 
    return model;
end;

## Utility functions
# These functions are used to update the EMS data dictionary with the results of the EMS model.
# They are used in the rolling horizon simulation and future works.
function update_storage_asset(stgAsset, results::Dict, key::String; typeOpt::String="MPC")
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    typeOpt == "MPC" ? shift = 1 : nothing;
    # take the length from results because its the one handled by handleInfeasible()
    typeOpt == "day-ahead" ? shift = Int(ceil(length(results["t"])/2))-1 : nothing;

    # Basic update
    stgAsset.GenInfo.SoC0 = results["SoC$key"][shift+1]
    haskey(results, "Q$key") ? stgAsset.GenInfo.initQ = results["Q$key"][shift+1] : nothing # Q is not updated when Wloss=0
    # Type-specific update
    if stgAsset.PerfParameters.type == "ECM"
        stgAsset.PerfParameters.iRn0[1] = results["iR1$key"][shift+1] # states
    elseif stgAsset.PerfParameters.type == "PBROM"
        stgAsset.PerfParameters.x0 = results["x$key"][shift+1] # states
    end
    return stgAsset;
end;

function update_measurements(results::Dict, s::modelSettings, data::Dict; typeOpt::String="MPC")
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    typeOpt == "MPC" ? shift = 1 : nothing;
    # take the length from results because its the one handled by handleInfeasible()
    typeOpt == "day-ahead" ? shift = Int(ceil(length(results["t"])/2))-1 : nothing;

    # Update TESS
    data["TESS"].SoC0 = results["SoCtess"][shift+1]

    # Update storage assets (BESS and EV)
    data["BESS"]=update_storage_asset(data["BESS"], results, "bess"; typeOpt=typeOpt);
    for n in s.nEV
        data["EV"][n].carBatteryPack=update_storage_asset(data["EV"][n].carBatteryPack, results, "ev[$n]"; typeOpt=typeOpt)
    end
    return data
end;

function set_warm_start(model::InfiniteModel, preRes::Dict)
# Warm starts the model variables with the results from the previous optimization.
    t = model[:t];
    # [set_start_value_function(x[i], t -> x_opt_interp[i](t)) for i ∈ 1:2]
    # x = all_variables(model)
    for x ∈ all_variables(model)
        if name(x) != "" # skip point variables
            #= # take the previous day and fill the rest with 0s
                x_opt = preRes[name(x)][97:end]
                println(length(x_opt))
                append!(x_opt, zeros(96))
            =#
            haskey(preRes, name(x)) ? x_opt = preRes[name(x)] : x_opt = zeros(length(supports(t)))
            length(x_opt) == length(supports(t)) ? nothing : append!(x_opt, zeros(length(supports(t))-length(x_opt)))
            # maybe change to just repeating the same value instead of 0s
            x_opt_interp = linear_interpolation(value.(t), x_opt, extrapolation_bc = Line())
            set_start_value_function(x, t -> x_opt_interp(t))
        end
    end
    return model
end

"""
    getResults(model::InfiniteModel)

Get the results of the optimized InfiniteModel.

# Arguments
- `model::InfiniteModel`: The InfiniteModel object representing the simulation model.

# Returns
- `xdict::Dict{Any, Any}`: A dictionary containing the names and values of the simulation results, including the simulation time, EV availability, and solver status.
"""
function getResults(model::InfiniteModel)
    x = all_variables(model)
    xdf=DataFrame(
        name = name.(x),
        value = value.(x),
    )
    idx=findall(x->x=="",name.(x)) # find point variables
    delete!(xdf, idx) # delete point variables
    xdict=Dict{Any,Any}(Pair.(xdf.name, xdf.value))
    # add simulation time, EV availability and status of the solver
    merge!(xdict, Dict("t"=>supports(model[:t]),
                        # Forecasts - simulated for now
                        "γf"=>value.(model[:γf]),   
                        "PpvMPPT"=>value.(model[:PpvMPPT]),
                        "Ple"=>value.(model[:Ple]),
                        "λsell"=>value.(model[:λsell]),
                        "λbuy"=>value.(model[:λbuy]),
                        # check if Plt is in the model
                        # "Plt"=>value.(model[:Plt]),
                        # s.season = "winter" ? "Plt"=>value.(model[:Plt]) : nothing,
                        haskey(xdict, "Phpe") ? "Plt"=>value.(model[:Plt]) : "Plt"=>zeros(size(supports(model[:t]))),
                        "status"=>termination_status(model)))
    return xdict
end
function getDuals(model::InfiniteModel)
    cs = all_constraints(model, MOI.EqualTo)
    ddict = Dict{Any,Any}(Pair.(name.(cs), dual.(cs)))
    sdpdict = Dict{Any,Any}(Pair.(name.(cs), shadow_price.(cs)))
    return ddict, sdpdict
end

"""
    concatResultsRH(results::Vector{Dict}; typeOpt::String="MPC")

Concatenates the results from a simulation into a single dictionary, with the option to choose the type of optimization.

## Arguments
- `results::Vector{Dict}`: A vector of dictionaries containing the simulation results.
- `typeOpt::String="MPC"`: The type of optimization. Valid options are "MPC" and "day-ahead".

## Returns
- `RHdict::Dict`: A dictionary containing the concatenated results.

"""

function concatResultsRH(results::Vector{Dict}; typeOpt::String="MPC")
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    
    keys_list = keys(results[1])
    steps = length(results)
    RHdict = Dict()
 
    if typeOpt == "MPC"
        shift = 1
        for key in keys_list
            # Check if the key requires special handling and skip the "status" key
            if key == :"status" || key == :"compTime" 
                continue
            elseif key == :"γf"
                # check nEV (number of EVs) to see if we need to handle the γf differently
                if size(results[1]["γf"],1) == length(results[1]["t"]) 
                    # only one EV
                    # γf = [results[st][key][shift+1] for st in 1:steps];
                    γf = [results[st][key][shift] for st in 1:steps];
                else
                    nEV = size(results[1]["γf"],1)
                    # γf = [[results[st][key][n][shift+1] for n ∈ 1:nEV] for st in 1:steps];
                    γf = [[results[st][key][n][shift] for n ∈ 1:nEV] for st in 1:steps];
                    # γf = [[results[st][key][1][2], results[st][key][2][2]] for st in 1:steps];
                end
                RHdict[key] = hcat(γf...)
            else
                # For other keys, use the original approach
                # hardcoding the shift for now just in case
                if length(results[1][key]) != 1
                    RHdict[key] = [results[st][key][shift+1] for st in 1:steps]
                else
                    RHdict[key] = [results[st][key][shift] for st in 1:steps]
                end
                # RHdict[key] = [results[st][key][shift+1] for st in 1:steps]
                # RHdict[key] = [results[st][key][shift] for st in 1:steps]
            end
        end
    elseif typeOpt == "day-ahead"
        # take the length from results because its the one handled by handleInfeasible()
        # shift = Int(ceil(length(results[1]["t"])/2))-1
        Δt = (results[1]["t"][2] - results[1]["t"][1])/3600;
        shift = Int((24 / Δt)-1);
        for key in keys_list
            # Check if the key requires special handling and skip the "status" key
            if key == :"status" || key == :"compTime" 
                continue
            elseif key == :"γf"
                #=
                # γf = [[results[st][key][1][1:shift+1], results[st][key][2][1:shift+1]] for st in 1:steps];
                # RHdict[key] = hcat(γf...)
                =#
                if size(results[1]["γf"],1) == length(results[1]["t"]) 
                    # only one EV
                    # γf = [results[st][key][shift+1] for st in 1:steps];
                    γf = vcat([results[st][key][1:shift+1] for st in 1:steps]...);
                else
                    nEV = size(results[1]["γf"],1)
                    γf = [vcat([results[st][key][n][1:shift+1] for st ∈ 1:steps]...) for n ∈ 1:nEV]
                end
                RHdict[key] = γf
            else
                # For other keys, use the original approach
                RHdict[key] = vcat([results[st][key][1:shift+1] for st in 1:steps]...);
            end
        end
    end
    return RHdict
end