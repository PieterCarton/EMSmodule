### EMS functions
# This file contains the functions used for building a EMS model object.
# The functions are mainly the device models (solar pv, bess, tess, etc.) and other elements (grid balances) of the EMS problem.

# By: Darío Slaifstein, PhD-student @TU Delft, DCES.
# Branch: main
# Version: 1.6
# Date: 16/09/2025

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
function spv!(model::InfiniteModel, data::Dict; add_noise::Bool = true) # solar pv panel
    # MPPT measurement
    t = model[:t];
    Dt = supports(t);
    t0 = supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(supports(t))-1;
    MPPTmeas = data["SPV"].MPPTData[it0:itend];
    MPPTnoisy = copy(MPPTmeas);
    if add_noise
        # simulated forecast
        MPPTnoisy[MPPTmeas .> 0] = [rand(Uniform(0.5*MPPTmeas[tt],1.1*MPPTmeas[tt])) for tt ∈ findall(MPPTmeas .> 0)]
    end
    # Add forecast
    @parameter_function(model, PpvMPPT == (t) -> MPPTnoisy[zero_order_time_index(Dt, t)])
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
        0 ≤ PgNeg ≤ -PgMin, Infinite(t)
        0 ≤ PgPos ≤ PgMax, Infinite(t) # Pg^+ out/buy power
    end)
    
    # Dummy variables for bidirectional power flow, ensuring only export or import
    # PgNeg + PgPos .== Pg
    @expression(model, Pg, PgPos * ηg .- PgNeg * (1/ηg))
    @constraint(model, PgNeg ⟂ PgPos) # MPEC with ⟂
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
function pei!(model::InfiniteModel, sets::modelSettings, data::Dict; add_noise::Bool = true) # power electronic interface
    # Extract data
    Dt = sets.dTime; nEV=sets.nEV;
    t = model[:t];
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # Load data
    loadElec = data["grid"].loadE[it0:itend];
    MPPTmeas = data["SPV"].MPPTData[it0:itend];
    if add_noise
        # white noise to simulate a forecast
        εₗᵉ= randn(Int(length(loadElec)*Δt/3600)) .* 0.2 # 200W
        εₗᵉ = repeat(εₗᵉ, inner = Int(3600/Δt))
        loadElec = loadElec .+ εₗᵉ; # add noise
        # loadElec can't be negative
        loadElec[loadElec .< 0] .= 0
    end
    @parameter_function(model, Ple == (t) -> loadElec[zero_order_time_index(Dt, t)])
    
    # Power balance DC busbar
    # If we have EVs
    if any(name.(all_variables(model)) .== "PevPos[1]")
        Pev_sum = nEV != 1 ? sum(model[:γf][n].*model[:Pev][n] for n in 1:nEV) : sum(model[:γf].*model[:Pev][n] for n in 1:nEV)
    else
        Pev_sum = 0
    end

    # If we have HP
    Phpe = any(name.(all_variables(model)) .== "Phpe") ? model[:Phpe] : 0
    # If we have a BESS
    Pbess = any(name.(all_variables(model)) .== "PbessPos") ? model[:Pbess] : 0

    model[:powerBalance]=@constraint(model, model[:PpvMPPT] + Pbess + Pev_sum + model[:Pg] .== Ple + Phpe)
    resLoad = loadElec .- MPPTmeas;
    resLoad⁺ = copy(resLoad); resLoad⁻ = copy(resLoad);
    resLoad⁺[resLoad .< 0] .= 0; resLoad⁻[resLoad .> 0] .= 0;
    set_start_value_function(model[:PgPos], t -> resLoad⁺[zero_order_time_index(Dt, t)])
    set_start_value_function(model[:PgNeg], t -> -resLoad⁻[zero_order_time_index(Dt, t)])
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
function costFunction!(model, sets::modelSettings, data::Dict; add_noise::Bool = true) # Objective function
    W = sets.costWeights;
    Dt = sets.dTime;
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,t0/Δt);
    itend = it0+length(Dt)-1;

    # Note: We need to use hcat() to form the axis cause otherwise λ[:,c] broadcasts into a vector.
    # In order to stay consistent with the other DenseAxisArray we need matrices, so here hcat does that for us.
    priceBuy=hcat(data["grid"].λ[it0:itend, 1]) # in [c€/kWs]
    priceSell=hcat(data["grid"].λ[it0:itend, 2]) # in [c€/kWs]
    if add_noise
        # simulated forecast
        ελ = randn(Int(length(priceBuy)*Δt/3600)) .* 20*1e-3/3600 # 20 €/MWh noise
        ελ = repeat(ελ, inner = Int(3600/Δt))
        priceBuy = priceBuy .+ ελ; # add noise
        priceSell = priceSell .+ ελ; # add noise
    end
    @parameter_function(model, λbuy == (t) -> priceBuy[zero_order_time_index(Dt, t)])
    @parameter_function(model, λsell == (t) -> priceSell[zero_order_time_index(Dt, t)])
    # Def cost
    # totCAPEX = [sum(data[k].capex .* data[k].RatedPower' for k ∈ keys)]
    # totCAPEX = totCAPEX+data["TESS"].capex*data["TESS"].Q;

    # Grid costs
    Wgrid = W[1]; # regularization factor for grid cost. max(λ)*max(P)
    cgrid = Wgrid .* (model[:PgPos]*λbuy - model[:PgNeg]*λsell);
    # cgridᴰᴬ = Wgrid .* ((λbuyᴰᴬ-λsellᴰᴬ)/2 * (-Pgᴰᴬ*bPg + Pgᴰᴬ*(1 - bPg)) +
    #                  (λbuyᴰᴬ+λsellᴰᴬ)/2 * Pgᴰᴬ)
    # cgridᶜᵀ = Wgrid .* ((λbuyᶜᵀ-λsellᶜᵀ)/2 * abs(Pgᶜᵀ - Pgfᴰᴬ)) +
    #                  (λbuyᶜᵀ+λsellᶜᵀ)/2 * (Pgᶜᵀ-Pgfᴰᴬ))

    # Aging costs CHECK
    clossbess = 1.2; # cost of lost capacity EUR/Ah
    clossbess /= 3600; # cost of lost capacity EUR/As
    Wloss=W[3]; # regularization factor for lost capacity
    if Wloss != 0.
        # cQloss = Wloss != 0 ? ( Wloss*(model[:ilossbess]+sum(model[:ilossev][n] for n ∈ 1:sets.nEV))/3600) : 0.0;
        lossBESS = data["BESS"].GenInfo.Ns * data["BESS"].GenInfo.Np * model[:ilossbess];
        if any(name.(all_variables(model)) .== "ilossev[1]")
            lossEV = [data["EV"][n].carBatteryPack.GenInfo.Ns * 
                    data["EV"][n].carBatteryPack.GenInfo.Np *
                    model[:ilossev][n] for n ∈ 1:sets.nEV]
        else
            lossEV = zeros(sets.nEV);
        end
        cQloss = ∫(( Wloss*(lossBESS+sum(lossEV[n] for n ∈ 1:sets.nEV))),t);
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
    @variable(model, auxTess ≥ 0., Infinite(t));
    @constraint(model, auxTess ≥ SoCtess - SoCtessMax);
    # @variable(model, auxBess ≥ 0, Infinite(t));
    # @constraint(model, auxBess ≥ SoCbessMin - SoCbess);

    # Define penalty for not charging
    WSoCDep = W[2]
    pDep = (any(name.(all_variables(model)) .== "PevPos[1]") ?
            WSoCDep*sum(model[:ϵSoC][n]^2 for n ∈ eachindex(model[:ϵSoC])) : 0);

    # Define objective function
    if any(name.(all_variables(model)) .== "ilossbess")
        @objective(model, Min, ∫(cgrid,t) +
            Wlims*∫(auxTess,t) +
            cQloss + 
            pDep)
        # @objective(model, Min, ∫(cgrid,t)/sum(Dt) + pDep + Wlims*∫(auxTess,t)/sum(Dt) +
        #     clossbess*cQloss/sum(Dt))
        # @objective(model, Max, -∫(cgrid,t)/sum(Dt) - pDep - Wlims*∫(auxTess + auxBess,t)/sum(Dt) -
        #     clossbess*cQloss/sum(Dt))
    else
        @objective(model, Min, ∫(cgrid,t) +
            pDep +
            Wlims*∫(auxTess,t))
        # @objective(model, Min, ∫(cgrid,t)/sum(Dt) + pDep + Wlims*∫(auxTess,t)/sum(Dt))
        # @objective(model, Min, ∫(cgrid,t)/sum(Dt) + pDep + Wlims*∫(auxTess + auxBess,t)/sum(Dt) + cπ/sum(Dt))
    end 
    return model;
end;

"""
    costFunctionDA!(model, sets::modelSettings, data::Dict)

Adds Objective function to the day-ahead model. The model may have three components:
``C_{\textrm{grid}}``, ``C_{\textrm{loss}}`` and a penalty for not charging the EVs.
The function includes weights and picks how to build the objective depending on the settings.

# Arguments
- `model`: The model object representing the FLEXINet simulation.
- `sets::modelSettings`: The settings object containing the model settings.
- `data`: The data object containing the simulation data.

# Returns
- `model`: The updated model object.

"""
function costFunctionDA!(model, sets::modelSettings, data::Dict; add_noise::Bool = true) # Objective function
    W = sets.costWeights; t = model[:t];
    Dt = supports(t);
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,t0/Δt);
    itend = it0+length(Dt)-1;

    # Note: We need to use hcat() to form the axis cause otherwise λ[:,c] broadcasts into a vector.
    # In order to stay consistent with the other DenseAxisArray we need matrices, so here hcat does that for us.
    priceBuyDA=hcat(data["grid"].λ[it0:itend, 1]) # convert from €/MWh to €/kWs
    priceSellDA=hcat(data["grid"].λ[it0:itend, 2]) # convert from €/MWh to €/kWs
    
    if add_noise
        # simulated forecast
        ελᴰᴬ = randn(Int(length(priceBuyDA)*Δt/3600)) .* 20*1e-3/3600 # 20 €/MWh noise
        ελᴰᴬ = repeat(ελᴰᴬ, inner = Int(3600/Δt))
        priceBuyDA = priceBuyDA .+ ελᴰᴬ; # add noise
        priceSellDA = priceSellDA .+ ελᴰᴬ; # add noise
        # ελᶜᵀ = randn(Int(length(priceBuyCT)*Δt/3600)) .* 20*1e-3/3600 # 20 €/MWh noise
        # ελᶜᵀ = repeat(ελᶜᵀ, inner = Int(3600/Δt))
        # priceBuyCT = priceBuyCT .+ ελᶜᵀ; # add noise
        # priceSellCT = priceSellCT .+ ελᶜᵀ; # add noise
    end
    # @parameter_function(model, λbuyᴰᴬ == (t) -> priceBuyDA[zero_order_time_index(Dt, t)])
    # @parameter_function(model, λsellᴰᴬ == (t) -> priceSellDA[zero_order_time_index(Dt, t)])
    @parameter_function(model, λbuy == (t) -> priceBuyDA[zero_order_time_index(Dt, t)])
    @parameter_function(model, λsell == (t) -> priceSellDA[zero_order_time_index(Dt, t)])
    # Grid costs
    Wgrid = W[1]; # regularization factor for grid cost. max(λ)*max(P)
    cgridᴰᴬ = Wgrid .* (λbuy - λsell)/2 * (model[:PgPos] .+ model[:PgNeg]) + (λbuy + λsell)/2 * model[:Pg];
    # cgridᴰᴬ = Wgrid .* (model[:PgPos]*λbuy .- model[:PgNeg]*λsell);
    
    # Aging costs CHECK
    closs = 1.2; # cost of lost capacity EUR/Ah
    closs /= 3600 # EUR/As
    Wloss=W[3]; # regularization factor for lost capacity
    if Wloss != 0.
        lossBESS = data["BESS"].GenInfo.Ns * data["BESS"].GenInfo.Np * model[:ilossbess];
        if any(name.(all_variables(model)) .== "ilossev[1]")
            lossEV = [data["EV"][n].carBatteryPack.GenInfo.Ns * 
                    data["EV"][n].carBatteryPack.GenInfo.Np *
                    model[:ilossev][n] for n ∈ 1:sets.nEV]
        else
            lossEV = zeros(sets.nEV);
        end
        cQloss = ∫(Wloss*(lossBESS+sum(lossEV[n] for n ∈ 1:sets.nEV)),t);
    else
        cQloss = 0.0;
    end

    # Penalty/soft constraint for TESS overcharging
    Wlims = W[4]; # penalty for TESS overcharging
    Tin = model[:Tin]; TinMax = 24 + 273; TinMin = 17 + 273;
    occ = data["Building"].occupancy[it0:itend]; 
    @parameter_function(model, 𝒪 == (t) -> occ[zero_order_time_index(Dt, t)])
    @variable(model, sTin ≥ 0., Infinite(t));
    @constraint(model, sTin ≥ Tin - TinMax);
    @constraint(model, sTin ≥ TinMin - Tin);
    # Define penalty for not charging
    WSoCDep = W[2]
    pDep = (any(name.(all_variables(model)) .== "PevPos[1]") ?
            WSoCDep*sum(model[:ϵSoC][n]^2 for n ∈ eachindex(model[:ϵSoC])) : 0);

    # Define objective function
    if any(name.(all_variables(model)) .== "ilossbess")
        # @objective(model, Min, ∫(cgridᴰᴬ,t)/Dt[end] +
        #     pDep +
        #     Wlims*∫(sTin .* 𝒪,t)/Dt[end] +
        #     closs*cQloss/Dt[end])
            @objective(model, Min, ∫(cgridᴰᴬ,t) +
            pDep +
            Wlims*∫(sTin .* 𝒪,t) +
            closs*cQloss)
    else
        @objective(model, Min, ∫(cgridᴰᴬ,t)/Dt[end] +
            pDep +
            Wlims*∫(sTin .* 𝒪,t)/Dt[end])
    end 
    return model;
end;

"""
    costFunctionCT!(model, sets::modelSettings, data::Dict)

Adds Objective function to the continous-time/intra-day model. The model may have three components:
``C_{\textrm{grid}}``, ``C_{\textrm{loss}}`` and a penalty for not charging the EVs.
The function includes weights and picks how to build the objective depending on the settings.

# Arguments
- `model`: The model object representing the FLEXINet simulation.
- `sets::modelSettings`: The settings object containing the model settings.
- `data`: The data object containing the simulation data.

# Returns
- `model`: The updated model object.

"""
function costFunctionCT!(model, sets::modelSettings, data::Dict, sol_DA; add_noise::Bool = true) # Objective function
    W = sets.costWeights;
    Dt = sets.dTime;
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,t0/Δt);
    itend = it0+length(Dt)-1;

    # Note: We need to use hcat() to form the axis cause otherwise λ[:,c] broadcasts into a vector.
    # In order to stay consistent with the other DenseAxisArray we need matrices, so here hcat does that for us.
    priceBuyCT=hcat(data["grid"].λ[it0:itend, 1]) # convert from €/MWh to €/kWs
    priceSellCT=hcat(data["grid"].λ[it0:itend, 2]) # convert from €/MWh to €/kWs
    if add_noise
        # simulated forecast
        ελ = randn(Int(length(priceBuyCT)*Δt/3600)) .* 20*1e-3/3600 # 20 €/MWh noise
        ελ = repeat(ελ, inner = Int(3600/Δt))
        priceBuyCT = priceBuyCT .+ ελ; # add noise
        priceSellCT = priceSellCT .+ ελ; # add noise
    end
    @parameter_function(model, λbuyᶜᵀ == (t) -> priceBuyCT[zero_order_time_index(Dt, t)])
    @parameter_function(model, λsellᶜᵀ == (t) -> priceSellCT[zero_order_time_index(Dt, t)])
    # reference grid power
    refPg = sol_DA["Pg"][1:Int(length(sol_DA["Pg"])/2)];
    Dt_DA = sol_DA["t"]; # discrete time of the day-ahead solution
    Δt_DA = Dt_DA[2]-Dt_DA[1];
    refPg=repeat(refPg, inner = Int(Δt_DA/Δt))
    @parameter_function(model, Pgfᴰᴬ == (t) -> refPg[zero_order_time_index(Dt, t)])
    # Pgfᴰᴬ = model[:Pgfᴰᴬ];

    # Grid costs
    Wgrid = W[1]; # regularization factor for grid cost. max(λ)*max(P)
    # cgridᶜᵀ = Wgrid .* ((λbuyᶜᵀ-λsellᶜᵀ)/2 * abs(Pgᶜᵀ - Pgfᴰᴬ) +
    #                  (λbuyᶜᵀ+λsellᶜᵀ)/2 * (Pgᶜᵀ-Pgfᴰᴬ))
    Pg = model[:Pg];
    cgridᶜᵀ = Wgrid .* ((λbuyᶜᵀ-λsellᶜᵀ)/2 * abs(Pg - Pgfᴰᴬ) +
                     (λbuyᶜᵀ+λsellᶜᵀ)/2 * (Pg-Pgfᴰᴬ))
    # cgridᶜᵀ = Wgrid .* ((λbuyᶜᵀ-λsellᶜᵀ)/2 * (model[:Pgᶜᵀ⁺]+model[:Pgᶜᵀ⁻]) +
    #                  (λbuyᶜᵀ+λsellᶜᵀ)/2 * (model[:Pgᶜᵀ]))
    # for the predictive DA+CT model
    # cgridᶜᵀ = Wgrid .* (λbuyᶜᵀ*model[:Pgᶜᵀ⁺] .- λsellᶜᵀ*model[:Pgᶜᵀ⁻]);
        
    # Aging costs CHECK
    closs = 1.2; # cost of lost capacity EUR/Ah
    closs /= 3600 # EUR/As
    Wloss=W[3]; # regularization factor for lost capacity
    if Wloss != 0.
        # cQloss = Wloss != 0 ? ( Wloss*(model[:ilossbess]+sum(model[:ilossev][n] for n ∈ 1:sets.nEV))/3600) : 0.0;
        lossBESS = data["BESS"].GenInfo.Ns * data["BESS"].GenInfo.Np * model[:ilossbess];
        if any(name.(all_variables(model)) .== "ilossev[1]")
            lossEV = [data["EV"][n].carBatteryPack.GenInfo.Ns * 
                    data["EV"][n].carBatteryPack.GenInfo.Np *
                    model[:ilossev][n] for n ∈ 1:sets.nEV]
        else
            lossEV = zeros(sets.nEV);
        end
        cQloss = ∫(Wloss*(lossBESS+sum(lossEV[n] for n ∈ 1:sets.nEV)),t);
    else
        cQloss = 0.0;
    end

    # Penalty/soft constraint for TESS overcharging
    Wlims = W[4]; # penalty for thermal comfort
    Tin = model[:Tin]; TinMax = 24 + 273; TinMin = 17 + 273;
    occ = data["Building"].occupancy[it0:itend]; 
    @parameter_function(model, 𝒪 == (t) -> occ[zero_order_time_index(Dt, t)])
    @variable(model, sTin ≥ 0., Infinite(t));
    @constraint(model, sTin ≥ Tin - TinMax);
    @constraint(model, sTin ≥ TinMin - Tin);

    # Define penalty for not charging
    WSoCDep = W[2]
    pDep = (any(name.(all_variables(model)) .== "PevPos[1]") ?
            WSoCDep*sum(model[:ϵSoC][n]^2 for n ∈ eachindex(model[:ϵSoC])) : 0);

    # Define objective function
    if any(name.(all_variables(model)) .== "ilossbess")
        @objective(model, Min, ∫(cgridᶜᵀ,t) +
            closs*cQloss +
            pDep +
            Wlims*∫(sTin .* 𝒪,t)
            )
    else
        @objective(model, Min, ∫(cgridᶜᵀ,t)/Dt[end] +
            pDep +
            Wlims*∫(sTin .* 𝒪,t)/Dt[end]
            )
    end 
    return model;
end;

## Utility functions
# These functions are used to update the EMS data dictionary with the results of the EMS model.
# They are used in the rolling horizon simulation and future works.
function update_storage_asset(stgAsset::BESSData, results::Dict, key::String; typeOpt::String="MPC")
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
            # length(x_opt) == length(supports(t)) ? nothing : append!(x_opt, zeros(length(supports(t))-length(x_opt)))
            if length(x_opt) != length(supports(t))
                if length(x_opt) > length(supports(t))
                    x_opt = x_opt[1:length(supports(t))]
                else # length(x_opt) < length(supports(t))
                    append!(x_opt, zeros(length(supports(t))-length(x_opt)))
                end
            end
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
                        haskey(xdict, "Pev[1]") ? "γf"=>value.(model[:γf]) : "γf"=>zeros(size(supports(model[:t]))),   
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

# Dispatch based on object type to fill out the data dictionary
function add_to_dict(data, sym, obj::GeneralVariableRef)
    if typeof(dispatch_variable_ref(obj)) <: Union{InfiniteVariableRef, 
        SemiInfiniteVariableRef, 
        PointVariableRef, 
        FiniteVariableRef,
        IndependentParameterRef,
        ParameterFunctionRef,
        }
        data[string(sym)] = value(obj)
    else
        add_to_dict(data, sym, dispatch_variable_ref(obj))  # Ensure correct dispatch
    end
end

function add_to_dict(data, sym, obj::AbstractJuMPScalar)
    data[string(sym)] = value(obj)  # Expressions don't have names
end

function add_to_dict(data, sym, obj::AbstractArray)
    for i in eachindex(obj)
        index_name = string(sym, "[", join(Tuple(i), ", "), "]")  # Ensure proper indexing
        add_to_dict(data, index_name, obj[i])
    end
end

function add_to_dict(data, sym, obj)
    return  # Fallback: do nothing
end

function getResultsFull(model)
    xdict = Dict()
    for (sym, obj) in object_dictionary(model)
        add_to_dict(xdict, sym, obj)
    end
    merge!(xdict, Dict("status"=>termination_status(model),
        # Empty vars for extra states
        "COPhatᴰ" => zeros(size(supports(model[:t]))),
        "COPhatᵗᵉˢˢ" => zeros(size(supports(model[:t]))),
        "T2" => zeros(size(supports(model[:t]))),
        "T3" => zeros(size(supports(model[:t]))),
        "T4" => zeros(size(supports(model[:t]))),
    ))
    return xdict
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
                    γf = [results[st][key][shift+1] for st in 1:steps];
                    # γf = [results[st][key][shift] for st in 1:steps];
                else
                    nEV = size(results[1]["γf"],1)
                    γf = [[results[st][key][n][shift+1] for n ∈ 1:nEV] for st in 1:steps];
                    # γf = [[results[st][key][n][shift] for n ∈ 1:nEV] for st in 1:steps];
                    # γf = [[results[st][key][1][2], results[st][key][2][2]] for st in 1:steps];
                end
                RHdict[key] = hcat(γf...)
            else
                # For other keys, use the original approach
                # hardcoding the shift for now just in case
                RHdict[key] = zeros(steps)
                for st ∈ 1:steps
                    if length(results[st][key]) != 1
                        RHdict[key][st] = results[st][key][shift+1]
                    else
                        RHdict[key][st] = results[st][key][shift]
                    end
                end
                # if length(results[1][key]) != 1
                #     RHdict[key] = [results[st][key][shift+1] for st in 1:steps]
                # else
                #     RHdict[key] = [results[st][key][shift] for st in 1:steps]
                # end
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
                try
                    RHdict[key] = vcat([results[st][key][1:shift+1] for st in 1:steps]...);
                catch e
                    continue
                end
            end
        end
    end
    return RHdict
end

function  concatResultsFlex(results::Vector{Dict}, typeOpt::CT_MPC)
    keys_list_DA = keys(results[1]);
    keys_list_CT = keys(results[2]);
    # outer join the two keys lists
    keys_list = Set(keys_list_DA) ∪ Set(keys_list_CT)

    steps = length(results)
    RHdict = Dict()
 
    shift = 1
    for key in keys_list
        # Check if the key requires special handling and skip the "status" key
        if key == :"status" || key == :"compTime"
            try
                RHdict[key] = [results[st][key] for st in 1:steps]
            catch e
                continue
            end
        elseif key == :"γf"
            # check nEV (number of EVs) to see if we need to handle the γf differently
            if size(results[1]["γf"],1) == length(results[1]["t"]) 
                # only one EV
                γf = [results[st][key][shift+1] for st in 1:steps];
            else
                nEV = size(results[1]["γf"],1)
                γf = [[results[st][key][n][shift+1] for n ∈ 1:nEV] for st in 1:steps];
            end
            RHdict[key] = hcat(γf...)
        else
            # For other keys, use the original approach
            # hardcoding the shift for now just in case
            RHdict[key] = zeros(steps)
            for st ∈ 1:steps
                try
                    if length(results[st][key]) != 1
                        RHdict[key][st] = results[st][key][shift+1]
                    else
                        RHdict[key][st] = results[st][key][shift]
                    end    
                catch e
                    continue
                end
            end
        end
    end
    return RHdict
end

function  concatResultsFlex(results::Vector{Dict}, typeOpt::DayAhead)
    keys_list_DA = keys(results[1]);
    keys_list_CT = keys(results[2]);
    # outer join the two keys lists
    keys_list = Set(keys_list_DA) ∪ Set(keys_list_CT)

    steps = length(results)
    RHdict = Dict()
 
    # take the length from results because its the one handled by handleInfeasible()
    # shift = Int(ceil(length(results[1]["t"])/2))-1
    Δt = (results[1]["t"][2] - results[1]["t"][1])/3600;
    shift = Int((24 / Δt)-1);
    for key in keys_list
        # Check if the key requires special handling and skip the "status" key
        if key == :"status" || key == :"compTime" 
            RHdict[key] = [results[st][key] for st in 1:steps]
        elseif key == :"γf"
            if size(results[1]["γf"],1) == length(results[1]["t"]) 
                γf = vcat([results[st][key][1:shift+1] for st in 1:steps]...);
            else
                nEV = size(results[1]["γf"],1)
                γf = [vcat([results[st][key][n][1:shift+1] for st ∈ 1:steps]...) for n ∈ 1:nEV]
            end
            RHdict[key] = γf
        else
            # For other keys, use the original approach
            try
                RHdict[key] = vcat([results[st][key][1:shift+1] for st in 1:steps]...);
            catch e
                continue
            end
        end
    end
    return RHdict
end