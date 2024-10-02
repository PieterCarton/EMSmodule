### sim EMS
# This file contains the functions used for simulating a Multicarrier Energy System for a given power setpoint.

# By: Darío Slaifstein, PhD-student @TU Delft, DCES.
# Branch: ageingModelling_v2
# Version: 1.5
# Date: 26/09/2024

function perfModel_matching(stgAsset::BESSData)
    if typeof(stgAsset.PerfParameters) == CIDRAPBROMPerfParams
        perfModel = stgAsset.PerfParameters
    else
        stgAsset.cellID == "SYNSANYO" ? perfModel = CIDRAPBROMPerfParams() : perfModel = CIDRAPBROMPerfParams(Cell=Construct(stgAsset.cellID));
        # perfModel = CIDRAPBROMPerfParams(Cell=Construct(stgAsset.cellID));
        perfModel.Cell.RA.Fs = 1 # Modify transfer function sampling frequency
        perfModel.Cell.RA.SamplingT = 1 # Modify final system sampling period
        perfModel.Cell.Const.T = 298.15 # Modify cell temperature
        # perfModel.Cell.Const.Debug = true # Modify cell printing level
        # modify the default parameters to match the ones in the stgAsset
        # SoHQ = Qₛₐₜ/Qₛₐ₀
        # θ_100 = θ_100 * SoHQ
        perfModel.Cell.Neg.θ_100 = perfModel.Cell.Neg.θ_100 * stgAsset.GenInfo.SoHQ;
        # RFilm = RFilm + ΔR0
        perfModel.Cell.Neg.RFilm = perfModel.Cell.Neg.RFilm + copy(stgAsset.GenInfo.SoHR0);
    end
    return perfModel;
end

function agingModel_matching(stgAsset::BESSData, perfModel::CIDRAPBROMPerfParams)
    if typeof(stgAsset.AgingParameters) == JinAgingParams
        agingModel = stgAsset.AgingParameters
    else
        # agingModel = JinAgingParams() # create submodel
        agingModel = JinAgingParams(
                                    # OCVn = perfModel.Cell.Const.Uocp, # <--- we have to extract the right function from the struct
                                    # OCVn = θ -> perfModel.Cell.Const.Uocp("Neg", θ); CHECK
                                    # λ = λ, check
                                    # εAM = perfModel.Cell.Neg.ϵ_s, check
                                    Rs = perfModel.Cell.Neg.Rs,
                                    An = perfModel.Cell.Const.CC_A,
                                    Ln = perfModel.Cell.Neg.L,
                                    z100p = perfModel.Cell.Neg.θ_100,
                                    z0p = perfModel.Cell.Neg.θ_0,
                                    # εₑ0 = perfModel.Cell.Neg.ϵ_e,
                                    t⁺₀ = perfModel.Cell.Const.t_plus,
                                    # DeRef=perfModel.Cell.Neg.De,
                                    ce_avg = perfModel.Cell.Const.ce0,
                                    # ce_max=perfModel.Cell.Const.ce0,
                                    σn = perfModel.Cell.Neg.σ,
                                    εₛ = perfModel.Cell.Neg.ϵ_s, # check
                                    )
        # adjust default parameters:
        # Initial conditions δSEI0, εₑ0, z100p
        # agingModel.z100p = agingModel.z100p * stgAsset.GenInfo.SoHQ; already did it in `perfModel_matching()`
        agingModel.δSEI0 = agingModel.δSEI0 * (1 + stgAsset.GenInfo.SoHR0);
        
        # parameters
        # agingModel.εAM = perfModel.Cell.Neg.ϵ_s;
        # agingModel.Rs = perfModel.Cell.Neg.Rs;
    end
    return agingModel;
end

function update_stgAsset_deg!(stgAsset::BESSData, results::Dict, key::String; typeOpt::String="MPC")
    # This function updates the storage asset with the degradation results
    
    # extract from the results dictionary the degradation results
    Qsa = copy(results["Q$key"])
    R0 = copy(results["R0$key"])
    Qsan = copy(stgAsset.GenInfo.initQ); # rated capacity
    SoHQ = copy(stgAsset.GenInfo.SoHQ); # initial SoHQ
    Qsa0 = Qsan * SoHQ; # initial capacity
    Qloss = Qsa0 .- Qsa
    δSEI = copy(results["δSEI$key"])
    εₑ = copy(results["εₑ$key"])
    if typeOpt == "MPC"
        shift = 1;
        stgAsset.GenInfo.SoHQ = copy(Qsa[shift+1]/Qsan);
        stgAsset.GenInfo.SoHR0 = copy(R0);
        if typeof(stgAsset.PerfParameters) == CIDRAPBROMPerfParams
            stgAsset.PerfParameters.Cell.Neg.θ_100 = copy(stgAsset.PerfParameters.Cell.Neg.θ_100 .- Qloss[shift+1] / Qsa0)
            stgAsset.PerfParameters.Cell.Neg.RFilm = copy(R0);
        elseif typeof(stgAsset.PerfParameters) == ECMPerfParams
            stgAsset.PerfParameters.R0Param = copy([R0]);
        end
        if typeof(stgAsset.AgingParameters) == JinAgingParams
            stgAsset.AgingParameters.z100p = copy(stgAsset.AgingParameters.z100p .- Qloss[shift+1] / Qsa0);
            stgAsset.AgingParameters.δSEI0 = copy(δSEI);
            stgAsset.AgingParameters.εₑ0 = copy(εₑ);
        end
    else # typeOpt == "day-ahead"
        shift = Int(24*3600/Δt)-1;
        stgAsset.GenInfo.SoHQ = copy(Qsa[end]/Qsan);
        stgAsset.GenInfo.SoHR0 = copy(R0[end]);
        if typeof(stgAsset.PerfParameters) == CIDRAPBROMPerfParams
            stgAsset.PerfParameters.Cell.Neg.θ_100 = copy(stgAsset.PerfParameters.Cell.Neg.θ_100 .- Qloss[end] / Qsa0)
            stgAsset.PerfParameters.Cell.Neg.RFilm = copy(R0[end]);
        elseif typeof(stgAsset.PerfParameters) == ECMPerfParams
            stgAsset.PerfParameters.R0Param = copy([R0[end]]);
        end
        if typeof(stgAsset.AgingParameters) == JinAgingParams
            stgAsset.AgingParameters.z100p = copy(stgAsset.AgingParameters.z100p .- Qloss[end] / Qsa0);
            stgAsset.AgingParameters.δSEI0 = copy(δSEI[end]);
            stgAsset.AgingParameters.εₑ0 = copy(εₑ[end]);
        end
    end
    return stgAsset;
end

# Now we calculate the aging independently for each of the storage devices
function simulate_storage_asset_deg!(stgAsset::BESSData, perfModel::CIDRAPBROMPerfParams, results::Dict, key::String; typeOpt::String="MPC")
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    typeOpt == "MPC" ? shift = 1 : nothing;
    # take the length from results because its the one handled by handleInfeasible()
    Δt = copy(results[:"t"][2]-results[:"t"][1]); # setpoint timestep
    typeOpt == "day-ahead" ? shift = Int(24*3600/Δt)-1 : nothing;

    t=results["t"][1:shift+1]; # time
    Δt=t[2]-t[1]; # time step
    # make keys for the results dictionary    
    Qsan = copy(stgAsset.GenInfo.initQ); # rated capacity
    SoHQ = copy(stgAsset.GenInfo.SoHQ); # initial SoHQ
    Qsa0 = Qsan * SoHQ; # initial capacity
    i_key = "i$key"; SoC_key = "SoC$key"; Q_key = "Q$key";
    if typeOpt == "MPC"
        # isa = copy(results[i_key][end]); # current
        # SoCsa = copy(results[SoC_key][end]); # SoC
        isa = copy(results[i_key][shift+1]); # current
        SoCsa = copy(results[SoC_key][shift+1]); # SoC
    else # typeOpt == "day-ahead"
        isa = copy(results[i_key]); # current
        SoCsa = copy(results[SoC_key]); # SoC
    end
    
    # Battery
    Tk = copy(perfModel.Cell.Const.T); # temperature in Kelvin
    Tref = copy(perfModel.Cell.Const.T_ref)
    R = 8.314; # gas constant

    # Pick ageing model
    agingModel = agingModel_matching(stgAsset, perfModel)
    @unpack_JinAgingParams agingModel # unpack submodel

    # the initial series resistance comes from perfModel
    R00 = copy(perfModel.Cell.Neg.RFilm);

    # # Empirical ageing submodel, Wang (2014)
    # c=[0.0008/3600, 0.39, 1.035, 50, 14.876/√(24*3600)]; # aging parameters
    # ilossCycle = c[1]*c[3]/c[4].*ℯ.^(c[2].*abs.(isa)).*(1 .-SoCsa).*abs.(isa).Δt; # cyclic aging
    # ilossCal = c[5]*.√Dt* ℯ^(-24e3/R/T); # calendar aging
    # iloss = ilossCycle .+ ilossCal # total aging
    
    # implement Jin (2022) for SEI and AM.
    # SEI
    ηk=2*R*Tk/F*asinh.(isa/nSEI/as/An/Ln/i0) # kinetic overpotential
    z = SoCsa*(z100p-z0p) .+ z0p
    # OCVn = 0.6379 .+ 0.5416*ℯ .^ (-305.5309 .* z) .+
    #     0.044*tanh.(-(z.-0.1958)/0.1088) .-
    #     0.1978*tanh.((z.-1.0571)/0.0854) .-
    #     0.6875*tanh.((z.+0.0117)/0.0529) .-
    #     0.0175*tanh.((z.-0.5692)/0.0875)
    κ = κref*exp(Eκ/R*(1/Tref-1/Tk))
    # De = DeRef*exp(EDe/R*(1/Tref-1/T))
    
    # DeEff = De*εₑ .^ brug
    # κDeff = 2*R*T*κeff*(t⁺₀-1)/F
    # σeff = σn*εₛ^brug 
    # Et = -(isa/κeff/An/Ln)*((β*εₑ-(1-t⁺₀))*κDeff/ce_avg/DeEff/F+1)
    # θ = ℯ .^ (nSEI*F/R/Tk*(ηk .+ OCVn("Neg", z) .- OCVs))
    θ = ℯ .^ (nSEI*F/R/Tk*(ηk .+ OCVn(z) .- OCVs)) # fitting param
    if typeOpt == "MPC"
        iSEI = (kSEI*ℯ^(-ESEI/R/Tk)) ./ (nSEI*(1 .+λ .* θ) .* .√ t[1]);
    else # typeOpt == "day-ahead"
        iSEI = (kSEI*ℯ^(-ESEI/R/Tk)) ./ (nSEI*(1 .+λ .* θ) .* .√ t);
    end
    # AM
    # iAM = kAM*ℯ^(-EAM/R/Tk) * SoCsa .* abs.(isa)*Qsa0*3600;
    iAM = kAM*ℯ^(-EAM/R/Tk) * SoCsa .* abs.(isa)*Qsa0;
    # # implement Xavier (2021) for Li-plating.
    # # ηoc,k(x) = φs,k(x) − φe,k(x) −Uref,oc − Fjoc,k(x)Rﬁlm,
    # Uref=0;
    # update Q and Rfilm
    iloss = iSEI .+ iAM; # total current loss
    Qloss = cumsum(iloss).*Δt./3600; # # Cell capacity, remember unit transf As <-> Ah
    δSEI = δSEI0 .+ cumsum(iSEI) * Δt * MSEI/nSEI/F/ρSEI/An # SEI layer thickness [m]
    # The S leaves the electrolyte to form the SEI layer.
    εₑ = εₑ0 .- as*Δt*cumsum(iSEI)*MSEI/nSEI/F/ρSEI/An; # electrolyte vol. fraction
    κeff = κ*εₑ .^ brug;
    R0 = R00 .+ εₛ ./ κeff * Δt .* cumsum(iSEI)*MSEI/nSEI/F/ρSEI/An # Series resistance update [Ω]

    if typeOpt == "MPC"
        # ALL OF THIS SHOULD BE IN results[ts+1]
        # # Update results dictionary
        haskey(results,"Q$key") ? results[Q_key][shift+1] = copy(Qsa0.-Qloss) : merge!(results, Dict("Q$key"=>copy(Qsa0.-Qloss)));
        haskey(results,"R0$key") ? results["R0$key"][shift+2] = copy(R0) : merge!(results,Dict("R0$key"=>copy(R0)));
        haskey(results,"δSEI$key") ? results["δSEI$key"][shift+2] = copy(δSEI) : merge!(results,Dict("δSEI$key"=>copy(δSEI))); 
        haskey(results,"εₑ$key") ? results["εₑ$key"][shift+2] = copy(εₑ) : merge!(results,Dict("εₑ$key"=>copy(εₑ)));
    else # typeOpt == "day-ahead"
        haskey(results,"Q$key") ? results[Q_key] = copy(Qsa0.-Qloss) : merge!(results, Dict("Q$key"=>copy(Qsa0.-Qloss)));
        haskey(results,"R0$key") ? results["R0$key"] = copy(R0) : merge!(results,Dict("R0$key"=>copy(R0)));
        haskey(results,"δSEI$key") ? results["δSEI$key"] = copy(δSEI) : merge!(results,Dict("δSEI$key"=>copy(δSEI))); 
        haskey(results,"εₑ$key") ? results["εₑ$key"] = copy(εₑ) : merge!(results,Dict("εₑ$key"=>copy(εₑ)));
    end

    # final state for the stgAsset
    update_stgAsset_deg!(stgAsset, results, key; typeOpt=typeOpt)
    return stgAsset, results
end

function simulate_storage_asset!(stgAsset::BESSData, results::Dict, key::String; typeOpt::String="MPC")
    # This function simulates a storage asset for a given PsaOpt using LiiBRA
    # The inputs are:
    # - stgAsset: a BESSData object containing the storage asset parameters
    # - results: a dictionary containing the setpoints of the controller
    # - key: the storage asset device name in the results dictionary
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    typeOpt == "MPC" ? shift = 1 : nothing;
    # take the length from results because its the one handled by handleInfeasible()
    Δt = copy(results[:"t"][2]-results[:"t"][1]); # setpoint timestep
    typeOpt == "day-ahead" ? shift = Int(24*3600/Δt)-1 : nothing;

    # The model to be simulated is a PBROM, but such model might not always be available in the EMS controller.
    # Thus we first check if the model is available in the EMS controller, if not we build one here.
    # If the model is available, we just simulate it.
    perfModel = perfModel_matching(stgAsset)

    if typeOpt == "MPC"
        # tk = copy(results[:"t"][shift]); # time vector
        tk = copy(results[:"t"][shift+1]); # time vector
    else
        tk = copy(results[:"t"][1:shift+1]); # time vector
    end
    # unfortunately the model is not realised in the same timesteps as the setpoints, 
    # ie. Δt ≠ perfModel.Cell.RA.SamplingT
    # Thus we need to upsample the time ref. of the setpoints to simulate the model
    # We do this by repeating the setpoints and generating a longer time ref.
    upSampRatio = Int(Δt / perfModel.Cell.RA.SamplingT) # ratio cell model sampling time
    tk = 0:1:((length(tk)*upSampRatio)-1);

    # Power setpoint
    if typeOpt == "MPC"
        PsaOpt = copy(results["P$key"][shift+1]); # [kW]
        PsaOpt = repeat([PsaOpt], inner=upSampRatio)
    else # typeOpt == "day-ahead"
        PsaOpt = copy(results["P$key"][1:shift+1]); # [kW]
        PsaOpt = repeat(PsaOpt, inner=upSampRatio)
    end
    # we adapt for the packs series and parallel cells and the units
    PsaOpt = 1e3*PsaOpt / stgAsset.GenInfo.Ns / stgAsset.GenInfo.Np;
    
    # The idea is to simulate the storage asset 
    # Tk = stgAsset.GenInfo.Tk; # Time constant
    Tk = ones(size(PsaOpt))*perfModel.Cell.Const.T # Cell Temperature;

    # Update the PBROM matrices
    SoC0 = copy(stgAsset.GenInfo.SoC0); # Starting SOC
    # SList = collect(1:-0.1:stgAsset.GenInfo.SoCLim[1]) # List of SOC points for model generation
    SList = [collect(1.0:-0.1:0.2); 0.15] # List of SOC points for model generation
    Sₑ = 4 # Spatial points in electrolyte
    Sₛ = 2 # Spatial point in solid
    # Spatial!(perfModel.Cell, Sₑ, Sₛ) # ideally we should be using this instead of Base.invokelatest(), but there is an issue with LiiBRA.jl
    Base.invokelatest(Spatial!, perfModel.Cell, Sₑ, Sₛ)
    # A, B, C, D = Realise(perfModel.Cell, SList)
    A, B, C, D = Base.invokelatest(Realise,perfModel.Cell, SList)
    # stgAsset.PerfParameters.A = A; stgAsset.PerfParameters.B = B;
    # stgAsset.PerfParameters.C = C; stgAsset.PerfParameters.D = D;
    
    # Calculate the transition Sₛₐ,ₜ₊₁==Sᴹ(Sₛₐ,ₜ , Pₛₐₜ*, Wₜ₊₁)
    # stgVars, ~ = Simulate(perfModel.Cell, PsaOpt, "Power", Tk, SList, SoC0, A, B, C, D, tk)
    stgVars, ~ = Base.invokelatest(Simulate, perfModel.Cell, PsaOpt, "Power", Tk, SList, SoC0, A, B, C, D, tk)

    if any(isnan.(stgVars.Cell_SOC))
        println("NaNs in the solution, $key is out of bounds.")
        # For the state Sₛₐ,ₜ₊₁, replace NaNs with the last valid value 
        iNaN = findall(isnan.(stgVars.Cell_SOC)); # indeces of the SoC NaNs
        # stgVars.Cell_SOC[iNaN] .= stgAsset.GenInfo.SoCLim[1]; # this assumes its NaN only in the lowerlimit
        # if the first NaN is in the first index, then we have to replace it with the initial SoC
        iNaN[1] == CartesianIndex(1,1) ? SoCNaN = SoC0 : SoCNaN = stgVars.Cell_SOC[iNaN[1]-CartesianIndex(1,0)];
        # check if SoCNaN is out of bounds
        SoCNaN < stgAsset.GenInfo.SoCLim[1] ? SoCNaN = (stgAsset.GenInfo.SoCLim[1] .+ 1e-4) : nothing;
        SoCNaN > stgAsset.GenInfo.SoCLim[2] ? SoCNaN = (stgAsset.GenInfo.SoCLim[2] .- 1e-4) : nothing;
        stgVars.Cell_SOC[iNaN] .= SoCNaN
        # indeces of the vₜ NaNs
        iNaN = findall(isnan.(stgVars.Cell_V));
        # stgVars.Cell_V[iNaN] .= stgAsset.GenInfo.vLim[1];
        iNaN[1] == CartesianIndex(1,1) ? vtNaN = stgAsset.GenInfo.vLim[1] : vtNaN = stgVars.Cell_V[iNaN[1]-CartesianIndex(1,0)];
        # check if vtNaN is out of bounds
        vtNaN < stgAsset.GenInfo.vLim[1] ? vtNaN = stgAsset.GenInfo.vLim[1] : nothing;
        vtNaN > stgAsset.GenInfo.vLim[2] ? vtNaN = stgAsset.GenInfo.vLim[2] : nothing;
        stgVars.Cell_V[iNaN] .= vtNaN
        # for the actions replace NaNs with 0
        iNaN = isnan.(stgVars.Iapp); # indeces of the NaNs
        stgVars.Iapp[iNaN] .= 0.;
        PsaOpt[iNaN[1:end-1]] .= 0.;
        # update the results dictionary with PsaOpt
        PsaOpt = 1e-3*PsaOpt * stgAsset.GenInfo.Ns * stgAsset.GenInfo.Np; # pack power in kW
        results["P$key"] = copy(PsaOpt[1:upSampRatio:end]); # save the downsampled array
    end
    # ensure that the SoC is within the limits
    if any(stgVars.Cell_SOC .< stgAsset.GenInfo.SoCLim[1]) || any(stgVars.Cell_SOC .> stgAsset.GenInfo.SoCLim[2])
        iMin = stgVars.Cell_SOC .< stgAsset.GenInfo.SoCLim[1]
        iMax = stgVars.Cell_SOC .> stgAsset.GenInfo.SoCLim[2]
        # replace the SoC values that are out of bounds
        stgVars.Cell_SOC[iMin] .= (stgAsset.GenInfo.SoCLim[1] .+ 1e-4);
        stgVars.Cell_SOC[iMax] .= (stgAsset.GenInfo.SoCLim[2] .- 1e-4);
        # same for the voltages
        stgVars.Cell_V[stgVars.Cell_V .< stgAsset.GenInfo.vLim[1]] .= (stgAsset.GenInfo.vLim[1] .+ 1e-4);
        stgVars.Cell_V[stgVars.Cell_V .> stgAsset.GenInfo.vLim[2]] .= (stgAsset.GenInfo.vLim[2] .- 1e-4);
        # for the actions replace NaNs with 0
        stgVars.Iapp[findall(iMin .| iMax)] .= 0.;
        PsaOpt[findall(iMin[1:end-1])] .= 0.;
        PsaOpt[findall(iMax[1:end-1])] .= 0.;
        # update the results dictionary with PsaOpt
        PsaOpt = 1e-3*PsaOpt * stgAsset.GenInfo.Ns * stgAsset.GenInfo.Np; # pack power in kW
        results["P$key"] = copy(PsaOpt[1:upSampRatio:end]); # save the downsampled array
    end

    # Update the results dictionary with the Performance Vars
    # modify the key for the EV case. Check if it is "evTot" or "evTot[$n]"
    # if it is "evTot[$n]" then we need to change it to "ev[$n]"
    if occursin("evTot", key)
        key = replace(key, "evTot" => "ev")
    end

    if typeOpt == "MPC"
        results["SoC$key"][shift+2] = copy(stgVars.Cell_SOC[end]);
        # stgAsset.GenInfo.SoC0 = copy(results["SoC$key"][2])
        stgAsset.GenInfo.SoC0 = copy(stgVars.Cell_SOC[end])
        results["vt$key"][shift+1] = copy(stgVars.Cell_V[end]); # this one is up for debate CHECK
        results["i$key"][shift+1] = copy(stgVars.Iapp[end]);
    else # typeOpt == "day-ahead"
        results["SoC$key"] = copy(stgVars.Cell_SOC[1:upSampRatio:end]);
        stgAsset.GenInfo.SoC0 = copy(stgVars.Cell_SOC[end])
        results["vt$key"] = copy(stgVars.Cell_V[1:upSampRatio:end]);
        results["i$key"] = copy(stgVars.Iapp[1:upSampRatio:end-1]);
    end

    # Update the rest of the states depending on the type of storage device
    
    # For ECMs
    # if PerfParams ≠ CIDRAPBROMPerfParams then update straight from the optimization results
    # not the most accurate but it is the best estimate
    typeof(stgAsset.PerfParameters) == ECMPerfParams ? stgAsset.PerfParameters.iRn0[1] = copy(stgVars.ϕ_ẽ1[2,end]+stgVars.ϕ_ẽ2[2,end]) : nothing; # should I recalculate results["iR1$key"]??

    # now use the current i to calculate the ageing of the storage asset
    # Update θsa, the ageing of the storage asset
    simulate_storage_asset_deg!(stgAsset, perfModel, results, key; typeOpt=typeOpt)

    return stgAsset, results
    # return stgAsset, stgVars
end

function simulate_storage_asset!(stgAsset::TESSData, results::Dict, key::String; typeOpt::String="MPC")
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    typeOpt == "MPC" ? shift = 1 : nothing;
    # take the length from results because its the one handled by handleInfeasible()
    Δt = copy(results[:"t"][2]-results[:"t"][1]); # setpoint timestep
    typeOpt == "day-ahead" ? shift = Int(24*3600/Δt)-1 : nothing;
    # typeOpt == "day-ahead" ? shift = Int(ceil(length(results["t"])/2))-1 : nothing;
    
    # Inputs:
    # time vector
    tk = results[:"t"]; 
    Δt = tk[2]-tk[1]; # setpoint timestep

    if typeOpt == "MPC"
        # easy simulation of a first order bucket model
        PsaOpt = results["P$key"][shift+1];
        SoC0 = copy(stgAsset.SoC0); # Starting SOC
        # easy simulation of a first order bucket model
        if PsaOpt < 0 # when the TESS is being charged
            # check if the TESS is being overcharged
            if (SoC0 .- PsaOpt .* Δt .* stgAsset.η ./ stgAsset.Q / 3600) > stgAsset.SoCLim[2]
                # check if the TESS approaching from below or above the SoCmax
                if SoCsa0 < stgAsset.SoCLim[2]
                    SoCsa = copy(stgAsset.SoCLim[2]); # SoCₜ₊₁ --> SoCmax
                    # reduce the power to avoid overcharging and reach SoCmax
                    PsaOpt = (stgAsset.SoCLim[2] - SoC0) / Δt / stgAsset.η * stgAsset.Q * 3600;
                else
                    PsaOpt = 0; # reject charging
                    SoCsa = copy(SoC0);
                end
            else
                SoCsa = copy(SoC0 .- PsaOpt .* Δt .* stgAsset.η ./ stgAsset.Q / 3600);
            end
        else
            SoCsa = copy(SoC0 .- PsaOpt .* Δt .* stgAsset.η ./ stgAsset.Q / 3600);
        end
        # previous implementation
        # results["SoC$key"] = copy(SoC0 .- cumsum(PsaOpt) .* Δt .* stgAsset.η ./ stgAsset.Q / 3600);
        # Save in the results Dict
        results["SoC$key"][shift+2] = copy(SoCsa);
        results["P$key"][shift+1] = copy(PsaOpt);
        # save the last SoC for the next initial SoC
        stgAsset.SoC0 = copy(SoCsa);
    else # typeOpt == "day-ahead"
        # easy simulation of a first order bucket model
        PsaOpt = results["P$key"][1:shift+1];
        SoCsa = zeros(shift+1); # initialize the SoC
        SoC0 = copy(stgAsset.SoC0); # Starting SOC
        SoCsa[1] = copy(SoC0); # initial SoC
        # easy simulation of a first order bucket model
        # with a check for overcharging
        for i ∈ 2:shift+1
            if PsaOpt[i] < 0 # when the TESS is being charged
                # check if the TESS is being overcharged
                if (SoCsa[i-1] .- PsaOpt[i] .* Δt .* stgAsset.η ./ stgAsset.Q / 3600) > stgAsset.SoCLim[2]
                    # check if the TESS approaching from below or above the SoCmax
                    if SoCsa[i-1] < stgAsset.SoCLim[2]
                        SoCsa[i] = copy(stgAsset.SoCLim[2]); # SoCₜ₊₁ --> SoCmax
                        # reduce the power to avoid overcharging and reach SoCmax
                        PsaOpt[i] = (stgAsset.SoCLim[2] - SoCsa[i-1]) / Δt / stgAsset.η * stgAsset.Q * 3600;
                    else
                        PsaOpt[i] = 0; # reject charging
                        SoCsa[i] = copy(SoCsa[i-1]);
                    end
                else
                    SoCsa[i] = copy(SoCsa[i-1] .- PsaOpt[i] .* Δt .* stgAsset.η ./ stgAsset.Q / 3600);
                end
            else
                SoCsa[i] = copy(SoCsa[i-1] .- PsaOpt[i] .* Δt .* stgAsset.η ./ stgAsset.Q / 3600);
            end
            # previous implementation
            # results["SoC$key"] = copy(SoC0 .- cumsum(PsaOpt) .* Δt .* stgAsset.η ./ stgAsset.Q / 3600);
            # Save in the results Dict
            results["SoC$key"] = copy(SoCsa);
            results["P$key"] = copy(PsaOpt);
        end
        # save the last SoC for the next initial SoC
        stgAsset.SoC0 = copy(results["SoC$key"])[end];
    end
    return stgAsset, results
end

function simTransitionFun!(results::Dict, data::Dict, s::modelSettings; typeOpt::String="MPC")
    # This function simulates the transition function Sₐ,ₜ₊₁ = Sₐ,ₜᴹ(Sₐ,ₜ , Pₐ,ₜ* , Wₐ,ₜ₊₁).
    # Sₐ,ₜ: state in timestep t
    # Pₐ,ₜ*: optimal setpoint in timestep t
    # Wₐ,ₜ₊₁: exogenous information in timestep t+1
    
    @assert typeOpt ∈ ["MPC", "day-ahead"];
    typeOpt == "MPC" ? shift = 1 : nothing;
    # take the length from results because its the one handled by handleInfeasible()
    Δt = copy(results[:"t"][2]-results[:"t"][1]); # setpoint timestep
    typeOpt == "day-ahead" ? shift = Int(24*3600/Δt)-1 : nothing;
    # typeOpt == "day-ahead" ? shift = Int(ceil(length(results["t"])/2))-1 : nothing;

    # Inputs:
    # These are the power setpoints of the storage devices (BESS, EVs, TESS) and the Heat Pump.
    nEV=s.nEV;
    
    # Pick the time window
    t = results["t"];
    t0=t[1]; Δt=t[2]-t[1];
    it0 = round(Int, (t0/Δt))
    itend = it0 + shift;

    # To simulate the transition we have to distinguish between the MPC and the day-ahead case.
    # In the MPC case, we have to simulate the transition for 1 timestep.
    # In the day-ahead case, we have to simulate the transition for 1 day.
    if typeOpt == "MPC"
        # First, we get the optimal decisions from our policy function.
        # PevKey = ["Pev[$n]" for n ∈ nEV];
        # keyOpt = ["Pbess", PevKey, "Ptess", "Phpe"];
        # Phpe = results[:"Phpe"];
        # PaOpt = Dict();
        Pbess = results["Pbess"][shift+1];
        Phpe = results["Phpe"][shift+1];
        Pev = [results["Pev[$n]"][shift+1] for n ∈ 1:nEV];
        
        # Second, we get the exogenous information.
        # P = P̂ + ϵ
        Ple = data["grid"].loadE[it0]
        Plt = data["grid"].loadTh[it0]
        PpvMPPT = data["SPV"].MPPTData[it0]
        Pst = PpvMPPT.*data["ST"].η

        # First of all we compensate for the EMS numerical errors by adjusting the HP and grid powers.
        # from the Thermal balance we adjust the HP
        # results["Phpe"][shift+1] =  copy((Plt - Pst - Ptess) ./ data["HP"].η) # check the case were Phpe<0
        # adjust TESS power instead of HP (to avoid negative HP power)
        results["Ptess"][shift+1] =  copy(Plt - Pst - Phpe .* data["HP"].η)
        # and from the electrical one the grid
        if nEV != 1
            γf = [results["γf"][n][shift+1] for n ∈ 1:nEV];
            results["Pg"][shift+1] = copy(Ple + copy(results["Phpe"][shift+1]) - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf[n] for n ∈ 1:nEV]));
        else
            γf = results["γf"][shift+1];
            results["Pg"][shift+1] = copy(Ple + copy(results["Phpe"][shift+1]) - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf for n ∈ 1:nEV]));
        end
    else # typeOpt == "day-ahead"
        # First, we get the optimal decisions from our policy function.
        Pbess = results["Pbess"][1:shift+1];
        Pev = [results["Pev[$n]"][1:shift+1] for n ∈ 1:nEV];
        # check if we have a heat pump
        haskey(results, "Phpe") ? Phpe = results["Phpe"][1:shift+1] : Phpe = zeros(shift+1);
        
        # Second, we get the exogenous information.
        # # P = P̂ + ϵ
        # Ple = data["grid"].loadE + rand(data["grid"].loadE);
        # Plt = data["grid"].loadTh + rand(data["grid"].loadTh);
        # PpvMPPT = data["SPV"].MPPTData + rand(data["SPV"].MPPTData);
        Ple = copy(data["grid"].loadE[it0:itend])
        Plt = copy(data["grid"].loadTh[it0:itend])
        PpvMPPT = copy(data["SPV"].MPPTData[it0:itend])
        Pst = copy(PpvMPPT.*data["ST"].η)
        # Third of all we compensate for the EMS numerical errors by adjusting the HP and grid powers.
        # from the Thermal balance we adjust the HP
        results["Ptess"] =  copy(Plt - Pst - Phpe .* data["HP"].η)
        # and from the electrical one the grid
        if nEV != 1
            γf = [results["γf"][n][1:shift+1] for n ∈ 1:nEV];
            results["Pg"][1:shift+1] = copy(Ple + Phpe - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf[n] for n ∈ 1:nEV]));
        else
            γf = results["γf"][1:shift+1];
            results["Pg"][1:shift+1] = copy(Ple + Phpe - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf for n ∈ 1:nEV]));
        end
        
    end
    
    # get all the SoC
    for sa ∈ ["BESS", "EV", "TESS"]
        if typeof(data[sa]) <: Vector
            for n ∈ 1:length(data[sa])
                # pick the right key
                if typeof(data[sa][n]) == BESSData
                    key="bess[$n]"
                    stgAsset = data[sa][n];
                elseif typeof(data[sa][n]) == EVData
                    key="evTot[$n]"
                    stgAsset = data[sa][n].carBatteryPack;
                elseif typeof(data[sa][n]) == TESSData
                    key="tess[$n]"
                    stgAsset = data[sa][n];
                end
                # Sₛₐ,ₜ₊₁ = Sₛₐ,ₜᴹ(Sₛₐ,ₜ , Pₛₐ,ₜ* , Wₛₐ,ₜ₊₁)
                simulate_storage_asset!(stgAsset, results, key; typeOpt=typeOpt)
            end
        else
            # pick the right key
            if typeof(data[sa]) == BESSData
                key="bess"
                stgAsset = data[sa];
            elseif typeof(data[sa]) == EVData
                key="evTot"
                stgAsset = data[sa].carBatteryPack;
            elseif typeof(data[sa]) == TESSData
                key="tess"
                stgAsset = data[sa];
            end
            # Sₛₐ,ₜ₊₁ = Sₛₐ,ₜᴹ(Sₛₐ,ₜ , Pₛₐ,ₜ* , Wₛₐ,ₜ₊₁)
            simulate_storage_asset!(stgAsset, results, key; typeOpt=typeOpt)
        end
    end
    # re-balance the power
    if typeOpt == "MPC"
        # from the Thermal balance we adjust the HP
        # since the TESS overcharge might have come from the ST or the HP this new HP power might be negative,
        # thus we have to check if the HP power goes negative and dump it in the house
        # excess heat, rejected from the TESS
        if (Plt - Pst - results["Ptess"][shift+1]) < 0.
            Qex = copy(-Plt + Pst + results["Ptess"][shift+1]) .* data["HP"].η
            results["Phpe"][shift+1] = 0.; # set the HP power to 0
            results["Plt"][shift+1] = copy(Plt .+ Qex); # dump the rejected heat in the house
        else
            results["Phpe"][shift+1] =  copy(Plt - Pst - results["Ptess"][shift+1]) / data["HP"].η
        end
        # electrical re-balance with the new HP power
        # re-extract the power of the sa (just in case we have hit the SoC limits)
        Pev = [results["Pev[$n]"][shift+1] for n ∈ 1:nEV];
        Pbess = results["Pbess"][shift+1];
        if nEV != 1
            γf = [results["γf"][n][shift+1] for n ∈ 1:nEV];
            Phpe = results["Phpe"][shift+1]
            results["Pg"][shift+1] = copy(Ple + Phpe - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf[n] for n ∈ 1:nEV]));
        else
            γf = results["γf"][shift+1];
            Phpe = results["Phpe"][shift+1]
            results["Pg"][shift+1] = copy(Ple + Phpe - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf for n ∈ 1:nEV]));
        end
    else # typeOpt == "day-ahead"
        # from the Thermal balance we adjust the HP
        results["Phpe"] =  copy(Plt - Pst - results["Ptess"]) / data["HP"].η
        # since the TESS overcharge might have come from the ST or the HP this new HP power might be negative,
        # thus we have to check if the HP power goes negative and dump it in the house
        # excess heat, rejected from the TESS
        Qex = zeros(shift+1); # excess heat, rejected from the TESS
        Qex[results["Phpe"] .< 0] = copy(-results["Phpe"][results["Phpe"] .< 0]) .* data["HP"].η;
        results["Phpe"][results["Phpe"] .< 0] .= 0; # set the HP power to 0
        results["Plt"] = copy(Plt .+ Qex); # dump the rejected heat in the house
        # electrical re-balance with the new HP power
        # re-extract the power of the sa (just in case we have hit the SoC limits)
        Pev = [results["Pev[$n]"][1:shift+1] for n ∈ 1:nEV];
        Pbess = results["Pbess"][1:shift+1];
        if nEV != 1
            γf = [results["γf"][n][1:shift+1] for n ∈ 1:nEV];
            results["Pg"][1:shift+1] = copy(Ple + results["Phpe"] - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf[n] for n ∈ 1:nEV]));
        else
            γf = results["γf"][1:shift+1];
            results["Pg"][1:shift+1] = copy(Ple + results["Phpe"] - PpvMPPT - # data
                Pbess - sum([Pev[n] .* γf for n ∈ 1:nEV]));
        end
    end
    return results, data
end
