## ESS functions
# This file contains the functions used for the electrical ESS (BESS and EVs) of the EMS.
# The functions include different models for battery performance and aging.
#
# By: Darío Slaifstein, PhD-student @TU Delft, DCES.
# Branch: new_mpec_1
# Version: 1.0
# Date: 10/01/2025


function add_battPerf(model::InfiniteModel, sets::modelSettings, data::BESSData)
# battPerf: Battery performance modeling function
# This function adds variables and constraints to the model obj following the different
# types of battery models (Bucket, ECM or PB)
    t=model[:t];
    t0=supports(t)[1];
    # bPbess = model[:bPbess];
    bbess⁻ = model[:bbess⁻];
    bbess⁺ = model[:bbess⁺];
    @unpack GenInfo, PerfParameters, AgingParameters=data
    @unpack type=PerfParameters   
    @unpack_Generic GenInfo
    @unpack ocvLine = OCVParam;
    Npbess = Np; Nsbess = Ns; 
    Qbess0 = initQ*SoHQ; ηbess = η;
    PbessMax = PowerLim[2]; # Max power [kW]
    SoCbess0 = SoC0;
    aOCV=ocvLine[1]; bOCV=ocvLine[2];
    vmin = vLim[1]; # Min voltage [V]
    vmax = vLim[2]; # Max voltage [V]
    imax = 1e3*PbessMax/Npbess/Nsbess/vmin; # Max current [A]

    if type == "bucket"
    # Model variables
        OCVbess0 = aOCV+bOCV*SoCbess0;
        @variables(model, begin
            # aOCV ≤ OCVbess ≤ aOCV+bOCV, Infinite(t), (start=OCVbess0) # open circuit voltage of the cell
            vmin ≤ OCVbess ≤ vmax, Infinite(t), (start=OCVbess0) # open circuit voltage of the cell
            -imax ≤ ibess ≤ imax, Infinite(t), (start=0) # current per branch
        end); 
    
    # Initial conditions
        #=
        @constraints(model, begin
                OCVbess(t0) ==  OCVbess0 
                ibess(t0) ==  1e3*Pbess0/Npbess/Nsbess/OCVbess0
            end); 
        =#
    # Model constraints
        @constraints(model, begin
            OCVbess == aOCV+bOCV*model[:SoCbess] # Linear voltage model
            # OCVbess == OCVfromSOCtemp(SoCbess, T, data["BESS"]) # Lookup table voltage model
            ibess * OCVbess == 1e3*model[:Pbess]/Npbess/Nsbess # current per branch. 1e3 to convert kW->W
        end);
        
        if sets.costWeights[3] != 0 # Aging check
            @variables(model, begin
                -imax ≤ ilossbess ≤ imax, Infinite(t), (start=0.0) # total aging
                0.8*Qbess0 ≤ Qbess ≤ Qbess0, Infinite(t), (start=Qbess0) # cell capacity
            end);
            @constraints(model, begin
                # ∂.(model[:SoCbess], t) * 1e3 .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess/3600 * 1e3 # Aging Qbess
                # ∂.(model[:SoCbess], t) * 1e3 .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess0/3600 * 1e3 # No-Aging Qbess
                ∂.(model[:SoCbess], t) * 1e3 .== -(ηbess * bbess⁻ + bbess⁺)*ibess/Qbess0/3600 * 1e3 # Aging Qbess
                Qbess(t0) .== Qbess0;
            end);
        else
            # @constraint(model, ∂.(model[:SoCbess], t) * 1e3 .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess0/3600 * 1e3) # Static Qbess
            @constraint(model, ∂.(model[:SoCbess], t) * 1e3 .== -(ηbess * bbess⁻ + bbess⁺)*ibess/Qbess0/3600 * 1e3) # Static Qbess
        end
    elseif type == "ECM"
    # Model variables
        @unpack R0Param, RParam, RCParam, iRn0, vt0 = PerfParameters
        R0bess0=R0Param[1]; # not very elegant, but it works
        vtbess0 = vt0;
        iR1bess0=iRn0[1]; # not very elegant, but it works
        R1bess=RParam; 
        taubess=RCParam;
        @unpack iloss0 = AgingParameters
        
        OCVbess0 = aOCV+bOCV*SoCbess0;
        @variables(model, begin
            # aOCV ≤ OCVbess ≤ aOCV+bOCV, Infinite(t),(start=OCVbess0) # open circuit voltage of the cell
            # aOCV ≤ vtbess ≤ aOCV+bOCV, Infinite(t),(start=OCVbess0)  # terminal voltage of the cell
            vmin ≤ OCVbess ≤ vmax, Infinite(t),(start=OCVbess0) # open circuit voltage of the cell
            vmin ≤ vtbess ≤ vmax, Infinite(t),(start=OCVbess0)  # terminal voltage of the cell
            -imax ≤ ibess ≤ imax, Infinite(t), (start=0.0) # total current per branch
            -imax ≤ iR1bess ≤ imax, Infinite(t), (start=0.0) # pole current
            -imax ≤ ilossbess ≤ imax, Infinite(t), (start=0.0) # total aging
            0.8*Qbess0 ≤ Qbess ≤ Qbess0, Infinite(t), (start=Qbess0) # cell capacity
        end);

    # Initial conditions
        @constraints(model, begin
            # OCVbess(t0) == OCVbess0;
            # ibess(t0) == 1e3*Pbess0/Npbess/Nsbess/vtbess0;
            iR1bess(t0) == iR1bess0; # will be an estimation
            # vtbess(t0) == vtbess0; # will be measurement
            # ilossbess(t0) == iloss0;
        end); 

    # Model constraints
        @constraints(model, begin
            ibess * vtbess == 1e3*model[:Pbess]/Npbess/Nsbess # current per branch. 1e3 to convert kW->W
            # Transition function
            # ∂.(model[:SoCbess], t)* 1e3 .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess/3600 * 1e3 # Aging Qbess
            # ∂.(model[:SoCbess], t)* 1e3 .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess0/3600 * 1e3 # No-Aging Qbess
            ∂.(model[:SoCbess], t)* 1e3 .== -(ηbess * bbess⁻ + bbess⁺)*ibess/Qbess0/3600 * 1e3 # Aging Qbess
            ∂.(iR1bess, t) * 1e3 .== (-1 ./taubess*iR1bess+1 ./taubess*ibess) * 1e3 # RC resistor currents
            #Hysterisis?
            # Output equations
            # OCVbess == OCVfromSOCtemp(SoCbess, T, data["BESS"]) # Lookup table voltage model
            OCVbess .== aOCV+bOCV*model[:SoCbess] # Linear OCV model
            vtbess .== OCVbess .- R1bess*iR1bess .- R0bess0*ibess # no power fade
        end); 
    elseif type == "PBROM"
        # PBROM-xRA parameter unpacking
        @unpack_CIDRAPBROMPerfParams PerfParameters
        # Get ROM orders
        # states
        nX=1:size(A)[1];
        # outputs
        nY=1:size(C)[1];
        @unpack iloss0 = AgingParameters

    # Model variables
        @variables(model, begin
            # Perfromance variables
            # aOCV ≤ OCVbess ≤ aOCV+bOCV, Infinite(t),(start=OCVbess0) # open circuit voltage of the cell
            # aOCV ≤ vtbess ≤ aOCV+bOCV, Infinite(t),(start=OCVbess0)  # terminal voltage of the cell
            vmin ≤ OCVbess ≤ vmax, Infinite(t),(start=OCVbess0) # open circuit voltage of the cell
            vmin ≤ vtbess ≤ vmax, Infinite(t),(start=OCVbess0)  # terminal voltage of the cell
            # total current per branch
            ibess, Infinite(t), (start=0.0)
            # internal cell states, no physical meaning
            xbess[xx ∈ nX], Infinite(t), (start=0.0)
            # internal cell variables, physical meaning 
            ybess[yy ∈ nY], Infinite(t), (start=0.0)
            # Aging variables
            ilossbess, Infinite(t), (start=0.0) # total aging
            0.8*Qbess0 ≤ Qbess ≤ Qbess0, Infinite(t), (start=Qbess0) # cell capacity
        end);
    # Initial conditions
        @constraints(model, begin
            # OCVbess(t0) == OCVbess0;
            # vtbess(t0) == OCVbess0; # check
            # ibess(t0) == 1e3*Pbess0/Npbess/Nsbess/OCVbess0;
            # ilossbess(t0) == iloss0;
            xbess[xx ∈ nX, t0] == 0.0; # check
            # ybess[yy ∈ nY, t0] == 0.0; # check
        end);
    # Model constraints
    # Parameter list:
    # - A, B, C, D state space matrices
    # - Rfilmₖ, k ∈ [pos,neg] film resistance
    # - z0pₖ, k ∈ [pos,neg] z @ 0% Lithium Concentration
    # - z100pₖ, k ∈ [pos,neg] z @ 100% Lithium Concentration
    # - cs_maxₖ, k ∈ [pos,neg] max electrode/solid concentration
    # - ce0 electrolyte concentration
    # - αₖ k ∈ [pos,neg] charge transfer coefficient
    # - kₖ k ∈ [pos,neg] reaction rate constant
    # - t₀⁺ initial transference number
        # Concentrations
        SOC_Pos0 = SoCbess0*(z100p_Pos-z0p_Pos)+z0p_Pos
        SOC_Neg0 = SoCbess0*(z100p_Neg-z0p_Neg)+z0p_Neg
        Cse_Pos = SOC_Pos0.*cs_max_Pos .+ ybess[CsePosInd](t);
        Cse_Neg = SOC_Neg0.*cs_max_Neg .+ ybess[CseNegInd](t);
        Ce = ce0 .+ ybess[CeInd](t);
        z_p = Cse_Pos[1]./cs_max_Pos;
        z_n = Cse_Neg[1]./cs_max_Neg;

        OCVn = 1.97938*2.7182818284*exp(-39.3631*z_n) + 0.2482 - 
                0.0909*tanh(29.8538*(z_n - 0.1234)) - 
                0.04478*tanh(14.9159*(z_n - 0.2769)) - 
                0.0205*tanh(30.4444*(z_n - 0.6103))
        OCVp = -0.8090*z_p + 4.4875 - 0.0428*tanh(18.5138*(z_p - 0.5542))-
                17.7326*tanh(15.7890*(z_p - 0.3117)) +
                17.5842*tanh(15.9308*(z_p - 0.3120))

        # flux of positive electrode at x=L
        jL = ybess[FluxPosInd[1]]
        m_Pos = k_Pos.*(Cse_Pos.^α_Pos).*(Ce[1].^(1-α_Pos)) # reaction rate
        j0_Pos = m_Pos.*(cs_max_Pos .- Cse_Pos).^(1-α_Pos)[1] # exchange current density
        # overpotential at x=L
        ηL = (2*R*T)/F*asinh(jL/(2*j0_Pos[1]))
        # flux of negative electrode at x=0
        j0 = ybess[FluxNegInd[1]]
        m_Neg = k_Neg.*(Cse_Neg.^α_Neg).*(Ce[1].^(1-α_Neg)) # reaction rate
        j0_Neg = m_Neg.*(cs_max_Neg .- Cse_Neg).^(1-α_Neg)[1] # exchange current density
        # overpotential at x=0
        η0 = (2*R*T)/F*asinh(j0/(2*j0_Neg[1]))
        
        # From Planden (2023), Eq. (23)
        # ̃ϕₑ/Iₐₚₚ=[̃ϕₑ(z,s)]₁+[̃ϕₑ(z,s)]₂ 
        ϕ_ẽ1 = ybess[ϕ_ẽInd](t)
        # check
        ϕ_ẽ2 = 2*R*T*(1-t₀⁺)/F*(log(Ce ./ Ce[1]))
        @constraints(model, begin
            # CHECK
            ibess == 1e3*model[:Pbess]/Npbess/Nsbess/vtbess # current per branch. 1e3 to convert kW->W
            # Transition function
            # ∂.(model[:SoCbess], t) .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess/3600 # Aging Qbess
            ∂.(model[:SoCbess], t) .== -(ηbess * bPbess + (1-bPbess))*ibess/Qbess0/3600 # Aging Qbess
            ∂.(xbess[o ∈ nO], t) .== A*xbess[o ∈ nO]+B*ibess # ROM state transition
            # Output equations
            ybess .== C*xbess[o ∈ nO] + D*ibess
            vtbess .== (OCVp-OCVn) + (ηL-η0) + (ϕ_ẽ1[end]+ϕ_ẽ2[end]) +
                     (RFilm_Pos*jL-RFilm_Neg*j0)*F  # ROM output equation
        end);
    else
        error("Battery performance model not recognized")
    end
    return model;
end

function add_battPerf(model::InfiniteModel, sets::modelSettings, data::EVData)
# battPerf: Battery performance modeling function
    t=model[:t]; # continuous time
    Dt=sets.dTime; # discrete time
    t0=supports(t)[1]; tend=supports(t)[end];
    Δt=supports(t)[2]-supports(t)[1];
    nEV=sets.nEV;
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    day = ceil(Int, it0/(24*3600/Δt));
    # bPev = model[:bPev];
    bev⁻ = model[:bev⁻];
    bev⁺ = model[:bev⁺];

    @unpack GenInfo, PerfParameters, AgingParameters=data.carBatteryPack
    @unpack_Generic GenInfo
    @unpack ocvLine = OCVParam;
    performanceParams=PerfParameters;
    agingParams=AgingParameters;
    Npev = Np; Nsev = Ns; ηev = η;
    Qev0 = initQ*SoHQ; SoCev0 = SoC0; Pev0=P0;
    aOCV = ocvLine[1]; bOCV = ocvLine[2];
    vmin = vLim[1]; # Min voltage [V]
    vmax = vLim[2]; # Max voltage [V]
    imax = 1e3*PowerLim[2]/Npev/Nsev/vmin; # Max current [A]
    @unpack type=data.carBatteryPack.PerfParameters;
    
    # Driving consumption
    # μDrive=data.driveInfo.μDrive;
    # σDrive=data.driveInfo.σDrive;
    # # depLims[n,:]=data[n].driveInfo.depLims;
    # # arrLims[n,:]=data[n].driveInfo.arrLims;
    # tDep = data.driveInfo.tDep;
    SoCdep = data.driveInfo.SoCdep;
    γ = data.driveInfo.γ[it0:itend];
    Pdrive = data.driveInfo.Pdrive[day]
    # Now we need to project it into the cont t-domain.
    γ_interp = linear_interpolation(Dt, γ)
    @parameter_function(model, γf == (t) -> γ_interp(t)) # make InfiniteOpt compatible
    # γ_interp = linear_interpolation(Dt, γ[:]) 
    # @parameter_function(model, γf == (t) -> γ_interp(t))
    # Pdrive = rand(truncated(Normal(μDrive, σDrive); lower = 0.01)); # Gaussian distribution
    # # Pdrive .< 0 ? Pdrive = 0.01 : nothing; # safe lock for negative driving power
    # Ereq=sum(Pdrive.*(1 .-γ)*Δt)./3600;
    # # check if the driving power is greater than the energy in the battery pack.
    # while  Ereq .> Qev0.*Npev.*Nsev.*(aOCV.+bOCV)./1000*0.8
    #     Pdrive = rand(truncated(Normal(μDrive, σDrive); lower = 0.01)); # Gaussian distribution
    #     # Pdrive .< 0 ? Pdrive = 0.01 : nothing; # safe lock for negative driving power
    #     Ereq=sum(Pdrive.*(1 .-γ)*Δt)./3600;
    # end

    # User requirement at departure time
    # check if tDep is inside Dt
    SoCev=model[:SoCev];
    depIdx = findfirst(diff(γ) .== -1);
    @expression(model, ϵSoC[n ∈ 1:nEV], (SoCev[n] .- SoCdep) * 0.)
    if !isnothing(depIdx) # if there's a departure time
        tDep = Dt[depIdx] # departure time
        [ϵSoC[n] = SoCev[n](tDep) .- SoCdep for n ∈ 1:nEV]
    end
    # model[:ϵSoC]=[]; # assign name in the model
    # model[:ϵSoC] = SoCev(tDep) - SoCdep;
    # SoCev=model[:SoCev];
    # model[:ϵSoC]=[]; # assign name in the model
    # for day in eachindex(tDep[:]) # if there's more than one day loop over them
    #     td = tDep[day] # pick value
    #     td = td * 3600 + (24 * 3600 * (day - 1)) # change to secs and add days
    #     # find the nearest td inside supports(t)
    #     td = findmin(abs.(td .- supports(t)))[1]
    #     # Check if the time index is within the bounds of SoCev
    #     if td >= t0 && td <= tend
    #         ϵSoCexpr = SoCev[1](td) .- SoCdep
    #         push!(model[:ϵSoC], ϵSoCexpr) # push to the expression vector
    #     end
    # end
    
    if type == "bucket"
        # Model variables
        OCVev0 = [aOCV[n]+bOCV[n]*SoCev0[n] for n in 1:nEV];    
        @variables(model, begin
            # aOCV[n] ≤ OCVev[n ∈ 1:nEV] ≤ aOCV[n]+bOCV[n], Infinite(t),(start=OCVev0[n]) # open circuit voltage of the cell
            vmin[n] ≤ OCVev[n ∈ 1:nEV] ≤ vmax[n], Infinite(t),(start=OCVev0[n]) # open circuit voltage of the cell
            -imax[n] ≤ iev[n ∈ 1:nEV] ≤ imax[n], Infinite(t)  # current per branch 
        end);
        # Initial conditions
            #=
                @constraints(model, begin
                    [n ∈ 1:nEV], OCVev[n](t0) == OCVev0[n];
                    [n ∈ 1:nEV], iev[n](t0) == 1e3*Pev0[n]/Npev[n]/Nsev[n]/OCVev0[n]; # the 1e3 to kW -> W
                end);
            =#
        # Model constraints
        @constraints(model, begin
            [n ∈ 1:nEV], bev⁻[n] .+ γf .≥ 1
            availability[n ∈ 1:nEV], model[:γf].*model[:Pev][n] + (1-model[:γf]).*Pdrive[n] - model[:PevTot][n] .== 0 # power balance
            [n ∈ 1:nEV], OCVev[n] .== aOCV[n]+bOCV[n]*model[:SoCev][n] # linear voltage model
            # [n ∈ 1:nEV], OCVev[n] == OCVfromSoC(SoCev[n]) # Lookup table voltage model
            [n ∈ 1:nEV], iev[n] .* OCVev[n] .== 1e3*model[:PevTot][n]/Npev[n]/Nsev[n] # current per branch     
        end);
        if sets.costWeights[3] != 0 # Aging check
                @variables(model, begin
                    -imax[n] ≤ ilossev[n ∈ 1:nEV] ≤ imax[n], Infinite(t) # total aging
                    0.8*Qev0[n] .≤ Qev[n ∈ 1:nEV] .≤ Qev0[n], Infinite(t) # cell capacity
                end);
                @constraints(model, begin
                    # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev[n]/3600 * 1e3 # Aging Qev
                    # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev0[n]/3600 * 1e3 # No-Aging Qev
                    [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bev⁻[n] + bev⁺[n])*iev[n]/Qev0[n]/3600 * 1e3
                    [n ∈ 1:nEV], Qev[n](t0) .== Qev0;
                end);
            else
                # Static Qbess
                # @constraint(model, [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev0[n]/3600 * 1e3)
                @constraint(model, [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bev⁻[n] + bev⁺[n])*iev[n]/Qev0[n]/3600 * 1e3)
            end
    elseif type == "ECM"
        # Model variables
            @unpack R0Param, RParam, RCParam, iRn0, vt0 = performanceParams
            R0ev0=R0Param[]; R1ev=RParam[]; tauev=RCParam[];
            vtev0= vt0; iR1ev0=iRn0[1];
            @unpack iloss0 = agingParams
            ilossev0=iloss0;
            OCVev0 = aOCV+bOCV*SoCev0;
            @variables(model, begin
                # aOCV[n] ≤ OCVev[n ∈ 1:nEV] ≤ aOCV[n]+bOCV[n], Infinite(t), (start=OCVev0[n]) # open circuit voltage of the cell
                # aOCV[n] ≤ vtev[n ∈ 1:nEV] ≤ aOCV[n]+bOCV[n], Infinite(t), (start=OCVev0[n]) # terminal voltage of the cell
                vmin[n] ≤ OCVev[n ∈ 1:nEV] ≤ vmax[n], Infinite(t), (start=OCVev0[n]) # open circuit voltage of the cell
                vmin[n] ≤ vtev[n ∈ 1:nEV] ≤ vmax[n], Infinite(t), (start=OCVev0[n]) # terminal voltage of the cell
                -imax[n] ≤ iev[n ∈ 1:nEV] ≤ imax[n], Infinite(t), (start=0.0) # total current per branch
                -imax[n] ≤ ilossev[n ∈ 1:nEV] ≤ imax[n], Infinite(t), (start=0.0) # total aging
                -imax[n] ≤ iR1ev[n ∈ 1:nEV] ≤ imax[n], Infinite(t), (start=0.0) # pole current
                0.8*Qev0[n] .≤ Qev[n ∈ 1:nEV] .≤ Qev0[n], Infinite(t), (start=Qev0[n]) # cell capacity
                # R0ev0[n] .≤ R0ev[n ∈ 1:nEV] .≤ R0ev0[n]*1.2, Infinite(t), (start=R0ev0[n]) # cell capacity
            end);
        # Initial conditions
            @constraints(model, begin
                # [n ∈ 1:nEV], OCVev[n](t0) == OCVev0[n];
                # [n ∈ 1:nEV], iev[n](t0) == 1e3*Pev0[n]/Npev[n]/Nsev[n]/vtev0[n]; # the 1e3 to kW -> W
                [n ∈ 1:nEV], iR1ev[n](t0) == iR1ev0[n];
                # [n ∈ 1:nEV], vtev[n](t0) == vtev0[n];
                # [n ∈ 1:nEV], ilossev[n](t0) == ilossev0[n];
            end); 
        # Model constraints
            @constraints(model, begin
                [n ∈ 1:nEV], bev⁻[n] .+ γf .≥ 1
                availability[n ∈ 1:nEV], model[:γf].*model[:Pev][n] + (1-model[:γf]).*Pdrive[n] - model[:PevTot][n] .== 0 # power balance
                [n ∈ 1:nEV], iev[n] * vtev[n] .== 1e3*model[:PevTot][n]/Npev[n]/Nsev[n] # current per branch. 1e3 to convert kW->W
                # Transition function
                # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev[n]/3600 * 1e3 # Aging Qev
                # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev0[n]/3600 * 1e3 # No-Aging Qev
                [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bev⁻[n] + bev⁺[n])*iev[n]/Qev0[n]/3600 * 1e3
                [n ∈ 1:nEV], ∂.(iR1ev[n], t) * 1e3 .== (-1/tauev[n]*iR1ev[n]+1/tauev[n]*iev[n]) * 1e3 # RC resistor currents
                #Hysterisis?
                # Output equations
                #[n ∈ 1:nEV], OCVev[n] == OCVfromSOCtemp(SoCev, T, data["EV"]) # Lookup table voltage model
                [n ∈ 1:nEV], OCVev[n] .== aOCV[n]+bOCV[n]*model[:SoCev][n] # Linear OCV model
                [n ∈ 1:nEV], vtev[n] .== OCVev[n] - R1ev[n]*iR1ev[n] - R0ev0[n]*iev[n]
                # [n ∈ 1:nEV], vtev[n] .== OCVev[n] - R1ev[n]*iR1ev[n] - R0ev[n]*iev[n]
            end); 
    # elseif type == "PBROM"
    end

    return model;
end

function add_battPerf(model::InfiniteModel, sets::modelSettings, data::Vector{EVData})
# battPerf: Battery performance modeling function
    t=model[:t]; # continuous time
    Dt=sets.dTime; # discrete time
    t0=supports(t)[1]; tend=supports(t)[end];
    Δt=supports(t)[2]-supports(t)[1];
    nEV=sets.nEV;
    # make new time window
    # it0 = round(Int,(t0/Δt) + 1);
    it0 = round(Int,(t0/Δt));
    itend = it0+length(Dt)-1;
    # calculate in which day is it0
    day = ceil(Int, it0/(24*3600/Δt));
    # bPev = model[:bPev];
    bev⁻ = model[:bev⁻];
    bev⁺ = model[:bev⁺];

    Npev = zeros(nEV,1); Nsev = zeros(nEV,1);
    ηev = zeros(nEV,1); Qev0 = zeros(nEV,1);
    SoCev0 = zeros(nEV,1); Pev0 = zeros(nEV,1);
    aOCV = zeros(nEV,1); bOCV = zeros(nEV,1);
    vmin = zeros(nEV,1); vmax = zeros(nEV,1);
    performanceParams=Vector{PerfParams}(undef, nEV);
    agingParams=Vector{AgingParams}(undef, nEV);
    for n in 1:nEV
        @unpack GenInfo, PerfParameters, AgingParameters=data[n].carBatteryPack
        @unpack_Generic GenInfo
        @unpack ocvLine = OCVParam;
        performanceParams[n]=PerfParameters;
        agingParams[n]=AgingParameters;
        Npev[n] = Np; Nsev[n] = Ns; ηev[n] = η;
        Qev0[n] = initQ*SoHQ; SoCev0[n] = SoC0; Pev0[n]=P0;

        aOCV[n]=ocvLine[1]; bOCV[n]=ocvLine[2];
        vmin[n]=vLim[1]; vmax[n]=vLim[2];
        imax[n]=1e3*PowerLim[2]/Np/Ns/vLim[1]; # Max current [A]
    end
    @unpack type=data[1].carBatteryPack.PerfParameters;

    # Driving consumption
    # μDrive=zeros(nEV,1);
    # σDrive=zeros(nEV,1);
    γ=zeros(nEV,length(Dt));
    # tDep=[];
    SoCdep = zeros(nEV)
    Pdrive = zeros(nEV)

    for n in 1:nEV
        # μDrive[n]=data[n].driveInfo.μDrive;
        # σDrive[n]=data[n].driveInfo.σDrive;
        SoCdep[n]=data[n].driveInfo.SoCdep;
        γ[n,:] = data[n].driveInfo.γ[it0:itend];
        Pdrive[n] = data[n].driveInfo.Pdrive[day]
        # push!(tDep, data[n].driveInfo.tDep);
    end
    # tDep=vcat(tDep'...); # reorganize in a matrix
    
    # Now we need to project it into the cont t-domain.
    # γ_interp = linear_interpolation((nEV, Dt), γ)
    γ_interp = linear_interpolation((1:nEV, Dt), γ)
    @parameter_function(model, γf[n ∈ 1:nEV] == (t) -> γ_interp(n, t)) # make InfiniteOpt compatible
    # Pdrive = [rand(truncated(Normal(μDrive[n], σDrive[n]); lower = 0.01)) for n in 1:nEV]; # Gaussian distribution
    # Pdrive[Pdrive .< 0] .= 0.01; # safe lock for negative driving power
    # Ereq=[sum(Pdrive[n].*(1 .-γ[n,:])*Δt)./3600 for n in 1:nEV];
    # # check if the driving power is greater than the energy in the battery pack.
    # while  any(Ereq .> Qev0.*Npev.*Nsev.*(aOCV.+bOCV)./1000*0.8)   
    #     Pdrive = [rand(truncated(Normal(μDrive[n], σDrive[n])); lower = 0.01) for n in 1:nEV]; # Gaussian distribution
    #     # Pdrive[Pdrive .< 0] .= 0.01; # safe lock for negative driving power
    #     Ereq=[sum(Pdrive[n].*(1 .-γ[n,:])*Δt)./3600 for n in 1:nEV];
    # end
    
    # User requirement at departure time
    # check if tDep is inside Dt
    SoCev=model[:SoCev];
    depIdx = [findfirst(diff(γ[n,:]) .== -1) for n ∈ 1:nEV];
    # initialize the ϵSoC in 0.0
    @expression(model, ϵSoC[n ∈ 1:nEV], (SoCev[n] .- SoCdep[n]) .* 0.)
    if .!isnothing(depIdx) # if there's a departure time
        # an element might be nothing, so we need to filter it out
        # for the non-nothing elements the εSoC is 0.
        for n in 1:nEV
            if !isnothing(depIdx[n])
                tDep = Dt[depIdx[n]] # departure time
                ϵSoC[n] = SoCev[n](tDep) .- SoCdep[n]
            end
        end
    end
    # model[:ϵSoC]=[]; # assign name in the model
    # tDep = [Dt[depIdx[n]] for n ∈ 1:nEV] # departure time
    # model[:ϵSoC] = [SoCev[n](tDep[n]) - SoCdep[n] for n ∈ 1:nEV];
    # for day in eachindex(tDep[1,:]) # if there's more than one day loop over them
    #     for n in 1:nEV # loop over the EVs
    #         td = tDep[n, day] # pick value
    #         td = td * 3600 + (24 * 3600 * (day - 1)) # change to secs and add days
    #         # find the nearest td inside supports(t)
    #         td = findmin(abs.(td .- supports(t)))[1]
    #         # Check if the time index is within the bounds of SoCev
    #         if td >= t0 && td <= tend
    #             ϵSoCexpr = SoCev[n](td) - SoCdep[n]
    #             push!(model[:ϵSoC], ϵSoCexpr) # push to the expression vector
    #         end
    #     end
    # end

    if type == "bucket"
    # Model variables
        OCVev0 = [aOCV[n]+bOCV[n]*SoCev0[n] for n in 1:nEV];    
        @variables(model, begin
            # aOCV[n] ≤ OCVev[n ∈ 1:nEV] ≤ aOCV[n]+bOCV[n], Infinite(t),(start=OCVev0[n]) # open circuit voltage of the cell
            vmin[n] ≤ OCVev[n ∈ 1:nEV] ≤ vmax[n], Infinite(t),(start=OCVev0[n]) # open circuit voltage of the cell
            -imax[n] ≤ iev[n ∈ 1:nEV] ≤ imax[n], Infinite(t)  # current per branch 
        end);
    # Initial conditions   
    # Model constraints
       @constraints(model, begin
            [n ∈ 1:nEV], bev⁻[n] .+ γf .≥ 1
            availability[n ∈ 1:nEV], model[:γf][n].*model[:Pev][n] + (1-model[:γf][n]).*Pdrive[n] - model[:PevTot][n] .== 0 # power balance
            [n ∈ 1:nEV], OCVev[n] .== aOCV[n]+bOCV[n]*model[:SoCev][n] # linear voltage model
            # [n ∈ 1:nEV], OCVev[n] == OCVfromSoC(SoCev[n]) # Lookup table voltage model
            [n ∈ 1:nEV], iev[n] .* OCVev[n] .== 1e3*model[:PevTot][n]/Npev[n]/Nsev[n] # current per branch
       end);
    if sets.costWeights[3] != 0 # Aging check
            @variables(model, begin
                -imax[n] ≤ ilossev[n ∈ 1:nEV] ≤ imax[n], Infinite(t) # total aging
                0.8*Qev0[n] .≤ Qev[n ∈ 1:nEV] .≤ Qev0[n], Infinite(t) # cell capacity
            end);
            @constraints(model, begin
                # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev[n]/3600 * 1e3 # Aging Qev
                # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev0[n]/3600 * 1e3 # No-Aging Qev
                [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bev⁻[n] + bev⁺[n])*iev[n]/Qev0[n]/3600 * 1e3
                [n ∈ 1:nEV], Qev[n](t0) .== Qev0;
            end);
        else
            # Static Qbess
            # @constraint(model, [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev0[n]/3600 * 1e3)
            @constraint(model, [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bev⁻[n] + bev⁺[n])*iev[n]/Qev0[n]/3600 * 1e3)
        end
    elseif type == "ECM"
    # Model variables
        R0ev0 = zeros(nEV); R1ev = zeros(nEV); tauev = zeros(nEV);
        vtev0 = zeros(nEV); iR1ev0 = zeros(nEV); ilossev0 = zeros(nEV);
        for n in 1:nEV
            @unpack R0Param, RParam, RCParam, iRn0, vt0 = performanceParams[n]
            R0ev0[n]=R0Param[]; R1ev[n]=RParam[]; tauev[n]=RCParam[];
            vtev0[n]= vt0; iR1ev0[n]=iRn0[1];
            @unpack iloss0 = agingParams[n]
            ilossev0[n]=iloss0;
        end
        OCVev0 = [aOCV[n]+bOCV[n]*SoCev0[n] for n in 1:nEV];
        @variables(model, begin
            # aOCV[n] ≤ OCVev[n ∈ 1:nEV] ≤ aOCV[n]+bOCV[n], Infinite(t), (start=OCVev0[n]) # open circuit voltage of the cell
            # aOCV[n] ≤ vtev[n ∈ 1:nEV] ≤ aOCV[n]+bOCV[n], Infinite(t), (start=OCVev0[n]) # terminal voltage of the cell
            vmin[n] ≤ OCVev[n ∈ 1:nEV] ≤ vmax[n], Infinite(t), (start=OCVev0[n]) # open circuit voltage of the cell
            vmin[n] ≤ vtev[n ∈ 1:nEV] ≤ vmax[n], Infinite(t), (start=OCVev0[n]) # terminal voltage of the cell
            -imax[n] ≤ iev[n ∈ 1:nEV] ≤ imax[n], Infinite(t), (start=0.0) # total current per branch
            -imax[n] ≤ ilossev[n ∈ 1:nEV] ≤ imax[n], Infinite(t), (start=0.0) # total aging
            -imax[n] ≤ iR1ev[n ∈ 1:nEV] ≤ imax[n], Infinite(t), (start=0.0) # pole current
            0.8*Qev0[n] .≤ Qev[n ∈ 1:nEV] .≤ Qev0[n], Infinite(t), (start=Qev0[n]) # cell capacity
            # R0ev0[n] .≤ R0ev[n ∈ 1:nEV] .≤ R0ev0[n]*1.2, Infinite(t), (start=R0ev0[n]) # cell capacity
        end);
    # Initial conditions
        @constraints(model, begin
            [n ∈ 1:nEV], iR1ev[n](t0) == iR1ev0[n];
        end); 
    # Model constraints
        @constraints(model, begin
            [n ∈ 1:nEV], bev⁻[n] .+ γf .≥ 1
            availability[n ∈ 1:nEV], model[:γf][n].*model[:Pev][n] + (1-model[:γf][n]).*Pdrive[n] - model[:PevTot][n] .== 0 # power balance
            [n ∈ 1:nEV], iev[n] .* vtev[n] .== 1e3*model[:PevTot][n]/Npev[n]/Nsev[n] # current per branch. 1e3 to convert kW->W
            # Transition function
            # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev[n]/3600 * 1e3 # Aging Qev
            # [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bPev[n] + (1-bPev[n]))*iev[n]/Qev0[n]/3600 * 1e3 # No-Aging Qev
            [n ∈ 1:nEV], ∂.(model[:SoCev][n], t) * 1e3 .== -(ηev[n] * bev⁻[n] + bev⁺[n])*iev[n]/Qev0[n]/3600 * 1e3
            [n ∈ 1:nEV], ∂.(iR1ev[n], t) * 1e3 .== (-1/tauev[n]*iR1ev[n]+1/tauev[n]*iev[n]) * 1e3 # RC resistor currents
            #Hysterisis?
            # Output equations
            #[n ∈ 1:nEV], OCVev[n] == OCVfromSOCtemp(SoCev, T, data["EV"]) # Lookup table voltage model
            [n ∈ 1:nEV], OCVev[n] .== aOCV[n]+bOCV[n]*model[:SoCev][n] # Linear OCV model
            [n ∈ 1:nEV], vtev[n] .== OCVev[n] - R1ev[n]*iR1ev[n] - R0ev0[n]*iev[n]
            # [n ∈ 1:nEV], vtev[n] .== OCVev[n] - R1ev[n]*iR1ev[n] - R0ev[n]*iev[n]
        end); 
        # elseif type == "PBROM"

    end

    return model;
end

function add_battDeg(model::InfiniteModel, data::BESSData)
# battDeg: Battery degradation modeling function
# This function adds variables and constraints to the model obj following the different
    # Sub-models available:
    # 1. empirical:
    # Wang et al (2014) doi: 10.1016/j.jpowsour.2014.07.030
    # 2. PB Jin: 
    # Jin (2022) doi: 10.1016/j.electacta.2021.139651
    # 3. PB Reniers:
    # Reniers et al (2018) 10.1016/j.jpowsour.2018.01.004
    # 4. PB Plett: 
    # G.L. Plett, Ch.  "Battery Management Systems, Volume I, Battery Modeling," Artech House, 2015.
    # G.L. Plett, Ch. 7 "Battery Management Systems, Volume II, Equivalent Circuit Methods" Artech House, 2016.
    #
    # Some useful refs.:
    # SEI: Solid Electrolyte Interface
    # AM: Active Material
    # c^*s,p: Bulk concentration of solvent reactant/reduction product at equilibrium state
    t=model[:t];
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    ilossbess=model[:ilossbess];
    Qbess=model[:Qbess];
    SoCbess=model[:SoCbess];
    ibess=model[:ibess];
    # bPbess=model[:bPbess];
    bess⁺=model[:bess⁺];
    bess⁻=model[:bess⁻];

    R = 8.314 # Gas constant [J/K/mol]
    T = 25+273 # pack temperature [K]
    F = 96485 # Faraday constant [C/mol]
    
    @unpack GenInfo, PerfParameters, AgingParameters=data
    @unpack initQ, SoHQ = GenInfo
    @unpack type=AgingParameters;
    Qbess0 = initQ*SoHQ;
    if type == "empirical"
    # Empirical from Wang et al (2014)
        @unpack c, initT = AgingParameters
        # ilossCyclebess = c[1]*c[3]/c[4]*ℯ^(c[2]*abs(ibess))*(1-SoCbess)*abs(ibess); # cyclic aging
        ilossCyclebess = c[1]*c[3]/c[4]*ℯ^(c[2]*((bess⁺-bess⁻)*ibess))*(1-SoCbess)*((bess⁺-bess⁻)*ibess);
        ilossCalbess = c[5]*√(initT + t)* ℯ^(-24e3/R/T); # calendar aging
        @constraints(model, begin
            ilossbess .== ilossCyclebess + ilossCalbess # total aging
            ∂.(Qbess, t) .== -ilossbess/3600 # Parameter update, remember unit transf As <-> Ah
            Qbess(t0) == Qbess0 # Initial cell capacity
        end);
    elseif type == "PB Jin"
    # Physics-based from Jin (2022)
    # Most equations and values come from Jin (2022), a small piece comes from Jin (2017) the original modeling paper.
        @unpack_JinAgingParams AgingParameters
        Tref=T; # Reference temperature 
    # SEI layer.
        # Parameters list
        # nSEI: number of electrons transferred in the SEI side-reaction;
        # λ=c^*s.√Ds/c^*p.√Dp=constant
        # OCVs: Open circuit potential of the side reaction [V];
        # OCVn: open circuit voltage of the anode [V] # object coming from Cell.
        # as: Specific surface area of the Anode [m2/m3]
        # A: Active surface area of the Anode [m2]
        # Ln: length of the anode [m]
        # i0: exchange current [A/m2]
        # kSEI: kinetic rate constant [1/√sec]
        # ESEI: activation energy of SEI side-reaction [J/mol]
        # MSEI: molar weight of the SEI layer [kg/mol]
        # ρSEI: density of the SEI layer [kg/m3]
        ηk = 2*R*T/F*asinh(ibess/nSEI/as/An/Ln/i0) # kinetic overpotential    
        z = SoCbess*(z100p-z0p)+z0p
        OCVn = 0.6379+0.5416*ℯ^(-305.5309*z) +
            0.044*tanh(-(z-0.1958)/0.1088) -
            0.1978*tanh((z-1.0571)/0.0854) -
            0.6875*tanh((z+0.0117)/0.0529) -
            0.0175*tanh((z-0.5692)/0.0875)
        θ = ℯ^(nSEI*F/R/T*(ηk+OCVn-OCVs)) # fitting param
        iSEI = (kSEI*ℯ^(-ESEI/R/T))/(nSEI*(1+λ*θ)*√(initT+t))*Qbess0;

        # @expressions(model, begin
        #     ηk, 2*R*T/F*asinh(ibess/nSEI/as/An/Ln/i0) # kinetic overpotential    
        #     z, SoCbess*(z100p-z0p)+z0p
        #     OCVn, 0.6379+0.5416*ℯ^(-305.5309*z) +
        #     0.044*tanh(-(z-0.1958)/0.1088) -
        #     0.1978*tanh((z-1.0571)/0.0854) -
        #     0.6875*tanh((z+0.0117)/0.0529) -
        #     0.0175*tanh((z-0.5692)/0.0875)
        #     θ, ℯ^(nSEI*F/R/T*(ηk+OCVn-OCVs)) # fitting param
        #     iSEI, (kSEI*ℯ^(-ESEI/R/T))/(nSEI*(1+λ*θ)*√(initT+t));    
        # end)

    # Loss of Active Material (AM)
        # Parameter list
        # kAM = kAM⁰/εAM⁰, [1/Ah]
        # EAM: activation energy [J/mol]
        # iAM = kAM*ℯ^(-EAM/R/T)*(SoCbess*100)*(-ibess*bPbess + ibess*(1 - bPbess))*Qbess0;
        iAM = kAM*ℯ^(-EAM/R/T)*(SoCbess*100)*(bbess⁺-bbess⁻)*ibess*Qbess0;

    # Lithium Plating
        # Parameter list:
        # αLi: cathodic transf. coeff. for Li plating
        # εₑ0: initial volume fraction of electrolyte
        # ηLiMin: user define limit/safety margin for no plating [V]
        # t⁺₀: transport/transference number
        # β: fitting parameter
        
        #= All of this should come from the performance model
            Eκ: Activation energy for κ [J/mol]
            EDe: Activation energy for De [J/mol]
            κref: Reference value for κ ionic conductivity at reference temperature [S/m]
            DeRef: Reference value for De at reference temperature [m2/s]
            brug: Bruggeman exponent
            ce_avg: volume-average concentration of Li in the electrolyte [mol/m3]
            σn: electronic conductivity of the electrode [S/m]
            εₛ: volume fraction of solid in the electrolyte
            # particular variables for aging-submodel
            # @variables(model, begin
            #     ηLibess ≥ ηLiMin, Infinite(t) # to avoid Li plating 
            #     1e-4 ≤ εₑbess ≤ 1, Infinite(t) # electrolyte volume fraction
            #     1e-6 ≤ δSEIbess, Infinite(t) # SEI layer thickness [m]
            # end)
            # κ = κref*exp(Eκ/R*(1/Tref-1/T))
            # De = DeRef*exp(EDe/R*(1/Tref-1/T))
            # κeff = κ*εₑbess^brug
            # DeEff = De*εₑbess^brug
            # κDeff = 2*R*T*κeff*(t⁺₀-1)/F
            # σeff = σn*εₛ^brug 
            # Et = -(ibess/κeff/An/Ln)*((β*εₑbess-(1-t⁺₀))*κDeff/ce_avg/DeEff/F+1)
        =#
        @constraints(model, begin
        #=
            # # Lithium overpotential at x=Ln
            # # when charging use ηLi ≥ ηLiMin from Jin (2022)
            # (ηLibess + ibess*Ln/2/σeff/An+Et/3*Ln^2-ηk-OCVn) * bPbess == 0;
            # # when discharging the jLi=0 --> ηLi=OCPLi the eq. potential
            # (ηLibess - 0.2) * (1-bPbess) == 0;
        =#
            ilossbess * 1e5 .== (iSEI + iAM) * 1e5 # total aging. Since ηLi  ≥ ηLiMin --> iLi=0
            ∂.(Qbess, t) * 1e5.== -ilossbess/3600 * 1e5 # Cell capacity, remember unit transf As <-> Ah
            # ∂.(R0bess, t) .== εₛ/κeff*∂.(δSEIbess, t) # Series resistance [Ω]
            # The S leaves the electrolyte to form the SEI layer.
            # ∂.(δSEIbess, t) .== -iSEI*MSEI/nSEI/F/ρSEI/An; # SEI layer thickness [m]
            # ∂.(εₑbess, t).== -as*∂.(δSEIbess, t); 
            # Initial condition
            Qbess(t0) == Qbess0; # cell capacity
            # εₑbess(t0) == εₑ0; # Initial electrolyte volume fraction
            # δSEIbess(t0) == 1e-6; 
            # ηLibess(t0) == 0.2; # Initial overpotential CHECK
        end);
    elseif type == "PB Reniers"
    # Physics-based from Reniers et al (2018)
        @unpack_ReniersAgingParams AgingParameters
        iSEI=An*ℯ^(-nSEI*F/R/T*η_neg)/((nSEI*F*kSEI*ℯ^(-nSEI*F/R/T*(OCVn-OCVs)))^-1+δ/nSEI/F/DSEI)
        # particular variables for aging-submodel
        @variable(model, 0 ≤ δSEIbess, Infinite(t)); # SEI layer thickness [m]
        @constraints(model, begin
            ∂.(δSEIbess, t) .== iSEI*MSEI/nSEI/F/ρSEI
            ilossbess .== iSEIbess # total aging. Since ηLi  ≥ ηLiMin --> iLi=0
            ∂.(Qbess, t) .== -ilossbess/3600 # Cell capacity, remember unit transf As <-> Ah
        end)
    elseif type == "PB Plett"
    # PBROM for SEI, from Plett (2016) Ch.7, BMS ECM Vol II.
        # i0SEI = 1.5e-6; # exchange current of SEI growth [A/m2]
        # θ=θ_min+z_cell*(θ_max-θ_min)
        # A =-i0SEI/F*exp⁡(F*(OCVn-OCVs)/2/R/T);
        # B =-ibess/(2*as*i0*A*Ln)
        # C=F/(2*i0)
        # jSEI=(A(θ)*B(ibess)+A(θ)*√(B(ibess)^2+(1-2*C*A(θ))))/(1-2CA(θ))
        # ∂.(Rfilm,t)=-(M_P*Δt)/(ρ_P.k_P )*jSEI
        # ∂.(Qbess,t)=as*A*F*Ln*jSEI/3600
    
    end
    return model;
end

function add_battDeg(model::InfiniteModel, sets::modelSettings, data::EVData)
    t=model[:t];
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    ilossev=model[:ilossev]; Qev=model[:Qev];
    SoCev=model[:SoCev]; iev=model[:iev]
    # bPev=model[:bPev];
    bev⁻=model[:bev⁻]; bev⁺=model[:bev⁺];
    R = 8.314 # Gas constant [J/K/mol]
    T = 25+273 # pack temperature [K]
    F = 96485.0 # Faraday constant [C/mol]
    nEV=sets.nEV;
    
    @unpack GenInfo, PerfParameters, AgingParameters=data.carBatteryPack
    agingParams=AgingParameters;
    @unpack initQ, SoHQ=GenInfo
    Qev0=initQ*SoHQ;
    @unpack type=agingParams; # change later

    if type == "empirical"
    # Empirical from Wang et al (2014)
        @unpack c, initT = agingParams # should change in the future.
        # Empirical from Wang et al (2014) doi: 10.1016/j.jpowsour.2014.07.030
        # ilossCycleev = [c[1]*c[3]/c[4]*ℯ^(c[2]*abs(iev[n]))*(1-SoCev[n])*abs(iev[n]) for n in 1:nEV]; # cyclic aging
        ilossCycleev = [c[1]*c[3]/c[4]*ℯ^(c[2]*((bev⁺[n]-bev⁻[n])*iev[n]))*(1-SoCev[n])*((bev⁺[n]-bev⁻[n])*iev[n]) for n in 1:nEV];
        ilossCalev = c[5]*√(initT+t)* ℯ^(-24e3/R/T); # calendar aging
        @constraints(model, begin
            [n ∈ 1:nEV], ilossev[n] .== ilossCycleev[n] + ilossCalev # total aging
            [n ∈ 1:nEV], ∂.(Qev[n], t) .== -ilossev[n]/3600 # Parameter update, remember unit transf As <-> Ah
            [n ∈ 1:nEV], Qev[n](t0) .== Qev0[n] # Initial cell capacity
        end);
    elseif type == "PB Jin"
    # Physics-based from Jin (2022)
    # Most equations and values come from Jin (2022), a small piece comes from Jin (2017) the original modeling paper.
        # R0ev=model[:R0ev];
        @unpack_JinAgingParams agingParams # should change in the future.
        Tref=T; # Reference temperature [K]
    # SEI layer.
        # Parameters list
        # nSEI: number of electrons transferred in the SEI side-reaction;
        # λ=c^*s.√Ds/c^*p.√Dp=constant
        # OCVs: Open circuit potential of the side reaction [V];
        # OCVn: open circuit voltage of the anode [V] # object coming from Cell.
        # as: Specific surface area of the Anode [m2/m3]
        # A: Active surface area of the Anode [m2]
        # Ln: length of the anode [m]
        # i0: exchange current [A/m2]
        # kSEI: kinetic rate constant [1/√sec]
        # ESEI: activation energy of SEI side-reaction [J/mol]
        # MSEI: molar weight of the SEI layer [kg/mol]
        # ρSEI: density of the SEI layer [kg/m3]
        ηk=[2*R*T/F*asinh(iev[n]/nSEI/as/An/Ln/i0) for n in 1:nEV] # kinetic overpotential
        z = [SoCev[n]*(z100p-z0p)+z0p for n in 1:nEV]
        OCVn = [0.6379+0.5416*ℯ^(-305.5309*z[n]) +
            0.044*tanh(-(z[n]-0.1958)/0.1088) -
            0.1978*tanh((z[n]-1.0571)/0.0854) -
            0.6875*tanh((z[n]+0.0117)/0.0529) -
            0.0175*tanh((z[n]-0.5692)/0.0875) for n in 1:nEV] 
        θ =[ℯ^(nSEI*F/R/T*(ηk[n]+OCVn[n]-OCVs)) for n in 1:nEV]  # fitting param
        iSEI = [(kSEI*ℯ^(-ESEI/R/T))/(nSEI*(1+λ*θ[n])*√(initT+t))*Qev0[n] for n in 1:nEV];
        
    # Loss of Active Material (AM)
        # Parameter list
        # kAM = kAM⁰/εAM⁰, [1/Ah]
        # EAM: activation energy [J/mol]
        # iAM =[kAM*ℯ^(-EAM/R/T)*(SoCev[n]*100)*(- iev[n]*bPev[n] + iev[n]*(1 .-bPev[n]))*Qev0[n] for n in 1:nEV];
        iAM =[kAM*ℯ^(-EAM/R/T)*(SoCev[n]*100)*(bev⁺[n]-bev⁻[n])*iev[n]*Qev0[n] for n in 1:nEV];

    # Lithium Plating
        # Parameter list:
        # αLi: cathodic transf. coeff. for Li plating
        # εₑ0: initial volume fraction of electrolyte
        # ηLiMin: user define limit/safety margin for no plating [V]
        # t⁺₀: transport/transference number
        # β: fitting parameter
        
        #= All of this should come from the performance model
        Eκ: Activation energy for κ [J/mol]
        EDe: Activation energy for De [J/mol]
        κref: Reference value for κ ionic conductivity at reference temperature [S/m]
        DeRef: Reference value for De at reference temperature [m2/s]
        brug: Bruggeman exponent
        ce_avg: volume-average concentration of Li in the electrolyte [mol/m3]
        σn: electronic conductivity of the electrode [S/m]
        εₛ: volume fraction of solid in the electrolyte
        =#
        # particular variables for aging-submodel
        # @variables(model, begin
        #     ηLiev[1:nEV] ≥ ηLiMin, Infinite(t) # to avoid Li plating 
        #     1e-4 ≤ εₑev[1:nEV] ≤ 1.0, Infinite(t) # electrolyte volume fraction
        #     1e-6 ≤ δSEIev[1:nEV], Infinite(t) # SEI layer thickness [m]
        # end)
        
        # κ = κref*exp(Eκ/R*(1/Tref-1/T))
        # De = DeRef*exp(EDe/R*(1/Tref-1/T))
        # κeff = [κ*εₑev[n]^brug for n in 1:nEV]
        # DeEff = [De*εₑev[n]^brug for n in 1:nEV]
        # κDeff = [2*R*T*κeff[n]*(t⁺₀-1)/F for n in 1:nEV]
        # σeff = σn*εₛ^brug
        
        # Et = [-(iev[n]/κeff[n]/An/Ln)*((κDeff[n]*β*εₑev[n]-(1-t⁺₀))/ce_avg/DeEff[n]/F+1) for n in 1:nEV]
        @constraints(model, begin
            #=
                # Lithium overpotential at x=Ln
                # when charging use ηLi ≥ ηLiMin from Jin (2022)
                [n ∈ 1:nEV], (ηLiev[n] +iev[n]*Ln/2/σeff/An+Et[n]/3*Ln^2-ηk[n]-OCVn[n]) .* bPev[n] .== 0
                # when discharging the jLi=0 --> ηLi=OCPLi the eq. potential
                [n ∈ 1:nEV], (ηLiev[n] .- 0.2) .* (1 .- bPev[n]) .== 0
            =#
            [n ∈ 1:nEV], ilossev[n] * 1e5 .== (iSEI[n] + iAM[n]) * 1e5 # total aging. Since ηLi  ≥ ηLiMin --> iLi=0
            [n ∈ 1:nEV], ∂.(Qev[n], t) * 1e5 .== -ilossev[n]/3600*1e5 # Cell capacity, remember unit transf As <-> Ah
            # [n ∈ 1:nEV], ∂.(R0ev[n], t) .== εₛ/κeff[n]*∂.(δSEIev[n], t) # Series resistance [Ω]
            # [n ∈ 1:nEV], ∂.(δSEIev[n], t) .== -iSEI[n]*MSEI/nSEI/F/ρSEI/An # SEI layer thickness [m]
            # [n ∈ 1:nEV], ∂.(εₑev[n], t).== -as*∂.(δSEIev[n], t); # The S leaves the electrolyte to form the SEI layer.
            # Initial condition
            [n ∈ 1:nEV], Qev[n](t0) .== Qev0[n]; # cell capacity
            # [n ∈ 1:nEV], εₑev[n](t0) .== εₑ0; # Initial electrolyte volume fraction
            # [n ∈ 1:nEV], δSEIev[n](t0) .== 1e-6; 
            # [n ∈ 1:nEV], ηLiev[n](t0) .== 0.2; # Initial overpotential CHECK
        end);
    elseif type == "PB Reniers"
    # Physics-based from Reniers et al (2018)
    @unpack_ReniersAgingParams AgingParameters
        iSEI=[An*ℯ^(-nSEI*F/R/T*η_neg)/((nSEI*F*kSEI*ℯ^(-nSEI*F/R/T*(OCVn[n]-OCVs)))^-1+δSEIev[n]/nSEI/F/DSEI) for n in 1:nEV]
        # particular variables for aging-submodel
        @variable(model, 0 ≤ δSEIev[1:nEV], Infinite(t)); # SEI layer thickness [m]
        @constraints(model, begin
            [n ∈ 1:nEV], ∂.(δSEIev[n], t) .== iSEI[n]*MSEI/nSEI/F/ρSEI
            [n ∈ 1:nEV], ilossev[n] .== iSEI[n] # total aging. Since ηLi  ≥ ηLiMin --> iLi=0
            [n ∈ 1:nEV], ∂.(Qev[n], t) .== -ilossev[n]/3600 # Cell capacity, remember unit transf As <-> Ah
        end)
    elseif type == "PB Plett"
    # PBROM for SEI, from Plett (2016) Ch.7, BMS ECM Vol II.
        # i0SEI = 1.5e-6; # exchange current of SEI growth [A/m2]
        # θ = θ_min+z_cell*(θ_max-θ_min)
        # A = -i0SEI/F*exp⁡(F*(OCVn-OCVs)/2/R/T);
        # B = -ibess/(2*as*i0*A*Ln)
        # C = F/(2*i0)
        # jSEI=(A(θ)*B(ibess)+A(θ)*√(B(ibess)^2+(1-2*C*A(θ))))/(1-2CA(θ))
        # ∂.(Rfilm,t)=-(M_P*Δt)/(ρ_P.k_P )*jSEI
        # ∂.(Qbess,t)=as*A*F*Ln*jSEI
    
    end
    return model;    
end

function add_battDeg(model::InfiniteModel, sets::modelSettings, data::Vector{EVData})
    t=model[:t];
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    ilossev=model[:ilossev]; Qev=model[:Qev];
    SoCev=model[:SoCev]; iev=model[:iev]
    # bPev=model[:bPev];
    bev⁻=model[:bev⁻]; bev⁺=model[:bev⁺];
    R = 8.314 # Gas constant [J/K/mol]
    T = 25+273 # pack temperature [K]
    F = 96485.0 # Faraday constant [C/mol]
    nEV=sets.nEV;
    Qev0=zeros(nEV,1);
    agingParams=Vector{AgingParams}(undef, nEV);
    for n in 1:nEV
        @unpack GenInfo, PerfParameters, AgingParameters=data[n].carBatteryPack
        agingParams[n]=AgingParameters;
        @unpack initQ, SoHQ=GenInfo
        Qev0[n]=initQ*SoHQ;
    end
    @unpack type=agingParams[1]; # change later

    if type == "empirical"
    # Empirical from Wang et al (2014)
        @expression(model, ilossCycleev[n ∈ 1:nEV], 0.); 
        @expression(model, ilossCalev[n ∈ 1:nEV], 0.)
        cev=[]; initTev=zeros(nEV);
        for n ∈ 1:nEV
            @unpack c, initT = agingParams[n] # should change in the future.
            # Empirical from Wang et al (2014) doi: 10.1016/j.jpowsour.2014.07.030
            push!(cev, c); initTev[n] = initT;
        end
        # ilossCycleev = [cev[n][1]*cev[n][3]/cev[n][4]*ℯ^(cev[n][2]*abs(model[:iev][n]))*(1-model[:SoCev][n])*abs(model[:iev][n]) for n ∈ 1:nEV]; # cyclic aging
        ilossCycleev = [cev[n][1]*cev[n][3]/cev[n][4]*ℯ^(cev[n][2]*((bev⁺-bev⁻)*iev[n]))*(1-SoCev[n])*((bev⁺-bev⁻)*iev[n]) for n ∈ 1:nEV];
        ilossCalev = [cev[n][5]*√(initTev[n]+t)* ℯ^(-24e3/R/T) for n ∈ 1:nEV]; # calendar aging
        @constraints(model, begin
            [n ∈ 1:nEV], ilossev[n] .== ilossCycleev[n] + ilossCalev[n] # total aging
            [n ∈ 1:nEV], ∂.(Qev[n], t) .== -ilossev[n]/3600 # Parameter update, remember unit transf As <-> Ah
            [n ∈ 1:nEV], Qev[n](t0) .== Qev0[n] # Initial cell capacity
        end);
    elseif type == "PB Jin"
    # Physics-based from Jin (2022)
    # Most equations and values come from Jin (2022), a small piece comes from Jin (2017) the original modeling paper.
        # R0ev=model[:R0ev];
        initTev = zeros(nEV);
        z100pev = zeros(nEV); z0pev = zeros(nEV)
        nSEIev = zeros(nEV); λev = zeros(nEV)
        OCVsev = zeros(nEV); OCVnev = zeros(nEV)
        asev = zeros(nEV); Anev = zeros(nEV)
        Lnev = zeros(nEV); i0ev = zeros(nEV)
        kSEIev = zeros(nEV); ESEIev = zeros(nEV)
        MSEIev = zeros(nEV); ρSEIev = zeros(nEV)
        kAMev = zeros(nEV); EAMev = zeros(nEV)
        for n ∈ 1:nEV
            @unpack_JinAgingParams agingParams[n] # should change in the future.
            initTev[n] = initT
            z100pev[n] = z100p; z0pev[n] = z0p;
            nSEIev[n] = nSEI; λev[n] = λ;
            OCVsev[n] = OCVs; OCVnev[n] = OCVn;
            asev[n] = as; Anev[n] = An;
            Lnev[n] = Ln; i0ev[n] = i0;
            kSEIev[n] = kSEI; ESEIev[n] = ESEI;
            MSEIev[n] = MSEI; ρSEIev[n] = ρSEI;
            kAMev[n] = kAM; EAMev[n] = EAM;
        end
    # SEI layer.
        # Parameters list
        # nSEI: number of electrons transferred in the SEI side-reaction;
        # λ=c^*s.√Ds/c^*p.√Dp=constant
        # OCVs: Open circuit potential of the side reaction [V];
        # OCVn: open circuit voltage of the anode [V] # object coming from Cell.
        # as: Specific surface area of the Anode [m2/m3]
        # A: Active surface area of the Anode [m2]
        # Ln: length of the anode [m]
        # i0: exchange current [A/m2]
        # kSEI: kinetic rate constant [1/√sec]
        # ESEI: activation energy of SEI side-reaction [J/mol]
        # MSEI: molar weight of the SEI layer [kg/mol]
        # ρSEI: density of the SEI layer [kg/m3]
        ηk=[2*R*T/F*asinh(iev[n]/nSEIev[n]/asev[n]/Anev[n]/Lnev[n]/i0ev[n]) for n ∈ 1:nEV] # kinetic overpotential
        z = [SoCev[n]*(z100pev[n]-z0pev[n])+z0pev[n] for n ∈ 1:nEV]
        OCVn = [0.6379+0.5416*ℯ^(-305.5309*z[n]) +
            0.044*tanh(-(z[n]-0.1958)/0.1088) -
            0.1978*tanh((z[n]-1.0571)/0.0854) -
            0.6875*tanh((z[n]+0.0117)/0.0529) -
            0.0175*tanh((z[n]-0.5692)/0.0875) for n ∈ 1:nEV] 
        θ =[ℯ^(nSEIev[n]*F/R/T*(ηk[n]+OCVn[n]-OCVsev[n])) for n ∈ 1:nEV]  # fitting param
        iSEI = [(kSEIev[n]*ℯ^(-ESEIev[n]/R/T))/(nSEIev[n]*(1+λev[n]*θ[n])*√(initTev[n]+t))*Qev0[n] for n ∈ 1:nEV];
        
    # Loss of Active Material (AM)
        # Parameter list
        # kAM = kAM⁰/εAM⁰, [1/Ah]
        # EAM: activation energy [J/mol]
        # iAM =[kAMev[n]*ℯ^(-EAMev[n]/R/T)*(SoCev[n]*100)*(- iev[n]*bPev[n] + iev[n]*(1 .-bPev[n]))*Qev0[n] for n in 1:nEV];
        iAM =[kAMev[n]*ℯ^(-EAMev[n]/R/T)*(SoCev[n]*100)*(bev⁺[n]-bev⁻[n])*iev[n]*Qev0[n] for n in 1:nEV];

    # Lithium Plating
        # Parameter list:
        # αLi: cathodic transf. coeff. for Li plating
        # εₑ0: initial volume fraction of electrolyte
        # ηLiMin: user define limit/safety margin for no plating [V]
        # t⁺₀: transport/transference number
        # β: fitting parameter
        
        #= All of this should come from the performance model
        Eκ: Activation energy for κ [J/mol]
        EDe: Activation energy for De [J/mol]
        κref: Reference value for κ ionic conductivity at reference temperature [S/m]
        DeRef: Reference value for De at reference temperature [m2/s]
        brug: Bruggeman exponent
        ce_avg: volume-average concentration of Li in the electrolyte [mol/m3]
        σn: electronic conductivity of the electrode [S/m]
        εₛ: volume fraction of solid in the electrolyte
        =#
        # particular variables for aging-submodel
        # @variables(model, begin
        #     ηLiev[1:nEV] ≥ ηLiMin, Infinite(t) # to avoid Li plating 
        #     1e-4 ≤ εₑev[1:nEV] ≤ 1.0, Infinite(t) # electrolyte volume fraction
        #     1e-6 ≤ δSEIev[1:nEV], Infinite(t) # SEI layer thickness [m]
        # end)
        
        # κ = κref*exp(Eκ/R*(1/Tref-1/T))
        # De = DeRef*exp(EDe/R*(1/Tref-1/T))
        # κeff = [κ*εₑev[n]^brug for n in 1:nEV]
        # DeEff = [De*εₑev[n]^brug for n in 1:nEV]
        # κDeff = [2*R*T*κeff[n]*(t⁺₀-1)/F for n in 1:nEV]
        # σeff = σn*εₛ^brug
        
        # Et = [-(iev[n]/κeff[n]/An/Ln)*((κDeff[n]*β*εₑev[n]-(1-t⁺₀))/ce_avg/DeEff[n]/F+1) for n in 1:nEV]
        @constraints(model, begin
            #=
                # Lithium overpotential at x=Ln
                # when charging use ηLi ≥ ηLiMin from Jin (2022)
                [n ∈ 1:nEV], (ηLiev[n] +iev[n]*Ln/2/σeff/An+Et[n]/3*Ln^2-ηk[n]-OCVn[n]) .* bPev[n] .== 0
                # when discharging the jLi=0 --> ηLi=OCPLi the eq. potential
                [n ∈ 1:nEV], (ηLiev[n] .- 0.2) .* (1 .- bPev[n]) .== 0
            =#
            [n ∈ 1:nEV], ilossev[n] * 1e5 .== (iSEI[n] + iAM[n]) * 1e5 # total aging. Since ηLi  ≥ ηLiMin --> iLi=0
            [n ∈ 1:nEV], ∂.(Qev[n], t) * 1e5 .== -ilossev[n]/3600 * 1e5 # Cell capacity, remember unit transf As <-> Ah
            # [n ∈ 1:nEV], ∂.(R0ev[n], t) .== εₛ/κeff[n]*∂.(δSEIev[n], t) # Series resistance [Ω]
            # [n ∈ 1:nEV], ∂.(δSEIev[n], t) .== -iSEI[n]*MSEI/nSEI/F/ρSEI/An # SEI layer thickness [m]
            # [n ∈ 1:nEV], ∂.(εₑev[n], t).== -as*∂.(δSEIev[n], t); # The S leaves the electrolyte to form the SEI layer.
            # Initial condition
            [n ∈ 1:nEV], Qev[n](t0) .== Qev0[n]; # cell capacity
            # [n ∈ 1:nEV], εₑev[n](t0) .== εₑ0; # Initial electrolyte volume fraction
            # [n ∈ 1:nEV], δSEIev[n](t0) .== 1e-6; 
            # [n ∈ 1:nEV], ηLiev[n](t0) .== 0.2; # Initial overpotential CHECK
        end);
    elseif type == "PB Reniers"
    # Physics-based from Reniers et al (2018)
    @unpack_ReniersAgingParams AgingParameters[1]
        iSEI=[An*ℯ^(-nSEI*F/R/T*η_neg)/((nSEI*F*kSEI*ℯ^(-nSEI*F/R/T*(OCVn[n]-OCVs)))^-1+δSEIev[n]/nSEI/F/DSEI) for n in 1:nEV]
        # particular variables for aging-submodel
        @variable(model, 0 ≤ δSEIev[1:nEV], Infinite(t)); # SEI layer thickness [m]
        @constraints(model, begin
            [n ∈ 1:nEV], ∂.(δSEIev[n], t) .== iSEI[n]*MSEI/nSEI/F/ρSEI
            [n ∈ 1:nEV], ilossev[n] .== iSEI[n] # total aging. Since ηLi  ≥ ηLiMin --> iLi=0
            [n ∈ 1:nEV], ∂.(Qev[n], t) .== -ilossev[n]/3600 # Cell capacity, remember unit transf As <-> Ah
        end)
    elseif type == "PB Plett"
    # PBROM for SEI, from Plett (2016) Ch.7, BMS ECM Vol II.
        # i0SEI = 1.5e-6; # exchange current of SEI growth [A/m2]
        # θ=θ_min+z_cell*(θ_max-θ_min)
        # A =-i0SEI/F*exp⁡(F*(OCVn-OCVs)/2/R/T);
        # B =-ibess/(2*as*i0*A*Ln)
        # C=F/(2*i0)
        # jSEI=(A(θ)*B(ibess)+A(θ)*√(B(ibess)^2+(1-2*C*A(θ))))/(1-2CA(θ))
        # ∂.(Rfilm,t)=-(M_P*Δt)/(ρ_P.k_P )*jSEI
        # ∂.(Qbess,t)=as*A*F*Ln*jSEI
    
    end
    return model;    
end

function bess!(model::InfiniteModel, sets::modelSettings, data::Dict) # stationary battery pack
    t=model[:t];
    t0=supports(t)[1]; tend = supports(t)[end];
    Δt = supports(t)[2]-supports(t)[1];
    fs = 3600/Δt; # sampling frequency [1/hr]

    # Extract data
    # General Info
    @unpack PowerLim, P0, SoCLim, SoC0, ηC, termCond = data["BESS"].GenInfo
    PbessMax = PowerLim[2]; # Max power [kW]
    PbessMin = PowerLim[1]; # Min power [kW]
    SoCbessMin = SoCLim[1]; # Min State of Charge [p.u.]
    SoCbessMax = SoCLim[2]; # Max State of Charge [p.u.]
    SoCbess0 = SoC0; # Initial SoC [p.u.]
    ηbess = ηC; # Charger/converter efficiency
    
    # Add variables
    @variables(model, begin
        SoCbessMin ≤ SoCbess ≤ SoCbessMax, Infinite(t) # State of Charge    
        # Pbess, Infinite(t) # output power
        # bPbess, Infinite(t), Bin # Binary variable for output power    
        # PbessNeg ≤ 0, Infinite(t) # Pbess^- in power
        bbess⁺, Infinite(t), Bin
        bbess⁻, Infinite(t), Bin
        0 ≤ PbessPos, Infinite(t) # Pbess^+ out power
        0 ≤ PbessNeg, Infinite(t)
    end);
    
    # Bidirectional power flow, ensuring only export or import
    # PbessNeg + PbessPos .== Pbess
    @expression(model, Pbess, PbessPos * ηbess .- PbessNeg * (1/ηbess))
    @constraints(model, begin
        # Base MPEC 1 Bin
        # PbessNeg * (1/ηbess) + PbessPos * ηbess .== Pbess
        # bPbess*PbessMin ≤ PbessNeg
        # PbessPos ≤ (1-bPbess)*PbessMax
        # Alt 1: MPEC 2 Bin
        bbess⁺ + bbess⁻ .≤ 1
        PbessNeg ≤ - bbess⁻ * PbessMin
        PbessPos ≤ bbess⁺ * PbessMax
        # Alt 2: with ⟂
        # PbessPos ⟂ PbessNeg
        # Initial conditions
        SoCbess(t0) ==  SoCbess0
    end);

    if termCond ≥ 0.
        t1 = t0 + termCond*3600
        @constraint(model, termC, SoCbess(t1) ==  SoCbess(t1+24*3600-Δt)) # periodic condition
    end

    model=add_battPerf(model, sets, data["BESS"]) # Operation model
    # check if aging model is needed
    if sets.costWeights[3] != 0
        model=add_battDeg(model, data["BESS"]) # Aging model
    end
    
    return model;
end

function ev!(model::InfiniteModel, sets::modelSettings, data::Dict) # electric vehicle
    t=model[:t];
    t0=supports(t)[1];
    nEV=sets.nEV;
    
    # Extract data
    # General Info
    # Battery pack data
    PevMax=zeros(nEV,1);
    PevMin=zeros(nEV,1);
    Pev0=zeros(nEV,1);
    SoCevMax=zeros(nEV,1);
    SoCevMin=zeros(nEV,1);
    SoCev0=zeros(nEV,1);
    ηev=zeros(nEV,1);
    
    # for n in nEV
    for n in 1:nEV
        @unpack PowerLim, P0, SoCLim, SoC0, ηC = data["EV"][n].carBatteryPack.GenInfo
        PevMax[n] = PowerLim[2]; # Max power [kW]
        PevMin[n] = PowerLim[1]; # Min power [kW]
        Pev0[n] = P0; # Initial power [kW]
        SoCevMin[n] = SoCLim[1]; # Min State of Charge [p.u.]
        SoCevMax[n] = SoCLim[2]; # Max State of Charge [p.u.]
        SoCev0[n] = SoC0; # Initial SoC [p.u.]
        ηev[n] = ηC; # Charger/converter efficiency
    end
        
    # Add variables
    # Pev > 0 -> out power and Pev < 0 -> in power
    @variables(model, begin
        # Pev[n ∈ 1:nEV], Infinite(t)  # EV charger power 
        # bPev[n ∈ 1:nEV], Infinite(t), Bin  # EV charger power
        PevTot[n ∈ 1:nEV], Infinite(t)  # total power of each EV, driving+V2G
        SoCevMin[n] .≤ SoCev[n ∈ 1:nEV] .≤ SoCevMax[n], Infinite(t) # State of Charge
        # Dummy variables for bidirectional flow
        bev⁺[n ∈ 1:nEV], Infinite(t), Bin
        bev⁻[n ∈ 1:nEV], Infinite(t), Bin
        # PevNeg[n ∈ 1:nEV] ≤ 0, Infinite(t) # Pev^- in power
        0 ≤ PevNeg[n ∈ 1:nEV], Infinite(t)
        0 ≤ PevPos[n ∈ 1:nEV], Infinite(t) # Pev^+ out power        
    end);
    
    # Bidirectional power flow, ensuring only export or import
    # [n ∈ 1:nEV], PevNeg[n] + PevPos[n] == Pev[n]
    @expression(model, Pev[n ∈ 1:nEV], PevPos[n] * ηev[n] .- PevNeg[n] * (1/ηev[n]))
    @constraints(model, begin
        # Base MPEC 1 Bin
        # [n ∈ 1:nEV], PevPos[n] * ηev[n] .- PevNeg[n] * (1/ηev[n]) .== Pev[n]
        # [n ∈ 1:nEV], bPev[n]*PevMin[n] ≤ PevNeg[n]
        # [n ∈ 1:nEV], PevPos[n] ≤ (1-bPev[n])*PevMax[n]
        # Alt 1: MPEC 2 Bin
        [n ∈ 1:nEV], bev⁺[n] .+ bev⁻[n] .≤ 1
        [n ∈ 1:nEV], PevNeg[n] .≤ - bev⁻[n] .* PevMin[n]
        [n ∈ 1:nEV], PevPos[n] .≤ bev⁺[n] .* PevMax[n]
        # Alt 2: with ⟂
        # [n ∈ 1:nEV], PevPos[n] ⟂ PevNeg[n]
        # Initial conditions
        [n ∈ 1:nEV], SoCev[n](t0) == SoCev0[n] 
    end);
    # Operation model
    length(data["EV"]) == 1 ? model=add_battPerf(model, sets, data["EV"][1]) : model=add_battPerf(model, sets, data["EV"])
    if sets.costWeights[3] != 0
        # Aging model
        length(data["EV"]) == 1 ? model=add_battDeg(model, sets, data["EV"][1]) : model=add_battDeg(model, sets, data["EV"])
    end;
    return model;
end