### makeEMSobjs
#  This scripts creates all the device structures necessary for the main EMS model.
#  The basic idea is that each device is a structure with its parameters in it. So first we declare the structure obj.
#  Then we load the params (rated power, fit coeffs., measurement data, etc.). Finally we define the current object.
#  Besides the definition of the objects themselves, more useful functions are defined. These are:
#  - fromSecTo15min: downsamples a time series from seconds to 15 minutes.
#  - processPrices: processes the prices data from EPEX in several formats.
#  - getSeason: gets the season from a month.
#  - getSeasonalProfiles: gets the seasonal profiles for each device. The profiles can be daily or weekly.
#  - availabilityEV: creates availability vectors for each EV in the disc t-domain.
# By: Darío Slaifstein, PhD-student @TU Delft, DCES.
# Branch: agingModeling
# Version: 1.4
# Date: 13/09/2023

# # example from the knapsack problem of the docs
# https://jump.dev/JuMP.jl/stable/tutorials/getting_started/design_patterns_for_larger_models/
# struct EMSObject
#     profit::Float64
#     weight::Float64
#     function EMSObject(profit::Float64, weight::Float64)
#         if weight < 0
#             throw(DomainError("Weight of object cannot be negative"))
#         end
#         return new(profit, weight)
#     end
# end

################ FUNCTION DEFINITIONS ################
function fromSecTo15min(x::Array)
    # this function downsamples a time series from seconds to 15 minutes.
    # Define the length of the original time series (in seconds)
    T_sec = length(x)
    # Define the timestep of the new time series (in minutes)
    Δt=15*60;
    # Create a new dictionary with downsampled values
    # Reshape the time series into a matrix with T_sec rows and N columns
    X = reshape(x, T_sec, :)
    # Compute the average of the values in each 15-minute interval
    Y = [mean(X[i:i+(Δt-1), :], dims=1) for i in 1:Δt:T_sec-(Δt-1)]
    # Reshape the downsampled values into a vector with T_min*N elements
    new_value = vec(hcat(Y...))
    return new_value=[0 ; new_value]
end

function processPrices(df;
    type::String="raw",
    upSampRatio::Int=4, # samples per hour
    profType::String="daily",
    season::String="summer",
    )
    # This function processes the prices data from EPEX in several formats.
    # type:
    # - "raw" raw data from EPEX
    # - "summary". Summary data coming from "makePriceForecast.ipynb" dashboard.
    # The price data is in [€/MWh]
    @assert type ∈ ["raw", "summary"] "Invalid type of data"
    @assert profType ∈ ["daily", "weekly", "biweekly","monthly","yearly"] "Invalid profile type"
    profType == "yearly" ? nothing : @assert season ∈ ["summer", "winter"] "Invalid season";
    
    if type=="raw"        
        select!(df, Not(:"Hour 3B")) # eliminate that weird column
        df=coalesce.(df, 0) # replace missing values with 0
        
        # Now we need to reshape the DataFrame into a timeseries
        priceData = Vector{Float64}();
        for r in 1:nrow(df)
            row=Vector{Float64}(df[r,2:25])
            append!(priceData, row)
        end
    elseif type=="summary"
        priceData=Vector{Float64}(df[:,:mean])
    end
    # Change resolution of the prices. From 1h/sample to 15min/sample
    priceData = repeat(priceData, inner=upSampRatio);
    priceData=[0; priceData]
    if profType == "yearly"
        return priceData = [priceData priceData.*0.95]
    end
    # for the rest of the profiles you continue with the seasonal profiles
    priceData=getSeasonalProfiles(priceData; type=profType, n_samples_per_hour = upSampRatio)[season];
    # for biweekly profiles we need to repeat the weekly profile twice and append the first day to the end
    if profType == "biweekly"
        priceData=repeat(priceData[1:(end-upSampRatio*24)], outer=2); append!(priceData, priceData[1:upSampRatio*24])
    end
    # Buy < sell prices 
    priceData = [priceData priceData.*0.95]
    # priceData=CSV.read("energy_prices.csv", DataFrame);  # old data from Wil
    return priceData;
end;

# define a function to get the season from a month
function getSeason(month::Int)
    if month in 3:5
        return "spring"
    elseif month in 6:8
        return "summer"
    elseif month in 9:11
        return "fall"
    else
        return "winter"
    end
end

function getSeasonalProfiles(data::Vector;
    n_samples_per_hour::Int=4, # number of samples per hour
    type::String="daily", # type of profile to return
    )
    # check the type of profile
    @assert type ∈ ["monthly", "biweekly", "weekly", "daily"] "Invalid profile type"    

    type == "daily" ? n_days=1 : nothing;
    type == "weekly" ? n_days=7 : nothing;
    type == "biweekly" ? n_days=7 : nothing;
    type == "monthly" ? n_days=30 : nothing;
    # Get seasonal profiles for each device
    start_date = DateTime("2023-01-01T00:00:00")
    end_date = DateTime("2023-12-31T24:00:00")
    step = Dates.Minute(60/n_samples_per_hour)    
    date_range = start_date:step:end_date

    seasons = ["winter", "spring", "summer", "fall"]

    profiles = Dict()
    # syncronize the data with the date range using a DataFrame
    data_df = DataFrame(date=date_range, values=data)
    # create a new column with the date only yyyy-mm-dd
    data_df[!, :date_only] = Dates.Date.(data_df.date);
    # create a new column with the time only hh:mm:ss
    data_df[!, :hour_minutes] = string.(Dates.Time.(data_df.date));
    # create a new column with the week number
    data_df[!, :week] = Dates.Week.(data_df.date)
    # create a new column with the day of the week
    data_df[!, :dayofweek] = Dates.dayofweek.(data_df.date)

    for season ∈ seasons
        prof_df=Vector();
        for nd ∈ 1:n_days
            # take the data for a season and the nd-th day of the week
            if type == "monthly"
                red_df = data_df[(getSeason.(month.(data_df.date)) .== season) .& (Dates.dayofmonth.(data_df.date) .== nd), :]
            else
                red_df = data_df[(getSeason.(month.(data_df.date)) .== season) .& (Dates.dayofweek.(data_df.date) .== nd), :]
            end
            # reshape in day x hour
            season_df = unstack(red_df,
                                :hour_minutes, :date_only, :values);
            # drop the first and last day of the season
            season_df = season_df[:,2:end-1];
            # Store the daily profile for the season
            push!(prof_df, mean.(eachrow(season_df[:,2:end]))');
        end
        # first concatenate all the profiles
        profiles[season] = hcat(prof_df...)';
        # copy the first day and push to the end of the profile
        # this is necessary because you always need an extra day to not deplete the ESS
        append!(profiles[season], profiles[season][1:n_samples_per_hour*24]);
    end
    return profiles
end

################ TYPES DEFINITION ################
## Make process and load input data
# Define model settings object
@with_kw mutable struct modelSettings # structure of sets
    nEV::Int64 # number of EVs
    t0::Float64 # initial time [hr]
    Tw::Float64 # time window [hr]
    Δt::Float64 # time length of each step [hr]
    tend::Float64 = t0+Tw # final time [hr]
    dTime::Array = collect(t0:Δt:tend)*3600 # discrete time array [s]
    steps::Int64 # number of times to move the window
    costWeights::Array=[1, 1000, 1e4]; # [Wgrid WSoC Wloss] array with the corresponding weigths for each cost
    # inputs for build_data
    season::String = "summer" # "summer" or "winter"
    profType::String # "daily" or "weekly" or "biweekly"
    loadType::String # "GV" or "mffbas" or "base_models"
    year::Int64 = 2023
    cellID::String = "SYNSANYO" # cell ID for the battery packs
end

#--------------- DEVICE DEFINITIONS ---------------#
# Solar PV data
@with_kw mutable struct SPVData
    RatedPower::Float64 = 270e-3; # rated power of the PV panel [kWp]
    capex::Float64 = 1100; # capital expenditure [USD/kW]
    MPPTData::Array # array of MPPT measurement
end

# Connection def
abstract type ConnectionData end

# Storage asset def
abstract type StorageAssetData end

# BESS data
abstract type OCVParams end

@with_kw mutable struct OCVlinearPerfParams <: OCVParams
# Open Circuit Voltage (OCV) linear model
# OCV=aOCV+bOCV*SOC
    ocvLine::Array{Float64}=[3.2, 0.7]; # params for OCV linear model
end

@with_kw mutable struct OCVnonlinPerfParams <: OCVParams
# Non linear OCV=f(SoC) 
# OCV[Tk]=OCV0[Tk]+OCVrel[Tk]*SOC for Tk in Ts
    # From SOC to OCV
    SOC::Array{Float64} # [p.u.]
    OCV0::Array{Float64} # Open Circuit Voltage offset [V]
    OCVrel::Array{Float64}  # Change of OCV with T [V/K]
    # From OCV to SOC
    # SOC[Tk]=SOC0[Tk]+SOCrel[Tk]*OCV for Tk in Ts
    OCV::Array{Float64} # [p.u.]
    SOC0::Array{Float64} # SOC offset [V]
    SOCrel::Array{Float64}  # Change of SOC with T [V/K]
end

@with_kw mutable struct Generic
# General Info for a batteryPack
    PowerLim::Array{Float64} # Min-Max power [kW]
    P0::Float64 # Initial P [kW]
    SoCLim::Array{Float64} # Min-Max SoC [p.u.]
    SoC0::Float64 # Initial SoC [p.u.]
    termCond::Float64=6. # Termination condition [h]
    # Base parameters for performance sub-models.
    initQ::Float64 # initial cell Capacity in [Ah]
    SoHQ::Float64=1 # State of Health [p.u.]
    SoHR0::Float64=0 # State of Health [Ω]
    Np::Float64 # pack parallel branches
    Ns::Float64 # pack series cells
    η::Float64 # coulombic efficiency [p.u.]
    OCVParam::OCVParams
    vLim::Array{Float64} # Min-Max voltage [V]
    ηC::Float64 # charger/converter efficiency [p.u.]
    # Cost info
    initValue::Float64 # Initial value of the BESS [USD/kWh]
end

abstract type PerfParams
# Performance sub-models available:
    # - Bucket
    # - ECM from Plett (2016).
    # - Blended SS-PBROM from Plett (2015).
    # - SPM-PBROM from Reniers (2018-2022).
end

@with_kw mutable struct bucketPerfParams <: PerfParams
# Bucket model.
    type="bucket";
end

@with_kw mutable struct ECMPerfParams <: PerfParams
# Equivalent Circuit Model (ECM) from Plett (2016).
    type="ECM";
    R0Param::Array{Float64} # Series resistance R0 [Ω]
    RParam::Array{Float64} # Resistances for the poles R1...Rn [Ω]
    RCParam::Array{Float64} # Cont time poles τ=RC0, RC1...RCn [s]
    GParam::Float64=0.0; # Hysteresis parameter for fac=exp(-abs(G*etaik*deltaT/(3600*Q)))
    MParam::Float64=0.0; # Hysteresis parameter module for the h state
    M0Param::Float64=0.0; # Hysteresis parameter module for the si state
    # initial conditions
    vt0::Float64 = 3.2+0.7*0.2 # terminal voltage [V]. Default comes from OCVlinearPerfParams at default SoC0 from GenInfo
    iRn0::Array{Float64}=zeros(size(R0Param)); # current for the poles R1...Rn [Ω]
end

@with_kw mutable struct CIDRAPBROMPerfParams <: PerfParams
# Physics-Based Reduced Order Model made with CIDRA (CIDRA-PBROM) from LiiBRA.jl, Planden (2022).
# Default values from Chen (2020), LG M50 cells.
    # [1] Planden (2022) doi: 10.1016/j.est.2022.105637
    # [2] Chen (2020) doi: 10.1016/j.est.2022.105637
    type="PBROM";
    Cell = Construct("LG M50") # Cell type from LiiBRA.jl
    #=
        # State space matrices
        # Transition matrices
        # xₖ₊₁=A*xₖ+B*uₖ
        A::Matrix{Float64}
        B::Matrix{Float64} 
        # Output matrices
        # yₖ =C*xₖ+D*uₖ
        C::Matrix{Float64}
        D::Matrix{Float64}
        # Electrochemical parameters
        # - z0pₖ, k ∈ [pos,neg] z @ 0% Lithium Concentration
        # - z100pₖ, k ∈ [pos,neg] z @ 100% Lithium Concentration
        # - cs_maxₖ, k ∈ [pos,neg] max electrode/solid concentration
        # - ce0 electrolyte concentration
        # - αₖ k ∈ [pos,neg] charge transfer coefficient
        # - kₖ k ∈ [pos,neg] reaction rate constant
        # - t₀⁺ initial transference number

        # Positive electrode
        z100p_Pos::Float64 = 0.910612 # Theta @ 100% Lithium Concentration
        z0p_Pos::Float64 = 0.0263473 # Theta @ 0% Lithium Concentration
        cs_max_Pos::Float64 = 63104 # Max Electrode Concentration
        α_Pos::Float64 = 0.5 # Alpha Factor
        RFilm_Pos::Float64 = 0 # Film Resistance - Ωm²
        k_Pos::Float64 = 3.54458e-11 # Reaction Rate

        # Negative electrode
        z100p_Neg::Float64 = 0.263849 #0.27 # Theta @ 100% Lithium Concentration
        z0p_Neg::Float64 = 0.853974 #0.9084 # Theta @ 0% Lithium Concentration
        cs_max_Neg::Float64 = 33133 # Max Electrode Concentration
        α_Neg::Float64 = 0.5 # Alpha Factor
        RFilm_Neg::Float64 = 0 # Film Resistance - Ωm²
        k_Neg::Float64 = 6.716047e-12 # Reaction Rate
        
        # Initial Electrolyte Concentration (mol m⁻³)
        ce0::Float64 = 1000
        t₀⁺::Float64 = 0.363 # Initial transference number
    =#
end

abstract type AgingParams
# Parameters for aging sub-models.
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
end

@with_kw mutable struct empAgingParams <: AgingParams
# Parameters for empirical aging sub-model.
# Wang (2014)
    # doi: 10.1016/j.jpowsour.2014.07.030
    # 1. empirical:
    type="empirical";
    c::Array{Float64}= [0.0008, 0.39, 1.035, 50, 14.876/√(24*3600)]; # ref values, cathd: NMC+LMO anode: graphite
    initT::Float64=900.; # initial life time
    iloss0::Float64=0.0; # initial value of the capacity fade current
end

@with_kw mutable struct JinAgingParams <: AgingParams
# Parameters for Jin (2022) aging sub-model.
    # [1] Jin (2022) doi: 10.1016/j.electacta.2021.139651
    # Cell: 18650 NMC/graphite Sanyo
    # [2] Jin (2017) doi: 10.1016/j.jpowsour.2016.12.099
    # Cell: 2.3Ah 26650 LFP cathode
    # [3] Yang (2017) doi: 10.1016/j.jpowsour.2017.05.110
    # graphite anode w/ Ethylene carbonate electrolyte -> SEI made of Ethylene dicarbonate
    # [4] Safari (2009) doi: 10.1149/1.3043429
    # ref. used by Yang (2017)
    # [5] Perkins (2012) doi: 10.1016/j.jpowsour.2012.03.003

    type="PB Jin";
    iloss0::Float64=0.0; # initial value of the capacity fade current
    # SEI layer
    nSEI::Float64=2.0; # number of e⁻ transfered in SEI side reaction
    λ::Float64=5.51e-5; # λ=c*ₛ.√Dₛ/c*ₚ.√Dₚ=constant
    OCVs::Float64=0.; # OCV of the side reaction, cause OCV = 0.4V is included in λ
    # OCV of the anode electrode [V]
    OCVn::Function = (z) ->
        OCVn = @. 0.6379+0.5416*ℯ^(-305.5309*z) +
            0.044*tanh(-(z-0.1958)/0.1088) -
            0.1978*tanh((z-1.0571)/0.0854) -
            0.6875*tanh((z+0.0117)/0.0529) -
            0.0175*tanh((z-0.5692)/0.0875)
    εAM::Float64=0.552; # Active Material Volume Fraction
    Rs::Float64=7.5e-6; # Particle radius [m]
    as::Float64=3*εAM/Rs; # Specific surface area of the Anode
    An::Float64=0.105; # Active surface area of the Anode
    Ln::Float64=50e-6; # Thickness of anode [m]
    i0::Float64=1.5; # exchange current of the intercalation current [A/m2]
    kSEI::Float64=66.85; # kinetic rate [1/√sec] CHECK
    ESEI::Float64=39146.0;# Activation energy [J/mol]
    initT::Float64=900.; # initial life time [s]
    δSEI0::Float64=2.0-9; # initial value of the SEI layer thickness [m]
    # initδSEI::Float64=2.0-4; # initial value of the SEI layer thickness [m]
    
    # Yang (2017) & Safari (2009) [3,4]
    MSEI::Float64=0.162; # molecular weight of the SEI layer [kg/mol]
    ρSEI::Float64=1690.0; # density of the SEI layer [kg/mol]

    # Active Material loss
    z100p=0.9; # full electrode stochiometry
    z0p=0; # empty electrode stochiometry
    kAM::Float64=0.0137; # k⁰AM/εAM₀, kinetic rate [1/Ah] CHECK
    EAM::Float64=39500.0; # Activation energy [J/mol]
    
    # Li plating
    αLi::Float64=0.5;# cathodic transf. coeff. for Li plating CHECK
    εₑ0::Float64=0.34; # initial volume fraction of electrolyte
    ηLiMin::Float64=5e-3; # [V] user defined limit
    β::Float64=1.7; # tuning parameter
    t⁺₀::Float64=0.363; # transport/transference number

    Eκ::Float64= 34700.0; # Activation energy for κ [J/mol]
    EDe::Float64= 34700.0; # Activation energy for De [J/mol]
    κref::Float64=0.174; # Reference ionic conductivity for κ at reference temperature [S/m]
    DeRef::Float64=7.5e-11; # Reference value for De at reference temperature [m2/s]
    brug::Float64=3/2;
    ce_avg::Float64=1000.0; # avg volume concentration of Li in the electrolyte
    ce_max::Float64=1000.0; # max volume concentration of Li in the electrolyte
    # this has to be constant because the plating doesn't have to occure ηLiMin>0
    σn::Float64=100.0; # Electronic conductivity [S/m] [5]
    εₛ::Float64=0.59; # Volume fraction of solid electrolyte [5]
end

@with_kw mutable struct ReniersAgingParams <: AgingParams
# Parameters for Reniers (2018, 19, 22) aging sub-model.
    # [1] Jin (2022) doi: 10.1016/j.electacta.2021.139651
    # Cell: 18650 NMC/graphite Sanyo
    # [2] Jin (2017) doi: 10.1016/j.jpowsour.2016.12.099
    # Cell: 2.3Ah 26650 LFP cathode
    type="PB Reniers"
    iloss0::Float64=0.0; # initial value of the capacity fade current
    # SEI layer
    nSEI::Float64=2; # number of e⁻ transfered in SEI side reaction
    λ::Float64=5.51e-5; # λ=c^*s.√Ds/c^*p.√Dp=constant
    OCVs::Float64=0.4; # OCV of the side reaction  [V]
    # OCV of the anode electrode [V]
    OCVn::Function = (Electrode, z) ->
        Uocp = @. 0.6379+0.5416*ℯ^(-305.5309*z) +
            0.044*tanh(-(z-0.1958)/0.1088) -
            0.1978*tanh((z-1.0571)/0.0854) -
            0.6875*tanh((z+0.0117)/0.0529) -
            0.0175*tanh((z-0.5692)/0.0875)
    εAM::Float64=0.552; # Active Material Volume Fraction
    Rs::Float64=7.5e-6; # Particle radius [m]
    as::Float64=3*εAM/Rs; # Specific surface area of the Anode
    A::Float64=0.105; # Active surface area of the Anode
    Ln::Float64=50e-6; # Thickness of anode [m]
    i0::Float64=1.5; # exchange current of the intercalation current [A/m2]
    kSEI::Float64=66.85; # kinetic rate [1/√sec]
    ESEI::Float64=39146; # Activation energy [J/mol]
    initT::Float64=900.; # initial life time [s]

    # Yang (2017) & Safari (2009) [3,4]
    MSEI::Float64=0.162; # molecular weight of the SEI layer [kg/mol]
    ρSEI::Float64=1690; # density of the SEI layer [kg/mol]

    # Active Material loss
    z100p=0.9; # full electrode stochiometry
    z0p=0; # empty electrode stochiometry
    kAM::Float64=0.0137; # kAM⁰/εAM⁰, kinetic rate [1/Ah]
    EAM::Float64=39500; # Activation energy [J/mol]
    
    # Li plating
    αLi::Float64=0.5;# cathodic transf. coeff. for Li plating CHECK
    εₑ0::Float64=0.34; # initial volume fraction of electrolyte
    ηLiMin::Float64=5e-3; # [V] user defined limit
    β::Float64=1.7; # tuning parameter
    t⁺₀::Float64=0.363; # transport/transference number

    Eκ::Float64= 34700; # Activation energy for κ [J/mol]
    EDe::Float64= 34700; # Activation energy for De [J/mol]
    κref::Float64=0.174; # Reference ionic conductivity for κ at reference temperature [S/m]
    DeRef::Float64=7.5e-11; # Reference value for De at reference temperature [m2/s]
    brug::Float64=3/2;
    ce_avg::Float64=1000; # volume-average concentration of Li in the electrolyte
    ce_max::Float64=1000; # volume-average concentration of Li in the electrolyte
    # this has to be constant because the plating doesn't have to occure ηLiMin>0
    σn::Float64=100; # Electronic conductivity [S/m] [5]
    εₛ::Float64=0.59; # Volume fraction of solid electrolyte [5]
end

@with_kw mutable struct BESSData <: StorageAssetData
    GenInfo::Generic
    PerfParameters::PerfParams
    AgingParameters::AgingParams
    
    cellID::String = "SYNSANYO" # cell ID for the battery packs
    # # Cost info
    # initValue::Float64 # Initial value of the BESS [USD/kWh]
end

# EV data
function rejCriteriaEnergy(μDrive, σDrive, tDep, tArr, Tconn, Ns, Np, vmax, Q0, chgPmax)
    Pdrive = rand(truncated(Normal(μDrive, σDrive); lower = 0.01), length(tDep));
    Ereq = Pdrive .* (tArr .- tDep) # energy required for each session [kWh]
    Emax = chgPmax .* Tconn # max energy that can be charged [kWh]
    Eev = Ns .* Np .* vmax .* Q0 .* 1e-3 * 0.8 # 80% energy of the EV battery [kWh]
    # Check if the energy required is less than the max energy that can be charged
    # if not adjust the demanded Pdrive
    for i ∈ eachindex(Ereq)
        if Ereq[i] > Emax[i] || Ereq[i] > Eev
            # pick the lowest bound between charger and EV battery
            println("Energy required for session $i is too high")
            if Emax[i] < Eev
                # adjust the power to the max that can be charged
                println("Adjusting power to the max that can be charged")
                println("From $(Ereq[i]) kWh to $(Emax[i]) kWh")
                Pdrive[i] = Emax[i] / (tArr[i] - tDep[i])

            else
                println("Adjusting power to the max that can be charged")
                println("From $(Ereq[i]) kWh to $(Eev) kWh")
                Pdrive[i] = Eev / (tArr[i] - tDep[i])
            end
        end
    end
    return Pdrive, Ereq, Emax
end
# EV availability modeling function
function availabilityEV(ns, # number of samples
    fs::Int64 = 4, # samples per hour
    μDrive::Float64 = 3.5, # mean of Pdrive [kW]
    σDrive::Float64 = 1.5, # standard dev. of Pdrive [kW]
    Ns::Int64 = 100, # number of series cells
    Np::Int64 = 25, # number of parallel Branches
    vmax::Float64 = 4.2, # max voltage [V]
    Q0::Float64 = 5.2, # initial capacity [Ah]
    chgPmax::Float64 = 17.5, # EV charger limits [kW]
    type::String = "det", # type of availability
    )
    ndays=Int(floor(ns/fs/24)); # number of days

    # First, we create availability vectors for each EV in the disc t-domain.
    # Get n samples from a uniform distribution with n the number of days
    # tDep = rand(depLims[1]:depLims[2], ndays)
    # tArr = tDep + rand(arrLims[1]:arrLims[2], ndays)
    γ = ones(ns)

    # using the data from Elaadusing Serialization
    # Deserialize the mixture model
    open("../data/gmmElaadFit.dat", "r") do f 
        global gmm = deserialize(f) # Gaussian Mixture Model
    end
    # load the lookup table of the connection times
    μtCon_tarr_df = CSV.read("../data/input/Elaad Data/Data downloaded/mean-session-length-per.csv",
                        DataFrame,silencewarnings=true);
    # sort following the arrival times
    sort!(μtCon_tarr_df, "Arrival Time")
    # add a column with the arrival time in hs
    μtCon_tarr_df.tArr = collect(0:0.5:23.5)
    gmmt = truncated(gmm, 0, 23.5) # truncate the GMM to the limits
    # get the median of the GMM
    if type == "mean"
        tArr  = [mean(gmm) for i ∈ 1:(ndays+1)]
    else
        tArr = rand(gmmt, ndays+1)
    end

    # Interpolate mean session length for arrival times
    tCon_interp = linear_interpolation(μtCon_tarr_df.tArr, μtCon_tarr_df.home)
    tDep = tArr .+ tCon_interp.(tArr)

    # Adjust departure times if outside the limits
    tDep = [td > 23.5 ? td - 23.5 : td for td ∈ tDep]

    # Adjust first arrival and last departure times
    tDep = tDep[1:end-1]
    tArr = tArr[2:end]

    # Ensure departure time is before arrival time and not 0
    tDep = [t == 0 ? 0.5 : t for t in tDep]
    tArr = [tDep[i] > t ? tDep[i] + 1.0 : t for (i, t) in enumerate(tArr)]

    # Create the availability signal
    for day in 1:ndays
        # Determine the time indices corresponding to arrival and departure for this day
        t = collect(0:1/fs:(24-1/fs))
        # get the index of the departure and arrival times
        depIdx = findmin(abs.(tDep[day] .- t))[2]
        arrIdx = findmin(abs.(tArr[day] .- t))[2]
        # modify the index to be in the range of the time series, to avoid modifying supports
        tDep[day] = t[depIdx]
        tArr[day] = t[arrIdx]
        # Correct for the day
        depIdx = depIdx + (day-1)*24*fs
        arrIdx = arrIdx + (day-1)*24*fs
        # Mark the time series as parked during the parked interval for this day
        if depIdx <= arrIdx
            γ[depIdx:arrIdx] .= 0
        end
    end
    # Since we removed the the first element of tArr, 
    # the session length is the same as the first tDep append
    # the first element of tDep in the first position of Tconn
    # Tconn = [tDep[1]; Tconn]; Tconn = Tconn[1:end-1];
    Tconn = [24. + tDep[i] - tArr[i-1] for i ∈ 2:ndays]
    Tconn = [tDep[1]; Tconn]

    # Create the driving signal
    Pdrive, Ereq, Emax = rejCriteriaEnergy(μDrive, σDrive, tDep, tArr, Tconn, Ns, Np, vmax, Q0, chgPmax)

    return γ, tDep, tArr, Pdrive, Tconn, Ereq, Emax;
end

@with_kw mutable struct driveData
    # Driving info
    # Consumed power
    Pdrive::Vector{Float64} # matrix of power consumption ℝ^{days}
    # arrival and departure times
    SoCdep::Float64 # desired SoC for departure [p.u.]
    γ::Vector{Float64} # matrix of availability of the EV ℝ^{t}.
    tDep::Vector{Float64} # matrix of departure times ℝ^{days}
    tArr::Vector{Float64} # matrix of arrival times ℝ^{days}
end

@with_kw mutable struct EVData <: StorageAssetData
    carBatteryPack::BESSData
    driveInfo::driveData
end

# Electro-thermal Converter (ST & HP) data
# no default values for this one
@with_kw mutable struct ElectroThermData
    RatedPower::Float64 # rated thermal power [kWt]
    capex::Float64 # capital expenditure [USD/kW]
    η::Float64 # conversion factor from Electric PV to thermal [p.u.]
end

# Thermal Energy Storage System from @Borg 
@with_kw mutable struct TESSData <: StorageAssetData
    Q::Float64 = 200; # Capacity [kWh]
    PowerLim::Array{Float64}= [-5, 5]; # Min-Max power [kW]
    SoCLim::Array{Float64} = [0.1, 0.95]; # Min-Max SoC [p.u.]
    SoC0::Float64=0.8; # Initial SoC [p.u.]
    η::Float64=0.95; # thermal efficiency [p.u.]
    capex::Float64= 15000; # capital expenditure [USD/kWh]
end

# Power Electronic Interface
@with_kw mutable struct peiData
    RatedPower::Float64=10; # rated power [kW]
    capex::Float64=300; # capital expenditure [USD/kW]
end;

@with_kw mutable struct gridData <: ConnectionData
    PowerLim::Array{Float64} # Max-Min power [kW]
    η::Float64 # multiport-converter efficiency [p.u.]
    λ::Array # energy prices. [buy; sell] 
    loadE::Array{Float64}; # electrical load measurement
    loadTh::Array{Float64}; # thermal load measurement
end

@with_kw mutable struct Costs
    objVal::Float64 # objective value [€]
    wObjVal::Float64 # weighted objective value [w€]
    wCgrid::Array{Float64} # weighted grid cost [w€]
    wpDep::Float64 # weighted EV departure penalty [w€]
    wCloss::Array{Float64} # weighted degradation cost [w€]
    weights::Array{Float64} # cost weights Wgrid, WSoC, Wloss
end

# ################ MODEL CREATION ################

# function build_data(; nEV::Int64=2, # number of EVs in the system
#     season::String = "winter", # "summer" or "winter"
#     profType::String, # "weekly" or "daily"
#     loadType::String, # "GV" or "mffbas" or "base_models"
#     year::Int64,
#     fs::Int64 = 4, # samples per hour
#     cellID::String = "SYNSANYO" # cell ID for the battery packs
#     )
#     # This function builds the data for the optimization problem.
#     # The inputs are:
#     # - season: "summer" or "winter"
#     # - profType: "weekly" or "daily"
#     # - loadType: "GV" or "mffbas" or "base_models"
#     # - year: the year of the simulation
#     # The output is a dictionary that contains the models for all the different devices in the Multi-carrier Energy System.
#     # Check inputs
#     # @assert typeof(nEV) == UnitRange{Int64} "nEV must be a range of integers"
#     @assert typeof(nEV) == Int64 "max nEV must be a integer"
#     @assert season ∈ ["summer", "winter"] "Invalid season"
#     @assert profType ∈ ["biweekly", "weekly", "daily", "yearly"] "Invalid profile type"
#     @assert loadType ∈ ["GV", "mffbas", "base_models"] "Invalid load type"
  
#     ## PV model
#     if loadType == "base_models" # From Joel's Base models
#         MPPT = CSV.read("../Base models/PV_15min.csv", DataFrame, header=false) # measurement of the max. power point tracking
#     else
#         # From Ibrahim's PV model for the NL
#         pvData = matread("../data/input/Arnhem_results_Saffirio_Diab.mat")
#         MPPT=fromSecTo15min(pvData["P_DC_sec"]);
#     end
#     MPPT=Array(MPPT);
#     Npv=17; # number of pv panels
#     if profType ≠ "yearly"
#         # get the daily seasonal profile for the PV
#         MPPT=getSeasonalProfiles(MPPT, type = profType)[season];
#         # for biweekly profiles we need to repeat the weekly profile twice and append the first day to the end
#         if profType == "biweekly" MPPT=repeat(MPPT[1:(end-fs*24)], outer=2); append!(MPPT, MPPT[1:fs*24]) end
#     end
    
#     spvModel = SPVData(MPPTData = MPPT/1e3*Npv);
    
#     ## BESS model
#     # A battery pack model is composed of three parts:
#     # - its general information, contained in the GenInfo type.
#     # - its performance submodel, contained in the PerfParameters type. This describes the SoC and terminal voltage of a cell.
#     # - its ageing submodel, contained in the AgingParameters type. This describes the evolution of the performance sub-model parameters.
    
#     # Read data from the E2 of the ESCtoolbox from Plett.
#         # OCV non linear data model
#         bessOCV=JSON3.read(open("../data/input/cell_models/$cellID-modelocv.json", "r"),
#                     Dict{String, Vector{Float64}});
#         # Dynamic ECM model, with np=1
#         bessDYN=JSON3.read(open("../data/input/cell_models/$cellID-modeldyn-no-hys.json", "r"),
#                     Dict{String, Union{Vector{Float64}, Vector}});
#         # bessOCVmat = matread("../data/input/cell_models/E2model-ocv.mat")
#         # # Dynamic ECM model
#         # bessDYNmat = matread("../data/input/cell_models/E2model.mat")
#         # # Unpack
#         # bessOCVmat=bessOCVmat["model"];
#         # bessDYNmat=bessDYNmat["model"];

#         # Bucket model info
#         Q0 = bessOCV["OCVQ"]; Q0 = mean(Q0); # [Ah/cell]
#         η = bessOCV["OCVeta"]; η = mean(η); # [p.u.]

#         # ECM info
#         T=bessDYN["temps"]; T=T[1,:]; # [C]
#         R0Param=bessDYN["R0Param"]; # [Ohm]
#         RParam=bessDYN["RParam"]; # [Ohm]
#         RCParam=bessDYN["RCParam"]; # τ [s]
#         Tind25=findall(x->x==25, T)
    
#     # Heliox 43kWh battery pack seems a little bit too much lets use a 20kWh pack.
#     PowerLim = [-17, 17]; P0 = 0; # Power limits and initial condition [kW]. 
#     SoCLim = [0.2, 0.95]; SoC0 = 0.9; # SoC limits and initial condition [p.u.].
#     # State of Health
#     SoHQ = 1; SoHR0 = R0Param[Tind25][1];
#     if cellID == "SYNSANYO"
#         Np = 10; Ns = 100; # Branches in parallel and series cells per branch.
#         ocv_params=OCVlinearPerfParams()
#         # Ageing submodel
#         # aging_params=empAgingParams();
#         aging_params=JinAgingParams();
#     elseif cellID == "A123"
#         Np = 25; Ns = 110; # Branches in parallel and series cells per branch.
#         # first option, 
#         # ocv_params.ocvLine = [2.5, 1.1]
#         # second option from 20% to 95% SoC
#         ocv_params=OCVlinearPerfParams(ocvLine = [3.2, 0.2105])
#         # Ageing submodel
#         # aging_params=empAgingParams();
#         Cell = Construct("A123");
#         aging_params=JinAgingParams(Rs = Cell.Neg.Rs,
#                             An = Cell.Const.CC_A,
#                             Ln = Cell.Neg.L,
#                             z100p = Cell.Neg.θ_100,
#                             z0p = Cell.Neg.θ_0,
#                             # εₑ0 = Cell.Neg.ϵ_e,
#                             t⁺₀ = Cell.Const.t_plus,
#                             # DeRef= Cell.Neg.De,
#                             ce_avg = Cell.Const.ce0,
#                             # ce_max=Cell.Const.ce0,
#                             σn = Cell.Neg.σ,
#                             εₛ = Cell.Neg.ϵ_s, # check
#                             );
#     end
#     initVal=494.246; # Cost info [USD/kWh]
#     # General info definition
#     gen_params=Generic(PowerLim, P0, SoCLim, SoC0, Q0, SoHQ, SoHR0, Np, Ns, η, ocv_params, initVal);
#     # Performance submodel
#     # perf_params=bucketPerfParams();
#     perf_params=ECMPerfParams(R0Param=R0Param[Tind25], RParam=RParam[Tind25], RCParam=RCParam[Tind25])
    
#     # wrap everything in a BESSData type
#     # battModel = BESSData(gen_params, perf_params, aging_params, initVal);
#     battModel = BESSData(gen_params, perf_params, aging_params, cellID);

#     ## EV MODEL
#     # An electric vehicle model is composed of two parts:
#     # - its battery pack, contained in a BESSData type.
#     # - its driving submodel, contained in the driveData type. This describes the availability, 
#     # times of departure and arrival, and the reference SoC of the EV.
#     # Battery pack definition
#     PowerLim = [-12.5, 12.5] #check
#     # The battery pack has to be around 400Vdc and 50kWh
#     if cellID == "SYNSANYO"
#         Ns = 100; Np = 25;
#     elseif cellID == "A123"
#         Ns = 110; Np = 61;
#     end
#     P0=[0, 0]; # initial
#     Q0 = Q0.*ones(nEV); # [Ah/cell]
#     SoC0 = [0.6, 0.8]; # initial
#     vt0=[gen_params.OCVParam.ocvLine[1]+gen_params.OCVParam.ocvLine[2].*SoC0[n] for n ∈ 1:nEV]
#     gen_params = [Generic(PowerLim, P0[n], SoCLim, SoC0[n], Q0[n], SoHQ, SoHR0, Np, Ns, η, ocv_params, initVal) for n ∈ 1:nEV]
#     perf_params=[ECMPerfParams(R0Param=R0Param[Tind25], RParam=RParam[Tind25],
#                 RCParam=RCParam[Tind25], vt0=vt0[n]) for n ∈ 1:nEV]
#     batteryPack = [BESSData(gen_params[n], perf_params[n], aging_params, cellID) for n ∈ 1:nEV]

#     # Driving information definition
#     μD = 3.5; σD=1.5; depParam=[6, 12]; arrParam=[1, 8]; # Parameters for the Gaussian distributions
#     # availability
#     av = [availabilityEV(length(MPPT), 4) for n ∈ 1:nEV];
#     refSoC=[0.85, 0.85]; # user requirement
#     drive_info = [driveData(μD, σD, depParam, arrParam, refSoC[n], av[n][1], av[n][2], av[n][3]) for n in 1:nEV]
#     # wrap everything in a EVData type
#     evModel = [EVData(batteryPack[n], drive_info[n]) for n in 1:nEV]

#     # Solar thermal model
#     Pn=0.5; capex= 1500; ηST = 0.6;
#     stModel = ElectroThermData(Pn, capex, ηST);

#     # TESS model
#     tessModel = TESSData();
#     # Power Electronic Interface model
#     peiModel = peiData();

#     ## Grid Model
#     # The grid is represented by:
#     # - λ energy prices. [buy; sell] 
#     # - loadE electrical load measurement
#     # - loadTh thermal load measurement
#     # for anual profiles the length is (365*24*fs)+1=35041.

#     ## pick file paths, read and convert to array
#     # for the prices
#     # Raw EPEX FTP server data
#     pricePath = "../data/input/EPEX/auction_spot_prices_netherlands_$year.csv"
#     spotPriceDF = CSV.read(pricePath, DataFrame, delim=',', header=2);

#     #= for IECON 2023 =#
#         # summary data
#         # pricePath = "../data/input/EPEX/summaryPrice_Winter.csv";
#         # spotPriceDF = CSV.read(pricePath, DataFrame, delim=',', header=1);
#         # pricePath = "../data/input/EPEX/summaryPrice_Summer.csv";
#         # spotPriceDF = CSV.read(pricePath, DataFrame, delim=',', header=1);

#     # for the electrical load
#     if loadType == "GV"
#         # processed synthezided load profile 1 year
#         loadEPath = "../data/input/GV/Load_1.csv";
#         loadE = CSV.read(loadEPath, DataFrame, header = false) # electric load
#         # loadE = CSV.read("../data/input/GV/Load_1.csv", DataFrame, header=false) # electric load
#         loadE= Vector(loadE[!,1]);
#     elseif loadType == "mffbas"
#         # From Market Facilitation Forum (MFF) and the Beheerder Afspraken Stelsel (BAS) i.e. mffbas
#         loadEPath = "../data/input/mffbas/summaryE1_$year.csv"
#         #= for IECON 2023 =#
#         # loadEPath = "../data/input/mffbas/summaryE1_mean_Winter.csv"
#         # loadEPath = "../data/input/mffbas/summaryE1_mean_Summer.csv"
#         loadE = CSV.read(loadEPath, DataFrame) # electric load
#         loadE= Vector(loadE[!,2]); 
        
#     elseif loadType == "base_models" # Joel's base models
#         loadEPath = "../data/input/Base models/Load_Profile_15min.csv";
#         loadE = CSV.read(loadEPath, DataFrame, transpose=true, header=false)
#         loadE= Vector(loadE[!,1]); 
#     end
#     # normalize loadE
#     loadE=loadE./maximum(loadE);
#     # peak of 2kW
#     loadE=loadE*2.5;
#     # for the thermal load
#     loadThPath= "../Base models/Thermal_load_15min.csv";    
#     loadTh = CSV.read(loadThPath, DataFrame, header = false) # thermal load
#     loadTh= Vector(loadTh[!,1]); loadTh=loadTh*1e-3;
    
#     # check length of loadE, loadTh
#     length(loadE) == (365*24*fs+1) ? nothing : loadE = [loadE[1]; loadE];

#     # Data processing (upsampling, seasonal patterns, etc.)
#     if profType ≠ "yearly"
#         priceData = processPrices(spotPriceDF; type="raw", profType = profType, season = season);
#         # get the profile for the loads
#         loadE=getSeasonalProfiles(loadE, type = profType)[season];
#         loadTh=getSeasonalProfiles(loadTh, type = profType)[season];
#         # for biweekly profiles we need to repeat the weekly profile twice and append the first day to the end
#         if profType == "biweekly"
#             loadE=repeat(loadE[1:(end-fs*24)], outer=2); append!(loadE, loadE[1:fs*24])
#             loadTh=repeat(loadTh[1:(end-fs*24)], outer=2); append!(loadTh, loadTh[1:fs*24])
#         end
#     else
#         priceData = processPrices(spotPriceDF; type="raw", profType = profType);
#     end
#     # Grid connection limits
#     gridModel = gridData([-17, 17], priceData, loadE, loadTh);
    
#     # Heat pump model
#     # only a uniderectional (heating) HP for now
#     Pn=4; capex= 500; ηHP = 3; 
#     # The initial condition follow the power balance
#     # Phpe0=(loadTh[1]-MPPT[1]*ηST)/ηHP;
#     hpModel = ElectroThermData(Pn, capex, ηHP);

#     data=Dict("SPV"=>spvModel, "BESS"=>battModel, "EV"=>evModel,
#         "ST"=>stModel, "HP"=>hpModel, "TESS"=>tessModel,
#         "grid"=>gridModel, "PEI"=>peiModel);
#     return data
# end