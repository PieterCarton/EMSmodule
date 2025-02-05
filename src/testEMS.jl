### Test Solution Feasibility
# This script tests the feasibility of the solution obtained from the EMS optimization problem.
using LinearAlgebra

function testPBalance(results::Dict, data::Dict)
    # Unwrap the result dict
    t=results[:"t"]; 
    t0=t[1]; Δt=t[2]-t[1];
    # make new time window
    it0 = round(Int,(t0/Δt));
    itend = it0+length(t)-1;

    Ppv=data[:"SPV"].MPPTData[it0:itend];
    Pst=Ppv*data[:"ST"].η;
    Phpe=results[:"Phpe"]; ηHP=data[:"HP"].η;
    Ptess=results[:"Ptess"];
    Plt=data[:"grid"].loadTh[it0:itend]

    Pg = results[:"Pg"];
    Ple = data[:"grid"].loadE[it0:itend];;
    # BESS
    Pbess=results[:"Pbess"];
    # EV
    nEV = length(data["EV"])
    if length(data["EV"]) !=1
        γ_cont = results[:"γ_cont"]
        Pev = [results["Pev[$n]"] for n ∈ 1:nEV];
    else
        γ_cont = results[:"γ_cont"][:]
        Pev = [results[:"Pev[1]"]]
    end
    # Electrical balance
    println("Electrical balance is:")
    if nEV !=1
        εE = Pg + Pbess + sum([Pev[n] .* γ_cont[n] for n ∈ 1:nEV]) + Ppv - Ple - Phpe;
        checkFeasibility(Pg + Pbess + sum([Pev[n] .* γ_cont[n] for n ∈ 1:nEV]) + Ppv, Ple + Phpe) # check
    else
        εE = Pg + Pbess + sum([Pev[n] .* γ_cont for n ∈ 1:nEV]) + Ppv - Ple - Phpe;
        checkFeasibility(Pg + Pbess + sum([Pev[n] .* γ_cont for n ∈ 1:nEV]) + Ppv, Ple + Phpe) # check
    end

    # Thermal balance
    εTh = Pst+Phpe.*ηHP+Ptess - Plt
    println("Thermal balance is:")
    checkFeasibility(Pst+Phpe.*ηHP+Ptess, Plt) # check

    # Plot balances
    CairoMakie.activate!()
    ff = Figure(size=(600, 300))
    ax = Axis(ff[1,1], xlabel = L"$t$ [h]", ylabel = L"$error$ [kW]",
            ytickformat = values -> [L"%$(value)" for value in values],
            xtickformat = values -> [L"%$(value)" for value in values])
    stairs!(ax, t./ 3600, εE, label = L"\varepsilon_{E}")
    stairs!(ax, t./ 3600, εTh, label = L"\varepsilon_{Th}")
    axislegend(ax, position = :rb)
    return ff
end

function testP(results::Dict, data::Dict)
    # Extract the optimal power time-series for the storage device
    t = results["t"]
    status = true

    # BESS
    ibess = results["ibess"];
    vtbess = results["vtbess"];
    Nsbess = data["BESS"].GenInfo.Ns;
    Npbess = data["BESS"].GenInfo.Np;
    Pbess_P = ibess.*vtbess.*Nsbess.*Npbess.*1e-3
    if !checkFeasibility(Pbess_P, results["Pbess"])
        status = false
        p1 = plot(t ./ 3600, Pbess_P, label = "P from i")
        plot!(t ./ 3600, results["Pbess"], label = "P^*")
    end

    # EV
    numEVs = length(data["EV"])
    for n in 1:numEVs
        iev = results["iev[$n]"]
        vtev = results["vtev[$n]"];
        Nsev = data["EV"][n].carBatteryPack.GenInfo.Ns;
        Npev = data["EV"][n].carBatteryPack.GenInfo.Np;
        Pev_P = iev.*vtev.*Nsev.*Npev.*1e-3
        if !checkFeasibility(Pbess_P, results["Pev[$n]"])
            status = false
            p1 = plot(t ./ 3600, Pev_P, label = "P from i")
            plot!(t ./ 3600, results["Pev[$n]"], label = "P^*")
        end
    end

end

function calculateSoCFromP(i, Δt, η, Q, SoC0)
    SoC_P = cumsum(-i) .* Δt .* η ./ Q / 3600 .+ SoC0
    return SoC_P
end

function checkFeasibility(ŝ, s)
    ε = norm(ŝ - s, 2) / norm(s, 2) # relative error
    if ε > 5e-2
        println("Not feasible, relative error is $(ε*100)%")
        return false
    else
        println("Feasible")
        return true
    end
end

function checkSoCFeasibility(ŝ, s)
    ε = norm(ŝ - s, 2) / norm(s, 2) # relative error
    if ε > 5e-2
        println("The Transition of SoC is not feasible")
        return false
    else
        println("The Transition of SoC is feasible")
        return true
    end
end

function checkBounds(ŝ)
    if any(ŝ .> 1.0) || any(ŝ .< 0.0)
        println("The SoC is out of bounds")
        return false
    else
        println("The SoC is within bounds")
        return true
    end
end

function testSoC(results::Dict, data::Dict)
    # Extract the optimal power time-series for the storage device
    t = results["t"]
    Δt = t[2] - t[1]
    status = true

    # BESS
    ibess = results["ibess"]; SoCbess = results["SoCbess"]; Qbess = results["Qbess"]
    ηbess = data["BESS"].GenInfo.η
    SoCbess_P = calculateSoCFromP(ibess, Δt, ηbess, Qbess, SoCbess[1])
    println("BESS")
    if !checkSoCFeasibility(SoCbess_P, SoCbess)
        status = false
        checkBounds(SoCbess_P)
    end

    # EV
    if haskey(data, "EV")
        numEVs = length(data["EV"])
        for n in 1:numEVs
            iev = results["iev[$n]"]; SoCev = results["SoCev[$n]"]; Qev = results["Qev[$n]"]
            ηev = data["EV"][n].carBatteryPack.GenInfo.η
            SoCev_P = calculateSoCFromP(iev, Δt, ηev, Qev, SoCev[1])
            println("EV $n")
            if !checkSoCFeasibility(SoCev_P, SoCev)
                status = false
                checkBounds(SoCev_P)
            end
        end
    end
    
    # TESS
    Ptess = results["Ptess"]; SoCtess = results["SoCtess"];
    Qtess = data["TESS"].Q; ηtess = data["TESS"].η;
    SoCtess_P = calculateSoCFromP(Ptess, Δt, ηtess, Qtess, SoCtess[1])
    println("TESS")
    if !checkSoCFeasibility(SoCtess_P, SoCtess)
        status = false
        checkBounds(SoCtess_P)
    end

    return status
    
end

function plotESSsummary(idx, data::Dict, i_key::String, SoC_key::String, Q_key::String, ff)
    if idx == "TESS"
        η = data["TESS"].η
        Q = data["TESS"].Q
        i_SoC = diff(-results[SoC_key]) ./ η ./ Δt * 3600 .* Q;
    elseif idx == "BESS"
        η = data["BESS"].GenInfo.η
        Q = results[Q_key]
        i_SoC = diff(-results[SoC_key]) ./ η ./ Δt * 3600 .* Q[2:end];
    else
        # check if the idx contains a number to indicate the EV number
        isnumeric(idx[4]) ? n = parse(Int64, idx[4]) : nothing
        η = data["EV"][n].carBatteryPack.GenInfo.η
        Q = results[Q_key]
        i_SoC = diff(-results[SoC_key]) ./ η ./ Δt * 3600 .* Q[2:end];
    end
    
    iopt = results[i_key];
    SoC_true = cumsum(-iopt .* η) .* Δt / 3600 ./ Q .+ results[SoC_key][1];
    ax1=Axis(ff, xlabel = L"$t$ [h]", ylabel = L"$i$ [A]")
    ax2=Axis(ff, xlabel = L"$t$ [h]", ylabel = L"$SoC$ [p.u]", yaxisposition=:right)
    colors=ColorSchemes.tab10.colors;
    # hide some things from the secondary axis
    hidespines!(ax2)
    hidexdecorations!(ax2)
    stairs!(ax1, t ./ 3600, iopt, label = L"i^*", step=:post, color=colors[1])
    stairs!(ax1, t[1:end-1] ./ 3600, i_SoC, label = L"i^{SoC}", step=:post, color=colors[2])
    axislegend(ax1, position = :lt)
    stairs!(ax2, t ./ 3600, SoC_true, label = L"Soc^{true}", step=:post, color=colors[3], ylims=(0,1.1))
    stairs!(ax2, t ./ 3600, results[SoC_key], label = L"SoC^*", step=:post, color=colors[4], ylims=(0,1.1))
    axislegend(ax2, position = :rb)
    return ff
end

function plotTestSoCeESS(results::Dict, data::Dict, key_list)
    t = results["t"]
    Δt = t[2] - t[1]
    CairoMakie.activate!()
    fig=Figure(size = (700, 900))
    i=0;
    # loop over given keys
    for key ∈ key_list
        i += 1;
        # check if key exists in data Dict
        if haskey(data, key)
            # check and make, if key is BESS, EV or TESS
            if key == "BESS"
                kk="bess"
            elseif key == "EV"
                numEVs=length(data["EV"])
                kk = ["ev[$n]" for n ∈ 1:numEVs];
            elseif key == "TESS"
                kk = "tess"
            end
            # use key to generate the subplots
            for k ∈ kk
                i_key = "i$k"; SoC_key = "SoC$k"; Q_key = "Q$k";
                # change for TESS
                i_key == "itess" ? i_key = "Ptess" : nothing
                plotESSsummary(key, data, i_key, SoC_key, Q_key, fig[i,1])
                i += 1;
            end
        end
    end
    return fig
end

function calcFullObj(result, data, s)
    # objective function calculation
    Wgrid, WSoC, Wloss=s.costWeights; # cost weight grid, SoC @tdep for the EVs, degradation
    nEV=s.nEV; # number of EVs
    t=result["t"]; # time
    Δt=t[2]-t[1]; # time step
    t0=t[1];
    # make new time window
    it0 = Int((t0/Δt) + 1);
    itend = it0+length(t)-1;
    # total cost
    # Wgrid*∫(cgrid,t)+pDep+Wloss*clossbess*∫(model[:ilossbess]/3600
    # Grid cost
    λbuy=data["grid"].λ[:,1]; # buy price [€/kWs]
    λsell=data["grid"].λ[:,2]; # sell price [€/kWs]
    Pgpos=result[:"PgPos"]; # positive grid power
    Pgneg=result[:"PgNeg"]; # negative grid power
    Cgrid=cumsum(Pgpos.*λbuy[it0:itend] - Pgneg.*λsell[it0:itend]).*Δt;
    
    # Degradation cost
    Qlossbess=result[:"Qbess"][1].-result[:"Qbess"];
    if length(data["EV"]) != 1
        Qlossev=[result[:"Qev[$n]"][1].-result[:"Qev[$n]"] for n ∈ 1:nEV];
        Qloss = data["BESS"].GenInfo.Ns * data["BESS"].GenInfo.Np * Qlossbess .+ # Ns * Np * Qbess
                    sum(data["EV"][n].carBatteryPack.GenInfo.Ns * # Ns * Np * Qev
                    data["EV"][n].carBatteryPack.GenInfo.Np * Qlossev[n] for n ∈ 1:nEV);
    else
        Qlossev=result[:"Qev[1]"][1].-result[:"Qev[1]"];
        Qloss = data["BESS"].GenInfo.Ns * data["BESS"].GenInfo.Np * Qlossbess .+
                data["EV"][n].carBatteryPack.GenInfo.Ns *
                data["EV"][n].carBatteryPack.GenInfo.Np * Qlossev;
    end    
    closs = 1.2; # [€/Ah]
    Closs = closs .* Qloss; # [€]

    # EV charging penalty
    SoCref=[data["EV"][n].driveInfo.SoCdep for n ∈ 1:nEV]; # desired SoC
    # EV SoC
    if length(data["EV"]) != 1
        SoCev=[result[:"SoCev[$n]"] for n ∈ 1:nEV];
        tdep=[data["EV"][n].driveInfo.tDep for n ∈ 1:nEV]; # departure times
        tdep = [tdep[n] .+ (0:length(tdep[n]) .- 1) .* 24 for n ∈ 1:nEV].*3600
        # find the day being simulated from t 
        day = Int(floor(t[1]/(24*3600)))+1; # this assumes that the time window is smaller than a day
        idtdep = [t.== tdep[n][day] for n ∈ 1:nEV]; # time index for departure
        ϵSoC=[SoCev[n][idtdep[n]].-SoCref[n] for n ∈ 1:nEV]; # SoC penalty
    else
        SoCev=result[:"SoCev[1]"];
        tdep=[data["EV"][n].driveInfo.tDep for n ∈ 1:nEV][1]; # departure times
        tdep = (tdep .+ (0:length(tdep) .- 1.0) .* 24) .*3600
        # find the day being simulated from t 
        day = Int(floor(t[1]/(24*3600)))+1; # this assumes that the time window is smaller than a day
        idtdep = t.== tdep[day]; # time index for departure
        ϵSoC=SoCev[idtdep] .- SoCref; # SoC penalty
    end
    
    ϵSoC=reduce(vcat,ϵSoC);
    pDep=WSoC*sum(ϵSoC.^2.0);
    # Total cost
    summary=Costs(Cgrid[end]+pDep+Closs[end],
                Wgrid*Cgrid[end]+pDep+Wloss*Closs[end],
                Wgrid*Cgrid, 
                pDep, # already weighted
                Wloss.*Closs,
                s.costWeights)  
    return summary
end
