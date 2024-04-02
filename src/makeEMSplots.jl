## make EMS Plots
# This script is used to make the plots for the EMS. It is called by the wrapper script, wrapperEMS.ipynb.
# The script is called with the JSON file name con containing the data to be plotted.


function makeInputsplot(gridModel::gridData, spvModel::SPVData)
        CairoMakie.activate!(type="svg")
        set_theme!(theme_latexfonts())
        f=Figure(size=(800, 400))
        colors=ColorSchemes.tab10.colors;
        ax1=Axis(f[1,1];
                xlabel=L"$t$ [hr]", ylabel=L"$P$ [kW]",
                )
        ax2=Axis(f[1,1];
                ylabel=L"\textrm{Prices}\ [€/MWh]",
                yaxisposition=:right,
                )

        # hide some things from the secondary axis
        hidespines!(ax2)
        hidexdecorations!(ax2)
        # Primary axis - Generation and demand
        stairs!(ax1,(1:length(gridModel.loadE))./4,gridModel.loadE, label=L"$P_{\text{load}^{\text{e}}}$ [kW]", step=:post, color=colors[1])
        stairs!(ax1,(1:length(gridModel.loadTh))./4,gridModel.loadTh, label=L"$P_{\text{load}^{\text{th}}}$ [kWth]", step=:post, color=colors[2])
        stairs!(ax1,(1:length(spvModel.MPPTData))./4, spvModel.MPPTData, label=L"$P_{PV}$ [kW]", step=:post, color=colors[3])
        axislegend(ax1; position=:lt)

        # Secondary axis - prices
        stairs!(ax2, (1:length(gridModel.λ[:,1]))./4, gridModel.λ[:,1],
                label=L"\textrm{Day-ahead Prices [€/MWh]}", step=:post, color=colors[4])
                # move 
        Makie.ylims!(ax2, 0, 1.1 .* maximum(gridModel.λ[:,1]))
        axislegend(ax2; position=:rt)
        return f
end

function makeEBplot(results::Dict, data::Dict)
        # Electricity balance plot for a single EMS run or a concatenated EMS run.
        CairoMakie.activate!()
        set_theme!(theme_latexfonts())
        fig = Figure(size=(1500, 800), fontsize=17)
        ax = Axis(fig[1, 1]; ylabel=L"$P$ [kW]", title=L"\text{Power Balance}",
                )
        colors=ColorSchemes.Dark2_7.colors;

        # Unwrap the result dict
        t=results[:"t"]; 
        t0=t[1]; Δt=t[2]-t[1];
        # make new time window
        # it0 = round(Int,(t0/Δt) + 1)
        it0 = round(Int,(t0/Δt));
        itend = it0+length(t)-1;
        # Inputs
        spvModel = data[:"SPV"];
        gridModel = data[:"grid"];
        # Results
        Pg=convert.(Float64, results[:"Pg"]);
        haskey(results, "Phpe") ? Phpe = convert.(Float64, results[:"Phpe"]) : nothing
        # BESSresults
        Pbess=convert.(Float64, results[:"Pbess"]);
        # EV
        nEV = length(data["EV"])
        if length(data["EV"]) !=1
                γ_cont = results[:"γ_cont"]
                Pev = convert.(Vector{Float64}, [results[:"Pev[1]"], results[:"Pev[2]"]]);
        else
                γ_cont = results[:"γ_cont"][:]
                Pev = convert.(Float64, results[:"Pev[1]"])
        end
        
        # Electric Power balance
        stairs!(ax, t/3600, gridModel.loadE[it0:itend], label=L"P_{\textrm{load}}^{e}", step=:post, color=colors[1], linewidth=1.5)
        stairs!(ax, t/3600, spvModel.MPPTData[it0:itend], label=L"P_{\textrm{PV}}", step=:post, color=colors[2], linewidth=1.5)
        stairs!(ax, t/3600, Pg, label=L"P_{\textrm{g}}", step=:post, color=colors[3], linewidth=1.5)
        if haskey(results, "Phpe")
                stairs!(ax, t/3600, Phpe, label=L"P_{\textrm{HP}}^{e}", step=:post, color=colors[4], linewidth=1.5)
        end
        stairs!(ax, t/3600, Pbess, label=L"P_{\textrm{BESS}}", step=:post, color=colors[5], linewidth=1.5)
        # if haskey(results, "Pev[1]")
        if length(data["EV"]) != 1
                [stairs!(ax, t/3600, Pev[n] .* γ_cont[n],color=colors[5+n], linewidth=1, label=L"P_{\textrm{EV}, %$n}") for n ∈ 1:nEV]
        else
                [stairs!(ax, t/3600, Pev .* γ_cont,color=colors[5+n], linewidth=1, label=L"P_{\textrm{EV}, %$n}") for n ∈ 1:nEV]
        end
        axislegend(ax; position=:rt); limits!(ax, t[1]/3600, t[end]/3600, nothing, nothing);
        # ylabel!(ax, L"P\ [kW]");
        # title!(L"\textbf{Power\ Balance}")
        # return display(GLMakie.Screen(), fig)
        return fig
end     

function makeEMSplots(results::Dict, data::Dict; backend::String="CairoMakie", filename::String="EMSplot.pdf")
        @assert backend ∈ ["CairoMakie", "GLMakie"] "Invalid backend. Choose CairoMakie or GLMakie."
        if backend == "CairoMakie"
                CairoMakie.activate!(type="svg")
        else
                GLMakie.activate!()
        end

        set_theme!(theme_latexfonts())
        fig = Figure(size=(1800, 800), fontsize=17)
        colors=ColorSchemes.Dark2_7.colors;
        
        t = results[:"t"]
        t0 = t[1]
        Δt = t[2] - t[1]
        it0 = round(Int, (t0 / Δt))
        itend = it0 + length(t) - 1
        # Inputs
        spvModel = data[:"SPV"];
        gridModel = data[:"grid"];
        # Results
        Pg = convert.(Float64, results[:"Pg"])
        haskey(results, "Phpe") ? Phpe = convert.(Float64, results[:"Phpe"]) : nothing
        Pbess = convert.(Float64, results[:"Pbess"])
        SoCbess = convert.(Float64, results[:"SoCbess"])
        Qbess = convert.(Float64,  results[:"Qbess"])
        R0bess = convert.(Float64,  results[:"R0bess"])
        nEV = length(data["EV"])
        if length(data["EV"]) !=1
                γ_cont = results[:"γ_cont"]
                Pev = [results[:"Pev[1]"], results[:"Pev[2]"]]
                SoCev = [results[:"SoCev[1]"], results[:"SoCev[2]"]]
                Qev = [convert.(Float64, results[:"Qev[1]"]), convert.(Float64, results[:"Qev[2]"])]
                R0ev = [convert.(Float64, results[:"R0ev[1]"]), convert.(Float64, results[:"R0ev[2]"])]
        else
                γ_cont = results[:"γ_cont"][:]
                Pev = convert.(Float64, results[:"Pev[1]"])
                SoCev = convert.(Float64, results[:"SoCev[1]"])
                Qev = convert.(Float64, results[:"Qev[1]"])
                R0ev = convert.(Float64, results[:"R0ev[1]"])
        end
        Ptess = convert.(Float64, results[:"Ptess"])
        SoCtess = convert.(Float64, results[:"SoCtess"])
        
        # Energy Balance plots
        eb_ax = Axis(fig[1, 1]; ylabel=L"$P$ [kW]", title=L"\text{Power Balance}",
                );
        stairs!(eb_ax, t/3600, gridModel.loadE[it0:itend], color=colors[1],linewidth=1.5, label=L"P_{\textrm{load}}^{e}", step=:post)
        stairs!(eb_ax, t/3600, spvModel.MPPTData[it0:itend], color=colors[2], linewidth=1.5, label=L"P_{\textrm{PV}}", step=:post)
        stairs!(eb_ax, t/3600, Pg, color=colors[3], linewidth=1.5, label=L"P_{\textrm{g}}", step=:post)
        if haskey(results, "Phpe")
                stairs!(eb_ax, t/3600, Phpe, color=colors[4], linewidth=1.5, label=L"P_{\textrm{HP}}^{e}", step=:post)
        end
        stairs!(eb_ax, t/3600, Pbess, color=colors[5], linewidth=1.5, label=L"P_{\textrm{BESS}}", step=:post)
        if length(data["EV"]) != 1
                [stairs!(eb_ax, t/3600, Pev[n] .* γ_cont[n], color=colors[5+n], linewidth=1.5, label=L"P_{\textrm{EV}, %$n}", step=:post) for n ∈ 1:nEV]
        else
                stairs!(eb_ax, t/3600, Pev .* γ_cont, color=colors[6], linewidth=1.5, label=L"P_{\textrm{EV}}", step=:post)
        end
        axislegend(eb_ax, position=:lt); limits!(eb_ax, t[1]/3600, t[end]/3600, nothing, nothing);
        
        # Thermal Balance plot
        tb_ax = Axis(fig[2,1]; ylabel=L"$P$ [kWt]", xlabel=L"$t$ [hr]", title=L"\text{Thermal Balance}",
                );
        stairs!(tb_ax, t/3600, gridModel.loadTh[it0:itend], color=colors[1], linewidth=1, label=L"P_{\textrm{load}}^{th}", step=:post)
        stairs!(tb_ax, t/3600, spvModel.MPPTData[it0:itend]*EMSData["ST"].η, color=colors[2], linewidth=1, label=L"P_{\textrm{ST}}", step=:post)
        stairs!(tb_ax, t/3600, Ptess, color=colors[3], linewidth=1, label=L"P_{\textrm{TESS}}", step=:post)
        if haskey(results, "Phpe")
                stairs!(tb_ax, t/3600, Phpe*EMSData["HP"].η, color=colors[4], linewidth=1, label=L"P_{\textrm{HP}}^{t}", step=:post)
        end
        axislegend(tb_ax); limits!(tb_ax, t[1]/3600, t[end]/3600, nothing, nothing);
        
        # Electric Storage plot
        eess_ax = Axis(fig[1, 2]; ylabel=L"[%]", title=L"\text{Electric Storage}",
                ytickformat = values -> [L"%$(value)" for value in values],
                xtickformat = values -> [L"%$(value)" for value in values]);
        stairs!(eess_ax, t/3600, SoCbess*100, color=colors[1], linewidth=1, label=L"SoC_{\textrm{BESS}}", step=:post)
        if length(data["EV"]) != 1
                [stairs!(eess_ax, t/3600, SoCev[n] .* γ_cont[n] * 100, 
                color=colors[1+n], linewidth=1, label=L"SoC_{\textrm{EV}, %$n}", step=:post) for n ∈ 1:nEV]
        else
                stairs!(eess_ax, t/3600, SoCev .* γ_cont .* 100, 
                color=colors[1+1], linewidth=1, label=L"SoC_{\textrm{EV}}", step=:post)
                # add a vspan for the EV availability
                # find the indeces were γ_cont changes from 0 to 1 and viceversa
                # arrival when γ_cont changes from 0 to 1
                indArr = findall(x -> x == 1, diff(γ_cont))
                # departure when γ_cont changes from 1 to 0
                indDep = findall(x -> x == -1, diff(γ_cont))
                # check if they are the same length.
                length(indDep) > length(indArr) ? append!(indArr, length(t)) : nothing
                vspan!(eess_ax, t[indDep]/3600,t[indArr]/3600, ymax=100, color = (:grey, 0.2))
        end
        axislegend(eess_ax); limits!(eess_ax, t[1]/3600, t[end]/3600, 0, 100);
        
        # Thermal Storage plot
        tess_ax = Axis(fig[2, 2]; ylabel=L"[%]", xlabel=L"$t$ [hr]", title=L"\text{Thermal Storage}",
                ytickformat = values -> [L"%$(value)" for value in values],
                xtickformat = values -> [L"%$(value)" for value in values]);
        stairs!(tess_ax, t/3600, SoCtess*100, color=colors[1], linewidth=1, label=L"SoC_{\textrm{TESS}}", step=:post)
        axislegend(tess_ax, position=:lt); limits!(tess_ax, t[1]/3600, t[end]/3600, nothing, nothing);
        
        # Ageing - ΔQₛₐ vs t, % of initial capacity
        qloss_ax = Axis(fig[1, 3]; ylabel=L"[%]", title=L"Ageing - $ΔQ_{sa}$ vs t",
        );
        stairs!(qloss_ax, t/3600, (Qbess .- Qbess[1]) ./ Qbess[1] * 100, color=colors[1], linewidth=1, label=L"Q_{loss , \textrm{BESS}}", step=:post)
        if length(data["EV"]) != 1
                [stairs!(qloss_ax, t/3600, (Qev[n] .- Qev[n][1]) ./ Qev[n][1] * 100,
                color=colors[1+n], linewidth=1, label=L"Q_{loss, \textrm{EV}, %$n}", step=:post) for n ∈ 1:nEV]
        else
                stairs!(qloss_ax, t/3600, (Qev .- Qev[1]) ./ Qev[1] * 100,
                color=colors[1+1], linewidth=1, label=L"Q_{loss, \textrm{EV}}", step=:post)
        end
        axislegend(qloss_ax); limits!(qloss_ax, t[1]/3600, t[end]/3600, nothing, nothing);
        
        # Ageing - ΔR₀,ₛₐ vs t
        rsa_ax = Axis(fig[2, 3]; ylabel=L"[%]", xlabel=L"$t$ [hr]", title=L"Ageing - $\Delta R_{0,sa}$ vs t");
        stairs!(rsa_ax, t/3600, (R0bess .- R0bess[1]) ./ R0bess[1] * 100, color=colors[1], linewidth=1, label=L"R_{0,\textrm{BESS}}", step=:post)
        if length(data["EV"]) != 1
                [stairs!(rsa_ax, t/3600, (Rev[n] .- R0ev[n][1]) ./ R0ev[n][1] * 100, color=colors[1+n], linewidth=1, label=L"R_{0,\textrm{EV}, %$n}", step=:post) for n ∈ 1:nEV]
        else
                stairs!(rsa_ax, t/3600, (R0ev .- R0ev[1]) ./ R0ev[1] * 100, color=colors[1+1], linewidth=1, label=L"R_{0,\textrm{EV}}", step=:post)
        end
        axislegend(rsa_ax, position=:rb); 
        limits!(rsa_ax, t[1]/3600, t[end]/3600, nothing, nothing);
        if backend == "CairoMakie"
                save(filename, fig)
                return fig
        else
                return display(GLMakie.Screen(), fig)
        end
end

function compareEB(results::Vector{Dict}, data::Dict)
        # This function overlaps the electric power balance of various EMS runs.
        # It is called by the wrapper script, wrapperEMS.ipynb.
        # The function takes the optimized model and the settings as input.
        GLMakie.activate!()
        set_theme!(theme_latexfonts())

        fig = Figure(size=(650,450))
        ax = Axis(fig[1, 1];
                ylabel=L"$P$ [kW]",
                xlabel=L"$t$ [hr]",
                );
        colors=ColorSchemes.Dark2_7.colors;
        # Inputs
        spvModel = data[:"SPV"];
        gridModel = data[:"grid"];

        for (idx, result) in enumerate(results)
                t=result[:"t"];
                t0=t[1]; Δt=t[2]-t[1];
                # make new time window
                # it0 = round(Int,(t0/Δt) + 1);
                it0 = round(Int,(t0/Δt));
                itend = it0+length(t)-1;
                Pg=result[:"Pg"];
                Phpe=result[:"Phpe"];
                # BESSresults
                Pbess=result[:"Pbess"];
                # EV
                nEV = length(data["EV"])
                if length(data["EV"]) !=1
                        γ_cont = result[:"γ_cont"]
                        Pev = convert.(Vector{Float64}, [result[:"Pev[1]"], result[:"Pev[2]"]]);
                else
                        γ_cont = result[:"γ_cont"][:]
                        Pev = convert.(Float64, result[:"Pev[1]"])
                end

                # Electric Power balance
                stairs!(ax,t/3600,gridModel.loadE[it0:itend], label=L"P_{\textrm{load}}^{e}",
                        step=:post, linestyle=:dash, linewidth=0.5, color=colors[1])
                stairs!(ax,t[1]/3600,gridModel.loadE[it0], label=L"P_{\textrm{load}}^{e}",
                        step=:post, linewidth=2, color=colors[1])
                stairs!(ax,t/3600,spvModel.MPPTData[it0:itend], label=L"P_{\textrm{PV}}",
                        step=:post, linestyle=:dash, linewidth=0.5, color=colors[2])
                stairs!(ax,t[1]/3600,spvModel.MPPTData[it0], label=L"P_{\textrm{PV}}",
                        step=:post, linewidth=2, color=colors[2])
                stairs!(ax,t/3600, convert.(Float64, Pg), label=L"P_{\textrm{g}}",
                        step=:post, linestyle=:dash, linewidth=0.5, color=colors[3])
                stairs!(ax,t[1]/3600,convert.(Float64, Pg[1]), step=:post, label=L"P_{\textrm{g}}",
                        linewidth=2, color=colors[3])
                stairs!(ax,t/3600,convert.(Float64, Phpe), step=:post, linestyle=:dash, label=L"P_{\textrm{HP}}^{e}",
                        linewidth=0.5, color=colors[4])
                stairs!(ax,t[1]/3600,convert.(Float64, Phpe[1]), step=:post, label=L"P_{\textrm{HP}}^{e}",
                        linewidth=2, color=colors[4])
                stairs!(ax,t/3600,convert.(Float64, Pbess), step=:post, linestyle=:dash, label=L"P_{\textrm{BESS}}",
                        linewidth=0.5, color=colors[5])
                stairs!(ax,t[1]/3600,convert.(Float64, Pbess[1]), step=:post, label=L"P_{\textrm{BESS}}",
                        linewidth=2, color=colors[5])
                if size(data["EV"],2) != 1
                        [stairs!(ax, t/3600, Pev[n].*γ_cont[n], label=L"P_{\textrm{EV}, %$n}", 
                                step=:post, linestyle=:dash, linewidth=0.5, color=colors[5+n]) for n ∈ 1:nEV]
                        [stairs!(ax, t[1]/3600, Pev[n][1].*γ_cont[n][1], label=L"P_{\textrm{EV}, %$n}",
                                step=:post, linewidth=2, color=colors[5+n]) for n ∈ 1:nEV]
                else
                        [stairs!(ax, t/3600, Pev.*γ_cont, label=L"P_{\textrm{EV}, %$n}", 
                                step=:post, linestyle=:dash, linewidth=0.5, color=colors[5+n]) for n ∈ 1:nEV]
                        [stairs!(ax, t[1]/3600, Pev[1].*γ_cont[1], label=L"P_{\textrm{EV}, %$n}",
                                step=:post, linewidth=2, color=colors[5+n]) for n ∈ 1:nEV]
                end

        end
        axislegend(merge=true)
        return fig
end

function compareEMSplots(results::Dict)
        # This function creates the plots for several EMS runs. It is called by the wrapper script, wrapperEMS.ipynb.
        # The plots are:
        #  - the carrier balances (electric and thermal),
        #  - the storage state of charge,
        #  - and the ageing of the storage.
        # The function takes the optimized model and the settings as input.
        CairoMakie.activate!(type="svg")
        set_theme!(theme_latexfonts())
        fig=Figure(size=(800,600),fontsize=17)
        ax=Axis(fig[1,1]; xlabel=L"$t$ [h]", ylabel=L"$SoC$ [%]",limits=(0,24,0,100),
                );
        colors=ColorSchemes.seaborn_deep.colors;
        # Unwrap the model
        CEmpDeg=results[:CEmpDeg];
        BNoDeg=results[:BNoDeg];
        
        # SoC comparison
        stairs!(ax,CEmpDeg["t"]/3600, CEmpDeg["SoCbess"].*100, label=L"\textrm{BESS\ - CEmpDeg}", 
            step=:post, color=colors[1], linewidth=4)
        stairs!(ax, BNoDeg["t"]/3600, BNoDeg["SoCbess"].*100, label=L"\textrm{BNoDeg}", 
            step=:post, color=colors[1], linestyle=:dash)
        stairs!(ax, CEmpDeg["t"]/3600, CEmpDeg["SoCev[1]"].*100 .*CEmpDeg["γ_cont"][1], label=L"\textrm{EV 1 - CEmpDeg}",
             step=:post, color=colors[2], linewidth=4)
        stairs!(ax, BNoDeg["t"]/3600, BNoDeg["SoCev[1]"].*100 .*BNoDeg["γ_cont"][1], label=L"\textrm{BNoDeg}",
             step=:post, color=colors[2], linestyle=:dash)
        stairs!(ax, CEmpDeg["t"]/3600, CEmpDeg["SoCev[2]"].*100 .*CEmpDeg["γ_cont"][2], label=L"\textrm{EV 2 - CEmpDeg}",
             step=:post, color=colors[3], linewidth=4)
        stairs!(ax, BNoDeg["t"]/3600, BNoDeg["SoCev[2]"].*100 .*BNoDeg["γ_cont"][2], label=L"\textrm{BNoDeg}",
             step=:post, color=colors[3], linestyle=:dash)
        axislegend(ax, position=:rt);
        # ylims!(ax, 0, 100);
        
        # # Add arrow annotation
        # quiver!([11], [15], quiver=([-1], [7]), color=:black)
        # annotate!([11.5], [12],L"\textrm{No\ deg.}")

        return display(GLMakie.Screen(), fig)
end

function compareTESS(results::Dict)
        # This function creates the plots for several EMS runs. It is called by the wrapper script, wrapperEMS.ipynb.
        # The plots are:
        #  - the carrier balances (electric and thermal),
        #  - the storage state of charge,
        #  - and the ageing of the storage.
        # The function takes the optimized model and the settings as input.
        CairoMakie.activate!()
        set_theme!(theme_latexfonts())
        fig=Figure(size=(800,400),fontsize=17)
        ax=Axis(fig[1,1]; xlabel=L"$t$ [h]", ylabel=L"$SoC$ [%]",limits=(0,24,75,95),
                );
        colors=ColorSchemes.tab10.colors;

        # Unwrap the model
        CEmpDeg=results[:CEmpDeg];
        BNoDeg=results[:BNoDeg];

        # SoC comparison
        stairs!(ax,CEmpDeg["t"]/3600, CEmpDeg["SoCtess"].*100, label=L"\textrm{TESS\ - CEmpDeg}", 
                step=:post, color=colors[1], linewidth=4)
        stairs!(ax,BNoDeg["t"]/3600, BNoDeg["SoCtess"].*100, label=L"\textrm{BNoDeg}",
                step=:post, color=colors[1], linestyle=:dash)
        axislegend(ax,position=:lt)
        # # Add arrow annotation
        # quiver!([11.5], [90], quiver=([1.5], [-4]), color=:black)
        # annotate!([11.5], [91],L"\textrm{No\ deg.}")

        return fig
end

function hist2d(x, y; nbins=10, normalize=false)
        # This function creates a 2D histogram of the input data.
        xedges = LinRange(minimum(x), maximum(x), nbins + 1)
        yedges = LinRange(minimum(y), maximum(y), nbins + 1)

        hist = zeros(Int, nbins, nbins)

        for i in 1:nbins, j in 1:nbins
                xmask = (x .>= xedges[i]) .& (x .< xedges[i + 1])
                ymask = (y .>= yedges[j]) .& (y .< yedges[j + 1])
                hist[i, j] = count(xmask .& ymask)
        end
        
        if normalize
                hist = hist ./ sum(hist)
        end

        return hist, xedges, yedges
end

function countourSP(results, keys, titles)
# This function creates the contour plots (SoC vs P) of the different ESS.
        i=0
        cmaps = [:algae, :reds, :blues, :balance]
        # Create a heatmap
        set_theme!(theme_latexfonts())
        fig = Figure(size = (1000, 350))
        axs=[];

        for kk in keys
            i=i+1
            Pkey = "P$kk"; SoCkey = "SoC$kk";
            # Create a 2D histogram
            histogram, xedges, yedges = hist2d(results[Pkey], results[SoCkey], nbins=60)
            ax = Axis(fig[1,i], xlabel=L"P_{sa}", ylabel=L"SoC_{sa}",
                    limits=(-10.0,10.0,0.1,0.9),
                    title = titles[i];
                    xticks = -10:2:10,
                    yticks = 0.:0.1:0.9,
                    )
            push!(axs,ax)
            # Makie.heatmap!(ax,xedges, yedges, histogram', colormap=:viridis, xlabel="Pa", ylabel="Pb")
            co = Makie.contourf!(ax, xedges[2:end], yedges[2:end], histogram,
                        levels=0.01:0.1:1, colormap=cmaps[i], mode = :relative, label = "$kk")
            Colorbar(fig[2,i], co, vertical = false)
        end
        hideydecorations!.(axs[2:end], grid = false)

        # return display(GLMakie.Screen(), fig)
        return fig
end

function compareEMSaging(results::Dict)
# This function creates the plots for several EMS runs. It is called by the wrapper script, wrapperEMS.ipynb.
        CairoMakie.activate!(type = "svg")
        set_theme!(theme_latexfonts())
        fig=Figure(size=(800,400),fontsize=11)        

        # extract the seasons from the results
            # loop over the season
        for (j,season) ∈ enumerate([:summer, :winter])
                # loop over the DAtype
                for (i, DAtype) ∈ enumerate(["BNoDeg_W10", "CEmpDeg_W11e4", "CPBDeg_W11e4"])
                        concRes = concatResultsRH(results_DA[season][Symbol(DAtype)]; typeOpt="day-ahead")
                        t = convert.(Float64, concRes["t"])
                        Psa = convert.(Float64, concRes[P_key])
                        
                        # stairs!(axes[j], t./3600, Psa, color=colors[i], label=L"%$DAtype", linewidth=2, step=:post)
                        # limits!(axes[j], t[1]/3600, t[end]/3600, nothing, nothing);
                end
        end

        return
end

function plotWeightedCost(costSumm::Dict{Symbol, Costs}, s::modelSettings)
        gr(dpi=500)
        t=s.dTime/3600;
        CEmpDeg=costSumm[:CEmpDeg];
        BNoDeg=costSumm[:BNoDeg];
        # ECM+Emp Degradation
        # weighted total cost
        p=plot(t,CEmpDeg.wCloss+CEmpDeg.wCgrid, label=L"\textrm{CEmpDeg} - W.C_{\textrm{total}}", color=1, width=4)
        # weighted degradation cost
        plot!(t,CEmpDeg.wCloss+CEmpDeg.wCgrid, label=L"W.C_{\textrm{loss}}", color=1, 
            fill=(CEmpDeg.wCgrid,0.5), fillstyle = :/) 
        # weighted grid cost
        plot!(t,CEmpDeg.wCgrid, label=L"W.C_{\textrm{grid}}", color=1, width=2,
            fill=(0,0.5), fillstyle = :\)
        # Bucket+ No Degradation
        # weighted total cost
        plot!(t,BNoDeg.wCloss+BNoDeg.wCgrid, label=L"\textrm{BNoDeg} - W.C_{\textrm{total}}", color=2, width=4)
        # weighted degradation cost
        plot!(t,BNoDeg.wCloss+BNoDeg.wCgrid, label=L"W.C_{\textrm{loss}}", color=2, 
            fill=(BNoDeg.wCgrid,0.5), fillstyle = :/)
        # weighted grid cost
        plot!(t,BNoDeg.wCgrid, label=L"W.C_{\textrm{grid}}", color=2, width=2,
            fill=(0,0.5), fillstyle = :\,
            xtickfontsize=13, ytickfontsize=13, legendfontsize=13, yguidefontsize=13, xguidefontsize=13,
            legend=:topleft, xformatter=:latex, yformatter=:latex)
        xlabel!(L"t\ [h]"); ylabel!(L"\textrm{Weighted\ Cost [€]}")
        return p
end

function compareCgrid(results::Dict,weights::Array)
        # plot winter
        t=0:0.25:24
        winter=results[:winter]
        summer=results[:summer]
        Wgrid=weights; # weight of grid prices
        gr(dpi=500)

        ax=plot(t, winter[:CEmpDeg].wCgrid ./ Wgrid[1], label=L"\textrm{CEmpDeg - Winter}", linewidth=2, color=1)
        plot!(t, winter[:BNoDeg].wCgrid ./ Wgrid[2], label=L"\textrm{BNoDeg}", linewidth=1.5, color=1, linestyle=:dash)
        plot!(t, summer[:CEmpDeg].wCgrid ./ Wgrid[3], label=L"\textrm{CEmpDeg - Summer}", linewidth=2, color=2)
        plot!(t, summer[:BNoDeg].wCgrid ./ Wgrid[4], label=L"\textrm{BNoDeg}", linewidth=1.5, color=2, linestyle=:dash,
             xformatter=:latex, yformatter=:latex, legend=:topleft, size(500,300),
             xtickfontsize=10, ytickfontsize=10, legendfontsize=9, yguidefontsize=10, xguidefontsize=10,)
        xlabel!(L"t\ [hr]"); ylabel!(L"\sum C_{\textrm{grid}}.\Delta t\ [€]")
        # Add arrow annotation
        quiver!([23], [1.17e-2], quiver=([0], [-0.55e-2]), color=:black)
        quiver!([24], [7.5e-3], quiver=([0], [-0.1e-2]), color=:black)
        annotate!([18], [0.8e-2],L"\textrm{Cost\ reduction}")
        return ax
end

function compareDiffCgrid(results::Dict)
        # plot winter
        t=0:0.25:23.75
        winter=results[:winter]
        summer=results[:summer]

        gr(dpi=500)

        ax=plot(t, diff(winter[:"CEmpDeg"][:"wCgrid"]).*3600/0.25, label=L"\textrm{CEmpDeg - Winter}", linewidth=2, color=1)
        plot!(t, diff(winter[:"BNoDeg"][:"wCgrid"]).*3600/0.25, label=L"\textrm{BNoDeg}", linewidth=1.5, color=1, linestyle=:dash)
        plot!(t, diff(summer[:"CEmpDeg"][:"wCgrid"]).*3600/0.25, label=L"\textrm{CEmpDeg - Summer}", linewidth=2, color=2)
        plot!(t, diff(summer[:"BNoDeg"][:"wCgrid"]).*3600/0.25, label=L"\textrm{BNoDeg}", linewidth=1.5, color=2, linestyle=:dash,
             xformatter=:latex, yformatter=:latex, legend=:topright, size(500,300),
             xtickfontsize=10, ytickfontsize=10, legendfontsize=9, yguidefontsize=10, xguidefontsize=10,)
        xlabel!(L"t\ [hr]"); ylabel!(L"W_{\textrm{grid}}.C_{\textrm{grid}}\ [€]")
        # # Add arrow annotation
        # quiver!([24], [0.255], quiver=([0], [-0.125]), color=:black)
        # quiver!([24], [0.095], quiver=([0], [-0.055]), color=:black)
        # annotate!([18], [0.12],L"\textrm{Cost\ reduction.}")
        return ax
end

function optStatusPlots(results::Vector{Dict})
        CairoMakie.activate!();
        set_theme!(theme_latexfonts())
        statusRes=[results[ts][:"status"] for ts in 1:length(results)];
        # what are the different statuses?
        statusTypes = unique(statusRes);
        # count many times each status appears
        barData = [count(x->x==st, statusRes) for st in statusTypes];
        println("Status types: ", statusTypes)
        println("Status counts: ", barData)
        # now we can plot the status of each simulation
        # barData=[sum(statusRes .== "LOCALLY_SOLVED"), sum(statusRes .== "LOCALLY_INFEASIBLE")]
        # # assign 0 to LOCALLY_SOLVED and 1 to LOCALLY_INFEASIBLE
        # statusRes=[statusRes[ts] == "LOCALLY_SOLVED" ? 0 : 1 for ts in 1:length(statusRes)]
        # assign integers from 0 to n to each status
        statusRes=[findfirst(x->x==st, statusTypes) for st in statusRes];
        # Define axes and figure
        fig=Figure(size=(1000,250))
        colors=ColorSchemes.tab10.colors;
        ga = fig[1, 1] = GridLayout()
        gb = fig[1, 2] = GridLayout()

        ax1=Axis(ga[1,1], width=200); 
        ax2=Axis(gb[1,1], ylabel=L"\text{Status}",
            yticks=1:length(statusTypes));
        ax3=Axis(gb[2,1], xlabel=L"$t$ [h]", ylabel="n");
        # piechart
        Makie.pie!(ax1, barData, color=colors[1:length(statusTypes)])
        hidedecorations!(ax1)
        hidespines!(ax1)
        # Legend(ga[1,2], [PolyElement(color=c) for c in colors],
        #         [L"\text{LOCALLY SOLVED}", L"\text{LOCALLY INFEASIBLE}"],framevisible=false)
        Legend(ga[1,2], [PolyElement(color=c) for c in colors[1:length(statusTypes)]],
                String.(statusTypes),framevisible=false)
        # status scatter plot
        Makie.scatter!(ax2,0.25:0.25:length(statusRes)/4, statusRes, label=L"\text{status}", strokewidth=0.5, strokecolor=:black);
        # consecutive INFEASIBLE runs
        xs = 0.25:0.25:length(results)/4;
        ys=[length(results[ts][:"t"]) for ts in 1:length(results)];
        barplot!(ax3, xs, -maximum(ys) .+ ys, color = colors[1], strokecolor = :black, strokewidth = 1)
        # fig
        return fig
end