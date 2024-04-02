### Make Forecasts
# This script makes forecasts for timeseries
# It uses the StateSpaceModels.jl package to fit a model to the data without feature engineering.
# It then uses the fitted model to make forecasts for a rolling horizon optimizer.

function makeForecast(y, modelType::String, horizon::Int64)
    # create a model
    if modelType == "LocalLevel"
        model = LocalLevel(y)
        StateSpaceModels.fit!(model)
        print_results(model) # print the results
    elseif modelType == "auto_SARIMA" # automatic fit a SARIMA model
        model = auto_arima(y, seasonal=24)
        print_results(model) # print the results
    elseif modelType == "ARIMA" # fit a ARIMA model
        model = SARIMA(y, order=(0,0,1), seasonal_order=(0,0,0,0))
        StateSpaceModels.fit!(model)
        print_results(model) # print the results
    elseif modelType == "auto_ETS" # fit a SARIMA model
        model = auto_ets(y, seasonal=23)
        print_results(model) # print the results
    end
    # make a forecast object
    forec = forecast(model, horizon) 
    # extract the expected value
    expected_value = forecast_expected_value(forec)

    return forec, model, expected_value
end

function updateForecast(model::InfiniteModel, data::Dict, sets::modelSettings)
    Dt = sets.dTime;
    t = model[:t]; 
    t0=supports(t)[1]; Δt=supports(t)[2]-supports(t)[1];
    # make new time window
    it0 = round(Int,t0/Δt);
    itend = it0+length(Dt)-1;

    # Exogenous information for the model
    # Load data
    loadElec = data["grid"].loadE[it0:itend];
    σLoadElec = std(loadElec); # standard deviation of the load
    εLoadElec = randn(length(Dt))*σLoadElec; # white noise
    loadElec = loadElec + εLoadElec; # add noise
    @parameter_function(model, Ple == (t) -> loadElec[zero_order_time_index(Dt, t)])
    
    # MPPT measurement
    # NOTE: ONLY ACTIVE WITH NO FCR
    MPPTmeas = data["SPV"].MPPTData[it0:itend];
    
    @parameter_function(model, PpvMPPT == (t) -> MPPTmeas[zero_order_time_index(Dt, t)])

    # Heat carrier
    # Thermal load
    loadTh = data["grid"].loadTh[it0:itend]; 
    @parameter_function(model, Plt == (t) -> loadTh[zero_order_time_index(Dt, t)])
    # Extract params
    ηST = data["ST"].η; # conversion factor from Electric PV to thermal
    @parameter_function(model, Pst == (t) -> ηST*MPPTmeas[zero_order_time_index(Dt, t)])
    
    # Prices
    # Note: We need to use hcat() to form the axis cause otherwise λ[:,c] broadcasts into a vector.
    # In order to stay consistent with the other DenseAxisArray we need matrices, so here hcat does that for us.
    priceBuy=hcat(data["grid"].λ[it0:itend, 1].*1e-3/3600) # convert from €/MWh to €/kWs
    @parameter_function(model, λbuy == (t) -> priceBuy[zero_order_time_index(Dt, t)])
    
    priceSell=hcat(data["grid"].λ[it0:itend, 2].*1e-3/3600) # convert from €/MWh to €/kWs
    @parameter_function(model, λsell == (t) -> priceSell[zero_order_time_index(Dt, t)])
end