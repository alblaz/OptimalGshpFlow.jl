###############################################
# extract heat exchanger temperature profiles #
###############################################

ϵ_counterflow(x,NTU,Cr) = (1 - exp(-NTU*(x)*(1-Cr)))/(1-Cr*exp(-NTU*(x)*(1-Cr)))

# the fluid with lower capacity Cc is flowing according to x (from 0 to 1), while the hot fluid with higher capacity is flowing on the contrary direction (from 1 to 0). This is important to keep in mind if we want to plot the inverse representation where the fluid with lower capacity enter the heat exchanger at x=1. In this case we first would perform the calculation above and then we would remap the the values of the

# function counterflow_heatexchanger(Ch,Thi,Tho,Cc,Tci,Tco,UA; x =0:0.1:1, defaultdirection = true) 
function counterflow_heatexchanger(heatexchanger_solution::HeatExchangerSolution; n=11 , defaultdirection = true) 

    x = range(0,1,n)
    @unpack Cc, Ch, NTU, Tch, Tcl, Thh, Thl, UA = heatexchanger_solution

    Thi, Tho = Thh, Thl
    Tci, Tco = Tcl, Tch

    Cmin = min(Cc,Ch)
    Cmax = max(Cc,Ch)
    # NTU = UA/Cmin
    Cr = Cmin/Cmax

    if Cc > Ch
        Th = @.  Thi - ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)
        Tc = @.  Tco - ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)*Ch/Cc
    else
        ϵx = ϵ_counterflow.(x,NTU,Cr)
        Tc = @.  Tci +  ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)
        Th = @.  Tho +  ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)*Cc/Ch
    end

    if defaultdirection == false
        reverse!(Tc) 
        reverse!(Th)
    end

    return x,Tc,Th
end