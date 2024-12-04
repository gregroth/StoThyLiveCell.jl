"""
    plotAll(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)

Export the plotly traces of all the model output
"""
function plotAll(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)
    (mnascentmrna_model, pburst_model, survivalspot_model,survivaldark_model, survivalnextburst_model, corr_interburst, intensity_model) = ModelOutput(model, parameters, maxrna,detectionlimitLC, detectionlimitNS,tmaxon,tmaxoff,tmaxnextburst,tmaxintensity)
    #traces
    offtimesurvival = plot(1:tmaxoff, survivaldark_model, linetype=:steppre, lc="black", lw=3, xguide="time (min)", yguide="OFF time survival probability", yscale=:log10);
    ontimesurvival = plot(1:tmaxon, survivalspot_model, linetype=:steppre, lc="black", lw=3, xguide="time (min)", yguide="ON time survival probability", yscale=:log10);
    nextburstsurvival = plot(1:tmaxnextburst, survivalnextburst_model, linetype=:steppre, lc="black", lw=3, xguide="time (min)", yguide="Next burst time survival probability", yscale=:log10); 
    pburst = plot(bar("model", pburst_model),title="Probability to detect a burst")
    mnasentrna = plot(bar("model", mnascentmrna_model),title="Mean number of nascent mRNA")
    interburstcorr = plot(bar("model", corr_interburst), title="Inter-burst correlation")
    avgintensity =  plot(1:tmaxintensity, intensity_model, linetype=:steppre, lc="black", lw=3, xguide="time (min)", yguide="Normalized intensity");
        
    return (offtimesurvival,ontimesurvival,nextburstsurvival,pburst,mnasentrna,interburstcorr,avgintensity)
end






"""
    plotAll(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)

Export the plotly traces of all the model output
"""
#= function plotAll(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)
    (mnascentmrna_model, pburst_model, survivalspot_model,survivaldark_model, survivalnextburst_model, corr_interburst, intensity_model) = ModelOutput(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64,tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)
    #traces
    offtimesurvival = scatter(x=1:tmaxoff, y=survivaldark_model, mode="lines", line=attr(color="black", width=3),showlegend=false);
    ontimesurvival = scatter(x=1:tmaxon, y=survivalspot_model, mode="lines", line=attr(color="black", width=3),showlegend=false);
    nextburstsurvival = scatter(x=1:tmaxnextburst, y=survivalnextburst_model, mode="lines",  line=attr(color="black", width=3),showlegend=false);
    pburst = bar(x="model", y=pburst_model)
    mnasentrna = bar(x="model", y=mnascentmrna_model)
    interburstcorr = bar(x="model", y=corr_interburst)
    avgintensity = scatter(x=1:tmaxintensity, y=intensity_model, mode="lines", line=attr(color="black", width=3),showlegend=false);

    #layouts
    lshape = "linear" #"hv"
    fontsizeaxis = 26
    fontsizeticks = 21
    figurehdim = 800
    figurevdim = 600
    fontfamily = "helvetica"
    graphgridcolor = "#a6a8a6"
    axislinecolor = "black"
    figurebgcolor = "white"
    layout_offsurvival = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "time (min)",
        title_font_size=fontsizeaxis,
        # range = [1,11],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "OFF time survival probability",
        title_font_size=fontsizeaxis,
        # range = [1,11],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor,yaxis_type="log")
        
    layout_onsurvival = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "time (min)",
        title_font_size=fontsizeaxis,
        # range = [1,11],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "ON time survival probability",
        title_font_size=fontsizeaxis,
        # range = [1,11],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor,yaxis_type="log")#,xaxis_type="log")

    layout_nextburstsurvial = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "time (min)",
        title_font_size=fontsizeaxis,
        #range = [0,300],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "Next burst time survival probability",
        title_font_size=fontsizeaxis,
        #range = [-1.6,0],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor, yaxis_type="log") #,xaxis_type="log") yaxis_type="log"
    
    layout_pburst = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "",
        title_font_size=fontsizeaxis,
        #range = [0,238],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "Probability to detect a burst",
        title_font_size=fontsizeaxis,
        #range = [-1.6,0],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"
    
    layout_interburstcorr = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "",
        title_font_size=fontsizeaxis,
        #range = [0,238],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "Inter-burst correlation",
        title_font_size=fontsizeaxis,
        #range = [-1.6,0],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"
    layout_mnascent = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "",
        title_font_size=fontsizeaxis,
        #range = [0,238],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "Mean number of nascent mRNA",
        title_font_size=fontsizeaxis,
        #range = [-1.6,0],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"

    layout_intensity = Layout(width=figurehdim, height=figurevdim, font_family=fontfamily, xaxis=attr(
        title_text = "Time (min)",
        title_font_size=fontsizeaxis,
        #range = [0,238],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        yaxis=attr(
        title_text = "Normalized intensity",
        title_font_size=fontsizeaxis,
        #range = [-1.6,0],
        gridcolor = graphgridcolor,
        linecolor= axislinecolor,
        tickfont=attr(family=fontfamily, color=axislinecolor, size=fontsizeticks)),
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"

    #plot

    plt_off = plot(offtimesurvival, layout_offsurvival);
    #display(plt_off)
    plt_on = plot(ontimesurvival, layout_onsurvival);
    #display(plt_on)
    plt_nextburst = plot(nextburstsurvival, layout_nextburstsurvial);
    #display(plt_nextburst)
    plt_pburst = plot(pburst, layout_pburst);
    #display(plt_pburst)
    plt_mnasentrna = plot(mnasentrna, layout_mnascent);
    #display(plt_mnasentrna)
    plt_interburstcorr = plot(interburstcorr, layout_interburstcorr);
    #display(plt_interburstcorr)
    plt_avgintensity = plot(avgintensity, layout_intensity);
    #display(plt_avgintensity)
    return (offtimesurvival,ontimesurvival,nextburstsurvival,pburst,mnasentrna,interburstcorr,avgintensity)
end


 =#



