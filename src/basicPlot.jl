"""
    plotAll(model::StandardStoModel, parameters::Vector{Float64},kini::Vector{Float64},delta::Float64, maxrna::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)

Export the plotly traces of all the model output
"""
function plotAll(model::StandardStoModel, parameters::Vector{Float64},kini::Vector{Float64},delta::Float64, maxrna::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)
    (mnascentmrna_model, pburst_model, survivalspot_model,survivaldark_model, survivalnextburst_model, corr_interburst, intensity_model) = ModelOutput(model::StandardStoModel, parameters::Vector{Float64},kini::Vector{Float64},delta::Float64, maxrna::Int64,tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)
    #traces
    offtimesurvival = scatter(x=1:tmaxoff, y=survivaldark_model, mode="lines", line=attr(color="black", width=3),showlegend=false);
    ontimesurvival = scatter(x=1:tmaxon, y=survivalspot_model, mode="lines", line=attr(color="black", width=3),showlegend=false);
    nextburstsurvival = scatter(x=1:tmaxnextburst, y=survivalnextburst_model, mode="lines",  line=attr(color="black", width=3),showlegend=false);
    pburst = bar(x="model", y=pburst_model)
    mnasentrna = bar(x="model", y=mnascentmrna_model)
    interburstcorr = bar(x="model", y=corr_interburst)
    avgintensity = scatter(x=1:tmaxintensity, y=intensity_model, mode="lines", line=attr(color="black", width=3),showlegend=false);

    #layouts
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
        plot_bgcolor  = figurebgcolor, yaxis_type="log") #,xaxis_type="log") yaxis_type="log"=#
    
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
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"=#
    
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
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"=#
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
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"=#

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
        plot_bgcolor  = figurebgcolor) #,xaxis_type="log") yaxis_type="log"=#

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




