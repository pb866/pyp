### Functions associated with plotting data

"""
    setup_plot(plot_list, twinax, ylab, logscale, logremove, xlims, ylims, major_yticks, minor_yticks, ptype, cs, lc, lt, pt, alpha, axcolour)

If `twinax` is an array split `plot_list` into 2 datasets based on the indices
in `twinax`. Assure all paramters concerning y data are arrays of length 2, if
a second axis is allowed or length 1 otherwise. Return `plot_list` as array of
2 or 1 distinct lists an the adjusted parameters.
"""
function setup_plot(plot_list, figsize, twinax, ylab, logscale, logremove, xlims, ylims,
                    major_yticks, minor_yticks, ptype, cs, lc, lt, pt, alpha, axcolour)

  # Start plot
  fig, ax = plt.subplots(figsize=figsize)
  ax = [ax]

  # Set flag for second axis, asume no second axis
  xlim = []; ylim = []

  [p.alpha = alpha for p in plot_list]

  # Set up second axis, if the twinax array is defined
  plt1 = PlotData[]; plt2 = PlotData[]
  cl1 = []; dt1 = []; mt1 = []
  cl2 = []; dt2 = []; mt2 = []
  if !isempty(twinax)
    # Set flag true and define 2nd axis in PyPlot
    push!(ax, plt.twinx())
    # Check correct input of twinax
    if length(twinax) ≠ length(plot_list)
      throw(ArgumentError(string("Array `twinax` must have the same length ",
        "as there are number of `PhotData` elements.")))
    end

    # Assign data to the axes
    for i = 1:length(twinax)
      if twinax[i] == 1
        push!(plt1,plot_list[i])
        # check whether an array of definitions or a single definition exists
        # and assign either the current array entry to the respective axis
        # or the default value for line colour, type and marker type
        if !isempty(lt) && (lt[i] isa Tuple || lt[i] isa Vector || lt[i] isa String)
          push!(dt1,lt[i])
        else
          push!(dt1,lt)
        end
        if pt isa Vector  push!(mt1,pt[i])  else  push!(mt1,pt)  end
        if lc isa Vector  push!(cl1,lc[i])  else  push!(cl1,lc)  end
      else
        push!(plt2,plot_list[i])
        # check whether an array of definitions or a single definition exists
        # and assign either the current array entry to the respective axis
        # or the default value for line colour, type and marker type
        if !isempty(lt) && (lt[i] isa Tuple || lt[i] isa Vector || lt[i] isa String)
          push!(dt2,lt[i])
        else
          push!(dt2,lt)
        end
        if pt isa Vector  push!(mt2,pt[i])  else  push!(mt2,pt)  end
        if lc isa Vector  push!(cl2,lc[i])  else  push!(cl2,lc)  end
      end
    end

    # Make sure, all parameters for both axes are arrays of length 2
    if logscale isa String logscale = String[logscale, logscale]  end
    if logremove isa String logremove = String[logremove, logremove]  end
    if axcolour isa String axcolour = String[axcolour, axcolour]  end
    if major_yticks isa Int64 major_yticks = Int64[major_yticks, major_yticks]  end
    if minor_yticks isa Int64 minor_yticks = Int64[minor_yticks, minor_yticks]  end
    if xlims == nothing
      xlimit = (nothing, nothing)
    elseif xlims isa Tuple
      xlimit = xlims
    elseif xlims isa Array
      throw(ArgumentError("xlims must be a tuple or nothing!"))
    end
    if ylims == nothing
      ylimit = [[nothing, nothing],[nothing, nothing]]
    elseif ylims isa Tuple
      yl = [ylims[1],ylims[2]]
      ylimit = [yl, yl]
    elseif ylims isa Array
      ylimit = [[ylims[1][1], ylims[1][2]], [ylims[2][1], ylims[2][2]]]
    end
    # if !isa(legpos, Array) legpos = [legpos, legpos]  end
    # if legcolumns isa Int64 legcolumns = [legcolumns, legcolumns]  end
  else
    # Assign all data to axis 1, if no second axis
    plt1 = plot_list
    # check whether an array of definitions or a single definition exists
    # and assign either the current array entry to the respective axis
    # or the default value for line colour, type and marker type
    if !isempty(lt) && (lt[1] isa Tuple || lt[1] isa Vector || lt[1] isa String)
      dt1 = lt
    else
      dt1 = [lt for i in plot_list]
    end
    if pt isa Vector  mt1 = pt  else  mt1 = [pt for i in plot_list]  end
    if lc isa Vector  cl1 = lc  else  cl1 = [lc for i in plot_list]  end
    # If no 2nd axis, make sure, all parameters for both axes are arrays of length 1
    # and not single parameters (numbers, strings...)
    if logscale isa String logscale = String[logscale]  end
    if logremove isa String logremove = String[logremove]  end
    if axcolour isa String axcolour = String[axcolour]  end
    if major_yticks isa Real major_yticks = Int64[major_yticks]  end
    if minor_yticks isa Real minor_yticks = Int64[minor_yticks]  end
    if xlims == nothing     xlimit = (nothing, nothing)
    elseif xlims isa Tuple  xlimit = xlims
    end
    if ylims == nothing     ylimit = [[nothing, nothing]]
    elseif ylims isa Tuple  ylimit = [[ylims[1],ylims[2]]]
    end
    # if !isa(legpos, Array) legpos = [legpos]  end
    # if legcolumns isa Int64  legcolumns = [legcolumns]  end
  end
  if ylab isa AbstractString  ylab = String[ylab]  end
  if length(ax) > 1
    pltdata = [plt1, plt2]
    if !isa(cs, Vector) cs = [cs, cs]  end
    lc = [cl1, cl2]
    lt = [dt1, dt2]
    pt = [mt1, mt2]
  else
    pltdata = [plt1]
    cs = [cs]
    lc = [cl1]
    lt = [dt1]
    pt = [mt1]
  end

  # Return adjusted data
  return pltdata, fig, ax, ylab, logscale, logremove, xlimit, ylimit,
         major_yticks, minor_yticks, cs, lc, lt, pt, axcolour
end #function setup_plot


"""
    set_style(pltdata, ptype, cs, lc, lt, pt, ax_2)

Set colour scheme for `pltdata` holding `PlotData` to `cs` from
function `sel_ls`. If colour scheme is set to `"own"` or `"mix"` for any array
element in `pltdata`, `lc`, `lt`, and `pt` must hold valid definitions of line colours,
line and point (marker) types at this array position. The boolean `ax_2` is needed
to specify, whether a second axis is used.
"""
function set_style(pltdata, ptype, cs, lc, lt, pt)

  # Reset line and marker types to solid lines and squares, respectively,
  # if plot_type is set to line, scatter or both (markers connected by lines)
  for i = 1:length(pltdata)
    if ptype == "line"
      [pltdata[i][j].dashes = Int64[] for j = 1:length(lt[i])]
      [pltdata[i][j].marker = "None" for j = 1:length(lt[i])]
    elseif ptype == "scatter"
      [pltdata[i][j].dashes = "None" for j = 1:length(lt[i])]
      [pltdata[i][j].marker = "s" for j = 1:length(lt[i])]
    elseif ptype == "both"
      [pltdata[i][j].dashes = Int64[] for j = 1:length(lt[i])]
      [pltdata[i][j].marker = "s" for j = 1:length(pt[i])]
    end
  end

  # Set index for sel_ls function, start at 1 for separate schemes
  # and use continuous indexing for the same scheme
  idx = []
  push!(idx, [i for i = 1:length(pltdata[1])])
  if length(cs) == 2
    cs[1] == cs[2] ? ind = length(idx[1]) : ind = 0
    push!(idx, [i+ind for i = 1:length(pltdata[2])])
  end
  # Set line/marker styles and colour for a chosen colour scheme
  # (overwrites default settings from plot_type)
  for i = 1:length(pltdata)  if cs[i]≠""
    for j = 1:length(pltdata[i])
      cl, dt, mt = sel_ls(cs[i], lc=idx[i][j], lt=idx[i][j], pt=idx[i][j])
      pltdata[i][j].colour = cl
      if ptype ≠ "scatter"  pltdata[i][j].dashes = dt  end
      if ptype ≠ "line"  pltdata[i][j].marker = mt  end
    end
  end  end

  # Set manually defined styles (overwrites previous settings)
  for i = 1:length(lc), j = 1:length(lc[i])
    if lc[i][j] ≠ "default"  pltdata[i][j].colour = lc[i][j]  end
  end
  for i = 1:length(lc), j = 1:length(lc[i])
    if lt[i][j] ≠ "default"  pltdata[i][j].dashes = lt[i][j]  end
  end
  for i = 1:length(lc), j = 1:length(lc[i])
    if pt[i][j] ≠ "default"  pltdata[i][j].marker = pt[i][j]  end
  end

  # Return adjusted PlotData
  return pltdata
end #function set_style


"""
    plt_DataWithErrors(pltdata, ax, offset)

For each `PlotData` in array `pltdata` (and `ax`), retrieve the errors in the `PlotData` and
plot according to specifications. For error bars of markers, the standard cap size is `3`,
which can be adjusted by an `offset` (positive or negative number to be added to cap size).

Returns `fig` and `ax` (array of axes) for further plot modifications or printing.
"""
function plt_DataWithErrors(pltdata, ax, offset)
  # Redefine errors for error bars
  xerr = redef_err(pltdata,:x,:xlerr,:xuerr)
  yerr = redef_err(pltdata,:y,:ylerr,:yuerr)
  # Loop over graph data and plot each data according to its type and format
  # defined by the struct PlotData
  for i = 1:length(pltdata)
    if pltdata[i].yuerr ≠ nothing && pltdata[i].marker == "None"
      p = ax.plot(pltdata[i].x, pltdata[i].y, lw = pltdata[i].lw,
          dashes=pltdata[i].dashes, color=pltdata[i].colour, label=pltdata[i].label)
      pltdata[i].colour = p[1].get_color()
      ax.fill_between(pltdata[i].x, pltdata[i].ylerr, pltdata[i].yuerr,
          color=pltdata[i].colour, alpha=0.2)
    elseif pltdata[i].xuerr ≠ nothing && pltdata[i].yuerr ≠ nothing
      if (!isempty(pltdata[i].dashes) && pltdata[i].dashes[1] == 0) || pltdata[i].dashes == "None"
        ax.errorbar(pltdata[i].x, pltdata[i].y, xerr=[xerr[i].lower, xerr[i].upper],
          yerr=[yerr[i].lower, yerr[i].upper], fmt=pltdata[i].marker, color=pltdata[i].colour,
          label=pltdata[i].label, capsize=3+offset, alpha=pltdata[i].alpha)
      else
        ax.errorbar(pltdata[i].x, pltdata[i].y, xerr=[xerr[i].lower, xerr[i].upper],
          yerr=[yerr[i].lower, yerr[i].upper], lw = pltdata[i].lw,
          marker=pltdata[i].marker, dashes=pltdata[i].dashes, color=pltdata[i].colour,
          label=pltdata[i].label, capsize=3+offset, alpha=pltdata[i].alpha)
      end
    elseif pltdata[i].yuerr ≠ nothing
      if (!isempty(pltdata[i].dashes) && pltdata[i].dashes[1] == 0) || pltdata[i].dashes == "None"
        ax.errorbar(pltdata[i].x, pltdata[i].y, yerr=[yerr[i].lower, yerr[i].upper],
          fmt=pltdata[i].marker, color=pltdata[i].colour, label=pltdata[i].label, capsize=3+offset,
          alpha=pltdata[i].alpha)
      else
        ax.errorbar(pltdata[i].x, pltdata[i].y, yerr=[yerr[i].lower, yerr[i].upper],
          lw = pltdata[i].lw, marker=pltdata[i].marker, dashes=pltdata[i].dashes,
          color=pltdata[i].colour, label=pltdata[i].label, capsize=3+offset, alpha=pltdata[i].alpha)
      end
    elseif pltdata[i].xuerr ≠ nothing
      if (!isempty(pltdata[i].dashes) && pltdata[i].dashes[1] == 0) || pltdata[i].dashes == "None"
        ax.errorbar(pltdata[i].x, pltdata[i].y, xerr=[xerr[i].lower, xerr[i].upper],
          lt= "None", marker=pltdata[i].marker,
          color=pltdata[i].colour, label=pltdata[i].label, capsize=3+offset, alpha=pltdata[i].alpha)
      else
        ax.errorbar(pltdata[i].x, pltdata[i].y, xerr=[xerr[i].lower, xerr[i].upper],
          lw = pltdata[i].lw, marker=pltdata[i].marker, dashes=pltdata[i].dashes,
          color=pltdata[i].colour, label=pltdata[i].label, capsize=3+offset, alpha=pltdata[i].alpha)
      end
    else
      if (!isempty(pltdata[i].dashes) && pltdata[i].dashes[1] == 0) || pltdata[i].dashes == "None"
        ax.scatter(pltdata[i].x, pltdata[i].y,lw = pltdata[i].lw, marker=pltdata[i].marker,
          color=pltdata[i].colour, label=pltdata[i].label, alpha=pltdata[i].alpha)
      else
        ax.plot(pltdata[i].x, pltdata[i].y,lw = pltdata[i].lw, marker=pltdata[i].marker,
          dashes=pltdata[i].dashes, color=pltdata[i].colour, label=pltdata[i].label,
          alpha=pltdata[i].alpha)
      end
    end
  end

  return pltdata, ax
end #function plt_DataWithErrors


"""
    set_axes(pltdata, ax, logscale, xlim, ylim)

Set x- and/or y-axis of `pltdata` to log scale for each `ax`, if the string
`logscale` consists of `"x"` or `"y"` or `"xy"` for the corresponding axis.
Adjust minimum/maximum values of the respective axis for logscale, if `xlim` or
`ylim` are `nothing`.
"""
function set_axes(pltdata, ax, logscale, xlim, ylim)

  xlimits = []; ylimits = []
  for i = 1:length(logscale)
    if occursin('x', lowercase(logscale[i]))
      ax[i].set_xlim(find_limits(pltdata[i], "x", xlim))
      ax[i].set_xscale("log")
    else
      ax[i].set_xlim(xlim)
    end
    if occursin('y', lowercase(logscale[i]))
      ax[i].set_ylim(find_limits(pltdata[i], "y", ylim[i]))
      ax[i].set_yscale("log")
    else
      ax[i].set_ylim(ylim[i])
    end
  end


  return ax
end #function set_axes


"""
    find_limits(pltdata, datacols, lims, limits)

Set the axes `limits` to the minimum and maximum values for log plots in `pltdata`
if `lims` contain `nothing`; `datacols` labels the axes that are logarithmic with
`"x"`, `"y"` or `"xy"`.
"""
function find_limits(pltdata, datacols, lims)

  # Find log axes and corresponding errors
  ctype = Symbol(datacols)
  low = Symbol("$(datacols)lerr"); high = Symbol("$(datacols)uerr")
  # Get data of log axes
  coldata = [getfield(p, ctype) for p in pltdata]
  # Revise minimum, if not pre-defined
  if lims[1] == nothing
    minerr = [getfield(p, low) for p in pltdata]
    minerr = minerr[minerr.≠nothing]
    isempty(minerr) ?
      xmin = minimum(minimum.(coldata)) : xmin = minimum(minimum.(minerr))
    xmin = 10^floor(log10(minimum(xmin)))
  else
    xmin = lims[1]
  end
  # Revise maximum, if not pre-defined
  if lims[2] == nothing
    maxerr = [getfield(p, high) for p in pltdata]
    maxerr = maxerr[maxerr.≠nothing]
    isempty(maxerr) ?
      xmax = maximum(maximum.(coldata)) : xmax = maximum(maximum.(maxerr))
    xmax = 10^ceil(log10(maximum(xmax)))
  else
    xmax = lims[2]
  end

  # Return revised limits
  return (xmin, xmax)
end #function find_limits


"""
    get_stack(plot_list...) -> xdata, ystack

From the `plot_list` with `PlotData`, return a Vector of combined `xdata`
and a vector with the `ydata` of each `PlotData`, where y data of missing `xdata`
values is filled with `NaN`.
"""
function get_stack(plot_list...)
  # Combine all xdata
  xdata = sort(union([p.x for p in plot_list]...))
  labels = [p.label for p in plot_list]
  ystack = []


  # Loop over plot data
  for data in plot_list
    # Find data different from combined xdata
    if data.x ≠ xdata
      # Warn about differences in x data
      @warn "x data in $(data.label) differs from common x data."
      if data.x ≠ sort(unique(data.x))
        # Warn about non-monotonic x data
        @warn string("x data in $(data.label) not strictly monotonic. ",
        "Only first (x, y) value pair considered.")
      end
      # Fill missing x data with zeros
      idx = findfirst.([xdata .== x for x in data.x])
      ydata = NaN.*zeros(length(xdata))
      ydata[idx] .= data.y
      push!(ystack, ydata)
    else
      # or use original y data, if x data doesn't differ from common x data
      push!(ystack, data.y)
    end
  end

  # Return combined x data and stacked y data
  return xdata, ystack, labels
end #function get_stack


"""
    format_stack(plot_list, cs, lc, alpha) -> clr, α

From the `PlotData` in `plot_list`, the definition of a colour scheme `cs` and,
if `cs = "own"`, a list of colours (`lc`), return a list of colours (`clr`).

Set `α` to `alpha`, if alpha is greater `0`, otherwise use the mean of the `alpha`
fields in the `PhotData` or (if `0`) use `1` (no transparency).
"""
function format_stack(cs, alpha, lc, lt, plot_list...)

  # Set transparency
  α = [a.alpha for a in plot_list]
  alpha < 1 ? α = alpha : α = stats.mean(α)
  # Set color scheme
  clr = []; ln = []
  for (i, plt) in enumerate(plot_list)
    if cs=="own"
      # Set own stack colours with lc list
      c = try lc[i]
      catch
        @warn "Colour not defined for data $i. Using default."
        sel_ls("default", lc=i)[1]
      end
      # Set own stack colours with lc list
      l = try lt[i]
      catch
        @warn "Line type not defined for data $i. Using default."
        sel_ls("default", lc=i)[2]
      end
    elseif cs==""
      # Use colours from PlotData
      c = plt.colour
      l = plt.dashes
    else
      # Use colour sccheme define by cs
      c, l = try sel_ls(cs,lc=i,lt=i)[1:2]
      catch
        @warn "Colour scheme $cs not defined. Using default."
        sel_ls("default", lc=i)[1:2]
      end
    end
    push!(clr,c); push!(ln,l)
  end

  return clr, ln, α
end #function format_stack


"""
    print_stack(alpha, xdata, ydata, colours, lt, ax) -> fig, ax

"""
function print_stack(xdata, ystack, ylines, border, labels, colours, lt, α, figsize)

  # Start plot
  fig, ax = plt.subplots(figsize=figsize)

  # Plot data
  ax.stackplot(xdata, ystack, labels=labels, colors=colours, alpha=α)

  # Resume, if optional border are skipped
  if border==0  return fig, ax  end
  [ax.plot(xdata, ylines[i], color=colours[i], dashes = lt[i], alpha=border)
    for i = 1:length(ylines)]

  return fig, ax
end #function print_stack


"""


"""
function format_axes_and_annotations(fig, ax, plot_list, ti, xlabel, ylabel,
  timeformat, timescale, major_interval, minor_interval, xlims, fontsize, legpos, legcolumns,
  axcolour, leg_offset, ti_offset, label_offset, ax_offset,
  major_xticks, major_yticks, minor_xticks, minor_yticks, mticks, ticksize, framewidth)

  # Set plot title
  ax[1].set_title(ti, fontsize=fontsize+ti_offset)

  # Generate axes labels
  ax[1].set_xlabel(xlabel,fontsize=fontsize+label_offset)
  for n = 1:length(ylabel)
    ax[n].set_ylabel(ylabel[n],fontsize=fontsize+label_offset, color=axcolour[n])
  end
  [plt.setp(ax[n].get_yticklabels(),color=axcolour[n]) for n = 1:length(axcolour)]

  # Set/Unset minor ticks
  if mticks == true
    plt.minorticks_on()
  else
    plt.minorticks_off()
  end

  # Format axes and ticks
  if typeof(plot_list[1][1].x) <: Vector{T} where T <: Real
    ax = format_xdata(plot_list, ax, major_xticks, minor_xticks, mticks,
      ticksize, framewidth, fontsize, ax_offset)
  else
    fig, ax = format_timeseries(fig, ax, plot_list, timeformat, timescale, major_interval, minor_interval,
      major_xticks, minor_xticks, mticks, xlims, ticksize, framewidth, fontsize, ax_offset)
  end

  # Format ticks and frame
  for i = 1:length(plot_list)
    ax[i].tick_params("both", which="both", direction="in", top=true, right=true,
      labelsize=fontsize+ax_offset, width=framewidth)
    ax[i].grid(linestyle=":", linewidth = framewidth)
    ax[i].spines["bottom"].set_linewidth(framewidth)
    ax[i].spines["top"].set_linewidth(framewidth)
    ax[i].spines["left"].set_linewidth(framewidth)
    ax[i].spines["right"].set_linewidth(framewidth)
  end

  # Set yticks and optional minor yticks
  for i = 1:length(major_yticks)  if major_yticks[i] > 0
    yint = collect(ax[i].get_ylim()[1]:major_yticks[i]:ax[i].get_ylim()[2])
    ax[i].set_yticks(yint)
  end  end
  # Set minor y ticks
  for i = 1:length(minor_yticks)  if minor_yticks[i] > 0
    my = plt.matplotlib.ticker.MultipleLocator(minor_yticks[i])
    ax[i].yaxis.set_minor_locator(my)
  end  end

  # Print legend
  if length(ax) > 1
    pleg = vcat(ax[1].get_legend_handles_labels()[1], ax[2].get_legend_handles_labels()[1])
    plab = vcat(ax[1].get_legend_handles_labels()[2], ax[2].get_legend_handles_labels()[2])
    if any(plab .≠ "") && legpos ≠ "None"
      ax[2].legend(pleg, plab, fontsize=fontsize+leg_offset, loc=legpos, ncol=legcolumns)
    end
  elseif legpos ≠ "None" && any([p.label≠"" for p in plot_list[1]])
    ax[1].legend(fontsize=fontsize+leg_offset, loc=legpos, ncol=legcolumns)
  end

  # Tight layout for plots with big labels
  if typeof(plot_list[1][1].x) <: Vector{Dates.TimeType}  fig.tight_layout()  end

  return fig, ax
end #function format_axes_and_annotations


function format_xdata(plot_list, ax, major_xticks, minor_xticks, mticks,
  ticksize, framewidth, fontsize, ax_offset)
  # Set ticks and optional minor ticks
  if major_xticks > 0
    xint = collect(ax[1].get_xlim()[1]:major_xticks:ax[1].get_xlim()[2])
    for i = 1:length(plot_list)  ax[i].set_xticks(xint)  end
  end
  # Set minor ticks
  if mticks && minor_xticks > 0
    mx = plt.matplotlib.ticker.MultipleLocator(minor_xticks)
    for i = 1:length(plot_list)
      ax[i].xaxis.set_minor_locator(mx)
    end
  end

  # Format ticks and frame
  for i = 1:length(plot_list)
    ax[i].tick_params("both", which="major", length=ticksize[1]⋅framewidth)
    ax[i].tick_params("both", which="minor", length=ticksize[2]⋅framewidth)
  end

  return ax
end #function format_xdata

function format_timeseries(fig, ax, plot_list, timeformat, timescale, major_interval, minor_interval,
  major_xticks, minor_xticks, mticks, xlims, ticksize, framewidth, fontsize, ax_offset)

  timeformat, majorlocator, minorlocator, majorformatter, minorformatter =
    set_locator_and_formatter(timeformat, timescale, major_xticks, minor_xticks,
    major_interval, minor_interval)

  # Format ticks and frame
  for i = 1:length(plot_list)
    ### Format xlim to beginning of period in timeseries
    ax[i] = set_timeperiod(ax[i], timescale, plot_list[1][1].x, xlims)
    # ax[i].set_xlim(left=plot_list[1][1].x[1], right=plot_list[1][1].x[end])
    ax[i].xaxis.set_major_formatter(majorformatter)
    ax[i].xaxis.set_minor_formatter(minorformatter)
    ax[i].xaxis.set_major_locator(majorlocator)
    ax[i].xaxis.set_minor_locator(minorlocator)
    fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
  end

  return fig, ax
end #function format_timeseries


function set_locator_and_formatter(timeformat, timescale, major_xticks, minor_xticks,
  major_interval, minor_interval)

  if timescale == timescale == "centuries"
    if timeformat == ""  timeformat = "%Y"  end
    if all(major_xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.YearLocator(major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.YearLocator(major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.YearLocator(25)
    end
    if all(minor_xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.YearLocator(minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.YearLocator(minor_interval)
    else #default
      minorlocator = plt.matplotlib.dates.YearLocator(5)
    end
  elseif timescale == "decades"
    if timeformat == ""  timeformat = "%Y"  end
    if all(major_xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.YearLocator(major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.YearLocator(major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.YearLocator(5)
    end
    if all(minor_xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.MonthLocator(interval=minor_interval)
    else #default
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=1)
    end
  elseif timescale == "years"
    if timeformat == ""  timeformat = "%Y"  end
    if all(major_xticks .> 0)  # Set major x ticks
      major_xticks isa Vector ?
        majorlocator = plt.matplotlib.dates.MonthLocator(bymonth=major_xticks) :
        majorlocator = plt.matplotlib.dates.YearLocator(major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.YearLocator(major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.YearLocator()
    end
    if all(minor_xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.MonthLocator(interval=minor_interval)
    else #default
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=[1,4,7,10])
    end
  elseif timescale == "months"
    if timeformat == ""  timeformat = "%B"  end
    if all(major_xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.MonthLocator(bymonth=major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.MonthLocator(interval=major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.MonthLocator(interval=1)
    end
    if all(minor_xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.DayLocator(bymonthday=minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.DayLocator(interval=minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.DayLocator(bymonthday=16)
    end
  elseif timescale == "weeks"
    if timeformat == ""  timeformat = "%-d. %b"  end
    if all(major_xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.DayLocator(bymonthday=major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.DayLocator(interval=major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.DayLocator(interval=1)
    end
    if all(minor_xticks .≥ 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.HourLocator(byhour=minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.HourLocator(interval=minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.HourLocator(byhour=[0,6,12,18])
    end
  elseif timescale == "days"
    if timeformat == ""
      major_xticks isa Vector ? timeformat = "%-d. %b, %H:%M" : timeformat = "%-d. %b"
    end
    if (major_xticks isa Vector && all(major_xticks .≥ 0)) || major_xticks > 0  # Set major x ticks
      major_xticks isa Vector ?
        majorlocator = plt.matplotlib.dates.HourLocator(byhour=major_xticks) :
        majorlocator = plt.matplotlib.dates.DayLocator(bymonthday=major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.DayLocator(interval=major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.DayLocator(interval=1)
    end
    if all(minor_xticks .≥ 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.HourLocator(byhour=minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.HourLocator(interval=minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.HourLocator(byhour=[0,6,12,18])
    end
  elseif timescale == "hours"
    if timeformat == ""  timeformat = "%H:%M"  end
    if all(major_xticks .≥ 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.HourLocator(byhour=major_xticks)
    elseif major_interval > 0
      majorlocator = plt.matplotlib.dates.HourLocator(interval=major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.HourLocator(interval=1)
    end
    if all(minor_xticks .≥ 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.MinuteLocator(byminute=minor_xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.MinuteLocator(interval=minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.MinuteLocator(byminute=[0,15,30,45])
    end
  end
  majorformatter = plt.matplotlib.dates.DateFormatter(timeformat)
  minorformatter = plt.matplotlib.dates.DateFormatter("")

  return timeformat, majorlocator, minorlocator, majorformatter, minorformatter
end

function set_timeperiod(ax, timescale, xdata, xlim)
  if xlim[1] == nothing
    if timescale == "centuries"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Year(25)))
    elseif timescale == "decades"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Year(5)))
    elseif timescale == "years"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Year))
    elseif timescale == "months"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Month))
    elseif timescale == "weeks" || timescale == "days"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Day))
    elseif timescale == "hours"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Hour))
    end
  end
  if xlim[2] == nothing
    if timescale == "centuries"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Year(25)))
    elseif timescale == "decades"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Year(5)))
    elseif timescale == "years"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Year))
    elseif timescale == "months"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Month))
    elseif timescale == "weeks" || timescale == "days"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Day))
    elseif timescale == "hours"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Hour))
    end
  end

  return ax
end #function set_timeperiod
