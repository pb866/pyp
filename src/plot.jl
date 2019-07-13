### Functions associated with plotting data

"""
    setup_plot(plot_list, twinax, ylab, logscale, logremove, xlims, ylims, Yticks, yticks, ptype, cs, lc, lt, pt, alpha, axcolour)

If `twinax` is an array split `plot_list` into 2 datasets based on the indices
in `twinax`. Assure all paramters concerning y data are arrays of length 2, if
a second axis is allowed or length 1 otherwise. Return `plot_list` as array of
2 or 1 distinct lists an the adjusted parameters.
"""
function setup_plot(plot_list, kw::kwargs)

  # Start plot
  fig, ax = plt.subplots(figsize=kw.figsize)
  ax = [ax]
  pltdata = deepcopy(plot_list)

  # Set alpha
  [p.alpha = kw.alpha for p in pltdata if kw.alpha > 0]
  # Initialise parameters
  plt1 = PlotData[]; plt2 = PlotData[]
  cl1 = []; dt1 = []; mt1 = []
  cl2 = []; dt2 = []; mt2 = []

  if !isempty(kw.twinax)
    # Set 2nd axis in PyPlot, if twinax is not empty
    push!(ax, plt.twinx())
    # Check correct input of twinax
    if length(kw.twinax) ≠ length(pltdata)
      throw(ArgumentError(string("Array `twinax` must have the same length ",
        "as there are number of `PlotData` elements.")))
    end

    # Assign data to the axes
    for i = 1:length(kw.twinax)
      if kw.twinax[i] == 1
        push!(plt1, pltdata[i])
        # check whether an array of definitions or a single definition exists
        # and assign either the current array entry to the respective axis
        # or the default value for line colour, type and marker type
        !isempty(kw.lt) && (kw.lt[i] isa Tuple || kw.lt[i] isa Vector || kw.lt[i] isa String) ?
          push!(dt1, kw.lt[i]) : push!(dt1, kw.lt)
        kw.pt isa Vector ? push!(mt1, kw.pt[i]) : push!(mt1, kw.pt)
        kw.lc isa Vector ? push!(cl1, kw.lc[i]) : push!(cl1, kw.lc)
      else
        push!(plt2, pltdata[i])
        # check whether an array of definitions or a single definition exists
        # and assign either the current array entry to the respective axis
        # or the default value for line colour, type and marker type
        !isempty(kw.lt) && (kw.lt[i] isa Tuple || kw.lt[i] isa Vector || kw.lt[i] isa String) ?
          push!(dt2, kw.lt[i]) : push!(dt2, kw.lt)
        kw.pt isa Vector ? push!(mt2, kw.pt[i]) : push!(mt2, kw.pt)
        kw.lc isa Vector ? push!(cl2, kw.lc[i]) : push!(cl2, kw.lc)
      end
    end

    # Make sure, all parameters for both axes are arrays of length 2
    if kw.logscale isa String kw.logscale = String[kw.logscale, kw.logscale]  end
    if kw.logremove isa String kw.logremove = String[kw.logremove, kw.logremove]  end
    if kw.axcolour isa String kw.axcolour = String[kw.axcolour, kw.axcolour]  end
    if kw.Yticks isa Int64 kw.Yticks = Int64[kw.Yticks, kw.Yticks]  end
    if kw.yticks isa Int64 kw.yticks = Int64[kw.yticks, kw.yticks]  end
    if kw.xlim == nothing
      kw.xlim = (nothing, nothing)
    elseif kw.xlim isa Array
      throw(ArgumentError("xlims must be a tuple or nothing!"))
    end
    if kw.ylim == nothing
      kw.ylim = [(nothing, nothing), (nothing, nothing)]
    elseif kw.ylim isa Tuple
      kw.ylim = [kw.ylim, kw.ylim]
    end
  else
    # Assign all data to axis 1, if no second axis
    plt1 = pltdata
    # check whether an array of definitions or a single definition exists
    # and assign either the current array entry to the respective axis
    # or the default value for line colour, type and marker type
    !isempty(kw.lt) && (kw.lt[1] isa Tuple || kw.lt[1] isa Vector || kw.lt[1] isa String) ?
      dt1 = kw.lt : dt1 = [kw.lt for i in pltdata]
    kw.pt isa Vector ? mt1 = kw.pt : mt1 = [kw.pt for i in pltdata]
    kw.lc isa Vector ? cl1 = kw.lc : cl1 = [kw.lc for i in pltdata]
    # If no 2nd axis, make sure, all parameters for both axes are arrays of length 1
    # and not single parameters (numbers, strings...)
    if kw.logscale isa String kw.logscale = String[kw.logscale]  end
    if kw.logremove isa String kw.logremove = String[kw.logremove]  end
    if kw.axcolour isa String kw.axcolour = String[kw.axcolour]  end
    if kw.Yticks isa Real Yticks = Int64[kw.Yticks]  end
    if kw.yticks isa Real yticks = Int64[kw.yticks]  end
    if kw.xlim == nothing     kw.xlim = (nothing, nothing)  end
    if kw.ylim == nothing     kw.ylim = [(nothing, nothing)]
    elseif kw.ylim isa Tuple  kw.ylim = [kw.ylim]
    end
  end
  if kw.ylabel isa AbstractString  kw.ylabel = String[kw.ylabel]  end
  if length(ax) > 1
    pltdata = [plt1, plt2]
    if !isa(kw.cs, Vector) kw.cs = [kw.cs, kw.cs]  end
    kw.lc = [cl1, cl2]
    kw.lt = [dt1, dt2]
    kw.pt = [mt1, mt2]
  else
    pltdata = [plt1]
    kw.cs = [kw.cs]
    kw.lc = [cl1]
    kw.lt = [dt1]
    kw.pt = [mt1]
  end

  # Return adjusted data
  return pltdata, fig, ax, kw
end #function setup_plot


"""
    set_style(pltdata, ptype, cs, lc, lt, pt, ax_2)

Set colour scheme for `pltdata` holding `PlotData` to `cs` from
function `sel_ls`. If colour scheme is set to `"own"` or `"mix"` for any array
element in `pltdata`, `lc`, `lt`, and `pt` must hold valid definitions of line colours,
line and point (marker) types at this array position. The boolean `ax_2` is needed
to specify, whether a second axis is used.
"""
function set_style(pltdata, kw::kwargs)

  # Reset line and marker types to solid lines and squares, respectively,
  # if plot_type is set to line, scatter or both (markers connected by lines)
  for i = 1:length(pltdata)
    if kw.plottype == "line"
      [pltdata[i][j].dashes = Real[] for j = 1:length(kw.lt[i]) if pltdata[i][j].dashes=="None"]
      [pltdata[i][j].marker = "None" for j = 1:length(kw.pt[i])]
    elseif kw.plottype == "scatter"
      [pltdata[i][j].dashes = "None" for j = 1:length(kw.lt[i])]
      [pltdata[i][j].marker = "s" for j = 1:length(kw.pt[i]) if pltdata[i][j].marker=="None"]
    elseif kw.plottype == "both"
      [pltdata[i][j].dashes = Real[] for j = 1:length(kw.lt[i]) if pltdata[i][j].dashes=="None"]
      [pltdata[i][j].marker = "s" for j = 1:length(kw.pt[i]) if pltdata[i][j].marker=="None"]
    end
  end

  # Set index for sel_ls function, start at 1 for separate schemes
  # and use continuous indexing for the same scheme
  idx = []
  push!(idx, [i for i = 1:length(pltdata[1])])
  if length(kw.cs) == 2
    ind = kw.cs[1] == kw.cs[2] ? length(idx[1]) : 0
    push!(idx, [i+ind for i = 1:length(pltdata[2])])
  end
  # Set line/marker styles and colour for a chosen colour scheme
  # (overwrites default settings from plot_type)
  for i = 1:length(pltdata)  if kw.cs[i]≠""
    for j = 1:length(pltdata[i])
      cl, dt, mt = sel_ls(cs=kw.cs[i], lc=idx[i][j], lt=idx[i][j], pt=idx[i][j])
      if kw.lc[i][j] ≠ "default"  pltdata[i][j].colour = cl  end
      if kw.plottype ≠ "scatter" && kw.lt[i][j] ≠ "default"  pltdata[i][j].dashes = dt  end
      if kw.plottype ≠ "line" && kw.pt[i][j] ≠ "default"  pltdata[i][j].marker = mt  end
    end
  end  end

  # Set manually defined styles (overwrites previous settings)
  for i = 1:length(kw.lc), j = 1:length(kw.lc[i])
    if !(kw.lc[i][j]=="" || kw.lc[i][j]=="default")
      pltdata[i][j].colour = kw.lc[i][j]
    elseif isnothing(pltdata[i][j].colour)
      pltdata[i][j].colour = sel_ls(cs="pyplot", lc=idx[i][j])[1]
    end
  end
  for i = 1:length(kw.lc), j = 1:length(kw.lc[i])
    if !(kw.lt[i][j]=="" || kw.lt[i][j]=="default") && pltdata[i][j].dashes ≠ "None"
      pltdata[i][j].dashes = kw.lt[i][j]
    end
  end
  for i = 1:length(kw.lc), j = 1:length(kw.lc[i])
    if !(kw.pt[i][j]=="" || kw.pt[i][j]=="default") && pltdata[i][j].dashes ≠ "None"
      pltdata[i][j].marker = kw.pt[i][j]
    end
  end

  # Return adjusted PlotData
  return pltdata
end #function set_style


"""
    plotDataWithErrors(pltdata, ax, offset)

For each `PlotData` in array `pltdata` (and `ax`), retrieve the errors in the `PlotData` and
plot according to specifications. For error bars of markers, the standard cap size is `3`,
which can be adjusted by an `offset` (positive or negative number to be added to cap size).

Returns `fig` and `ax` (array of axes) for further plot modifications or printing.
"""
function plotDataWithErrors(pltdata, ax, offset)
  # Redefine errors for error bars
  xerr = redef_err(pltdata,:x,:x_lower,:x_upper)
  yerr = redef_err(pltdata,:y,:y_lower,:y_upper)
  # Loop over graph data and plot each data according to its type and format
  # defined by the struct PlotData
  for i = 1:length(pltdata)
    if !isempty(pltdata[i].y_upper) && pltdata[i].marker == "None"
      p = ax.plot(pltdata[i].x, pltdata[i].y, lw = pltdata[i].lw,
          dashes=pltdata[i].dashes, color=pltdata[i].colour, label=pltdata[i].label)
      pltdata[i].colour = p[1].get_color()
      ax.fill_between(pltdata[i].x, pltdata[i].y_lower, pltdata[i].y_upper,
          color=pltdata[i].colour, alpha=0.2)
    elseif !isempty(pltdata[i].x_upper) && !isempty(pltdata[i].y_upper)
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
    elseif !isempty(pltdata[i].y_upper)
      if (!isempty(pltdata[i].dashes) && pltdata[i].dashes[1] == 0) || pltdata[i].dashes == "None"
        ax.errorbar(pltdata[i].x, pltdata[i].y, yerr=[yerr[i].lower, yerr[i].upper],
          fmt=pltdata[i].marker, color=pltdata[i].colour, label=pltdata[i].label, capsize=3+offset,
          alpha=pltdata[i].alpha)
      else
        ax.errorbar(pltdata[i].x, pltdata[i].y, yerr=[yerr[i].lower, yerr[i].upper],
          lw = pltdata[i].lw, marker=pltdata[i].marker, dashes=pltdata[i].dashes,
          color=pltdata[i].colour, label=pltdata[i].label, capsize=3+offset, alpha=pltdata[i].alpha)
      end
    elseif !isempty(pltdata[i].x_upper)
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
end #function plotDataWithErrors


"""
    set_axes(pltdata, ax, logscale, xlim, ylim)

Set x- and/or y-axis of `pltdata` to log scale for each `ax`, if the string
`logscale` consists of `"x"` or `"y"` or `"xy"` for the corresponding axis.
Adjust minimum/maximum values of the respective axis for logscale, if `xlim` or
`ylim` are `nothing`.
"""
function set_axes(pltdata, ax, kw::kwargs)

  for i = 1:length(kw.logscale)
    if occursin('x', lowercase(kw.logscale[i]))
      ax[i].set_xlim(find_limits(pltdata[i], "x", kw.xlim))
      ax[i].set_xscale("symlog")
    else
      ax[i].set_xlim(kw.xlim)
    end
    if occursin('y', lowercase(kw.logscale[i]))
      ax[i].set_ylim(find_limits(pltdata[i], "y", kw.ylim[i]))
      ax[i].set_yscale("symlog")
    else
      ax[i].set_ylim(kw.ylim[i])
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
  low = Symbol("$(datacols)_lower"); high = Symbol("$(datacols)_upper")
  # Get data of log axes
  coldata = [getfield(p, ctype)[isfinite.(getfield(p, ctype))] for p in pltdata]
  # Revise minimum, if not pre-defined
  if lims[1] == nothing
    minerr = [getfield(p, low) for p in pltdata]
    minerr = minerr[@. !isempty(minerr)]
    isempty(minerr) ?
      xmin = minimum(minimum.(coldata)) : xmin = minimum(minimum.(minerr))
    xmin = 10^floor(log10(minimum(xmin)))
  else
    xmin = lims[1]
  end
  # Revise maximum, if not pre-defined
  if lims[2] == nothing
    maxerr = [getfield(p, high) for p in pltdata]
    maxerr = maxerr[@. !isempty(maxerr)]
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
function format_stack(kw::kwargs, plot_list...)

  # Set opacity
  if kw.alpha == 0
    α = [a.alpha for a in plot_list]
    kw.alpha = stats.mean(α) > 0 ? stats.mean(α) : 1
  end
  # Set color scheme
  clr = []; ln = []
  for (i, plt) in enumerate(plot_list)
    if kw.cs=="own"
      # Set own stack colours with lc list
      c = try kw.lc[i]
      catch
        @warn "Colour not defined for data $i. Using default."
        sel_ls(cs="default", lc=i)[1]
      end
      # Set own stack colours with lc list
      l = try kw.lt[i]
      catch
        @warn "Line type not defined for data $i. Using default."
        sel_ls(cs="default", lt=i)[2]
      end
    elseif kw.cs==""
      # Use colours from PlotData
      c = plt.colour
      l = plt.dashes
    else
      # Use colour sccheme define by cs
      c, l = try sel_ls(cs=kw.cs,lc=i,lt=i)[1:2]
      catch
        @warn "Colour scheme $(kw.cs) not defined. Using default."
        sel_ls(cs="default", lc=i)[1:2]
      end
    end
    push!(clr,c); push!(ln,l)
  end
  kw.lc = clr; kw.lt = ln

  return kw
end #function format_stack


"""
    print_stack(alpha, xdata, ydata, colours, lt, ax) -> fig, ax

"""
function print_stack(xdata, ystack, ylines, labels, kw::kwargs)

  # Start plot
  fig, ax = plt.subplots(figsize=kw.figsize)

  # Plot data
  ax.stackplot(xdata, ystack, labels=labels, colors=kw.lc, alpha=kw.alpha)

  # Resume, if optional border are skipped
  if kw.border==0  return fig, ax  end
  [ax.plot(xdata, ylines[i], color=kw.lc[i], dashes = kw.lt[i], alpha=kw.border)
    for i = 1:length(ylines)]

  return fig, ax
end #function print_stack


"""


"""
function format_axes_and_annotations(fig, ax, plot_list, kw::kwargs)

  # Set plot title
  ax[1].set_title(kw.ti, fontsize=kw.fontsize+kw.ti_offset)

  # Generate axes labels
  ax[1].set_xlabel(kw.xlabel,fontsize=kw.fontsize+kw.lbl_offset)
  for n = 1:length(kw.ylabel)
    ax[n].set_ylabel(kw.ylabel[n],fontsize=kw.fontsize+kw.lbl_offset, color=kw.axcolour[n])
  end
  [plt.setp(ax[n].get_yticklabels(),color=kw.axcolour[n]) for n = 1:length(kw.axcolour)]

  # Set/Unset minor ticks
  if kw.mticks == true
    plt.minorticks_on()
  else
    plt.minorticks_off()
  end

  # Format axes and ticks
  if typeof(plot_list[1][1].x) <: Vector{<:Real}
    ax = format_xdata(plot_list, ax, kw)
  else
    fig, ax = format_timeseries(fig, ax, plot_list, kw)
  end

  # Format ticks and frame
  for i = 1:length(plot_list)
    ax[i].tick_params("both", which="both", direction="in", top=true, right=true,
      labelsize=kw.fontsize+kw.tick_offset, width=kw.framewidth)
    ax[i].grid(linestyle=":", linewidth=kw.framewidth)
    ax[i].spines["bottom"].set_linewidth(kw.framewidth)
    ax[i].spines["top"].set_linewidth(kw.framewidth)
    ax[i].spines["left"].set_linewidth(kw.framewidth)
    ax[i].spines["right"].set_linewidth(kw.framewidth)
  end

  # Set yticks and optional minor yticks
  for i = 1:length(kw.Yticks)  if kw.Yticks[i] > 0
    yint = collect(ax[i].get_ylim()[1]:kw.Yticks[i]:ax[i].get_ylim()[2])
    ax[i].set_yticks(yint)
  end  end
  # Set minor y ticks
  for i = 1:length(kw.yticks)  if kw.yticks[i] > 0
    my = plt.matplotlib.ticker.MultipleLocator(kw.yticks[i])
    ax[i].yaxis.set_minor_locator(my)
  end  end

  # Print legend
  if length(ax) > 1
    pleg = vcat(ax[1].get_legend_handles_labels()[1], ax[2].get_legend_handles_labels()[1])
    plab = vcat(ax[1].get_legend_handles_labels()[2], ax[2].get_legend_handles_labels()[2])
    if any(plab .≠ "") && kw.legpos ≠ "None"
      ax[2].legend(pleg, plab, fontsize=kw.fontsize+kw.leg_offset, loc=kw.legpos, ncol=kw.legcols)
    end
  elseif kw.legpos ≠ "None" && any([p.label≠"" for p in plot_list[1]])
    ax[1].legend(fontsize=kw.fontsize+kw.leg_offset, loc=kw.legpos, ncol=kw.legcols)
  end

  # Tight layout for plots with big labels
  if typeof(plot_list[1][1].x) <: Vector{Dates.TimeType}  fig.tight_layout()  end

  return fig, ax
end #function format_axes_and_annotations


function format_xdata(plot_list, ax, kw::kwargs)
  # Set ticks and optional minor ticks
  if kw.Xticks > 0
    xint = collect(ax[1].get_xlim()[1]:kw.Xticks:ax[1].get_xlim()[2])
    for i = 1:length(plot_list)  ax[i].set_xticks(xint)  end
  end
  # Set minor ticks
  if kw.mticks && kw.xticks > 0
    mx = plt.matplotlib.ticker.MultipleLocator(kw.xticks)
    for i = 1:length(plot_list)
      ax[i].xaxis.set_minor_locator(mx)
    end
  end

  # Format ticks and frame
  for i = 1:length(plot_list)
    ax[i].tick_params("both", which="major", length=kw.ticksize[1]⋅kw.framewidth)
    ax[i].tick_params("both", which="minor", length=kw.ticksize[2]⋅kw.framewidth)
  end

  return ax
end #function format_xdata

function format_timeseries(fig, ax, plot_list, kw::kwargs)

  kw.timeformat, majorlocator, minorlocator, majorformatter, minorformatter =
    set_locator_and_formatter(kw)

  # Format ticks and frame
  for i = 1:length(plot_list)
    ### Format xlim to beginning of period in timeseries
    ax[i] = set_timeperiod(ax[i], plot_list[1][1].x, kw)
    # ax[i].set_xlim(left=plot_list[1][1].x[1], right=plot_list[1][1].x[end])
    ax[i].xaxis.set_major_formatter(majorformatter)
    ax[i].xaxis.set_minor_formatter(minorformatter)
    ax[i].xaxis.set_major_locator(majorlocator)
    ax[i].xaxis.set_minor_locator(minorlocator)
    fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
  end

  return fig, ax
end #function format_timeseries


function set_locator_and_formatter(kw::kwargs)

  if kw.timescale == "centuries"
    if kw.timeformat == ""  kw.timeformat = "%Y"  end
    if all(kw.Xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.YearLocator(kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.YearLocator(kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.YearLocator(25)
    end
    if all(kw.xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.YearLocator(kw.xticks)
    elseif kw.minor_interval > 0
      minorlocator = plt.matplotlib.dates.YearLocator(kw.minor_interval)
    else #default
      minorlocator = plt.matplotlib.dates.YearLocator(5)
    end
  elseif kw.timescale == "decades"
    if kw.timeformat == ""  kw.timeformat = "%Y"  end
    if all(kw.Xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.YearLocator(kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.YearLocator(kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.YearLocator(5)
    end
    if all(kw.xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=kw.xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.MonthLocator(interval=kw.minor_interval)
    else #default
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=1)
    end
  elseif kw.timescale == "years"
    if kw.timeformat == ""  kw.timeformat = "%Y"  end
    if all(kw.Xticks .> 0)  # Set major x ticks
      kw.Xticks isa Vector ?
        majorlocator = plt.matplotlib.dates.MonthLocator(bymonth=kw.Xticks) :
        majorlocator = plt.matplotlib.dates.YearLocator(kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.YearLocator(kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.YearLocator()
    end
    if all(kw.xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=kw.xticks)
    elseif minor_interval > 0
      minorlocator = plt.matplotlib.dates.MonthLocator(interval=kw.minor_interval)
    else #default
      minorlocator = plt.matplotlib.dates.MonthLocator(bymonth=[1,4,7,10])
    end
  elseif kw.timescale == "months"
    if kw.timeformat == ""  kw.timeformat = "%B"  end
    if all(kw.Xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.MonthLocator(bymonth=kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.MonthLocator(interval=kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.MonthLocator(interval=1)
    end
    if all(kw.xticks .> 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.DayLocator(bymonthday=kw.xticks)
    elseif kw.minor_interval > 0
      minorlocator = plt.matplotlib.dates.DayLocator(interval=kw.minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.DayLocator(bymonthday=16)
    end
  elseif kw.timescale == "weeks"
    if kw.timeformat == ""  kw.timeformat = "%-d. %b"  end
    if all(kw.Xticks .> 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.DayLocator(bymonthday=kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.DayLocator(interval=kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.DayLocator(interval=1)
    end
    if all(kw.xticks .≥ 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.HourLocator(byhour=kw.xticks)
    elseif kw.minor_interval > 0
      minorlocator = plt.matplotlib.dates.HourLocator(interval=kw.minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.HourLocator(byhour=[0,6,12,18])
    end
  elseif kw.timescale == "days"
    if kw.timeformat == ""
      kw.Xticks isa Vector ? kw.timeformat = "%-d. %b, %H:%M" : kw.timeformat = "%-d. %b"
    end
    if (kw.Xticks isa Vector && all(kw.Xticks .≥ 0)) || kw.Xticks > 0  # Set major x ticks
      kw.Xticks isa Vector ?
        majorlocator = plt.matplotlib.dates.HourLocator(byhour=kw.Xticks) :
        majorlocator = plt.matplotlib.dates.DayLocator(bymonthday=kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.DayLocator(interval=kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.DayLocator(interval=1)
    end
    if all(kw.xticks .≥ 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.HourLocator(byhour=kw.xticks)
    elseif kw.minor_interval > 0
      minorlocator = plt.matplotlib.dates.HourLocator(interval=kw.minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.HourLocator(byhour=[0,6,12,18])
    end
  elseif kw.timescale == "hours"
    if kw.timeformat == ""  kw.timeformat = "%H:%M"  end
    if all(kw.Xticks .≥ 0)  # Set major x ticks
      majorlocator = plt.matplotlib.dates.HourLocator(byhour=kw.Xticks)
    elseif kw.major_interval > 0
      majorlocator = plt.matplotlib.dates.HourLocator(interval=kw.major_interval)
    else #default
      majorlocator = plt.matplotlib.dates.HourLocator(interval=1)
    end
    if all(kw.xticks .≥ 0)  # Set minor x ticks
      minorlocator = plt.matplotlib.dates.MinuteLocator(byminute=kw.xticks)
    elseif kw.minor_interval > 0
      minorlocator = plt.matplotlib.dates.MinuteLocator(interval=kw.minor_interval)
    else #default
        minorlocator = plt.matplotlib.dates.MinuteLocator(byminute=[0,15,30,45])
    end
  end
  majorformatter = plt.matplotlib.dates.DateFormatter(kw.timeformat)
  minorformatter = plt.matplotlib.dates.DateFormatter("")

  return kw.timeformat, majorlocator, minorlocator, majorformatter, minorformatter
end


function set_timeperiod(ax, xdata, kw::kwargs)
  if kw.xlim[1] == nothing
    if kw.timescale == "centuries"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Year(25)))
    elseif kw.timescale == "decades"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Year(5)))
    elseif kw.timescale == "years"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Year))
    elseif kw.timescale == "months"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Month))
    elseif kw.timescale == "weeks" || timescale == "days"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Day))
    elseif kw.timescale == "hours"
      ax.set_xlim(left=Dates.floor(xdata[1], Dates.Hour))
    end
  end
  if kw.xlim[2] == nothing
    if kw.timescale == "centuries"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Year(25)))
    elseif kw.timescale == "decades"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Year(5)))
    elseif kw.timescale == "years"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Year))
    elseif kw.timescale == "months"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Month))
    elseif kw.timescale == "weeks" || timescale == "days"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Day))
    elseif kw.timescale == "hours"
      ax.set_xlim(right=Dates.ceil(xdata[end], Dates.Hour))
    end
  end

  return ax
end #function set_timeperiod
