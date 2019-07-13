### Functions associated with load_PlotData

"""
    renameDF(pltdata, DFnames, err)

Rename columns in DataFrame `pltdata` with names specified in array `DFnames` with
entries for x, y, lower y error, upper y error, lower x error, and upper x error
based on the keyword given in `err`.

`err` (`String`):
- `"None"` (no errors, **default**)
- `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
- `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
- `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
- `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
- `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
- `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
- `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
"""
function renameDF(pltdata, DFnames, err)
  if err == "None"
    df.names!(pltdata,DFnames[1:2])
  elseif (startswith(err,"pm") && endswith(err,"x")) || err == "valuex"
    df.names!(pltdata,DFnames[[1,2,5,6]])
  elseif endswith(err,"x")
    df.names!(pltdata,DFnames[[1,2,5]])
  elseif (startswith(err,"pm") && endswith(err,"y")) || err == "valuey"
    df.names!(pltdata,DFnames[1:4])
  elseif endswith(err,"y")
    df.names!(pltdata,DFnames[1:3])
  elseif startswith(err,"pm") || err == "value"
    df.names!(pltdata,DFnames)
  else
    df.names!(pltdata,DFnames[[1,2,3,5]])
  end

  return pltdata
end #function renameDF


"""
    calc_errors(pltdata, col, err)

Calculate errors in DataFrame `pltdata` for columns with the specified column
names `col` based on the keyword `err`.

`err` (`String`):
- `"None"` (no errors, **default**)
- `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
- `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
- `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
- `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
- `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
- `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
- `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
"""
function calc_errors(pltdata, col, err)
  # ± errors (error range for x, y or both)
  # where `err` starts with `pm`, different values are
  # used for positive and negative errors
  if err == "rangex"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
  elseif err == "pmrangex"
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]]
  elseif err == "rangey"
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
  elseif err == "pmrangey"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]]
  elseif err == "range"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
  elseif err == "pmrange"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]]

  # ± (factor × value) (adding a percentage range)
  elseif err == "percentx"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]].⋅abs.(pltdata[col[1]])
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅abs.(pltdata[col[1]])
  elseif err == "pmpercentx"
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅abs.(pltdata[col[1]])
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]].⋅abs.(pltdata[col[1]])
  elseif err == "percenty"
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]].⋅abs.(pltdata[col[2]])
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅abs.(pltdata[col[2]])
  elseif err == "pmpercenty"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅abs.(pltdata[col[2]])
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]].⋅abs.(pltdata[col[2]])
  elseif err == "percent"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]].⋅abs.(pltdata[col[1]])
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅abs.(pltdata[col[1]])
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]].⋅abs.(pltdata[col[2]])
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅abs.(pltdata[col[2]])
  elseif err == "pmpercent"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅abs.(pltdata[col[2]])
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]].⋅abs.(pltdata[col[2]])
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅abs.(pltdata[col[1]])
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]].⋅abs.(pltdata[col[1]])

  # error × value (adding a factor error range)
  elseif err == "factorx"
    pltdata[col[5]][pltdata[col[5]].<1] .= inv.(pltdata[col[5]][pltdata[col[5]].<1])
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[5]].^sign.(pltdata[col[1]])
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]].^sign.(pltdata[col[1]])
  elseif err == "pmfactorx"
    pltdata[col[5]][pltdata[col[5]].<1] .= inv.(pltdata[col[5]][pltdata[col[5]].<1])
    pltdata[col[6]][pltdata[col[6]].<1] .= inv.(pltdata[col[6]][pltdata[col[6]].<1])
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]].^sign.(pltdata[col[1]])
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[6]].^sign.(pltdata[col[1]])
  elseif err == "factory"
    pltdata[col[3]][pltdata[col[3]].<1] .= inv.(pltdata[col[3]][pltdata[col[3]].<1])
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[3]].^sign.(pltdata[col[2]])
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]].^sign.(pltdata[col[2]])
  elseif err == "pmfactory"
    pltdata[col[3]][pltdata[col[3]].<1] .= inv.(pltdata[col[3]][pltdata[col[3]].<1])
    pltdata[col[4]][pltdata[col[4]].<1] .= inv.(pltdata[col[4]][pltdata[col[4]].<1])
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]].^sign.(pltdata[col[2]])
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[4]].^sign.(pltdata[col[2]])
  elseif err == "factor"
    pltdata[col[3]][pltdata[col[3]].<1] .= inv.(pltdata[col[3]][pltdata[col[3]].<1])
    pltdata[col[5]][pltdata[col[5]].<1] .= inv.(pltdata[col[5]][pltdata[col[5]].<1])
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[5]].^sign.(pltdata[col[1]])
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]].^sign.(pltdata[col[1]])
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[3]].^sign.(pltdata[col[2]])
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]].^sign.(pltdata[col[2]])
  elseif err == "pmfactor"
    pltdata[col[3]][pltdata[col[3]].<1] .= inv.(pltdata[col[3]][pltdata[col[3]].<1])
    pltdata[col[4]][pltdata[col[4]].<1] .= inv.(pltdata[col[4]][pltdata[col[4]].<1])
    pltdata[col[5]][pltdata[col[5]].<1] .= inv.(pltdata[col[5]][pltdata[col[5]].<1])
    pltdata[col[6]][pltdata[col[6]].<1] .= inv.(pltdata[col[6]][pltdata[col[6]].<1])
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]].^sign.(pltdata[col[2]])
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[4]].^sign.(pltdata[col[2]])
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]].^sign.(pltdata[col[1]])
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[6]].^sign.(pltdata[col[1]])
  end

  return pltdata
end #function calc_errors


"""
    redef_err(pltdata,val,low,high)

Recalculate errors relative to actual values (rather than absolute values in plot)
for use with error bars of markers. The function generates an array of DataFrames
with the revised errors in the columns `:upper` and `:lower`.
Data is generated from the array `pltdata` with the columns `val` with the measured/modelled
value and the columns with names `low` and `high` holding the
actual error values (rather than values relative to measured/modelled value.)
"""
function redef_err(pltdata,val,low,high)
  err = []
  for i = 1:length(pltdata)
    # Recalculate errors for error bars, if markers are plotted
    if pltdata[i].marker ≠ "None" && !isempty(getfield(pltdata[i],high))
      push!(err, DataFrame(upper = getfield(pltdata[i],high) .- getfield(pltdata[i],val),
            lower = getfield(pltdata[i],val) .- getfield(pltdata[i],low)))
    else
      # Otherwise, save errors as they are und column names `upper` and `lower`
      push!(err, DataFrame(upper = Float64[], lower = Float64[]))
    end
  end

  return err
end #function redef_err


"""
function setup_log(pltdata, logremove::Union{String, Vector{String}}, logscale::Union{String, Vector{String}})

Depending on the keyword of `logremove`, set all positive or negative values to
zero in PlotData `pltdata`, if `logscale` is set.
"""
function setup_log(pltdata, kw::kwargs)

  for i = 1:length(pltdata)
    if kw.logremove[i] == "pos" && occursin("x", kw.logscale[i])
      pltdata[i] = rm_log(pltdata[i], "x", >)
    end
    if kw.logremove[i] == "pos" && occursin("y", kw.logscale[i])
      pltdata[i] = rm_log(pltdata[i], "y", >)
    end
    if kw.logremove[i] == "neg" && occursin("x", kw.logscale[i])
      pltdata[i] = rm_log(pltdata[i], "x", <)
    end
    if kw.logremove[i] == "neg" && occursin("y", kw.logscale[i])
      pltdata[i] = rm_log(pltdata[i], "y", <)
    end
  end

  return pltdata
end #function setup_log


"""
    rm_log(p, x::String, rel)

In the PlotData of the current axis `p`, set x or y data (including errors)
defined by `x` to zero, if it is `>` or `<` `0` as defined by `rel`.
"""
function rm_log(p, x::String, rel)
  for i = 1:length(p)
    rv = findall(rel.(getfield(p[i], Symbol(x)), 0.))
    pv = getfield(p[i], Symbol(x)); pv[rv] .= 0.
    if !isempty(getfield(p[i], Symbol("$(x)_lower")))
      rl = findall(rel.(getfield(p[i], Symbol("$(x)_lower")), 0.))
      pl = getfield(p[i], Symbol("$(x)_lower")); pl[rl] .= 0.
    end
    if !isempty(getfield(p[i], Symbol("$(x)_upper")))
      ru = findall(rel.(getfield(p[i], Symbol("$(x)_upper")), 0.))
      pu = getfield(p[i], Symbol("$(x)_upper")); pu[ru] .= 0.
    end
  end

  return p
end


function interpolate_stack(xdata, ystack, kw::kwargs)
  ranges, extraranges = zip(index2range.(ystack)...)
  if kw.interpolate == "linreg" || kw.interpolate == "linear"
    kw.interpolate = "spline"
    kw.kspline = 1
  end
  if kw.interpolate == "ffill"
    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= ystack[i][j[1]-1]
    end
    if !(kw.extrapolate == false || kw.extrapolate == "None")
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][end]] .= ystack[i][extraranges[i][end][1]-1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][end]] .= 0 for i in 1:length(extraranges)]
    end
    if kw.extrapolate == "both" || kw.extrapolate == "nearest"
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][1]] .= ystack[i][extraranges[i][1][end]+1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][1]] .= 0 for i in 1:length(extraranges)]
    end

  elseif kw.interpolate == "bfill"

    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= ystack[i][j[end]+1]
    end
    if !(kw.extrapolate == false || kw.extrapolate == "None")
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][1]] .= ystack[i][extraranges[i][1][end]+1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][1]] .= 0 for i in 1:length(extraranges)]
    end
    if kw.extrapolate == "both" || kw.extrapolate == "nearest"
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][end]] .= ystack[i][extraranges[i][end][1]-1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][end]] .= 0 for i in 1:length(extraranges)]
    end

  elseif kw.interpolate == "mean"

    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= (ystack[i][j[1]-1] .+ ystack[i][j[end]+1]) ./ 2
    end
    if kw.extrapolate == true || kw.extrapolate == "nearest"
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][1]] .= ystack[i][extraranges[i][1][end]+1]
        catch; end
      end
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][end]] .= ystack[i][extraranges[i][end][1]-1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][1]] .= 0 for i in 1:length(extraranges)]
      [ystack[i][extraranges[i][end]] .= 0 for i in 1:length(extraranges)]
    end

  elseif kw.interpolate == "spline"

    xspl = [xdata[isfinite.(ystack[i])] for i = 1:length(ystack)]
    yspl = [ystack[i][isfinite.(ystack[i])] for i = 1:length(ystack)]
    if kw.extrapolate == false
      spline = [spl.Spline1D(xspl[i], yspl[i], bc="zero", k=kw.kspline)
                for i = 1:length(ystack)]
    elseif kw.extrapolate == true
      spline = [spl.Spline1D(xspl[i], yspl[i], bc="extrapolate", k=kw.kspline)
                for i = 1:length(ystack)]
    else
      spline = [spl.Spline1D(xspl[i], yspl[i], bc=extrapolate, k=kw.kspline)
                for i = 1:length(ystack)]
    end
    ystack = [spline[i](xdata) for i = 1:length(ystack)]

  else

    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= kw.interpolate
    end
    if kw.extrapolate == true
      [ystack[i][extraranges[i][1]] .= kw.interpolate for i in 1:length(extraranges)]
      [ystack[i][extraranges[i][end]] .= kw.interpolate for i in 1:length(extraranges)]
    else
      [ystack[i][extraranges[i][1]] .= 0 for i in 1:length(extraranges)]
      [ystack[i][extraranges[i][end]] .= 0 for i in 1:length(extraranges)]
    end
  end
  ylines = cumsum(ystack)

  return ystack, ylines
end #function interpolate_stack


function index2range(ydata)
  idx = findall(isnan.(ydata))
  global ptr = 0
  index = []
  while ptr < length(idx)
    global ptr += 1
    global istart = idx[ptr]
    global iend = idx[ptr] + 1
    while ptr < length(idx) && iend == idx[ptr+1]
      if ptr == length(idx)  break  end
      global iend += 1
      global ptr += 1
    end
    push!(index, istart:iend-1)
  end
  if isempty(index)
    firstextra, lastextra = 1:0, 1:0
  else
    firstextra = index[1][1] == 1 ? popfirst!(index) : 1:0
    lastextra = index[end][end] == length(ydata) ? pop!(index) : 1:0
  end

  return index, (firstextra, lastextra)
end #function index2range


"""


"""
function def_aliases(kw_aliases)
  # Define all sets of keyword argument aliases
  ti_aliases = (:ti, :title)
  lt_aliases = (:lt, :ls, :dt, :linetype, :linestyle, :line_type, :line_style,
    :dashes, :dashtype)
  pt_aliases = (:pt, :mt, :pointtype, :point_type, :marker_type, :marker_type, :marker)
  lc_aliases = (:lc, :mc, :linecolour, :markercolour, :line_colour, :marker_colour,
    :linecolor, :markercolor, :line_color, :marker_color, :colour, :color)
  lw_aliases = (:lw, :linewidth, :line_width)
  cs_aliases = (:cs, :colourscheme, :colour_scheme, :colorscheme, :color_scheme)
  plt_aliases = (:plottype, :plot_type)
  loc_aliases = (:legpos, :loc, :legloc)
  lcol_aliases = (:legcols, :leg_cols, :legcolumns, :leg_columns,
    :legendcols, :legend_cols, :legendcolumns, :legend_columns)
  xlim_aliases = (:xlim, :xlims)
  ylim_aliases = (:ylim, :ylims)
  mticks_aliases = (:mticks, :minor_ticks)
  tioff_aliases = (:ti_offset, :title_offset)
  lbloff_aliases = (:lbl_offset, :lab_offset, :label_offset,
    :axlbl_offset, :axlab_offset, :axlabel_offset,
    :ax_lbl_offset, :ax_lab_offset, :ax_label_offset)
  legoff_aliases = (:leg_offset, :legend_offset)
  axoff_aliases = (:tick_offset, :ax_offset, :axes_offset)
  axclr_aliases = (:axcolour, :ax_colour, :axcolor, :ax_color)
  mxt_aliases = (:xticks, :mxticks, :min_xticks, :minor_xticks)
  Mxt_aliases = (:Xticks, :Mxticks, :maj_xticks, :major_xticks)
  myt_aliases = (:yticks, :myticks, :min_yticks, :minor_yticks)
  Myt_aliases = (:Yticks, :Myticks, :maj_yticks, :major_yticks)
  α_aliases = (:alpha, :α)
  fig_aliases = (:figsize, :figuresize)

  # Combine aliases in an overall tuple
  aliases = [ti_aliases, lt_aliases, pt_aliases, lc_aliases, lw_aliases, cs_aliases,
    plt_aliases, loc_aliases, lcol_aliases, xlim_aliases, ylim_aliases, mticks_aliases,
    tioff_aliases, lbloff_aliases, legoff_aliases, axoff_aliases, axclr_aliases,
    mxt_aliases, Mxt_aliases, myt_aliases, Myt_aliases, α_aliases, fig_aliases,
    # Add remaining kwargs without aliases
    [:xlabel], [:ylabel], [:timeformat], [:timescale], [:major_interval], [:minor_interval],
    [:logscale], [:logremove], [:border], [:interpolate], [:extrapolate], [:kspline],
    [:twinax], [:fontsize], [:framewidth], [:ticksize], [:cap_offset]]

  # Find keywords used in the argument list
  kw_args = Dict(kw_aliases)
  input_kwargs = [intersect(alias, keys(kw_args)) for alias in aliases]
  # Define default parameters of type kwargs
  refined_kwargs = Dict()

  # Overwrite default values with values of argument list
  # and warn of duplicate definitions
  for (i, kw) in enumerate(input_kwargs)
    if length(kw) ≥ 1
      refined_kwargs[aliases[i][1]] = kw_args[kw[1]]
      if length(kw) > 1
        @warn "multiple aliases defined for $(kw[1]); others ignored", kw
      end
    end
  end

  # Return final kwargs
  return refined_kwargs
end #function def_aliases


function adjust_kwargs(kw::kwargs, pltdata)
  kw.ylabel = [kw.ylabel]
  if kw.xlim == nothing  kw.xlim = (nothing, nothing)  end
  kw.ylim == nothing ? kw.ylim = [(nothing, nothing)] : kw.ylim = [kw.ylim]
  kw.axcolour = ["black"]
  kw.Yticks = [kw.Yticks]
  kw.yticks = [kw.yticks]
  kw.logscale = [kw.logscale]
  kw.lt isa String || kw.lt isa Tuple{Real,Real} || typeof(kw.lt) <: Vector{<:Real} ?
    kw.lt = [[kw.lt for l in pltdata]] : kw.lt = [kw.lt]
  kw.lc isa String || kw.lc isa Symbol ?
    kw.lc = [[kw.lc for l in pltdata]] : kw.lc = [kw.lc]
  kw.pt = [["default" for p in pltdata]]
  if kw.cs isa String  kw.cs = [kw.cs]  end
  if kw.alpha == 0
    α = [a.alpha for a in pltdata]
    kw.alpha = stats.mean(α) > 0 ? stats.mean(α) : 1
  end

  for (i, p) in enumerate(pltdata)
    if p.colour == nothing
      p.colour = sel_ls(cs="pyplot", lc=i)[1]
    end
  end
  return kw
end #function adjust_kwargs


function create_PlotData_with_errors(plotdata, err, SF, select_cols)
  # Make copy of plotdata that can be altered
  pltdata = deepcopy(plotdata)
  # (Re-)define column names of DataFrame
  if isempty(select_cols)
    DFnames = Symbol[:x, :y, :y_lower, :y_upper, :x_lower, :x_upper]
    pltdata = renameDF(pltdata, DFnames, err)
  else
    if length(select_cols) ≠ 6
      throw(ArgumentError(string("for `select_cols`.\n",
        "Define all column names for `x`, `y`, `y_lower`, `y_upper`, ",
        "`x_lower`, and `x_upper`.")))
    end
    DFnames = deepcopy(select_cols)
  end

  # Calculate error columns depending on choice of `err`
  pltdata = calc_errors(pltdata, DFnames, err)

  # Scale data
  for s in DFnames[2:4]
    try pltdata[s] .*= SF
    catch; continue
    end
  end

  # Save errors
  errors = []
  if err ≠ "None"
    for s in DFnames[3:end]
      try push!(errors,pltdata[s])
      catch; push!(errors,Real[])
      end
    end
  else
    errors = [Real[], Real[], Real[], Real[]]
  end

  # Return PlotData type
  return PlotData(x = pltdata[DFnames[1]], y = pltdata[DFnames[2]],
    x_upper = errors[4], x_lower = errors[3], y_upper = errors[2], y_lower = errors[1])
end #function create_PlotData_with_errors


function def_PlotDataFormats(pltdata, kw::kwargs, kwdict::Dict, label, alpha)
  # Loop over fields as used in kwargs and PlotData
  for (k, p) in zip([:pt, :lt, :lc, :lw], [:marker, :dashes, :colour, :lw])
    if k in keys(kwdict)
      # Set fields in PlotData if defined in kwargs
      setfield!(pltdata, p, getfield(kw, k))
    end
  end
  pltdata.label = label
  pltdata.alpha = alpha
  return pltdata
end #function def_PlotDataFormats

"""


"""
function def_kwargs(kw::Dict; calledby=:other)
  # Refine type checks before instantiating kwargs with simple type tests
  if haskey(kw, :lt)  kw[:lt] = vectype(kw[:lt])  end
  if calledby ≠ :load_PlotData && haskey(kw, :lc) && kw[:lc] == nothing
    throw(ErrorException(
      "`nothing` in `lc` is only allowed in function `load_PlotData`"))
  end
  if calledby == :sel_ls
    if haskey(kw, :lt) && !(typeof(kw[:lt]) <: Int ||
      typeof(kw[:lt]) <: Vector{<:Int} || kw[:lt] isa UnitRange{Int})
      throw(ErrorException(
        "only `Int`, `Vector{<:Int}` or `UnitRange{Int}` allowed for `lt` in function `sel_ls`"))
    end
    if haskey(kw, :pt) && !(typeof(kw[:pt]) <: Int ||
      typeof(kw[:pt]) <: Vector{<:Int} || kw[:pt] isa UnitRange{Int})
      throw(ErrorException(
        "only `Int`, `Vector{<:Int}` or `UnitRange{Int}` allowed for `pt` in function `sel_ls`"))
    end
    if haskey(kw, :lc) && !(typeof(kw[:lc]) <: Int ||
      typeof(kw[:lc]) <: Vector{<:Int} || kw[:lc] isa UnitRange{Int})
      throw(ErrorException(
        "only `Int`, `Vector{<:Int}` or `UnitRange{Int}` allowed for `lc` in function `sel_ls`"))
    end
  else
    if haskey(kw, :lt) && (kw[:lt] isa UnitRange || typeof(kw[:lt]) <: Int)
      throw(ErrorException(
        "`UnitRange` or `Vector{<:Int}` for `lt` is only allowed in function `sel_ls`"))
    end
    if haskey(kw, :pt) && kw[:pt] isa UnitRange
      throw(ErrorException("`UnitRange` for `pt` is only allowed in function `sel_ls`"))
    end
    if haskey(kw, :lc) && (kw[:lc] isa UnitRange ||
      typeof(kw[:lc]) <: Vector{<:Int} || typeof(kw[:lc]) <: Int)
      throw(ErrorException(
        "`UnitRange`, `Vector{<:Int}` or `Int` for `lc` is only allowed in function `sel_ls`"))
    end
  end
  # Initialise default values that deviate from defaults in kwargs
  if calledby == :load_PlotData
    if !haskey(kw, :lc)  kw[:lc] = nothing  end
    if !haskey(kw, :lt)  kw[:lt] = Real[]  end
    if !haskey(kw, :pt)  kw[:pt] = "None"  end
  elseif calledby == :sel_ls
    if !haskey(kw, :cs)  kw[:cs] = "default"  end
    if !haskey(kw, :lc)  kw[:lc] = 1  end
    if !haskey(kw, :lt)  kw[:lt] = 1  end
    if !haskey(kw, :pt)  kw[:pt] = 0  end
  end

  # Instantiate and return kwargs
  kw_args = kwargs()
  for (key, val) in zip(keys(kw), values(kw))
    setfield!(kw_args, key, val)
  end
  kwargscheck(kw_args.ylabel, kw_args.lt, kw_args.pt, kw_args.lc, kw_args.cs,
    kw_args.ylim, kw_args.Yticks, kw_args.yticks, kw_args.logscale,
    kw_args.logremove, kw_args.axcolour, kw_args.twinax)

  return kw_args
end


function kwargscheck(ylabel::Union{AbstractString, Vector{<:AbstractString}},
  lt::lt_type, pt::Union{String,Int,UnitRange{Int},Vector}, lc::lc_type,
  cs::Union{String,Vector{<:String}}, ylim, Yticks::Union{Real,Vector{<:Real}},
  yticks::Union{Real,Vector{<:Real}}, logscale::Union{String,Vector{<:String}},
  logremove::Union{String,Vector{<:String}},
  axcolour::Union{String,Symbol,Vector}, twinax::Vector{<:Int})

  if ylabel isa Vector && (isempty(twinax) || length(ylabel) ≠ 2)
    throw(DefinitionError(length(ylabel)))
  end
  if cs isa Vector && (isempty(twinax) || length(cs) ≠ 2)
    throw(DefinitionError(length(cs)))
  end
  if ylim isa Vector
    if isempty(twinax) || length(ylim) ≠ 2
      throw(DefinitionError(length(ylim)))
    elseif !all([types <: limtype for types in typeof.(ylim)])
      throw(DefinitionError(
        "`ylim` must be a vector with `limtype` objects", typeof(ylim)))
    end
  end
  if Yticks isa Vector && (isempty(twinax) || length(Yticks) ≠ 2)
    throw(DefinitionError(length(Yticks)))
  end
  if yticks isa Vector && (isempty(twinax) || length(yticks) ≠ 2)
    throw(DefinitionError(length(yticks)))
  end
  if logscale isa Vector && (isempty(twinax) || length(logscale) ≠ 2)
    throw(DefinitionError(length(logscale)))
  end
  if logremove isa Vector && (isempty(twinax) || length(logremove) ≠ 2)
    throw(DefinitionError(length(logremove)))
  end
  if axcolour isa Vector
    if isempty(twinax) || length(axcolour) ≠ 2
      throw(DefinitionError(length(axcolour)))
    elseif !all([types in DataType[String, Symbol] for types in typeof.(axcolour)])
      throw(DefinitionError(
        "`axcolour` must be a vector with `String` and/or `Symbol` objects",
        typeof(axcolour)))
    end
  end
end #function kwargscheck

function vectorcheck(kw, len)
  kw.lt = vectype(kw.lt)
  if typeof(kw.lt)<:Vector{<:Real}
    if isodd(length(kw.lt)) || 0 < length(kw.lt) < 4
      throw(DefinitionError(
        "`lt` must be a vector with an even number of at least 4 elements; length(lt) =",
        length(kw.lt)))
    end
  elseif kw.lt isa Vector
    for (i, el) in enumerate(kw.lt)
      if !(typeof(el) <: Vector{<:Real} || el isa Tuple{Real,Real})
        throw(DefinitionError(
          string("`lt` must be a vector of `Vector{<:Real}` or `Tuple{Real,Real}` elements; ",
            "DataType of $i. element = "), typeof(kw.lt)))
      elseif el isa Vector && (isodd(length(el)) || 0 < length(el) < 4)
        throw(DefinitionError(
          string("`lt` must be a vector with an even number of at least 4 elements; ",
            "length(lt) in $i. element = "), length(el)))
      end
    end
    if length(kw.lt) ≠ len
      throw(DefinitionError("length of `lt` and `PlotData` must be the same; length(lt) = ",
        length(kw.lt)))
    end
  end
  if kw.pt isa Vector
    if length(kw.pt) ≠ len
      throw(DefinitionError("length of `pt` and `PlotData` must be the same; length(pt) = ",
        length(kw.pt)))
    elseif !all([types in DataType[String, Int] for types in typeof.(kw.pt)])
      throw(DefinitionError(
        "`pt` must be a vector with `String` and/or `Int` objects", typeof(kw.pt)))
    end
  end
  if kw.lc isa Vector
    if length(kw.lc) ≠ len
      throw(DefinitionError("length of `lc` and `PlotData` must be the same; length(lc) = ",
        length(kw.lc)))
    elseif !all([types in DataType[String, Symbol] for types in typeof.(kw.lc)])
      throw(DefinitionError(
        "`lc` must be a vector with `String` and/or `Symbol` objects", typeof(kw.lc)))
    end
  end
end


vectype(vec, dt::DataType=Real) = if vec isa Vector
  isempty(vec) ? Real[] :
    try [dt.(v) for v in vec]
    catch; vec
    end
  else
    return vec
end
