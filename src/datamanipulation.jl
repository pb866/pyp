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
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
  elseif err == "pmpercentx"
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]].⋅pltdata[col[1]]
  elseif err == "percenty"
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
  elseif err == "pmpercenty"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]].⋅pltdata[col[2]]
  elseif err == "percent"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
  elseif err == "pmpercent"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]].⋅pltdata[col[2]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]].⋅pltdata[col[1]]

  # error × value (adding a factor error range)
  elseif err == "factorx"
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
  elseif err == "pmfactorx"
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[6]]
  elseif err == "factory"
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
  elseif err == "pmfactory"
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[4]]
  elseif err == "factor"
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
  elseif err == "pmfactor"
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[4]]
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[6]]
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
function setup_log(pltdata, logremove::Union{String, Vector{String}},
  logscale::Union{String, Vector{String}})

  for i = 1:length(pltdata)
    if logremove[i] == "pos" && occursin("x", logscale[i])
      pltdata[i] = rm_log(pltdata[i], "x", >)
    end
    if logremove[i] == "pos" && occursin("y", logscale[i])
      pltdata[i] = rm_log(pltdata[i], "y", >)
    end
    if logremove[i] == "neg" && occursin("x", logscale[i])
      pltdata[i] = rm_log(pltdata[i], "x", <)
    end
    if logremove[i] == "neg" && occursin("y", logscale[i])
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


function interpolate_stack(xdata, ystack, interpolate, extrapolate, kspline)
  ranges, extraranges = zip(index2range.(ystack)...)
  if interpolate == "linreg" || interpolate == "linear"
    interpolate = "spline"
    kspline = 1
  end
  if interpolate == "ffill"
    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= ystack[i][j[1]-1]
    end
    if !(extrapolate == false || extrapolate == "None")
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][end]] .= ystack[i][extraranges[i][end][1]-1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][end]] .= 0 for i in 1:length(extraranges)]
    end
    if extrapolate == "both" || extrapolate == "nearest"
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][1]] .= ystack[i][extraranges[i][1][end]+1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][1]] .= 0 for i in 1:length(extraranges)]
    end

  elseif interpolate == "bfill"

    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= ystack[i][j[end]+1]
    end
    if !(extrapolate == false || extrapolate == "None")
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][1]] .= ystack[i][extraranges[i][1][end]+1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][1]] .= 0 for i in 1:length(extraranges)]
    end
    if extrapolate == "both" || extrapolate == "nearest"
      for i in 1:length(extraranges)
        try ystack[i][extraranges[i][end]] .= ystack[i][extraranges[i][end][1]-1]
        catch; end
      end
    else
      [ystack[i][extraranges[i][end]] .= 0 for i in 1:length(extraranges)]
    end

  elseif interpolate == "mean"

    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= (ystack[i][j[1]-1] .+ ystack[i][j[end]+1]) ./ 2
    end
    if extrapolate == true || extrapolate == "nearest"
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

  elseif interpolate == "spline"

    xspl = [xdata[isfinite.(ystack[i])] for i = 1:length(ystack)]
    yspl = [ystack[i][isfinite.(ystack[i])] for i = 1:length(ystack)]
    if extrapolate == false
      spline = [spl.Spline1D(xspl[i], yspl[i], bc="zero", k=kspline)
                for i = 1:length(ystack)]
    elseif extrapolate == true
      spline = [spl.Spline1D(xspl[i], yspl[i], bc="extrapolate", k=kspline)
                for i = 1:length(ystack)]
    else
      spline = [spl.Spline1D(xspl[i], yspl[i], bc=extrapolate, k=kspline)
                for i = 1:length(ystack)]
    end
    ystack = [spline[i](xdata) for i = 1:length(ystack)]

  else

    for i in 1:length(ranges), j in ranges[i]
      ystack[i][j] .= interpolate
    end
    if extrapolate == true
      [ystack[i][extraranges[i][1]] .= interpolate for i in 1:length(extraranges)]
      [ystack[i][extraranges[i][end]] .= interpolate for i in 1:length(extraranges)]
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
function def_aliases(kw_aliases...)
  # Define all sets of keyword argument aliases
  ti_aliases = (:ti, :title)
  lt_aliases = (:lt, :ls, :dt, :linetype, :linestyle, :line_type, :line_style,
    :dashes, :dashtype)
  pt_aliases = (:pt, :mt, :pointtype, :point_type, :marker_type, :marker_type)
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

  # Combine aliases in an overall tuple
  aliases = [ti_aliases, lt_aliases, pt_aliases, lc_aliases, lw_aliases, cs_aliases,
    plt_aliases, loc_aliases, lcol_aliases, xlim_aliases, ylim_aliases, mticks_aliases,
    tioff_aliases, lbloff_aliases, legoff_aliases, axoff_aliases, axclr_aliases]

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


function adjust_kwargs(kw, pltdata)
  if kw.xlim == nothing  kw.xlim = (nothing, nothing)  end
  kw.ylim == nothing ? kw.ylim = [(nothing, nothing)] : kw.ylim = [kw.ylim]
  kw.axcolour = ["black"]

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


function def_PlotDataFormats(pltdata, kw, kwdict, label, alpha)
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
function def_kwargs(kw; calledby=:other)
  # Refine type checks before instantiating kwargs with simple type tests
  if calledby ≠ :load_PlotData && haskey(kw, :lc) && kw[:lc] == nothing
    throw(ErrorException(
      "`nothing` in `lc` is only allowed in function `load_PlotData`"))
  end
  if calledby == :sel_ls
    if haskey(kw, :lt) && !(typeof(kw[:lt]) <: Int ||
      typeof(kw[:lt]) <: Vector{<:Int} || typeof(kw[:lt]) isa UnitRange{Int})
      throw(ErrorException(
        "only `Int`, `Vector{<:Int}` or `UnitRange{Int}` allowed for `lt` in function `sel_ls`"))
    end
    if haskey(kw, :pt) && !(typeof(kw[:pt]) <: Int ||
      typeof(kw[:pt]) <: Vector{<:Int} || typeof(kw[:pt]) isa UnitRange{Int})
      throw(ErrorException(
        "only `Int`, `Vector{<:Int}` or `UnitRange{Int}` allowed for `pt` in function `sel_ls`"))
    end
    if haskey(kw, :lc) && !(typeof(kw[:lc]) <: Int ||
      typeof(kw[:lc]) <: Vector{<:Int} || typeof(kw[:lc]) isa UnitRange{Int})
      throw(ErrorException(
        "only `Int`, `Vector{<:Int}` or `UnitRange{Int}` allowed for `lc` in function `sel_ls`"))
    end
  else
    if haskey(kw, :lt) && (typeof(kw[:lt]) isa UnitRange || typeof(kw[:lt]) <: Int)
      throw(ErrorException(
        "`UnitRange` or `Vector{<:Int}` for `lt` is only allowed in function `sel_ls`"))
    end
    if haskey(kw, :pt) && typeof(kw[:pt]) isa UnitRange
      throw(ErrorException("`UnitRange` for `pt` is only allowed in function `sel_ls`"))
    end
    if haskey(kw, :lc) && (typeof(kw[:lc]) isa UnitRange ||
      typeof(kw[:lc]) <: Vector{<:Int} || typeof(kw[:lc]) <: Int)
      throw(ErrorException(
        "`UnitRange`, `Vector{<:Int}` or `Int` for `lc` is only allowed in function `sel_ls`"))
    end
  end
  # Initialise default values that deviate from defaults in kwargs
  if calledby == :load_PlotData
    if !haskey(kw, :lc)  kw[:lc] = nothing  end
    if !haskey(kw, :lt)  kw[:lt] = []  end
    if !haskey(kw, :pt)  kw[:pt] = "None"  end
  elseif calledby == :sel_ls
    if !haskey(kw, :cs)  kw[:cs] = "default"  end
    if !haskey(kw, :lc)  kw[:lc] = 1  end
    if !haskey(kw, :lt)  kw[:lt] = 1  end
    if !haskey(kw, :pt)  kw[:pt] = 0  end
  end

  # Instantiate and return kwargs
  kw_args = kwargs()
  for key in keys(kw)
    setfield!(kw_args, key, kw[key])
  end

  return kw_args
end
