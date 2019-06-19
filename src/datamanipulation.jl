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
    if pltdata[i].marker ≠ "None" && getfield(pltdata[i],high) ≠ nothing
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
    if getfield(p[i], Symbol("$(x)lerr")) ≠ nothing
      rl = findall(rel.(getfield(p[i], Symbol("$(x)lerr")), 0.))
      pl = getfield(p[i], Symbol("$(x)lerr")); pl[rl] .= 0.
    end
    if getfield(p[i], Symbol("$(x)uerr")) ≠ nothing
      ru = findall(rel.(getfield(p[i], Symbol("$(x)uerr")), 0.))
      pu = getfield(p[i], Symbol("$(x)uerr")); pu[ru] .= 0.
    end
  end

  return p
end


function interpolate_stack(xdata, ystack, interpolate, extrapolate, kspline)
  # imiss = [findall(isnan.(ys)) for ys in ystack]
  # idata = [findall(isfinite.(ys)) for ys in ystack]
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
  firstextra = index[1][1] == 1 ? popfirst!(index) : 1:0
  lastextra = index[end][end] == length(ydata) ? pop!(index) : 1:0

  return index, (firstextra, lastextra)
end #function index2range


"""
    checkaliases(args)

Checks that aliases are used correctly without accidental multiple definitions
for the kw`args` of function `plot_data`.
"""
function checkaliases(ti, title, cs, colorscheme, colourscheme,
  lt, linestyle, linetype, dashes, pt, mt, marker, lc, linecolor, linecolour, color, colour,
  legpos, legloc, loc)

  ti = checkalias(kwargs = [ti, title], alias = ["ti", "title"], default = "")
  cs = checkalias(kwargs = [cs, colorscheme, colourscheme],
    alias = ["cs", "colorscheme", "colourscheme"], default = "")
  lt = checkalias(kwargs = [lt, linestyle, linetype, dashes],
    alias = ["lt", "linestyle", "linetype", "dashes"], default = "default")
  pt = checkalias(kwargs = [pt, mt, marker], alias = ["pt", "mt", "marker"],
    default = "default")
  lc = checkalias(kwargs = [lc, linecolor, linecolour, color, colour],
    alias = ["lc", "linecolor", "linecolour", "color", "colour"], default = "default")
  legpos = checkalias(kwargs = [legpos, legloc, loc],
    alias = ["legpos", "legloc", "loc"], default = "best")

  return ti, cs, lt, pt, lc, legpos
end

"""
    checkalias(;kwargs, alias=[""], default=[""])

Checks that `kwargs` of a function are used correctly without accidental
multiple definitions different from the `default` value. Otherwise warns about
the multiple definitions giving naming the `alias`es with multiple definitions
and chosing the value of the first kwarg for use in the function.
"""
function checkalias(;kwargs, alias=[""], default=[""])
  setkw = findall(kwargs.≠default)
  kw = kwargs[setkw]
  al = alias[setkw]
  if length(unique(kw)) > 1
    warning = string("\033[93mMultiple definitions of keyword aliases ",
      join(al, ", ", " and "), "!\n\33[0m`$(al[1]) = \"$(kw[1])\"` used.")
    @warn(warning)
    kwargs[1] = kw[1]
  elseif length(unique(kw)) == 1
    kwargs[1] = kw[1]
  end

  return kwargs[1]
end
