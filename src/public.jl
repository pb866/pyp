

##########################
###  PUBLIC FUNCTIONS  ###
##########################


"""
    load_PlotData(pltdata::DataFrames.DataFrame;  kwargs...) -> PlotData

Pack x and y data with possible assoiciated errors from DataFrame `pltdata`
as well as formatting parameters into a new DataType `PlotData` for `PyPlot` plotting.


### kwargs...

+ `err` (`String`):
  - `"None"` (no errors, **default**)
  - `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
  - `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
  - `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
  - `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
  - `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
  - `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
  - `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
+ `pt` (`Union{AbstractString,Int64}`): PyPlot marker type, see e.g. https://matplotlib.org/api/markers_api.html
  - `"None"` (**default**) to suppress use of markers
+ `lt` (`Union{String,Tuple{Int64,Int64},Array{Int64,1}}`), tuple or array with
  on/off PyPlot dash definitions
  - empty array means a solid line (**default**)
  - Tuple with 2 entries or array with 4 entries where odd entries are `0` suppresses lines
+ `lc` (`Union{Tuple{Real,Real},Vector{T}} where T=Real`) PyPlot line/marker colours
  (**default:** `nothing` – use PyPlot default colours), see e.g.
  https://matplotlib.org/examples/color/named_colors.html
+ `lw` (`Real`): linewidth (**default:** `1.4`)
+ `SF` (`Real`): scaling factor of y data and associated errors (**default:** `1`, no scaling)
+ `label` (`AbstractString`): Label for legend
  - `""` (empty string, **default**): no legend label
+ `select_cols` (`Vector{Symbol}`):
  - `Symbol[]` (**default**) (assume columns in the order: `x`, `y`, `ylerr`, `yuerr`, `xlerr`, `xuerr`)
  - `Vector{Symbol}`: give array with column names for `x`, `y`, `ylerr`, `yuerr`, `xlerr`, `xuerr`
   (give complete list, even if columns are incomplete due to the choice of `err`)
"""
function load_PlotData(plotdata::DataFrame; err::String="None",
         pt::Union{AbstractString,Int64}="None",
         lt::Union{Tuple{Real,Real},Vector{T} where T<:Real} = Int64[],
         lc::Union{Nothing,String,Symbol}=nothing, lw::Real=1.4, SF::Real=1,
         label::AbstractString="", alpha::Real=1, select_cols::Vector{Symbol}=Symbol[])

  # Make copy of plotdata that can be altered
  pltdata = deepcopy(plotdata)
  # (Re-)define column names of DataFrame
  if isempty(select_cols)
    DFnames = Symbol[:x, :y, :ylerr, :yuerr, :xlerr, :xuerr]
    pltdata = renameDF(pltdata, DFnames, err)
  else
    if length(select_cols) ≠ 6
      throw(ArgumentError(string("for `select_cols`.\n",
        "Define all column names for `x`, `y`, `ylerr`, `yuerr`, `xlerr`, and `xuerr`.")))
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
      catch; push!(errors,nothing)
      end
    end
  else
    errors = [nothing, nothing, nothing, nothing]
  end

  # Return PlotData type
  return PlotData(x = pltdata[DFnames[1]], y = pltdata[DFnames[2]],
         xuerr = errors[4], xlerr = errors[3], yuerr = errors[2], ylerr = errors[1],
         label = label, marker = pt, dashes = lt, colour = lc, lw = lw, alpha=alpha)
end #function load_PlotData


"""
    plot_data(plot_list::PlotData...; kwargs...) -> fig, ax1[, ax2]

Generate plots with scatter and/or line data from the varargs `plot_data` of Type `PlotData`
and return `PyPlot` figure and axes data for further plot formatting.


### kwargs...

For all arguments concerning y data the following applies:
- The dataset can be split into 2 subsets and a 2. y-axis with a different scale
  can be assigned (see optional parameter `twinax`)
- If a second axis is defined, parameters concerning y-data can be compiled in an
  array of length `2` with different values for the first and second y-axis,
  respectively
- If the second y-axis is defined and only a single parameter concerning y-data
  is defined, then this parameter applies to both axes

Keyword arguments for `function plot_data` are:
+ `xlabel` (`AbstractString`): x axis label
  - **default:** `"model time / hours"`
+ `ylabel` (`Union{AbstractString, Vector{T}} where T=AbstractString`):
  y axis label, `Array` can be used, if a second axis with a different label is introduced
  - **default:** `"concentration / mlc cm\$^{-3}\$ s\$^{-1}\$"`
+ `ti` (`Union{String, LaTeXString}`): Plot title
  - **default:** `""` (empty string) for no title
  - **keyword aliases:** `title`
+ `twinax` (`Array{Int64,1}`): optional array to devide data into 2 subsets and assign second dataset to a 2. y-axis; use array of `length` of the PlotData with integers `1` and `2` to assign each data to axis 1 and 2
  (**default:**: empty array `Int64[]`, no 2. y-axis)
+ `logscale` (`Union{String, Array{String, 1}}`): use `"x"`, `"y"` or `"xy"` to
  set either x-, y- or both axes to a logarithmic scale;
  if 2. y-axis is present, arrays may be used to have different scales on each y-axis;
  (**default:** `""` – linear scale)
+ `logremove` (`Union{String, Vector{String}}`; **default:** `"neg"`): removes positve (`"pos"`) or negative (`"neg"`) values in the plotting data for logarithmic plots; leaving data unmanipulated with `"none"` might result in PyPlot failure; use a vector of strings for different treatment of the second axis.
+ `plot_type` (`String`): Choose between the following plot types:
  - `"line"`: line plot (solid lines or as defined by `lt`)
  - `"scatter"`: scatter plot (squares or as defined by `pt`)
  - `"both"`: markers connected by lines (squares connected by solid lines or as defined by `lt` and/or `pt`)
  -  `"default"`: Use formats as predefined by the PlotData or defined by `kwargs` `ls`, `pt`, and `lc`/`cs`
+ `cs` (`Union{String, Array{String,1}}`): Define colour scheme `default`, `source`
  or `sink` from `function sel_ls`; array with 2 different colour schemes may be
  used in case of 2. y-axis
  - `""` (**default:**): colours as defined in the `PlotData` or default PyPlot colour scheme)
  - **keyword aliases:** `colorscheme`, `colourscheme`
+ `lt`: dash type (defined by PyPlot's `dashes`)
  - overwrites other colour scheme settings
  - Array of `length 1 or 2` with Arrays of `length(plot_data)` where inner arrays specify line type of each `PlotData`
  - use tuple or array of integers to define on/off dash pixels
  - `[]` (empty array, **default**): solid lines
  - tuple with 2 entries or array with at least 4 entries where odd entries are `0` for no lines
  - **keyword aliases:** `linestyle`, `linetype`, `dashes`
+ `lc`: line/marker colour (defined by PyPlot's `color` keywords, see e.g. https://matplotlib.org/examples/color/named_colors.html)
  - for use only, when `plot_type` is switched to `"mix"`
  - Array of `length 1 or 2` with Arrays of `length(plot_data)` where inner arrays specify the line/marker colour of each `PlotData`
  - **keyword aliases:** `linecolor`, `linecolour`, `color`, `colour`
+ `pt`: marker type (defined by PyPlot's `marker`)
  - overwrites other colour scheme settings
  - String or integer keyword for markers from PyPlot
  - **keyword aliases:** `mt`, `marker`
+ `alpha`: if alpha is set to a value `> 0` then all graphs will be displayed with this transparency value
+ `mticks` (`Bool`): Switch minor ticks on/off (**default:** `true`)
+ `minor_xticks` (`Union{Real,Vector{Int64},Vector{Float64}}`): Size of interval between minor x ticks (**default:** `0` – automatic determination by PyPlot)
+ `minor_yticks` (`Union{Int64, Array{Int64,1}}`):Size of interval between minor y ticks
  - `0` (**default**): automatic determination by PyPlot
  - If 2. y axis is defined, an array can be used to used different specifications for each axis
+ `major_xticks`/`major_yticks` (`Union{Real,Vector{Int64},Vector{Float64}}`/`Real`)
  interval size of major ticks in x- and y-axis, respectively; asign different intervals
  to 2nd y-axis using an array of integers
+ `xlims`/`ylims`: Tuple of minimum/maximum value for each axis
  - `nothing` may be used for default PyPlot values (**default**)
  - `nothing` may be used for one value within the tuple to use the default, e.g. `ylims = (0, nothing)`
  - if 2. y-axis is defined, `ylims` may be an array of tuples with different specifications
+ `figsiz` (`Tuple{Real,Real}`): figure size width × height in inches
  (**default:** `(6,4)`)
+ `fntsiz` (`Real`): default font size (**default:** `12`)
+ `framewidth` (`Real`): default line width of outer axes frame (**default:** `1`)
+ `ticksize` (`Tuple{Real,Real}`): Tuple with scaling factors for length of major
  (first tuple entry) and minor (second tuple entry) ticks in relation to their
  width (**default:** `(4.5,2.5)`)
+ `cap_offset` (`Real`): deviation from default cap size of `3` of marker error bars
  (**default:** `0`)
+ `ti_offset`, `label_offset`, `leg_offset`, `ax_offset`: Offsets for fontsizes of title, axes
  labels, and legend, respectively, from default sizes. Numbers
  (positive or negative) will be added to fontsizes
  (**defaults:** `4`, `2`, `0`)
+ `axcolour` (`Union{String,Array{String,1}}`): colour of y axis label and tick numbers;
  if 2. y axis is defined, different values may be defined in an array for each axis
  (**default:** `"black"`)
+ `legpos` (`Union{String, Int64}`):
  position of the legend; choose from the following options:
  - `"best"` or `0` (default)
  - `"upper right"` or `1`
  - `"upper left"` or `2`
  - `"lower left"` or `3`
  - `"lower right"` or `4`
  - `"right"` or `5`
  - `"center left"` or `6`
  - `"center right"` or `7`
  - `"lower center"` or `8`
  - `"upper center"` or `9`
  - `"center"` or `10`
  - **keyword aliases:** `legloc`, `loc`
+ `legcolumns` (`Int64`): number of legend columns
  - **default:** `1`
"""
function plot_data(plot_list::PlotData...;
                  xlabel::AbstractString="model time / hours",
                  ylabel::Union{AbstractString, Vector{T} where T<:AbstractString} =
                  "concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
                  ti::AbstractString="", twinax::Vector{Int}=Int64[],
                  logscale::Union{String,Vector{T} where T<:String}="",
                  logremove::Union{String,Vector{T} where T<:String}="neg",
                  plot_type::String="default", cs::Union{String,Vector{T} where T<:String}="",
                  lt="default", pt="default",
                  lc::Union{String,Symbol, Vector{T} where T<:String, Vector{T} where T<:Symbol}="default",
                  alpha::Real=1, xlims=nothing, ylims=nothing, mticks::Bool=true,
                  minor_xticks::Union{Real,Vector{Int}}=-1,
                  major_xticks::Union{Real,Vector{Int}}=-1,
                  minor_yticks::Union{Real,Vector{T} where T<:Real}=0,
                  major_yticks::Union{Real,Vector{T} where T<:Real}=0,
                  timeformat::String="", timescale::String="days",
                  major_interval::Int=0, minor_interval::Int=0,
                  figsize::Tuple{Real,Real}=(6,4), fontsize::Real=12,
                  framewidth::Real=1, ticksize::Tuple{Real,Real}=(4.5,2.5),
                  cap_offset::Real=0, ti_offset::Real=4, label_offset::Real=2,
                  leg_offset::Real=0, ax_offset::Real=0, legcolumns::Int64=1,
                  legpos::Union{String, Int64, Tuple{Real, Real}}="best",
                  axcolour::Union{String, Vector{T} where T<:String}="black",
                  # Aliases:
                  title::AbstractString="",
                  colorscheme::Union{String, Vector{T} where T<:String}="",
                  colourscheme::Union{String, Vector{T} where T<:String}="",
                  linestyle="default", linetype="default",
                  dashes="default", mt="default", marker="default",
                  linecolor::Union{String, Vector{T} where T<:String}="default",
                  linecolour::Union{String, Vector{T} where T<:String}="default",
                  color::Union{String, Vector{T} where T<:String}="default",
                  colour::Union{String, Vector{T} where T<:String}="default",
                  legloc::Union{String, Int64}="best", loc::Union{String, Int64}="best")

  # Check kwarg aliases as set all to the same value
  ti, cs, lt, pt, lc, legpos = checkaliases(ti, title, cs, colorscheme, colourscheme,
    lt, linestyle, linetype, pt, mt, marker, dashes, lc, linecolor, linecolour,
    color, colour, legpos, legloc, loc)

  # Check for twin axes and devide datasets
  pltdata, fig, ax, ylabel, logscale, logremove, xlimit, ylimit,
    major_yticks, minor_yticks, cs, lc, lt, pt, axcolour =
    setup_plot(plot_list, figsize, twinax, ylabel, logscale, logremove, xlims, ylims,
    major_yticks, minor_yticks, plot_type, cs, lc, lt, pt, alpha, axcolour)

  # Ensure strictly positive or negative values for log plots
  pltdata = setup_log(pltdata, logremove, logscale)

  # set colour scheme and line/marker types
  pltdata = set_style(pltdata, plot_type, cs, lc, lt, pt)

  # Plot data and associated errors
  for i = 1:length(pltdata)
    pltdata[i], ax[i] = plt_DataWithErrors(pltdata[i], ax[i], cap_offset)
  end

  # Define logscales and set axes limits
  ax = set_axes(pltdata, ax, logscale, xlimit, ylimit)


  # Format plot
  fig, ax = format_axes_and_annotations(fig, ax, pltdata, ti, xlabel, ylabel,
    timeformat, timescale, major_interval, minor_interval, xlimit, fontsize, legpos, legcolumns,
    axcolour, leg_offset, ti_offset, label_offset, ax_offset, major_xticks, major_yticks,
    minor_xticks, minor_yticks, mticks, ticksize, framewidth)

  # Add nothing to ax in case 2nd axis is missing
  if length(ax) < 2  push!(ax, nothing)  end

  # Return PyPlot data
  return fig, ax[1], ax[2]
end #function plot_data


"""
    plot_stack(plot_list::Union{PlotData,pyp.PlotData}...; \\*\\*kwargs)

Plot the `PlotData` listed in vararg `plot_list` as a stack plot, where the following
keyword arguments are possible.

## kwargs

+ `xlabel` (`Union{String, LaTeXString}`): x axis label
  - **default:** `"model time / hours"`
+ `ylabel` (`Union{String, LaTeXString, Array{String,1}, Array{LaTeXString,1}}`): y axis label, `Array` can be used, if a second axis with a different label is introduced
  - **default:** `"concentration / mlc cm\$^{-3}\$ s\$^{-1}\$"`
+ `ti` (`String`): Plot title
  - **default:** `""` (empty string) for no title
+ `border` (`Real`): transparency value for boundary lines of each stack data (`0` **default** for no border)
+ `alpha`: (`Real`): set transparency of areas of each stack data
  - `0` (**default**): Use transparency setting of the original `PlotData`
  - `1`: no transparency
+ `logscale` (`String`): use `"x"`, `"y"` or `"xy"` to set either x-, y- or both axes to a logarithmic scale
  (**default:** `""` – linear scale)
+ `cs` (`String`): Define colour scheme `default`, `source` or `sink` from `function sel_ls`
  - `""` (**default:**): colours as defined in the `PlotData`, where `nothing` is replaced by colours of the default colour scheme)
  - `"own"`: colours as define by array `lc`
+ `lt`: dash type (defined by PyPlot's `dashes`) for border
  - for use only, when `border` are switched on (with a value greater 0)
  - use tuple or array of integers/floats to define on/off dash pixels
  - `[]` (empty array, **default**): solid lines
  - tuple with 2 entries or array with at least 4 entries where odd entries are `0` for no lines
+ `lc`: stack/contour colour (defined by PyPlot's `color` keywords, see e.g. https://matplotlib.org/examples/color/named_colors.html)
  - for use only, when `cs` is switched to `"own"`
  - Array of `plot_list` with PyPlot colours or RGB code
+ `alpha`: if alpha is set to a value `> 0` then all graphs will be displayed with this transparency value
+ `xlims`/`ylims`: Tuple of minimum/maximum value for each axis
  - `nothing` may be used for default PyPlot values (**default**)
  - `nothing` may be used for one value within the tuple to use the default, e.g. `ylims = (0, nothing)`
+ `mticks` (`Bool`): Switch minor ticks on/off (**default:** true)
+ `minor_xticks` (`Union{Real,Vector{Int64},Vector{Float64}}`): Size of interval between minor x ticks
  (**default:** `0` – automatic determination by PyPlot)
+ `minor_yticks` (`Real`): Size of interval between minor y ticks
  - `0` (**default**): automatic determination by PyPlot
+ `major_xticks`/`major_yticks` (`Union{Real,Vector{Int64},Vector{Float64}}`/`Real`) interval size of major ticks in x- and y-axis, respectively;
+ `figsiz` (`Tuple{Real,Real}`): figure size width × height in inches
  (**default:** `(6,4)`)
+ `fntsiz` (`Real`): default font size (**default:** `12`)
+ `framewidth` (`Real`): default line width of outer axes frame (**default:** `1`)
+ `ticksize` (`Tuple{Real,Real}`): Tuple with scaling factors for length of major
  (first tuple entry) and minor (second tuple entry) ticks in relation to their
  width (**default:** `(4.5,2.5)`)
+ `ti_offset`, `label_offset`, `leg_offset`, `ax_offset`: Offsets for fontsizes of title, axes
  labels, and legend, respectively, from default sizes. Numbers
  (positive or negative) will be added to fontsizes
  (**defaults:** `4`, `2`, `0`)
+ `legpos` (`Union{String, Int64, Array{String,1}, Array{Int64,1}}`):
  position of the legend; choose from the following options:
  - `"best"` or `0` (default)
  - `"upper right"` or `1`
  - `"upper left"` or `2`
  - `"lower left"` or `3`
  - `"lower right"` or `4`
  - `"right"` or `5`
  - `"center left"` or `6`
  - `"center right"` or `7`
  - `"lower center"` or `8`
  - `"upper center"` or `9`
  - `"center"` or `10`
- `legcolumns` (`Union{Int64, Array{Int64,1}}`): number of legend columns
  (**default:** `1`)
"""
function plot_stack(plot_list::PlotData...;
         xlabel::AbstractString="model time / hours",
         ylabel::AbstractString="concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
         ti::AbstractString="", logscale::String="", logremove::String="neg",
         border=0, cs::String="", lt=[], lc=[], alpha::Real=1,
         xlims=nothing, ylims=nothing, mticks::Bool=true,
         minor_xticks::Union{Real,Vector{Int}} = -1, minor_yticks::Real = 0,
         major_xticks::Union{Real,Vector{Int}} = -1, major_yticks::Real = 0,
         timeformat::String="", timescale::String="days",
         major_interval::Int=0, minor_interval::Int=0,
         figsize::Tuple{Real,Real}=(6,4), fontsize::Real=12, framewidth::Real=1,
         ticksize::Tuple{Real,Real}=(4.5,2.5), ti_offset::Real=4,
         label_offset::Real=2, leg_offset::Real=0, ax_offset::Real=0, legcolumns::Int64=1,
         legpos::Union{String, Int64, Tuple{Real, Real}}="best",
         interpolate=0, extrapolate::Union{Bool,String}=false, kspline::Int64=3)

  # Get x and y data
  pltdata = deepcopy(plot_list)
  xdata, ystack, labels = get_stack(pltdata...)
  ystack, ylines = interpolate_stack(xdata, ystack, interpolate, extrapolate, kspline)

  # Format plot, set colour scheme
  colour, linetype, α = format_stack(cs, alpha, lc, lt, pltdata...)

  # Plot data as stack with optional boundary lines
  fig, ax = print_stack(xdata, ystack, ylines, border, labels, colour, linetype, α, figsize)

  # Set axis limits and log scales
  if xlims == nothing  xlims = (nothing, nothing)  end
  if ylims == nothing  ylims = (nothing, nothing)  end
  ax = set_axes([pltdata], [ax], [logscale], xlims, [ylims])

  # Format plot
  fig, ax = format_axes_and_annotations(fig, ax, [pltdata], ti, xlabel, [ylabel],
    timeformat, timescale, major_interval, minor_interval, xlims, fontsize, legpos, legcolumns,
    ["black"], leg_offset, ti_offset, label_offset, ax_offset, major_xticks, [major_yticks],
    minor_xticks, [minor_yticks], mticks, ticksize, framewidth)


  # Return PyPlot data
  return fig, ax[1]
end #function plot_stack


"""
    sel_ls(;cs::String="line",nc=1,nt=1)

Select colours with index `nc` from a colour scheme with the key word `cs`
(currently `"line"`, `"sink"`, and `"source"` available) and a line type
(from a range of solid, dashed, dash-dotted, and dotted lines) with the index `nt`.

Single colours/line types using integers or subsets using unit ranges `m:n` may be used.
"""
function sel_ls(cs::String="default";lc::Union{Int64,UnitRange{Int64},Vector{Int64}}=1,
    lt::Union{Int64,UnitRange{Int64},Vector{Int64}}=1,
    pt::Union{Int64,UnitRange{Int64},Vector{Int64}}=0)
  # Manually define sets of colour schemes
  def = String[
         "black", "red", "blue", "green", "#FF8533",
         "#D499ED", "#AC6012", "#FF00FF", "#A50029", "#808000",
         "#A0A0A0", "#FF7F7F", "#3D9EFF", "#99CC99", "#FBC96D",
         "#7F7FFF", "#DDBFA0", "#FFB2FF", "#C68696", "#D8D8B2"]
  src = String[
         "#330099", "#9900CC", "#003DF5", "#7547FF", "#339D9E",
         "#33FFCC", "#66FF33", "#009933", "#998000", "#DBA708",
         "#99CC66", "#00FF80", "#00FFFF", "#6699CC", "#006699",
         "#0033FF", "#000088", "#274E13", "#A7B763", "#48BD14"]
  snk = String[
         "#FF0000", "#FF8000", "#FFFF00", "#FF7A7A", "#FF0080",
         "#FF00FF", "#CC0099", "#CC0033", "#CC3300", "#FF470A",
         "#FFA347", "#FFCC99", "#FF9999", "#FF99CC", "#FFFF99",
         "#EEB702", "#AA8C2C", "#6D5A1C", "#B7AD8E", "#E69138"]
  # Combine colour schemes in a dataframe
  colourschemes = Dict{String,Vector{String}}("default"=>def, "source"=>src, "sink"=>snk)
  # Manually define line types
  dt = [Float64[],Float64[3,2,3,2],Float64[4,1.6,1,1.6,1,1.6],
        Float64[6,2,1,1,1,1,1,2],Float64[1,1.6,1,1.6],Float64[3,1.6,3,1.6,1,1.6],
        Float64[4,1.6,1,1.6,1,1.6,4,1.6,1,1.6,1,1.6,1,1.6],
        Float64[4,1.6,4,1.6,1,1.6,1,1.6,4,1.6,4,1.6,1,1.6],
        Float64[4,1.6,4,1.6,4,1.6,1,1.6,1,1.6],Float64[5,1.6,2,1.6],
        Float64[],Float64[3,2,3,2],Float64[4,1.6,1,1.6,1,1.6],
        Float64[6,2,1,1,1,1,1,2],Float64[1,1.6,1,1.6],Float64[3,1.6,3,1.6,1,1.6],
        Float64[4,1.6,1,1.6,1,1.6,4,1.6,1,1.6,1,1.6,1,1.6],
        Float64[4,1.6,4,1.6,1,1.6,1,1.6,4,1.6,4,1.6,1,1.6],
        Float64[4,1.6,4,1.6,4,1.6,1,1.6,1,1.6],Float64[5,1.6,2,1.6]]
  # Define marker types
  mrk=String["s", "o", "H", "^", "X", "D", "*", "p", "v", "P",
             "x", "2", "|", ".", "3", "+", "1", "_", ">", "4"]


  # Save selected colour or subset of colours
  colours = colourschemes[cs][lc]
  # Save selected line styles
  if lt == 0
    lines = "None"
  elseif lt == -1
    lines = String["None" for dt in colourschemes[cs][lc]]
  elseif lt == -2
    lines = String["None" for dt in mt[lc]]
  else
    lines = dt[lt]
  end
  # Save selected marker styles
  if pt == 0
    markers = "None"
  elseif pt == -1
    markers = String["None" for dt in colourschemes[cs][lc]]
  elseif pt == -2
    markers = String["None" for dt in mt[lc]]
  else
    markers = mrk[pt]
  end

  # Return a set of colours/styles depending on the input choice
  return colours, lines, markers
end #function sel_ls
