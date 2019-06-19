"""
# Module pyp

Plot data with PyPlot as line, scatter or mixed plots with optional errors
as well as stacked area plots.


# Data structures
- PlotData


# Functions

## Public
- load_PlotData
- plot_data
- plot_stack
- sel_ls

## Private
- renameDF
- calc_errors
- setup_plot
- set_style
- plt_DataWithErrors
- redef_err
- set_axes
- find_limits
- setup_log
- rm_log
- checkaliases
- checkalias
- get_stack
- format_stack
- print_stack
- format_axes_and_annotations
"""
module pyp

# Track changes during development
using Revise

###  Load Julia packages  ###
using LaTeXStrings
import Dates
import PyPlot; const plt = PyPlot
import Parameters; const par = Parameters
import DataFrames; const df = DataFrames
import DataFrames.DataFrame
import Statistics; const stats = Statistics
import Dierckx; const spl = Dierckx
import LinearAlgebra.⋅


# Export public functions
export load_PlotData,
       plot_data,
       plot_stack,
       sel_ls,
       PlotData

### NEW TYPES
"""
    par.@with_kw mutable struct PlotData

Mutable struct holding information of individual graphs in a plot
that can be addressed with keyword arguments for the following fields:
- `x`/`y`: x/y data for plotting, x/y must be vectors of `Real`s (or `DateTime`s for x)
  of same length
- `xlerr`/`xuerr`/`ylerr`/`yuerr`: upper and lower errors for x and y data,
  must either be `nothing` or vector of `Real`s of same length as `x`/`y`
  (default: `nothing`)
- `label`: `AbstractString` with label for legend (default: `""` – empty String, no entry)
- `marker`: `String` or `Int64` specifying `PyPlot` marker type (default: `"None"` – no markers)
- `dashes`: tuple or vector of `Real`s specifying an 'on-off' sequence for line types
  (default: `[]` – empty array = solid lines); use `"None"` for no lines
- `colour`: `String` or `Symbol` with `PyPlot` colour name or String with RGB code
  (default: `nothing` – use `PyPlot` default colours)
- `lw`: `Real` to specify linewidth
- `alpha`: `Real` to specify opaqueness
"""
par.@with_kw mutable struct PlotData
  x::Union{Vector{Dates.DateTime}, Vector{T} where T<:Real}
  y::Vector{T} where T<:Real
  xuerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  xlerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  yuerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  ylerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  label::AbstractString=""
  marker::Union{String,Int64}="None"
  dashes::Union{Tuple{Real,Real},Vector{T} where T<:Real}=[]
  colour::Union{Nothing,String,Symbol}=nothing
  lw::Real=1.4
  alpha::Real=1
end


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
+ `min_xticks` (`Union{Real,Vector{Int64},Vector{Float64}}`): Size of interval between minor x ticks (**default:** `0` – automatic determination by PyPlot)
+ `min_yticks` (`Union{Int64, Array{Int64,1}}`):Size of interval between minor y ticks
  - `0` (**default**): automatic determination by PyPlot
  - If 2. y axis is defined, an array can be used to used different specifications for each axis
+ `maj_xticks`/`maj_yticks` (`Union{Real,Vector{Int64},Vector{Float64}}`/`Real`)
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
+ `ti_offset`, `ax_offset`, `leg_offset`: Offsets for fontsizes of title, axes
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
                  date_format::String="",
                  xlabel::AbstractString="model time / hours",
                  ylabel::Union{AbstractString, Vector{T} where T<:AbstractString} =
                  "concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
                  ti::AbstractString="", twinax::Vector{T} where T<:Integer=Int64[],
                  logscale::Union{String,Vector{T} where T<:String}="",
                  logremove::Union{String,Vector{T} where T<:String}="neg",
                  plot_type::String="default", cs::Union{String,Vector{T} where T<:String}="",
                  lt="default", pt="default",
                  lc::Union{String,Symbol, Vector{T} where T<:String, Vector{T} where T<:Symbol}="default",
                  alpha::Real=-1, xlims=nothing, ylims=nothing, mticks::Bool=true,
                  min_xticks::Union{Real,Vector{T} where T<:Real}=0,
                  min_yticks::Union{Real,Vector{T} where T<:Real}=0,
                  maj_xticks::Union{Real,Vector{T} where T<:Real}=0,
                  maj_yticks::Union{Real,Vector{T} where T<:Real}=0,
                  figsize::Tuple{Real,Real}=(6,4), fontsize::Real=12,
                  framewidth::Real=1, ticksize::Tuple{Real,Real}=(4.5,2.5),
                  cap_offset::Real=0, ti_offset::Real=4, ax_offset::Real=2,
                  leg_offset::Real=0, legcolumns::Int64=1,
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
    maj_yticks, min_yticks, cs, lc, lt, pt, axcolour =
    setup_plot(plot_list, figsize, twinax, ylabel, logscale, logremove, xlims, ylims,
    maj_yticks, min_yticks, plot_type, cs, lc, lt, pt, alpha, axcolour)

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
  fig, ax = format_axes_and_annotations(fig, ax, pltdata, ti, xlabel, ylabel, date_format,
    fontsize, legpos, legcolumns, axcolour, leg_offset, ti_offset, ax_offset,
    maj_xticks, maj_yticks, min_xticks, min_yticks, mticks, ticksize, framewidth)

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
+ `boundaries` (`Real`): transparency value for boundary lines of each stack data (`0` **default** for no boundaries)
+ `alpha`: (`Real`): set transparency of areas of each stack data
  - `0` (**default**): Use transparency setting of the original `PlotData`
  - `1`: no transparency
+ `logscale` (`String`): use `"x"`, `"y"` or `"xy"` to set either x-, y- or both axes to a logarithmic scale
  (**default:** `""` – linear scale)
+ `cs` (`String`): Define colour scheme `default`, `source` or `sink` from `function sel_ls`
  - `""` (**default:**): colours as defined in the `PlotData`, where `nothing` is replaced by colours of the default colour scheme)
  - `"own"`: colours as define by array `lc`
+ `lt`: dash type (defined by PyPlot's `dashes`) for boundaries
  - for use only, when `boundaries` are switched on (with a value greater 0)
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
+ `min_xticks` (`Union{Real,Vector{Int64},Vector{Float64}}`): Size of interval between minor x ticks
  (**default:** `0` – automatic determination by PyPlot)
+ `min_yticks` (`Real`): Size of interval between minor y ticks
  - `0` (**default**): automatic determination by PyPlot
+ `maj_xticks`/`maj_yticks` (`Union{Real,Vector{Int64},Vector{Float64}}`/`Real`) interval size of major ticks in x- and y-axis, respectively;
+ `figsiz` (`Tuple{Real,Real}`): figure size width × height in inches
  (**default:** `(6,4)`)
+ `fntsiz` (`Real`): default font size (**default:** `12`)
+ `framewidth` (`Real`): default line width of outer axes frame (**default:** `1`)
+ `ticksize` (`Tuple{Real,Real}`): Tuple with scaling factors for length of major
  (first tuple entry) and minor (second tuple entry) ticks in relation to their
  width (**default:** `(4.5,2.5)`)
+ `ti_offset`, `ax_offset`, `leg_offset`: Offsets for fontsizes of title, axes
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
         date_format::String="",
         xlabel::AbstractString="model time / hours",
         ylabel::AbstractString="concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
         ti::AbstractString="", boundaries=0,
         logscale::String="", logremove::String="neg",
         cs::String="", lt=[], lc=[], alpha::Real=0,
         xlims=nothing, ylims=nothing, mticks::Bool=true,
         min_xticks::Real = 0,
         min_yticks::Real = 0,
         maj_xticks::Real = 0,
         maj_yticks::Real = 0,
         figsize::Tuple{Real,Real}=(6,4), fontsize::Real=12, framewidth::Real=1,
         ticksize::Tuple{Real,Real}=(4.5,2.5), ti_offset::Real=4,
         ax_offset::Real=2, leg_offset::Real=0, legcolumns::Int64=1,
         legpos::Union{String, Int64, Tuple{Real, Real}}="best",
         interpolate=0, extrapolate::Union{Bool,String}=false, kspline::Int64=3)

  # Get x and y data
  xdata, ystack, labels = get_stack(plot_list...)
  ystack, ylines = interpolate_stack(xdata, ystack, interpolate, extrapolate, kspline)

  # Format plot, set colour scheme
  colour, linetype, α = format_stack(cs, alpha, lc, lt, plot_list...)

  # Plot data as stack with optional boundary lines
  fig, ax = print_stack(xdata, ystack, ylines, boundaries, labels, colour, lt, α, figsize)

  # Set axis limits and log scales
  if xlims == nothing  xlims = (nothing, nothing)  end
  if ylims == nothing  ylims = (nothing, nothing)  end
  ax = set_axes([plot_list], [ax], [logscale], [xlims], [ylims])

  # Format plot
  fig, ax = format_axes_and_annotations(fig, ax, [plot_list], ti, xlabel, [ylabel], date_format,
    fontsize, legpos, legcolumns, ["black"], leg_offset, ti_offset, ax_offset,
    maj_xticks, [maj_yticks], min_xticks, [min_yticks], mticks, ticksize, framewidth)


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


###########################
###  PRIVATE FUNCTIONS  ###
###########################

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


### Functions associated with plot_data

"""
    setup_plot(plot_list, twinax, ylab, logscale, logremove, xlims, ylims, maj_yticks, min_yticks, ptype, cs, lc, lt, pt, alpha, axcolour)

If `twinax` is an array split `plot_list` into 2 datasets based on the indices
in `twinax`. Assure all paramters concerning y data are arrays of length 2, if
a second axis is allowed or length 1 otherwise. Return `plot_list` as array of
2 or 1 distinct lists an the adjusted parameters.
"""
function setup_plot(plot_list, figsize, twinax, ylab, logscale, logremove, xlims, ylims,
                    maj_yticks, min_yticks, ptype, cs, lc, lt, pt, alpha, axcolour)

  # Start plot
  fig, ax = plt.subplots(figsize=figsize)
  ax = [ax]

  # Set flag for second axis, asume no second axis
  xlim = []; ylim = []

  [p.alpha = alpha for p in plot_list if alpha ≥ 0]

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
    if maj_yticks isa Int64 maj_yticks = Int64[maj_yticks, maj_yticks]  end
    if min_yticks isa Int64 min_yticks = Int64[min_yticks, min_yticks]  end
    if xlims == nothing
      xlimit = [[nothing, nothing],[nothing, nothing]]
    elseif xlims isa Tuple
      xl = [xlims[1],xlims[2]]
      xlimit = [xl, xl]
    elseif xlims isa Array
      xlimit = [[xlims[1][1], xlims[1][2]], [xlims[2][1], xlims[2][2]]]
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
    if maj_yticks isa Real maj_yticks = Int64[maj_yticks]  end
    if min_yticks isa Real min_yticks = Int64[min_yticks]  end
    if xlims == nothing     xlimit = [[nothing, nothing]]
    elseif xlims isa Tuple  xlimit = [[xlims[1],xlims[2]]]
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
         maj_yticks, min_yticks, cs, lc, lt, pt, axcolour
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
      ax[i].set_xlim(find_limits(pltdata[i], "x", xlim[i]))
      ax[i].set_xscale("log")
    else
      ax[i].set_xlim(xlim[i])
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
  alpha > 0 ? α = alpha : stats.mean(α) > 0 ? α = stats.mean(α) : α = 1
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
function print_stack(xdata, ystack, ylines, boundaries, labels, colours, lt, α, figsize)

  # Start plot
  fig, ax = plt.subplots(figsize=figsize)

  # Plot data
  ax.stackplot(xdata, ystack, labels=labels, colors=colours, alpha=α)

  # Resume, if optional boundaries are skipped
  if boundaries==0  return fig, ax  end
  [ax.plot(xdata, ylines[i], color=colours[i], dashes = lt[i], alpha=boundaries)
    for i = 1:length(ylines)]

  return fig, ax
end #function print_stack


"""


"""
function format_axes_and_annotations(fig, ax, plot_list, ti, xlabel, ylabel, date_format,
  fontsize, legpos, legcolumns, axcolour, leg_offset, ti_offset, ax_offset,
  maj_xticks, maj_yticks, min_xticks, min_yticks, mticks, ticksize, framewidth)

  # Set plot title
  ax[1].set_title(ti, fontsize=fontsize+ti_offset)

  # Generate axes labels and legend, define axes label/tick colours
  ax[1].set_xlabel(xlabel,fontsize=fontsize+ax_offset)
  for n = 1:length(ylabel)
    ax[n].set_ylabel(ylabel[n],fontsize=fontsize+ax_offset, color=axcolour[n])
  end
  [plt.setp(ax[n].get_yticklabels(),color=axcolour[n]) for n = 1:length(axcolour)]

  # Set ticks and optional minor ticks
  if typeof(plot_list[1][1].x) ≠ Vector{Dates.DateTime}  if maj_xticks > 0
    xint = collect(ax[1].get_xlim()[1]:maj_xticks:ax[1].get_xlim()[2])
    for i = 1:length(plot_list)  ax[i].set_xticks(xint)  end
  end  end
  if maj_yticks[1] > 0  for i = 1:length(maj_yticks)
    yint = collect(ax[i].get_ylim()[1]:maj_yticks[i]:ax[i].get_ylim()[2])
    ax[i].set_yticks(yint)
  end  end
  if mticks == true
    plt.minorticks_on()
  else
    plt.minorticks_off()
  end
  # Set minor x ticks
  if typeof(plot_list[1][1].x) ≠ Vector{Dates.DateTime}
    if min_xticks > 0
      mx = plt.matplotlib.ticker.MultipleLocator(min_xticks)
      for i = 1:length(plot_list)
        ax[i].xaxis.set_minor_locator(mx)
      end
    end
  elseif min_xticks isa Real
    min_xticks = [6,12,18]
  end
  # Set minor y ticks
  for i = 1:length(min_yticks)
    if min_yticks[i] > 0
      my = plt.matplotlib.ticker.MultipleLocator(min_yticks[i])
      ax[i].yaxis.set_minor_locator(my)
    end
  end
  # Format ticks and frame
  Mtlen = ticksize[1]⋅framewidth
  mtlen = ticksize[2]⋅framewidth
  for i = 1:length(plot_list)
    ax[i].tick_params("both", which="both", direction="in", top=true, right=true,
      labelsize=fontsize, width=framewidth)
    ax[i].grid(linestyle=":", linewidth = framewidth)
    ax[i].spines["bottom"].set_linewidth(framewidth)
    ax[i].spines["top"].set_linewidth(framewidth)
    ax[i].spines["left"].set_linewidth(framewidth)
    ax[i].spines["right"].set_linewidth(framewidth)
    if typeof(plot_list[1][1].x) ≠ Vector{Dates.DateTime}
      ax[i].tick_params("both", which="major", length=Mtlen)
      ax[i].tick_params("both", which="minor", length=mtlen)
    else
      ax[i].set_xlim(left=plot_list[1][1].x[1], right=plot_list[1][1].x[end])
      if date_format == ""
        maj_xticks isa Vector ? date_format = "%d. %b, %H:%M" : date_format = "%d. %b"
      end
      majorformatter = plt.matplotlib.dates.DateFormatter(date_format)
      minorformatter = plt.matplotlib.dates.DateFormatter("")
      majorlocator = plt.matplotlib.dates.HourLocator(byhour=maj_xticks)
      minorlocator = plt.matplotlib.dates.HourLocator(byhour=min_xticks)
      ax[i].xaxis.set_major_formatter(majorformatter)
      ax[i].xaxis.set_minor_formatter(minorformatter)
      ax[i].xaxis.set_major_locator(majorlocator)
      ax[i].xaxis.set_minor_locator(minorlocator)
      fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right")
    end
  end

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
  if typeof(plot_list[1][1].x) ≠ Vector{Dates.DateTime}  fig.tight_layout()  end

  return fig, ax
end #function format_axes_and_annotations


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
end #module pyp
