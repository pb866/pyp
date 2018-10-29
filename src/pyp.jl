"""
# Module pyp

Load data and plot with PyPlot.


# Data structures

- PlotData


# Functions

## Public
- rd_data
- load_plotdata
- plot_data
- plot_stack
- sel_ls

## Private
- rename_DF
- calc_errors
- setup_axes
- set_cs
- plt_DataWithErrors
- redef_err
- set_log
- logBounds
- set_log2
"""
module pyp


###  Load Julia packages  ###
using DataFrames
using PyPlot
using LaTeXStrings
using Parameters
using Dates

# Load further modules
using filehandling: test_file


# Export public functions
export rd_data,
       lineplot,
       load_PlotData,
       plot_data,
       plot_stack,
       sel_ls,
       PlotData

### NEW TYPES
@with_kw mutable struct PlotData
  x::Union{Vector{Int64},Vector{Float64},Vector{DateTime}}
  y::Union{Vector{Int64},Vector{Float64}}
  xuerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing
  xlerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing
  yuerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing
  ylerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing
  label::Union{Nothing,AbstractString}=""
  marker::Union{AbstractString,Int64}="None"
  dashes::Union{Tuple{Number,Number},Vector{Float64},Vector{Int64},Vector{Any}}=Float64[]
  colour::Union{Nothing,AbstractString}=nothing
  lw::Number=1.4
  alpha::Number=1
end


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    rd_data(ifile; \*\*kwargs)


Read data from text file `ifile` in the following format:

    # Use '#' as first character for comments or data to be ignored.
    # Optional: specify column header names for DataFrame keywords
    # with keyword `jlheaders` as:
    # jlheaders: x1, y1, y2, ..., yn
    # (You may use whitespace, commas (,), semicolons (;), or pipes (|) as
    #  separators in the list above)
    <data separated by whitespace or any character/string>

If columns don't have the same length, it can be specified whether the first or last
columns will be filled with `NaN`s.

The function uses several keyword arguments (\*\*kwargs) for more freedom in the
file format or the selection of data.

### \*\*kwargs

- `dir` (`String`; default: `.`): Directory of the input file `ifile`
- `ix` (`Int64`; default: `1`): Column index for column in `ifile` holding the x data
  (default column name in output DataFrame: `x`). If `ix` is set to `0`, no x column
  is assigned and only y columns are used in the DataFrame.
- `iy` (default: `0`): Column index/indices of y data columns in `ifile`. If `iy`
  is set to `0`, all columns starting at column 2 are assigned as y columns
  (default column name(s) in output DataFrame: `y1` ... `yn`).
  Columns can be specified using an integer for the selection of a single column,
  ranges (`<n(low)>:<n(high)>`), or arrays with special selections (`[n1, n2, ..., nn]`)
  where column order can be rearranged.
- `headers` (`Bool`; default: `false`): If headers is set to `true`, you need to specify
  column header names for all columns for the output DataFrame as described above
  using the keyword `jlheaders`. All columns will be read in and saved the same
  order as in `ifile`.
- `SF` (default value: `1` (no scaling)): You can optionally apply scaling factors
  to y data. If `SF` is an integer, the scaling factor will be applied to all
  y columns. You can apply scaling to each y column individually by providing an
  array of scaling factors of `length number of columns - number of x column`. If
  you only want to scale certain column(s), set the scaling factors for these columns
  and use `1` otherwise.
- `sep` (default: `whitespace`): You can specify any column separator with the
  keyword charactar `sep`. Separators can be any unicode character (even special
  characters such as `≠` or `α`) or string series of unicode characters
  (including whitespace).
- `colfill` (`Int64`; default: `"last"`): If the column length of the input file varies,
  the `first` or `last` columns of the file are filled with `NaN`s according to
  the keyword. If you have a file with shorter columns to the right and the left,
  you either need to rearrange columns in the original data file or try to work
  with a specifically defined separator `sep`.
- `ncols` (`String`; default: `-1`): Defines the number of columns (x + y columns) in a file.
  If set to a negative number, the number of columns is derived from the `jlheaders`
  array or, if obsolete, from the first non-comment line of the file. You should
  only have to set the number of columns, if you have columns of different length
  with leading missing numbers.
- `skip_header` (`Int64`; default: `0`): Define how many lines to skip at the beginning
  of a file in case comments aren't used
- `skip_footer` (`Int64`; default: `0`): Define how many lines to skip at the end
  of a file in case comments aren't used
"""
function rd_data(ifile; dir::String=".", ix::Int64=1, iy::Union{Int64,Vector{Int64}}=0,
  headers::Bool=false, SF=1, sep::String="",
  colfill::String="last", ncols::Int64=-1, skip_header::Int64=0, skip_footer::Int64=0)
  # Read input file
  ifile = test_file(ifile, dir = dir) # check existence of file
  lines = String[]; colnames = String[] # initialise arrays
  open(ifile,"r") do f
    # Read file
    lines = readlines(f)
    # Extract column header names, if specified
    if headers == true
      icol = findfirst([occursin("jlheaders:", line) for line in lines])
      colnames = split(lines[icol],"jlheaders:")[2]
      colnames = replace(colnames,r",|;|\|" => " ")
      colnames = split(colnames)
      # Make sure, columns aren't rearranged
      if ix > 1  ix = 1  end
      iy = 0
    end
    # Skip first lines of a file, if skip_header is set to integer > 0
    deleteat!(lines,1:skip_header)
    # Skip last lines of a file, if skip_footer is set to integer > 0
    deleteat!(lines,1+length(lines)-skip_footer:length(lines))
    # Find and delete comment lines
    del=findall(startswith.(lines,"#"))
    deleteat!(lines,del)
    # Find and delete empty lines
    del=findall(lines.=="")
    deleteat!(lines,del)
  end

  # Set number of columns
  if ncols < 0 && length(colnames) > 0
    ncols = length(colnames)
  elseif ncols < 0
    ncols = length(split(lines[1]))
  end
  # Determine number of y columns for default case
  if iy == 0  && ix == 0
    iy = 1:ncols
  elseif iy == 0  && ix > 0
    iy = 2:ncols
  end

  # Initilise x and y data
  if ix > 0  x = Float64[]  end
  y = Matrix{Float64}(undef, 0, length(iy))
  # Loop over data lines
  for line in lines
    # Split into columns
    if sep == ""
      # Assume whitespace as default separator
      raw = split(line)
    else
      # Use separator, if specified
      raw = split(line,sep)
    end
    # Check number of current columns against maximum number of columns
    if length(raw) > ncols
      println("WARNING! Number of columns read in greater than defined number of columns.")
      println("The $colfill $(length(raw)-ncols) columns are ignored.")
      if lowercase(colfill[1]) == "l"
        raw = raw[1:ncols]
      else
        raw = raw[length(raw)-ncols+1:end]
      end
    end
    # Save current line to respective data arrays
    if ix > 0  push!(x,parse(Float64, raw[ix]))  end
    ydat = transpose(parse.(Float64, raw[1:end .!= ix]).*SF)
    ix == 0 ? nx = 0 : nx = 1
    if length(ydat)<ncols && colfill == "last"
      for i=length(ydat)+1:ncols-nx  ydat = hcat(ydat,NaN)  end
    elseif length(ydat)<ncols && colfill == "first"
      for i=length(ydat)+1:ncols-nx  ydat = hcat(NaN,ydat)  end
    end
    y = vcat(y,ydat)
  end

  # Generate output DataFrame
  if ix == 0  output = DataFrame()
  else        output = DataFrame(x = x)
  end
  for i = 1:length(iy)
    output[Symbol("y$i")] = y[:,i]
  end
  # Rename headers, if names were specified
  if headers == true
    if length(output) != length(colnames)
      println("Warning! Length of column names not equal to number of columns.")
      println("Using standard names x for first column and y1...yn for remaining columns.")
    else
      for i = 1:length(colnames)
        rename!(output,names(output)[i],Symbol(colnames[i]))
      end
    end
  end

  # Return file data as DataFrame
  return output
end #function rd_data


"""
    load_PlotData(pltdata::DataFrames.DataFrame;  \*\*kwargs)

Pack x and y data with possible assoiciated errors from DataFrame `pltdata`
as well as formatting parameters into a new DataType `PlotData`.


### \*\*kwargs

+ `err` (`String`):
  - `"None"` (no errors, **default**)
  - `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
  - `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
  - `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
  - `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
  - `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
  - `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
  - `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
+ `mt` (`Union{AbstractString,Int64}`): PyPlot marker type, see e.g. https://matplotlib.org/api/markers_api.html
  - `"None"` (**default**) to suppress use of markers
+ `lt` (`Union{String,Tuple{Int64,Int64},Array{Int64,1}}`), tuple or array with
  on/off PyPlot dash definitions
  - empty array means a solid line (**default**)
  - Tuple with 2 entries or array with 4 entries where odd entries are `0` suppresses lines
+ `lc` (`Union{Nothing,AbstractString}`) PyPlot line/marker colours (**default:** `nothing` –
  use PyPlot default colours), see e.g. https://matplotlib.org/examples/color/named_colors.html
+ `lw` (`Number`): linewidth (**default:** `1.4`)
+ `SF` (`Number`): scaling factor of y data and associated errors (**default:** `1`, no scaling)
+ `label` (`String`): Label for legend
  - `""` (empty string, **default**): no legend label
+ `renameDF` (`Union{Bool,Array{Symbol,1}}`):
  - `"true"` (assume columns in the order: `x`, `y`, `ylerr`, `yuerr`, `xlerr`, `xuerr`)
  - `"false"` (columns already in correct order with names as for `true`)
  - `Array{Symbol, 1}`: give array with column names for `x`, `y`, `ylerr`, `yuerr`, `xlerr`, `xuerr`
   (give complete list, even if columns are incomplete due to the choice of `err`)
"""
function load_PlotData(plotdata::DataFrames.DataFrame;  err::String="None",
         mt::Union{AbstractString,Int64}="None",
         lt::Union{String,Tuple{Number,Number},Vector{Int64},Vector{Float64}}=Float64[],
         lc::Union{Nothing,AbstractString}=nothing, lw::Number=1.4, SF::Number=1,
         label::String="", alpha::Number=1, renameDF::Union{Bool,Vector{Symbol}}=true)

  # Make copy of plotdata that can be altered
  pltdata = deepcopy(plotdata)
  # (Re-)define column names of DataFrame
  if renameDF == true
    DFnames = Symbol[:x, :y, :ylerr, :yuerr, :xlerr, :xuerr]
    pltdata = rename_DF(pltdata, DFnames, err)
  elseif renameDF isa Array{Symbol, 1}
    if length(renameDF) ≠ 6
      println("\'renameDF\' not correctly defined. Define all column names for")
      println("x, y, ylerr, yuerr, xlerr, and xuerr. Script stopped."); exit()
    end
    DFnames = deepcopy(renameDF)
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
         label = label, marker = mt, dashes = lt, colour = lc, lw = lw, alpha=alpha)
end #function load_PlotData


"""
    plot_data(plot_list::Union{Array{PlotData, 1}, Array{pyp.PlotData, 1}}, \*\*kwargs)

Generate scatter and/or line plots from the varargs `plot_data` of Type `PlotData`.


### \*\*kwargs

For all arguments concerning y data the following applies:
- The dataset can be split into 2 subsets and a 2. y-axis with a different scale
  can be assigned (see optional parameter `twinax`)
- If a second axis is defined, parameters concerning y-data can be compiled in an
  array of length `2` with different values for the first and second y-axis,
  respectively
- If the second y-axis is defined and only a single parameter concerning y-data
  is defined, then this parameter applies to both axes

Keyword arguments for `function plot_data` are:
+ `xlabel` (`Union{String, LaTeXString}`): x axis label
  - **default:** `"model time / hours"`
+ `ylabel` (`Union{String, LaTeXString, Array{String,1}, Array{LaTeXString,1}}`):
  y axis label, `Array` can be used, if a second axis with a different label is introduced
  - **default:** `"concentration / mlc cm\$^{-3}\$ s\$^{-1}\$"`
+ `ti` (`Union{String, LaTeXString}`): Plot title
  - **default:** `""` (empty string) for no title
+ `twinax` (`Array{Int64,1}`): optional array to devide data into 2 subsets and assign second dataset to a 2. y-axis; use array of `length` of the PlotData with integers `1` and `2` to assign each data to axis 1 and 2
  (**default:**: empty array `Int64[]`, no 2. y-axis)
+ `logscale` (`Union{String, Array{String, 1}}`): use `"x"`, `"y"` or `"xy"` to
  set either x-, y- or both axes to a logarithmic scale;
  if 2. y-axis is present, arrays may be used to have different scales on each y-axis;
  (**default:** `""` – linear scale)
+ `plot_type` (`String`): Choose between the following plot types:
  - `"line"`: line plot
  - `"scatter"`: scatter plot
  - `"both"`: markers connected by lines
  - `"mix"`: mix of line and scatter plots or both; for this plot type `lc`, `lt`, and `pt` need to be defined (see below)
  -  `"own"` (**default**): Use formats as predefined by the PlotData
+ `cs` (`Union{String, Array{String,1}}`): Define colour scheme `default`, `source`
  or `sink` from `function sel_ls`; array with 2 different colour schemes may be
  used in case of 2. y-axis
  - `""` (**default:**): colours as defined in the `PlotData` or default PyPlot colour scheme)
  - `"own"`: set your own colour schemes using parameters `lt` and `lc` (mandatory, if `cs` is set to `"own"`)
+ `lt`: dash type (defined by PyPlot's `dashes`)
  - for use only, when `plot_type` is switched to `"mix"`
  - Array of `length 1 or 2` with Arrays of `length(plot_data)` where inner arrays specify line type of each `PlotData`
  - use tuple or array of integers to define on/off dash pixels
  - `[]` (empty array, **default**): solid lines
  - tuple with 2 entries or array with at least 4 entries where odd entries are `0` for no lines
+ `lc`: line/marker colour (defined by PyPlot's `color` keywords, see e.g. https://matplotlib.org/examples/color/named_colors.html)
  - for use only, when `plot_type` is switched to `"mix"`
  - Array of `length 1 or 2` with Arrays of `length(plot_data)` where inner arrays specify the line/marker colour of each `PlotData`
+ `alpha`: if alpha is set to a value `> 0` then all graphs will be displayed with this transparency value
+ `mticks` (`String`): Switch minor ticks on/off (**default:** "on")
+ `nmxt` (`Union{Number,Vector{Int64},Vector{Float64}}`): Size of interval between minor x ticks (**default:** `0` – automatic determination by PyPlot)
+ `nmyt` (`Union{Int64, Array{Int64,1}}`):Size of interval between minor y ticks
  - `0` (**default**): automatic determination by PyPlot
  - If 2. y axis is defined, an array can be used to used different specifications for each axis
+ `Mxtint`/`Mytint` (`Union{Number,Vector{Int64},Vector{Float64}}`/`Number`)
  interval size of major ticks in x- and y-axis, respectively; asign different intervals
  to 2nd y-axis using an array of integers
+ `xlims`/`ylims`: Tuple of minimum/maximum value for each axis
  - `nothing` may be used for default PyPlot values (**default**)
  - `nothing` may be used for one value within the tuple to use the default, e.g. `ylims = (0, nothing)`
  - if 2. y-axis is defined, `ylims` may be an array of tuples with different specifications
+ `figsiz` (`Tuple{Number,Number}`): figure size width × height in inches
  (**default:** `(6,4)`)
+ `fntsiz` (`Number`): default font size (**default:** `12`)
+ `frw` (`Number`): default line width of outer axes frame (**default:** `1`)
+ `tsc` (`Tuple{Number,Number}`): Tuple with scaling factors for length of major
  (first tuple entry) and minor (second tuple entry) ticks in relation to their
  width (**default:** `(4.5,2.5)`)
+ `cap_offset` (`Number`): deviation from default cap size of `3` of marker error bars
  (**default:** `0`)
+ `ti_offset`, `ax_offset`, `leg_offset`: Offsets for fontsizes of title, axes
  labels, and legend, respectively, from default sizes. Numbers
  (positive or negative) will be added to fontsizes
  (**defaults:** `4`, `2`, `0`)
+ `axcol` (`Union{String,Array{String,1}}`): colour of y axis label and tick numbers;
  if 2. y axis is defined, different values may be defined in an array for each axis
  (**default:** `"black"`)
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
- `legcol` (`Union{Int64, Array{Int64,1}}`): number of legend columns
  (**default:** `1`)
"""
function plot_data(plot_list::Union{PlotData,pyp.PlotData}...;
                  xlabel::Union{String, LaTeXString}="model time / hours",
                  ylabel::Union{String, LaTeXString, Vector{String}, Vector{LaTeXString}}=
                  "concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
                  ti::Union{String, LaTeXString}="", twinax::Vector{Int64}=Int64[],
                  logscale::Union{String, Vector{String}}="",
                  plot_type::String="own", cs::Union{String, Vector{String}}="", lt=[], pt="s",
                  lc::Union{String, Vector{String}}="black", alpha::Number=-1,
                  xlims=nothing, ylims=nothing, mticks::String="on",
                  nmxt::Union{Int64,Vector{Int64},Vector{Float64}}=0, nmyt::Union{Int64, Vector{Int64}}=0,
                  Mxtint::Union{Number,Vector{Int64},Vector{Float64}}=0, Mytint::Union{Number,Vector{Int64},Vector{Float64}}=0,
                  figsiz::Tuple{Number,Number}=(6,4), fntsiz::Number=12,
                  frw::Number=1, tsc::Tuple{Number,Number}=(4.5,2.5), cap_offset::Number=0,
                  ti_offset::Number=4, ax_offset::Number=2,
                  axcol::Union{String,Vector{String}}="black", leg_offset::Number=0,
                  legpos::Union{String, Int64, Vector{String}, Vector{Int64}}="best",
                  legcol::Union{Int64, Vector{Int64}}=1)

  # Start plot
  fig, ax1 = subplots(figsize=figsiz)
  # Check for twin axes and devide datasets
  plt, ax2, ax_2, ylabel, logscale, xlim, ylim, Mytint, nmyt, lc, lt, pt, axcol, legpos, legcol =
    setup_axes(plot_list, twinax, ylabel, logscale, xlims, ylims, Mytint, nmyt,
    plot_type, lc, lt, pt, alpha, axcol, legpos, legcol)
  if ax_2  ax = [ax1, ax2]
  else     ax = [ax1]
  end

  # set colour scheme
  plt = set_cs(plt, plot_type, cs, lc, lt, pt, ax_2)

  # Plot data and associated errors
  plt[1], ax1 = plt_DataWithErrors(plt[1], ax1, cap_offset)
  if ax_2  plt[2], ax2 = plt_DataWithErrors(plt[2], ax2, cap_offset)  end

  # Define logscales
  xlim, ylim = logBounds(xlim, ylim, ax)
  for i = 1:length(logscale)
    ax[i], xlim[i] = set_log(plt[i], ax[i], xlim[i], logscale[i], "x",
                      :x, :xlerr, :xuerr, :set_xlim, :set_xscale)
    ax[i], ylim[i] = set_log(plt[i], ax[i], ylim[i], logscale[i], "y",
                      :y, :ylerr, :yuerr, :set_ylim, :set_yscale)
  end

  # Set axes limits
  ax1[:set_xlim](xlim[1]); ax1[:set_ylim](ylim[1])
  if ax_2  ax2[:set_xlim](xlim[2]); ax2[:set_ylim](ylim[2])  end

  # Set plot title
  ax1[:set_title](ti, fontsize=fntsiz+ti_offset)

  # Generate axes labels and legend, define axes label/tick colours
  ax1[:set_xlabel](xlabel,fontsize=fntsiz+ax_offset)
  for n = 1:length(ylabel)
    ax[n][:set_ylabel](ylabel[n],fontsize=fntsiz+ax_offset, color=axcol[n])
  end
  [setp(ax[n][:get_yticklabels](),color=axcol[n]) for n = 1:length(axcol)]

  if legpos[1] ≠ "None" && any([p.label≠"" for p in plot_list])
    ax1[:legend](fontsize=fntsiz+leg_offset, loc=legpos[1], ncols=legcol[1])
  end
  if ax_2 && legpos[2] ≠ "None" && any([p.label≠"" for p in plot_list])
    ax2[:legend](fontsize=fntsiz+leg_offset, loc=legpos[2], ncols=legcol[2])
  end

  # Set ticks and optional minor ticks
  if typeof(plot_list[1].x) ≠ Vector{DateTime}  if Mxtint > 0
    xint = collect(ax1[:get_xlim]()[1]:Mxtint:ax1[:get_xlim]()[2])
    ax1[:set_xticks](xint)
    if ax_2  ax2[:set_xticks](xint)  end
  end  end
  if Mytint[1] > 0  for i = 1:length(Mytint)
    yint = collect(ax[i][:get_ylim]()[1]:Mytint[i]:ax[i][:get_ylim]()[2])
    ax[i][:set_yticks](yint)
  end  end
  if mticks == "on"
    minorticks_on()
  else
    minorticks_off()
  end
  # Set minor x ticks
  if typeof(plot_list[1].x) ≠ Vector{DateTime}
    if nmxt > 0
      mx = matplotlib[:ticker][:MultipleLocator](nmxt)
      for a in ax
        a[:xaxis][:set_minor_locator](mx)
      end
    end
  elseif nmxt isa Number
    nmxt = [6,12,18]
  end
  # Set minor y ticks
  for i = 1:length(nmyt)
    if nmyt[i] > 0
      my = matplotlib[:ticker][:MultipleLocator](nmyt[i])
      ax[i][:yaxis][:set_minor_locator](my)
    end
  end
  # Format ticks and frame
  Mtlen = tsc[1]⋅frw
  mtlen = tsc[2]⋅frw
  for a in ax
    a[:tick_params]("both", which="both", direction="in", top="on", right="on",
      labelsize=fntsiz, width=frw)
    a[:grid](linestyle=":", linewidth = frw)
    a[:spines]["bottom"][:set_linewidth](frw)
    a[:spines]["top"][:set_linewidth](frw)
    a[:spines]["left"][:set_linewidth](frw)
    a[:spines]["right"][:set_linewidth](frw)
    if typeof(plot_list[1].x) ≠ Vector{DateTime}
      a[:tick_params]("both", which="major", length=Mtlen)
      a[:tick_params]("both", which="minor", length=mtlen)
    else
      a[:set_xlim](xmin=plot_list[1].x[1], xmax=plot_list[1].x[end])
      if Mxtint isa Vector
        majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b, %H:%M")
      else
        majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b")
      end
      minorformatter = matplotlib[:dates][:DateFormatter]("")
      majorlocator = matplotlib[:dates][:HourLocator](byhour=Mxtint)
      minorlocator = matplotlib[:dates][:HourLocator](byhour=nmxt)
      a[:xaxis][:set_major_formatter](majorformatter)
      a[:xaxis][:set_minor_formatter](minorformatter)
      a[:xaxis][:set_major_locator](majorlocator)
      a[:xaxis][:set_minor_locator](minorlocator)
      fig[:autofmt_xdate](bottom=0.2,rotation=-30,ha="left")
    end
  end

  tight_layout()

  # Return PyPlot data
  return fig, ax
end #function plot_data


"""
    plot_stack(plot_list::Union{PlotData,pyp.PlotData}...; \*\*kwargs)

Plot the `PlotData` listed in vararg `plot_list` as a stack plot, where the following
keyword arguments are possible.

## kwargs

+ `xlabel` (`Union{String, LaTeXString}`): x axis label
  - **default:** `"model time / hours"`
+ `ylabel` (`Union{String, LaTeXString, Array{String,1}, Array{LaTeXString,1}}`): y axis label, `Array` can be used, if a second axis with a different label is introduced
  - **default:** `"concentration / mlc cm\$^{-3}\$ s\$^{-1}\$"`
+ `ti` (`String`): Plot title
  - **default:** `""` (empty string) for no title
+ `boundaries` (`Number`): transparency value for boundary lines of each stack data (`0` **default** for no boundaries)
+ `alpha`: (`Number`): set transparency of areas of each stack data
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
+ `mticks` (`String`): Switch minor ticks on/off (**default:** "on")
+ `min_xticks` (`Union{Number,Vector{Int64},Vector{Float64}}`): Size of interval between minor x ticks
  (**default:** `0` – automatic determination by PyPlot)
+ `min_yticks` (`Number`): Size of interval between minor y ticks
  - `0` (**default**): automatic determination by PyPlot
+ `maj_xticks`/`maj_yticks` (`Union{Number,Vector{Int64},Vector{Float64}}`/`Number`) interval size of major ticks in x- and y-axis, respectively;
+ `figsiz` (`Tuple{Number,Number}`): figure size width × height in inches
  (**default:** `(6,4)`)
+ `fntsiz` (`Number`): default font size (**default:** `12`)
+ `frw` (`Number`): default line width of outer axes frame (**default:** `1`)
+ `tsc` (`Tuple{Number,Number}`): Tuple with scaling factors for length of major
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
- `legcol` (`Union{Int64, Array{Int64,1}}`): number of legend columns
  (**default:** `1`)
"""
function plot_stack(plot_list::Union{PlotData,pyp.PlotData}...;
         xlabel::Union{String, LaTeXString}="model time / hours",
         ylabel::Union{String, LaTeXString}=
         "concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
         ti::String="", boundaries::Number=0,
         logscale::String="", cs::String="", lt=[], lc=[],
         alpha::Number=0, xlims=nothing, ylims=nothing, mticks::String="on",
         min_xticks::Union{Number,Vector{Int64},Vector{Float64}}=0, min_yticks::Number=0,
         maj_xticks::Union{Number,Vector{Int64},Vector{Float64}}=0, maj_yticks::Number=0,
         figsiz::Tuple{Number,Number}=(6,4), fntsiz::Number=12, frw::Number=1,
         tsc::Tuple{Number,Number}=(4.5,2.5), ti_offset::Number=4,
         ax_offset::Number=2, leg_offset::Number=0,
         legpos::Union{String, Int64}="best", legcol::Int64=1)

  # Start plot
  fig, ax = subplots(figsize=figsiz)

  # set colour scheme

  # Plot data
  xdata = plot_list[1].x; ystack = []; α = []
  for data in plot_list
    if data.x ≠ xdata
      println("Warning! X data in $(data.label) differs from x data in first item of stack.")
      println("Data ignored in stack.")
      continue
    end
    push!(ystack,data.y); push!(α,data.alpha)
  end
  if alpha > 0
    α = alpha
  elseif mean(α) > 0
    α = mean(α)
  else
    α = 1
  end
  clr = []
  for (i, plt) in enumerate(plot_list)
    if cs=="own"
      try c = lc[i]
      catch
        println("Warning! Colour not defined for data $i. Using default.")
        c, lstyle = sel_ls("default",lc=i)
      end
    elseif cs==""
      if plt.colour == nothing
        println("Warning! Colour not defined for data $i. Using default.")
        c, lstyle = sel_ls("default",lc=i)
      else c = plt.colour
      end
    else
      c, lstyle = sel_ls(cs,lc=i)
    end
    push!(clr,c)
  end
  ax[:stackplot](xdata,ystack,colors=clr,alpha=α)
  ax=draw_boundaries(boundaries,xdata,ystack,clr,lt,ax)

  # Define logscales
  xlim = get_lim(xlims,ax,:get_xlim)
  ax[:set_xlim](xlim)
  ylim = get_lim(ylims,ax,:get_ylim)
  ax[:set_ylim](ylim)
  if contains(logscale,"x")
    xlim = set_log2(plot_list,xlim,:x)
    ax[:set_xlim](xlim)
    ax[:set_xscale]("log")
  end
  if contains(logscale,"y")
    ylim = set_log2(plot_list,ylim,:y)
    ax[:set_ylim](ylim)
    ax[:set_yscale]("log")
  end

  # Set plot title
  ax[:set_title](ti, fontsize=fntsiz+ti_offset)

  # Generate axes labels and legend, define axes label/tick colours
  ax[:set_xlabel](xlabel,fontsize=fntsiz+ax_offset)
  ax[:set_ylabel](ylabel,fontsize=fntsiz+ax_offset)

  if legpos ≠ "None" && any([p.label≠"" for p in plot_list])
    lbl = [p.label for p in plot_list]
    ax[:legend](lbl, fontsize=fntsiz+leg_offset, loc=legpos, ncols=legcol)
  end

  # Set ticks and optional minor ticks
  if typeof(plot_list[1].x) ≠ Vector{DateTime}  if maj_xticks > 0
    xint = collect(ax[:get_xlim]()[1]:maj_xticks:ax[:get_xlim]()[2])
    ax[:set_xticks](xint)
  end  end
  if maj_yticks > 0
    yint = collect(ax[:get_ylim]()[1]:maj_yticks:ax[:get_ylim]()[2])
    ax[:set_yticks](yint)
  end
  if mticks == "on"
    minorticks_on()
  else
    minorticks_off()
  end
  # Set minor x ticks
  if typeof(plot_list[1].x) ≠ Vector{DateTime}
    if min_xticks > 0
      mx = matplotlib[:ticker][:MultipleLocator](min_xticks)
      ax[:xaxis][:set_minor_locator](mx)
    end
  elseif min_xticks isa Number
    min_xticks = [6,12,18]
  end
  # Set minor y ticks
  for i = 1:length(min_yticks)
    if min_yticks > 0
      my = matplotlib[:ticker][:MultipleLocator](min_yticks)
      ax[:yaxis][:set_minor_locator](my)
    end
  end
  # Format ticks and frame
  Mtlen = tsc[1]⋅frw
  mtlen = tsc[2]⋅frw
  ax[:tick_params]("both", which="both", direction="in", top="on", right="on",
    labelsize=fntsiz, width=frw)
  ax[:grid](linestyle=":", linewidth = frw)
  ax[:spines]["bottom"][:set_linewidth](frw)
  ax[:spines]["top"][:set_linewidth](frw)
  ax[:spines]["left"][:set_linewidth](frw)
  ax[:spines]["right"][:set_linewidth](frw)
  if typeof(plot_list[1].x) ≠ Vector{DateTime}
    ax[:tick_params]("both", which="major", length=Mtlen)
    ax[:tick_params]("both", which="minor", length=mtlen)
  else
    ax[:set_xlim](xmin=plot_list[1].x[1], xmax=plot_list[1].x[end])
    if maj_xticks isa Vector
      majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b, %H:%M")
    else
      majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b")
    end
    minorformatter = matplotlib[:dates][:DateFormatter]("")
    majorlocator = matplotlib[:dates][:HourLocator](byhour=maj_xticks)
    minorlocator = matplotlib[:dates][:HourLocator](byhour=min_xticks)
    ax[:xaxis][:set_major_formatter](majorformatter)
    ax[:xaxis][:set_minor_formatter](minorformatter)
    ax[:xaxis][:set_major_locator](majorlocator)
    ax[:xaxis][:set_minor_locator](minorlocator)
    fig[:autofmt_xdate](bottom=0.2,rotation=-30,ha="left")
  end

  tight_layout()

  # Return PyPlot data
  return fig, ax
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
    rename_DF(pltdata, DFnames, err)

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
function rename_DF(pltdata, DFnames, err)
  if err == "None"
    names!(pltdata,DFnames[1:2])
  elseif startswith(err,"pm") && endswith(err,"x")
    names!(pltdata,DFnames[[1,2,5,6]])
  elseif endswith(err,"x")
    names!(pltdata,DFnames[[1,2,5]])
  elseif startswith(err,"pm") && endswith(err,"y")
    names!(pltdata,DFnames[1:4])
  elseif endswith(err,"y")
    names!(pltdata,DFnames[1:4])
  elseif startswith(err,"pm")
    names!(pltdata,DFnames)
  else
    names!(pltdata,DFnames[[1,2,3,5]])
  end

  return pltdata
end #function rename_DF


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
setup_axes(plot_list, twinax, ylab, logscale, xlims, ylims, Mytint, nmyt, cs, axcol, legpos, legcol)

If `twinax` is an array split `plot_list` into 2 datasets based on the indices
in `twinax`. Assure all paramters concerning y data are arrays of length 2, if
a second axis is allowed or length 1 otherwise. Return `plot_list` as array of
2 or 1 distinct lists an the adjusted parameters.
"""
function setup_axes(plot_list, twinax, ylab, logscale, xlims, ylims, Mytint, nmyt,
                    ptype, lc, lt, pt, alpha, axcol, legpos, legcol)

  # Set flag for second axis, asume no second axis
  ax_2 = false
  ax2 = nothing
  xlim = Any[]; ylim = Any[]

  [p.alpha = alpha for p in plot_list if alpha > 0]

  # Set up second axis, if the twinax array is defined
  plt1 = PlotData[]; plt2 = PlotData[]
  cl1 = String[]; dt1 = []; mt1 = []
  cl2 = String[]; dt2 = []; mt2 = []
  if !isempty(twinax)
    # Set flag true and define 2nd axis in PyPlot
    ax_2 = true
    ax2 = twinx()
    # Check correct input of twinax
    if length(twinax) ≠ length(plot_list)
      println("Array `twinax` must have the same length as array `phot_list`.")
      println("Script stopped."); exit()
    end

    # Assign data to the axes
    for i = 1:length(twinax)
      if twinax[i] == 1
        push!(plt1,plot_list[i])
        # check whether an array of definitions or a single definition exists
        # and assign either the current array entry to the respective axis
        # or the default value for line colour, type and marker type
        if ptype == "own" || ptype == "mix"
          if !isempty(lt) && (lt[i] isa Tuple || lt[i] isa Vector || lt[i] isa String)
            push!(dt1,lt[i])
          else
            push!(dt1,lt)
          end
          if pt isa Vector
            push!(mt1,pt[i])
          else
            push!(mt1,pt)
          end
          if lc isa Vector
            push!(cl1,lc[i])
          else
            push!(cl1,lc)
          end
        end
      else
        push!(plt2,plot_list[i])
        # check whether an array of definitions or a single definition exists
        # and assign either the current array entry to the respective axis
        # or the default value for line colour, type and marker type
        if ptype == "own" || ptype == "mix"
          if !isempty(lt) && (lt[i] isa Tuple || lt[i] isa Vector || lt[i] isa String)
            push!(dt2,lt[i])
          else
            push!(dt2,lt)
          end
          if pt isa Vector
            push!(mt2,pt[i])
          else
            push!(mt2,pt)
          end
          if lc isa Vector
            push!(cl2,lc[i])
          else
            push!(cl2,lc)
          end
        end
      end
    end

    # Make sure, all parameters for both axes are arrays of length 2
    if logscale isa String logscale = String[logscale, logscale]  end
    if axcol isa String axcol = String[axcol, axcol]  end
    if Mytint isa Int64 Mytint = Int64[Mytint, Mytint]  end
    if nmyt isa Int64 nmyt = Int64[nmyt, nmyt]  end
    if xlims == nothing
      xlim = [[nothing, nothing],[nothing, nothing]]
    elseif xlims isa Tuple
      xl = [xlims[1],xlims[2]]
      xlim = [xl, xl]
    elseif xlims isa Array
      xlim = [[xlims[1][1], xlims[1][2]], [xlims[2][1], xlims[2][2]]]
    end
    if ylims == nothing
      ylim = [[nothing, nothing],[nothing, nothing]]
    elseif ylims isa Tuple
      yl = [ylims[1],ylims[2]]
      ylim = [yl, yl]
    elseif ylims isa Array
      ylim = [[ylims[1][1], ylims[1][2]], [ylims[2][1], ylims[2][2]]]
    end
    if !isa(legpos, Array) legpos = [legpos, legpos]  end
    if legcol isa Int64 legcol = [legcol, legcol]  end
  else
    # Assign all data to axis 1, if no second axis
    plt1 = plot_list
    # check whether an array of definitions or a single definition exists
    # and assign either the current array entry to the respective axis
    # or the default value for line colour, type and marker type
    if ptype == "own" || ptype == "mix"
      if !isempty(lt) && (lt[1] isa Tuple || lt[1] isa Vector || lt[1] isa String)
        dt1 = lt
      else
        dt1 = [lt for i in plot_list]
      end
      if pt isa Vector
        mt1 = pt
      else
        mt1 = [pt for i in plot_list]
      end
      if lc isa Vector
        cl1 = lc
      else
        cl1 = [lc for i in plot_list]
      end
    end
    # If no 2nd axis, make sure, all parameters for both axes are arrays of length 1
    # and not single parameters (numbers, strings...)
    if logscale isa String logscale = String[logscale]  end
    if axcol isa String axcol = String[axcol]  end
    if Mytint isa Int64 Mytint = Int64[Mytint]  end
    if nmyt isa Int64 nmyt = Int64[nmyt]  end
    if xlims == nothing
      xlim = [[nothing, nothing]]
    elseif xlims isa Tuple
      xlim = [[xlims[1],xlims[2]]]
    end
    if ylims == nothing
      ylim = [[nothing, nothing]]
    elseif ylims isa Tuple
      ylim = [[ylims[1],ylims[2]]]
    end
    if !isa(legpos, Array) legpos = [legpos]  end
    if legcol isa Int64  legcol = [legcol]  end
  end
  if ylab isa String  ylab = String[ylab]  end
  if ax_2
    plt = [plt1, plt2]
    lc = [cl1, cl2]
    lt = [dt1, dt2]
    pt = [mt1, mt2]
  else
    plt = [plt1]
    lc = [cl1]
    lt = [dt1]
    pt = [mt1]
  end

  # Return adjusted data
  return plt, ax2, ax_2, ylab, logscale, xlim, ylim, Mytint, nmyt, lc, lt, pt, axcol, legpos, legcol
end #function setup_axes


"""
    set_cs(plt, ptype, cs, lc, lt, pt, ax_2)

Set colour scheme for `plt` holding `PlotData` to `cs` from
function `sel_ls`. If colour scheme is set to `"own"` or `"mix"` for any array
element in `plt`, `lc`, `lt`, and `pt` must hold valid definitions of line colours,
line and point (marker) types at this array position. The boolean `ax_2` is needed
to specify, whether a second axis is used.
"""
function set_cs(plt, ptype, cs, lc, lt, pt, ax_2)
  if ptype == "own"
    return plt
  elseif ptype == "mix"
    for i = 1:length(plt[1])
      plt[1][i].colour = lc[1][i]
      plt[1][i].dashes = lt[1][i]
      plt[1][i].marker = pt[1][i]
    end
    if ax_2  for i = 1:length(plt[2])
      plt[2][i].colour = lc[2][i]
      plt[2][i].dashes = lt[2][i]
      plt[2][i].marker = pt[2][i]
    end end
  elseif cs≠""
    for i = 1:length(plt[1])
      if cs isa String
        cl, dt, mt = sel_ls(cs, lc=i, lt=i, pt=i)
      else
        cl, dt, mt = sel_ls(cs[1], lc=i, lt=i, pt=i)
      end
      plt[1][i].colour = cl
      if ptype == "line"
        plt[1][i].dashes = dt
        plt[1][i].marker = "None"
      elseif ptype == "scatter"
        plt[1][i].dashes = (0,1)
        plt[1][i].marker = mt
      elseif ptype == "both"
        plt[1][i].dashes = dt
        plt[1][i].marker = mt
      end
    end
    if ax_2  for i = 1:length(plt[2])
      if cs isa String
        j = i + length(plt[1])
        cl, dt, mt = sel_ls(cs, lc=j, lt=j, pt=j)
      else
        cl, dt, mt = sel_ls(cs[2], lc=i, lt=i, pt=i)
      end
      plt[2][i].colour = cl
      if ptype == "line"
        plt[2][i].dashes = dt
        plt[2][i].marker = "None"
      elseif ptype == "scatter"
        plt[2][i].dashes = (0,1)
        plt[2][i].marker = mt
      elseif ptype == "both"
        plt[2][i].dashes = dt
        plt[2][i].marker = mt
      end
    end end
  else
    for i = 1:length(plt[1])
      cl, dt, mt = sel_ls(lt=i, pt=i)
      if ptype == "line"
        plt[1][i].dashes = dt
        plt[1][i].marker = "None"
      elseif ptype == "scatter"
        plt[1][i].dashes = (0,1)
        plt[1][i].marker = mt
      elseif ptype == "both"
        plt[1][i].dashes = dt
        plt[1][i].marker = mt
      end
    end
    if ax_2  for i = 1:length(plt[2])
      j = i + length(plt[1])
      cl, dt, mt = sel_ls(lt=j, pt=j)
      if ptype == "line"
        plt[2][i].dashes = dt
        plt[2][i].marker = "None"
      elseif ptype == "scatter"
        plt[2][i].dashes = (0,1)
        plt[2][i].marker = mt
      elseif ptype == "both"
        plt[2][i].dashes = dt
        plt[2][i].marker = mt
      end
    end end
  end

  return plt
end #function set_cs


"""
    plt_DataWithErrors(plt, ax, offset)

For each `PlotData` in array `plt` (and `ax`), retrieve the errors in the `PlotData` and
plot according to specifications. For error bars of markers, the standard cap size is `3`,
which can be adjusted by an `offset` (positive or negative number to be added to cap size).

Returns `fig` and `ax` (array of axes) for further plot modifications or printing.
"""
function plt_DataWithErrors(plt, ax, offset)
  # Redefine errors for error bars
  xerr = redef_err(plt,:x,:xlerr,:xuerr)
  yerr = redef_err(plt,:y,:ylerr,:yuerr)
  # Loop over graph data and plot each data according to its type and format
  # defined by the struct PlotData
  for i = 1:length(plt)
    if plt[i].yuerr ≠ nothing && plt[i].marker == "None"
      p=ax[:plot](plt[i].x, plt[i].y, lw = plt[i].lw,
          dashes=plt[i].dashes, color=plt[i].colour, label=plt[i].label)
      plt[i].colour = p[1][:get_color]()
      ax[:fill_between](plt[i].x, plt[i].ylerr, plt[i].yuerr,
          color=plt[i].colour, alpha=0.2)
    elseif plt[i].xuerr ≠ nothing && plt[i].yuerr ≠ nothing
      if !isempty(plt[i].dashes) && plt[i].dashes[1] == 0
        ax[:errorbar](plt[i].x, plt[i].y, xerr=[xerr[i][:lower], xerr[i][:upper]],
          yerr=[yerr[i][:lower], yerr[i][:upper]], ls = "None",
          marker=plt[i].marker, color=plt[i].colour,
          label=plt[i].label, capsize=3+offset, alpha=plt[i].alpha)
      else
        ax[:errorbar](plt[i].x, plt[i].y, xerr=[xerr[i][:lower], xerr[i][:upper]],
          yerr=[yerr[i][:lower], yerr[i][:upper]], lw = plt[i].lw,
          marker=plt[i].marker, dashes=plt[i].dashes, color=plt[i].colour,
          label=plt[i].label, capsize=3+offset, alpha=plt[i].alpha)
      end
    elseif plt[i].yuerr ≠ nothing
      if !isempty(plt[i].dashes) && plt[i].dashes[1] == 0
        ax[:errorbar](plt[i].x, plt[i].y, yerr=[yerr[i][:lower], yerr[i][:upper]],
          ls = "None", marker=plt[i].marker,
          color=plt[i].colour, label=plt[i].label, capsize=3+offset, alpha=plt[i].alpha)
      else
        ax[:errorbar](plt[i].x, plt[i].y, yerr=[yerr[i][:lower], yerr[i][:upper]],
          lw = plt[i].lw, marker=plt[i].marker, dashes=plt[i].dashes,
          color=plt[i].colour, label=plt[i].label, capsize=3+offset, alpha=plt[i].alpha)
      end
    elseif plt[i].xuerr ≠ nothing
      if !isempty(plt[i].dashes) && plt[i].dashes[1] == 0
        ax[:errorbar](plt[i].x, plt[i].y, xerr=[xerr[i][:lower], xerr[i][:upper]],
          ls = "None", marker=plt[i].marker,
          color=plt[i].colour, label=plt[i].label, capsize=3+offset, alpha=plt[i].alpha)
      else
        ax[:errorbar](plt[i].x, plt[i].y, xerr=[xerr[i][:lower], xerr[i][:upper]],
          lw = plt[i].lw, marker=plt[i].marker, dashes=plt[i].dashes,
          color=plt[i].colour, label=plt[i].label, capsize=3+offset, alpha=plt[i].alpha)
      end
    else
      if !isempty(plt[i].dashes) && plt[i].dashes[1] == 0
        ax[:scatter](plt[i].x, plt[i].y,lw = plt[i].lw, marker=plt[i].marker,
          color=plt[i].colour, label=plt[i].label, alpha=plt[i].alpha)
      else
        ax[:plot](plt[i].x, plt[i].y,lw = plt[i].lw, marker=plt[i].marker,
          dashes=plt[i].dashes, color=plt[i].colour, label=plt[i].label,
          alpha=plt[i].alpha)
      end
    end
  end

  return plt, ax
end #function plt_DataWithErrors


"""
    redef_err(plt,val,low,high)

Recalculate errors relative to actual values (rather than absolute values in plot)
for use with error bars of markers. The function generates an array of DataFrames
with the revised errors in the columns `:upper` and `:lower`.
Data is generated from the array `plt` with the columns `val` with the measured/modelled
value and the columns with names `low` and `high` holding the
actual error values (rather than values relative to measured/modelled value.)
"""
function redef_err(plt,val,low,high)
  err = []
  for i = 1:length(plt)
    # Recalculate errors for error bars, if markers are plotted
    if plt[i].marker ≠ "None" && getfield(plt[i],high) ≠ nothing
      push!(err, DataFrame(upper = getfield(plt[i],high) .- getfield(plt[i],val),
            lower = getfield(plt[i],val) .- getfield(plt[i],low)))
    else
      # Otherwise, save errors as they are und column names `upper` and `lower`
      push!(err, DataFrame(upper = Float64[], lower = Float64[]))
    end
  end

  return err
end #function redef_err


"""
    set_log(plt, ax, lims, logscale, x, val, low, high, cmd1, cmd2)

Set x- and/or y-axis to log scale, if the string `logscale` consists of `"x"` or
`"y"` or `"xy"`. Adjust minimum/maximum values of the respective axis for logscale.
"""
function set_log(plt, ax, lims, logscale, x, val, low, high, cmd1, cmd2)
  if contains(logscale,x)
    xmin = Float64[]; xmax = Float64[]
    for p in plt
      if getfield(p, high) == nothing
        push!(xmin, minimum(getfield(p, val)))
        push!(xmax, maximum(getfield(p, val)))
      else
        push!(xmin, minimum(getfield(p, low)))
        push!(xmax, maximum(getfield(p, high)))
      end
    end
    xmin = 10^floor(log10(minimum(xmin)))
    xmax = 10^ceil(log10(maximum(xmax)))
    lims = [xmin, xmax]
    ax[cmd1](lims)
    ax[cmd2]("log")
    return ax, lims
  end

  return ax, lims
end #function set_log


"""
    logBounds(xlims_orig, ylims_orig, ax)

Set x and y limits for current plot's axes layers in array `ax`. Replace nothings with limits
from current axis layer, otherwise use predefined values from `xlims_orig` and
`ylims_orig`.
"""
function logBounds(xlims_orig, ylims_orig, ax)
  xlims = []; ylims = []
  for l = 1:length(xlims_orig)
    xlims_orig[l][1] == nothing ? xl = ax[l][:get_xlim]()[1] : xl = xlims_orig[l][1]
    xlims_orig[l][2] == nothing ? xu = ax[l][:get_xlim]()[2] : xu = xlims_orig[l][2]
    push!(xlims,[xl, xu])
    ylims_orig[l][1] == nothing ? yl = ax[l][:get_ylim]()[1] : yl = ylims_orig[l][1]
    ylims_orig[l][2] == nothing ? yu = ax[l][:get_ylim]()[2] : yu = ylims_orig[l][2]
    push!(ylims,[yl, yu])
  end
  return xlims, ylims
end


function get_lim(lim,ax,cmd)
  if lim == nothing
    lim = ax[cmd]()
  else
    lims = []
    for (i,l) in enumerate(lim)
      if l == nothing
        push!(lims,ax[cmd]()[i])
      else
        push!(lims,l)
      end
    end
    lim = (lims[1],lims[2])
  end

  return lim
end


function set_log2(plt,lim,f)
  mn = lim[1]; mx = lim[2]
  if mn ≤ 0 && mx > 0
    println("Warning! Values must be strictly positive for logarithmic plotting.")
    println("Negative values and zeros ignored.")
    data = Float64[]
    for d in plt
      if f==:x      data = vcat(data,d.x)
      elseif f==:y  data = vcat(data,d.y)
      end
    end
    data=sort!(data)
    mn=data[findfirst(data.>0.)]
  elseif mx ≥ 0 && mn < 0
    println("Warning! Values must be strictly negative for logarithmic plotting.")
    println("Positive values and zeros ignored.")
    data = Float64[]
    for d in plt
      if f==:x      data = vcat(data,d.x)
      elseif f==:y  data = vcat(data,d.y)
      end
    end
    data=sort!(data)
    mx=data[findlast(data.<0.)]
  end
  try mn = 10^floor(log10(mn))
  catch
    mn = -10^floor(log10(abs(mn)))
  end
  try mx = 10^ceil(log10(mx))
  catch
    mx = -10^ceil(log10(abs(mx)))
  end
  lim = (mn, mx)
end


function draw_boundaries(alpha,xdata,ydata,colours,lt,ax)
  if alpha==0  return ax  end

  ylines=[ydata[1]]
  for i = 2:length(ydata)  push!(ylines,sum(ydata[1:i]))  end
  [plot(xdata,ylines[i], color=colours[i], dashes = lt) for i = 1:length(ylines)]

  return ax
end

end #module pyp
