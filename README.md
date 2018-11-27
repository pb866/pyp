Julia package pyp for plotting with PyPlot
==========================================

This package provides some functions for easy loading of data and plotting.

Installation
------------

Install into your preferred environment (`v1.0` for general Julia env or
to any project) with

```julia
using Pkg
Pkg.add("https://github.com/pb866/pyp.git")
Pkg.instantiate()
using pyp
```

or go to the package manager in your preferred environment by typing `]`
followed by

```julia
add https://github.com/pb866/pyp.git
instantiate
```

Leave the package manager with the backspace key and `import pyp` or
`using pyp`.


Available functions
-------------------

See also the documentation of each functions by typing `?<function name>`
or `?pyp.<function name>` in the julia console.


### load_PlotData

Loads data from a `DataFrame` (`plotdata`) into the data type `PlotData` (see [PlotData](#mutable-struct-plotdata)).
Each `PlotData` holds the properties of one dataset, which can be combined
in a graph. The below properties can be defined.


#### Error margins

Up to 6 columns can be specified for one dataset: x, y, lower y error,
upper y error, lower x error, upper x error with the column symbols `[:x, :y, :ylerr, :yuerr, :xlerr, :xuerr]`.

The keyword `err` can be used to load data with errors in a certain format.
In `PlotData`, the error are saved in fields with the actual values of the
bounds rather than spans.
For data with 3 columns (x, y, and error), `err` can be set as:
- `"rangex"` or `"rangey"`: ± error range
- `"percentx"` or `"percenty"`: ±err ⋅ x or y
- `"factorx"` or `"factory"`: x or y ⋅ 1/err and ⋅ err, respectively

For data with 4 columns (x, y, lower error, upper error), `err` can be set
to `"pmrangex"`, `"pmrangey"`, `"pmfactorx"`, `"pmfactory"`, `"pmpercentx"`, and `"pmpercenty"`
which is like the above, but with different upper and lower bounds.
Moreover, `err` can be set to `"rangexy"`, `"percentxy"`, and `"factorxy"`,
which is like the above, but with values set for x and y simultaneously.

For data with 6 columns (x, y, lower error, upper error), `err` can be set
to `"pmrangexy"`, `"pmfactorxy"`, and `"pmpercentxy"`, which is like the above,
but with different bounds for both x and y.

Furthermore, `"valuex"`, `"valuey"`, and `"value"` can be used for `err`
to directly set the values for the upper and lower bounds for x, y or x and y.

For the data, always the first n columns are used unless you redefine the
column names with the keyword `renameDF`. Default names are
`[:x, :y, :ylerr, :yuerr, :xlerr, :xuerr]`. If you redefine the names,
you have to define __the whole vector of symbols__ rather than a single
symbol. For values with equal errors, like `"valuey"`, redefine the lower
bound symbol.

Errors get calucalated and are added to the `PlotData` fields `xuerr`,
`xlerr`, `yuerr`, and `ylerr`; `nothing` is used for data with no errors.


#### Line and marker styles

You can set the marker style using `PyPlot` keywords with the keyword argument `pt`, e.g. `"s"` for squares. `"None"` is used for no markers,
i.e. line plots. You are allowed to enter strings or integers.

The line style is set with the keyword argument `lt` using PyPlot nomenclature
for `dashes`. Use Tuples of integers or floats to define the on/off point series or arrays for on/off series with different alternating values. Use
`"None"` or tuples/arrays with odd values set to `0` for no lines
(i.e., scatter plots) and `[]` (empty array). Omit keyword for solid lines or use `[]` (empty array).

The keyword argument `lc` sets the colour of the plot data. Use `PyPlot`
colour names or RGB code (`#XXXXXX`) to set the colours.

The keyword `lw` defines the linewidth in points. You can use integers and floats.

`alpha` set the opaqueness, from `0` (completely transparent/invisible) to
`1` (completely opaque).


#### Data scaling and labelling

You can scale y data with the keyword `SF` as shown for function [rd_data](#rd_data).
Data labels are defined by `String`s with the keyword `label`.


### Mutable struct PlotData

The `PlotData` data type holds fields to define datasets and plotting
parameters for easy plotting with `PyPlot`. It is mutable, so fields can
be redefined by `PlotData.var = ...`.

The following fields exist; `x` and `y` have to be defined. For the remaining
fields exist default values.
- `x::Union{Vector{Int64},Vector{Float64},Vector{DateTime}}`
- `y::Union{Vector{Int64},Vector{Float64}}`
- `xuerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing`
- `xlerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing`
- `yuerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing`
- `ylerr::Union{Nothing,Vector{Int64},Vector{Float64}}=nothing`
- `label::Union{Nothing,AbstractString}=""`
- `marker::Union{AbstractString,Int64}="None"`
- `dashes::Union{Tuple{Float64,Float64},Tuple{Int64,Int64},Vector{Float64},Vector{Int64},Vector{Any},String}=Float64[]`
- `colour::Union{Nothing,AbstractString}=nothing`
- `lw::Number=1.4`
- `alpha::Number=1`


### plot_data

Function to generate formatted `PyPlot` line and/or scatter plots from a
list of `PlotData`. Returns `fig, ax1[, ax2]` for further formatting of
the plots.

The following sections detail formatting options by keyword arguments.


#### General settings

- `ti` or `title`: sets the title (_default:_ `""`; no title)
- `xlabel`: sets x label (_default:_ `"model time / hours"`)
- `ylabel`: sets y label (_default:_ `"concentration / mlc cm\$^{-3}\$ s\$^{-1}\$"`)
- `figsize`: (x, y) tuple with size as used in `PyPlot` (_default:_ `1`)
- `fontsize`: font size in pt as used by `PyPlot` (_default:_ `12`)
- `framewidth`: linewidth of the frame (axes surrounding the plot) in pt (_default:_ `(6, 4)`)
- `ticksize`: tuple of scaling factors for the length of (major, minor) ticks in relation to their line width (_default:_ `(4.5,2.5)`)
- `cap_offset`, `ti_offset`, `ax_offset`, `leg_offset`: offsets for cap size of error bars in scatter plots, title font size, axes label font size, and legend font size, respectively, with values added or substracted from the regular font size (_defaults:_ `0`, `4`, `2`, and `0`)
- `xlims`/`ylims`: define tuples with minimum and maximum values the x and y axes limits (uses `PyPlot`'s limits as _default_)


#### Formatting graphs

General formats for graphs and labels for the data are defined by the
`PlotData`, which is handed over as a list (`vararg`).

With the keyword `plot_type`, the following plots can be chosen, which
overwrites settings in `PlotData`:
- `"line"`: line plots with solid lines and no markers
- `"scatter"`: scatter plots with no lines and all square markers
- `"both"`: square markers connected by solid lines

Leave blank to use settings from `PlotData`.

Predefined colour schemes from function `sel_ls` can be chosen
(`"default"` with a range of various colours, `source` with mostly blue
shades, `sink` with mostly red shades) with the keyword `cs`, `color`,
`colour`, `colorscheme` or `colourscheme`. Moreover, different dash and
marker types are provided. These settings overwrite any settings in `PlotData` and from `plot_type`.

Line types can be redefined by feeding `PyPlot` commands for `dashes` to
the keyword `lt`, `linetype`, `linestyle` or `dashes` (see also above).
To redefine `PyPlot`'s `marker`, use keyword `pt`, `mt`, or `marker`.
Both settings overwrite the previous ones.
Keyword `lc`, `linecolor`/`linecolour` or `color`/`colour` can be used to
define the colour of the graph with `PyPlot`'s named colours or an RGB code.
For all of the above settings you can either set a default value for all
data by using a single value or define each dataset in the `PlotData` list
by specifying a vector of length of the `PlotData` list with separate settings.

With `alpha` you can specify the opaqueness of the graphs. Only a default
value is possible and no separate settings in a vector.

#### Advanced settings

- `mticks`: switches minor ticks on/off (_default:_ `"on"`)
- `min_xticks`, `maj_xticks`, `min_yticks`, `maj_yticks`: sets the intervals of minor/major x/y ticks (uses `PyPlot`'s values as _default_)
- `logscale`: sets x, y or both axes to logarithmic scale if provided with `"x"`, `"y"` or `"xy"`
- `logremove`: sets all negative or positve values from a dataset to `0.0`, if set to `"neg"` or `"pos"` and `logscale` is switched on for that axis
  (It is important to have strictly positive or negative datasets for log plots.)
- `legpos`, `legloc` or `loc` defines the position of the legend using `PyPlot`'s keywords or integers (_default:_ `"best"`)
- `legcolumns`: sets the number of columns in the legend (_default:_ `1`)

You can define a secondary y axis, which is shown to the right.
Use the keyword argument `twinax` and define a vector of integers with `1`
and `2` of length of the `PlotData` list. `1`'s assign data in the `PlotData` list to the left axis, `2`'s to the right.
For the y-axis related properties `ylabel`, `logscale`,
`cs` (`colourscheme`), `ylimes`, `min_yticks`, and `maj_yticks`
you can either specify a global value for both axes or you can define an
array `[setting1, setting2]` with different values for the left and right
axis, respectively. For other settings like `lt` or `pt` you can still
define the array as for a single y axis and the settings get automatically
distributed to the correct axis.
Furthermore, you can specify the colour of the y axis labels with the
keyword `axcolour` as global value or array to better distinguish two
different datasets on two y axes.

### sel_ls

Function `sel_ls(cs; kwargs)` allows you to choose colours, line and marker
styles from predefined schemes. Currently `"default"`, `"source"`, and
`"sink"` are available with varying colours, mostly blue and red shades,
respectively.
Furthermore, different dash and marker types are available. Choose line
colours with keyword `lc`, dash types with `lt`, and marker types with `pt`.
You can select single values by handing over an integer, or ask for arrays
of settings by handing over arrays or ranges. Array sizes, do not have to
match, you can as for four colour setting, 2 line type settings, and a default marker type setting.  
Returns a single value or array of values in the order `lc`, `lt`, `pt`


Version history
===============

Version 0.1.1
-------------
- Julia dependence corrected to ≥0.7

Version 0.1.0
-------------
- initial version to read simple numeric text file with a certain format
- function to load new type PlotData
- function to use PyPlot to plot PlotData as scatter and/or line plots
