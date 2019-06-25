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
import filehandling; const fh = filehandling
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
  x::Union{Vector{T} where T<:Dates.TimeType, Vector{T} where T<:Real}
  y::Vector{T} where T<:Real
  xuerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  xlerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  yuerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  ylerr::Union{Nothing,Vector{T} where T<:Real}=nothing
  label::AbstractString=""
  marker::Union{String,Int64}="None"
  dashes::Union{Tuple{Real,Real},Vector{T} where T<:Real}=Int64[]
  colour::Union{Nothing,String,Symbol}=nothing
  lw::Real=1.4
  alpha::Real=1
end

include("public.jl")
include("datamanipulation.jl")
include("plot.jl")

end #module pyp
