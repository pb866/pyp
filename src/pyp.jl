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
- plotDataWithErrors
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

### Error handling

# Define Error type
struct DefinitionError <: Exception
  msg::String
  data
end
# Format Error message
Base.showerror(io::IO, e::DefinitionError) =
  print(io, typeof(e), "\n", e.msg, "\033[0m", e.data)
# Define default message
DefinitionError(data) =
  DefinitionError("`Vector` must be of length(2) and `twinax` must be set; length: ", data)


### NEW TYPES
const lt_type = Union{String,Tuple{Real,Real},Int,UnitRange{Int},Vector}
const lc_type = Union{Nothing,String,Symbol,Int,UnitRange{Int},Vector}
const limtype = Union{Nothing,Tuple{Union{Nothing,<:Real},Union{Nothing,<:Real}}}


"""
    par.@with_kw mutable struct PlotData

Mutable struct holding information of individual graphs in a plot
that can be addressed with keyword arguments for the following fields:
- `x`/`y`: x/y data for plotting, x/y must be vectors of `Real`s (or `DateTime`s for x)
  of same length
- `x_lower`/`x_upper`/`y_lower`/`y_upper`: upper and lower errors for x and y data,
  must either be `nothing` or vector of `Real`s of same length as `x`/`y`
  (default: `nothing`)
- `label`: `AbstractString` with label for legend (default: `""` – empty String, no entry)
- `marker`: `String` or `Int` specifying `PyPlot` marker type (default: `"None"` – no markers)
- `dashes`: tuple or vector of `Real`s specifying an 'on-off' sequence for line types
  (default: `[]` – empty array = solid lines); use `"None"` for no lines
- `colour`: `String` or `Symbol` with `PyPlot` colour name or String with RGB code
  (default: `nothing` – use `PyPlot` default colours)
- `lw`: `Real` to specify linewidth
- `alpha`: `Real` to specify opaqueness
"""
@par.with_kw mutable struct PlotData
  x::Union{Vector{<:Dates.TimeType}, Vector{<:Real}}
  y::Vector{<:Real}
  x_lower::Vector{<:Real}=Real[]
  x_upper::Vector{<:Real}=Real[]
  y_lower::Vector{<:Real}=Real[]
  y_upper::Vector{<:Real}=Real[]
  label::AbstractString=""
  marker::Union{String,Int}="None"
  dashes::lt_type=Real[]
  colour::Union{Nothing,String,Symbol}=nothing
  lw::Real=1.4
  alpha::Real=1

  function PlotData(x::Union{Vector{<:Dates.TimeType}, Vector{<:Real}},
    y::Vector{<:Real}, x_lower::Vector{<:Real}=Real[], x_upper::Vector{<:Real}=Real[],
    y_lower::Vector{<:Real}=Real[], y_upper::Vector{<:Real}=Real[],
    label::AbstractString="", marker::Union{String,Int}="None",
    dashes::lt_type=Real[], colour::Union{Nothing,String,Symbol}=nothing,
    lw::Real=1.4, alpha::Real=1)
    if length(x) == 0 || length(y) == 0
      throw(ErrorException("No empty `x` and `y` data allowed"))
    end
    if length(x) ≠ length(y)
      throw(ErrorException("`x` and `y` must be vectors of same length"))
    end
    if !isempty(x_lower) && length(x) ≠ length(x_lower)
      throw(ErrorException("`x` and `x_lower` must be vectors of same length"))
    end
    if !isempty(x_upper) && length(x) ≠ length(x_upper)
      throw(ErrorException("`x` and `x_upper` must be vectors of same length"))
    end
    if !isempty(y_lower) && length(x) ≠ length(y_lower)
      throw(ErrorException("`y` and `y_lower` must be vectors of same length"))
    end
    if !isempty(y_upper) && length(x) ≠ length(y_upper)
      throw(ErrorException("`y` and `y_upper` must be vectors of same length"))
    end
    if dashes isa Vector && isempty(dashes)
      dashes = Real[]
    end
    if dashes isa Vector && !(typeof(dashes) <: Vector{<:Real})
      throw(TypeError(:PlotData, "", Vector{T} where T<:Real, typeof(dashes)))
    end
    new(x,y,x_lower,x_upper,y_lower,y_upper,label,marker,dashes,colour,lw,alpha)
  end #function PlotData
end


"""


"""
@par.with_kw mutable struct kwargs
  ti::AbstractString=""
  xlabel::AbstractString="model time / hours"
  ylabel::Union{AbstractString, Vector{<:AbstractString}}=
    L"concentration / mlc cm$^{-3}$ s$^{-1}$"
  lt::lt_type=""
  pt::Union{String,Int,UnitRange{Int},Vector}=""
  lc::lc_type=""
  lw::Real=1.4
  cs::Union{String,Vector{<:String}}=""
  plottype::String="default"
  legpos::Union{String,Int,Tuple{Real,Real}}="best"
  legcols::Int=1
  xlim::limtype=nothing
  ylim::Union{limtype,Vector}=nothing
  mticks::Bool=true
  Xticks::Union{Real,Vector{<:Int}}=-1
  xticks::Union{Real,Vector{<:Int}}=-1
  Yticks::Union{Real,Vector{<:Real}}=0
  yticks::Union{Real,Vector{<:Real}}=0
  timeformat::String=""
  timescale::String="days"
  major_interval::Int=0
  minor_interval::Int=0
  logscale::Union{String,Vector{<:String}}=""
  logremove::Union{String,Vector{<:String}}="clip"
  border::Real=0
  alpha::Real=0
  interpolate=0
  extrapolate::Union{Bool,String}=false
  kspline::Int=3
  twinax::Vector{<:Int}=Int[]
  axcolour::Union{String,Symbol,Vector}="black"
  figsize::Tuple{Real,Real}=(6,4)
  fontsize::Real=12
  framewidth::Real=1
  ticksize::Tuple{Real,Real}=(4.5,2.5)
  cap_offset::Real=0
  ti_offset::Real=4
  lbl_offset::Real=2
  leg_offset::Real=0
  tick_offset::Real=0

  function kwargs(ti::AbstractString="", xlabel::AbstractString="model time / hours",
  ylabel::Union{AbstractString, Vector{<:AbstractString}} =
    L"concentration / mlc cm$^{-3}$ s$^{-1}$", lt::lt_type = "",
    pt::Union{String,Int,UnitRange{Int},Vector}="", lc::lc_type="", lw::Real=1.4,
    cs::Union{String,Vector{<:String}}="", plottype::String="",
    legpos::Union{String,Int,Tuple{Real,Real}}="best", legcols::Int=1,
    xlim::Union{Nothing,Tuple{Union{Nothing,<:Real},Union{Nothing,<:Real}}}=nothing,
    ylim::Union{limtype,Vector}=nothing, mticks::Bool=true,
    Xticks::Union{Real,Vector{<:Int}}=-1, xticks::Union{Real,Vector{<:Int}}=-1,
    Yticks::Union{Real,Vector{<:Real}}=0, yticks::Union{Real,Vector{<:Real}}=0,
    timeformat::String="", timescale::String="days", major_interval::Int=0,
    minor_interval::Int=0, logscale::Union{String,Vector{<:String}}="",
    logremove::Union{String,Vector{<:String}}="clip", border::Real=0,
    alpha::Real=0, interpolate=0, extrapolate::Union{Bool,String}=false,
    kspline::Int=3, twinax::Vector{<:Int}=Int[],
    axcolour::Union{String,Symbol,Vector{<:Union{String,Symbol}}}="black",
    figsize::Tuple{Real,Real}=(6,4), fontsize::Real=12,
    framewidth::Real=1, ticksize::Tuple{Real,Real}=(4.5,2.5), cap_offset::Real=0,
    ti_offset::Real=4, lbl_offset::Real=2, leg_offset::Real=0, tick_offset::Real=0)

    # Make sure empty Vectors are of type Real
    lt = emptyvec(lt)
    # Refined tests for types in vectors
    kwargscheck(ylabel, lt, pt, lc, cs, ylim, Yticks, yticks,
      logscale, logremove, axcolour, twinax)

    # Instantiate kwargs object
    new(ti,xlabel,ylabel,lt,pt,lc,lw,cs,plottype,legpos,legcols,xlim,ylim,
      mticks,Xticks,xticks,Yticks,yticks,timeformat,timescale,major_interval,
      minor_interval,logscale,logremove,border,alpha,interpolate,extrapolate,
      kspline,twinax,axcolour,figsize,fontsize,framewidth,ticksize,cap_offset,
      ti_offset,lbl_offset,leg_offset,tick_offset)
  end #function kwargs
end

include("public.jl")
include("datamanipulation.jl")
include("plot.jl")

end #module pyp
