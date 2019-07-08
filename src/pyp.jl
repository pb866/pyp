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
  dashes::lt_type=[]
  colour::Union{Nothing,String,Symbol}=nothing
  lw::Real=1.4
  alpha::Real=1

  function PlotData(x::Union{Vector{<:Dates.TimeType}, Vector{<:Real}},
    y::Vector{<:Real}, x_lower::Vector{<:Real}=Real[], x_upper::Vector{<:Real}=Real[],
    y_lower::Vector{<:Real}=Real[], y_upper::Vector{<:Real}=Real[],
    label::AbstractString="", marker::Union{String,Int}="None",
    dashes::lt_type=[], colour::Union{Nothing,String,Symbol}=nothing,
    lw::Real=1.4, alpha::Real=1)
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
    if dashes isa Vector && !isempty(dashes)
      if !(typeof(dashes) <: Vector{<:Real})
        throw(TypeError(:PlotData, "", Vector{T} where T<:Real, typeof(dashes)))
      end
    end
    new(x,y,x_lower,x_upper,y_lower,y_upper,label,marker,dashes,colour,lw,alpha)
  end #function PlotData
end


"""


"""
@par.with_kw mutable struct kwargs
  ti::AbstractString=""
  lt::lt_type = "default"
  pt::Union{String,Int,UnitRange{Int},Vector}="default"
  lc::lc_type="default"
  lw::Real=1.4
  cs::Union{String,Vector{<:String}}=""
  plottype::String="default"
  legpos::Union{String,Int,Tuple{Real,Real}}="best"
  legcols::Int=1
  xlim::limtype=nothing
  ylim::Union{limtype,Vector}=nothing
  mticks::Bool=true
  ti_offset::Real=4
  lbl_offset::Real=2
  leg_offset::Real=0
  tick_offset::Real=0
  axcolour::Union{String,Symbol,Vector}="black"

  function kwargs(ti::AbstractString="", lt::lt_type = [],
    pt::Union{String,Int,UnitRange{Int},Vector}="default",
    lc::lc_type="default", lw::Real=1.4,
    cs::Union{String,Vector{<:String}}="", plottype::String="default",
    legpos::Union{String,Int,Tuple{Real,Real}}="best", legcols::Int=1,
    xlim::Union{Nothing,Tuple{Union{Nothing,<:Real},Union{Nothing,<:Real}}}=nothing,
    ylim=nothing, mticks::Bool=true, ti_offset::Real=4, lbl_offset::Real=2,
    leg_offset::Real=0, tick_offset::Real=0,
    axcolour::Union{String,Symbol,Vector{<:Union{String,Symbol}}}="black")

    # Test types in Vector
    if typeof(lt) isa Vector && !isempty(lt)
      if !(typeof(lt) <: Vector{<:Real})
        throw(TypeError(:PlotData, "", Vector{T} where T<:Real, typeof(lt)))
      end
    end
    if typeof(pt) isa Vector && !all([types in DataType[String, Int] for types in typeof.(pt)])
      throw(ErrorException("`pt` must be a vector with `String` and/or `Int` objects"))
    end
    if typeof(lc) isa Vector{Int} ||
      (typeof(lc) isa Vector && !all([types in DataType[String, Symbol] for types in typeof.(lc)]))
      throw(ErrorException("`lc` must be a vector with `String` and/or `Symbol` objects"))
    end
    if typeof(axcolour) isa Vector &&
      !all([types in DataType[String, Symbol] for types in typeof.(axcolour)])
      throw(ErrorException("`axcolour` must be a vector with `String` and/or `Symbol` objects"))
    end
    if typeof(ylim) isa Vector && !all([types <: limtype for types in typeof.(ylim)])
      throw(ErrorException("`ylim` must be a vector with `limtype` objects"))
    end

    # Instantiate kwargs object
    new(ti,lt,pt,lc,lw,cs,plottype,legpos,legcols,xlim,ylim,mticks,ti_offset,
      lbl_offset,leg_offset,tick_offset,axcolour)
  end #function kwargs
end

include("public.jl")
include("datamanipulation.jl")
include("plot.jl")

end #module pyp
