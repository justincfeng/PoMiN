# Home

## Introduction

Relativistic positioning refers to the concept of establishing spacetime
positions from proper time broadcasts emitted by a system of satellites.
Central to relativistic positioning is the relativistic location
location problem, which is the problem of finding the intersection of
future pointing light cones from a collection of at least four emission
points. `cereal.jl` contains a collection of functions for the
relativistic location problem in flat spacetime.

The name `cereal` is derived from an attempted pronounciation of the 
acronym SRL for **S**pecial **R**elativistic **L**ocator.

## Short tutorial

### Setup

The `cereal.jl` code was written for and tested in Julia 1.6; we
recommend Julia 1.6 or newer---please refer to the installation
instructions [`here`](https://julialang.org/downloads/platform/). To add
`cereal.jl` as a package in Julia, open the Julia `REPL`, then open the
the package manager by typing `]`. In the package manager, run the
following command:

    pkg> add https://github.com/justincfeng/cereal.jl/

The package manager can be exited by pressing the backspace key. Once 
added, one may access the `cereal` module with the following command:

    julia> using cereal

### Relativistic locator

To run the cereal code, one begins by generating a set of emission
points with the following:

    julia> ( X , Xtar ) = cereal.ceval.pgen(Float64,5)

The quantity `Xtar` is a four component vector representing the true
intersection point, and `X` is a ``4×5`` matrix consisting set of `4`
column vectors representing the coordinates of the emission points. The
emission points are constructed by finding points on the past light cone
of the target point `Xtar`.

Three different methods for finding the intersection point have been
implemented, which are represented by the strings `CFM10`, `FHC22` and
`RTC21`. The method `RTC21` (see reference below) is recommended, but
requires at least five emission points. 

To select the locator function associated with the method `RTC21`, use
the `cereal.locatorselect` function, which outputs the appropriate
locator function:

    julia> locator = cereal.locatorselect(5,"RTC21")

The first argument is the number of emission points; for the method
`RTC21`, the value should be at least `5` (larger values yield functions
which take additional points into consideration). Once the locator
function is selected, one may feed the emission point matrix `X` into
the locator function to obtain the intersection point `Xc`:

    julia> Xc = locator(X)

The intersection point may then be compared with `Xtar`:

    julia> Xc - Xtar

In most cases, the differences in the components should be on the order
of the machine precision (``∼10^{-15}`` for the default floating point
type `Float64`).

### Other methods

The methods `CFM10`, `FHC22` are four-point methods, the first argument
of `cereal.locatorselect` can have a value of at least `4`. 

To try out the method `CFM10`, one may use the following command:

    julia> locator4a = cereal.locatorselect(4,"CFM10")

To try out the method `FHC22`, use:

    julia> locator4b = cereal.locatorselect(4,"FHC22")

Since four-point methods generally suffer from the bifurcation problem
(see Coll et al., Phys. Rev. D 86, 084036 (2012)), these locator
functions return a tuple of points. It should be mentioned that if one
feeds ``4×5`` matrix `X`, the functions `locator4a` and `locator4b` only
use the first four emission points.

    julia> Xca = locator4a(X)

    julia> Xcb = locator4b(X)

In the tuple `Xca`, either `Xca[1]` or `Xca[2]` should be close to the
point `Xtar`.

If one increases the number of emission points, then the resulting
functions take additional emission points into consideration for the
purpose of minimizing errors and avoiding the bifurcation problem:

    julia> locator5a = cereal.locatorselect(5,"CFM10")

    julia> locator5b = cereal.locatorselect(5,"FHC22")

### Evaluation

Routines have been written to evaluate the methods more comprehensively.
The function `cereal.ceval.main(locator,N,q,ne)` takes a locator
function `locator` generates `N` sets of `ne` emission points `X` on the past
light cone of target points `Xtar`, feeds each set into `locator`, and
checks that `locator` yields results `Xc` that differ from `Xtar` by a
factor less than a threshold value `q`. One may run the following:

    julia> cereal.ceval.main(cereal.locatorselect(4,"CFM10"),100000,1e-6,4)

    julia> cereal.ceval.main(cereal.locatorselect(4,"FHC22"),100000,1e-6,4)

    julia> cereal.ceval.main(cereal.locatorselect(5,"RTC21"),100000,1e-9,5)

    julia> cereal.ceval.main(cereal.locatorselect(6,"RTC21"),100000,5e-13,6)

One should encounter less than `10` failures in each case.

## References

`CFM10`: Coll, B. and Ferrando, J. J. and Morales-Lladosa, J. A., *Positioning Systems in Minkowski Space-Time: from Emission to Inertial Coordinates*, Class. Quant. Grav. **27**, 065013 (2010)  
doi:[10.1088/0264-9381/27/6/065013](https://doi.org/10.1088/0264-9381/27/6/065013) [\[arXiv:0910.2568\]](https://arxiv.org/abs/0910.2568)

`RTC21`: Ruggiero, M. L., Tartaglia, A., Casalino, L., *Geometric approach to the definition of emission coordinates*, (2021)  
[\[arXiv:2111.13423\]](https://arxiv.org/abs/2111.13423)

`FHC22`: Feng, J. C., Hejda, F., Carloni, S., *Relativistic location algorithm in curved spacetime*, (2022)
[\[arXiv:2201.01774\]](https://arxiv.org/abs/2201.01774)

### Citation

If you use `cereal.jl` in your work, please cite [the `FHC22` paper](https://arxiv.org/abs/2201.01774) (and the `CFM10` and `RTC21` papers above if you have employed those methods):

```bib
@article{Feng2022relloc,
    author = "Feng, Justin C. and Hejda, Filip and Carloni, Sante",
    title = "{Relativistic location algorithm in curved spacetime}",
    eprint = "2201.01774",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "1",
    year = "2022"
}
```