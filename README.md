# moominJL

This is a Julia implementation of MOOMIN, a method for analysing differential expression (DE) data with the help of a metabolic network, that was previously
released for MATLAB. If you prefer the MATLAB version, you can find it here:

https://github.com/htpusa/moomin

All the core functionality is the same, however, some eternal utilities for manipulating metabolic networks are currently lacking in Julia.
In addition, only .mat files can be used with this version.

The purpose of MOOMIN is to take DE results and metabolic network, and produce a hypothesis of a metabolic shift. This is represented by colours (coded as integers
in the code):

-'r.red'/2: increase in flux, reaction operates in reverse direction (with respect to the direction considered as "primary" in the network)

-'red'/1: increase in flux

-'grey'/0: no change in flux

-'blue'/-1: decrease in flux

-'r.blue'/-2: decrease in flux, reaction operates in reverse direction

-'yellow'/6: indeterminate change ie a change in flux took place but the direction of the reaction and hence the nature of the change are indeterminate

For full details of the methodology see the original publication:

https://doi.org/10.1093/bioinformatics/btz584

# Dependencies

How install and use Julia, see: https://julialang.org/

The CPLEX solver is needed. Installing CPLEX for Julia: https://github.com/jump-dev/CPLEX.jl

Note that academic licenses are available for CPLEX.

# Usage

The easiest way to run moominJL, once you have installed Julia and CPLEX, is to use the Julia REPL. Navigate to this directory, and activate the package by entering
Pkg by pressing `]` and running `activate .`. After that you can test that everything is working by running `test moomin`. Exit the Pkg with backspace or ^C, and
run `using moomin`. After that you can use all functions in the package.

The only functions you would normally need are `runMoomin` and `writeOutput`. The former reads the data and the network file and runs the algorithm. It returns a
custom structure that contains all the inputs and outputs of the method. This structure can then be passed to `writeOutput` to export results. Type `?` and then
the name of the function to see how to use them.

Note that MOOMIN requires DE results obtained using Bayesian methods. That is, the input should have posterior probabilities of differential expression (PPDE). For
Bayesian DE analysis, I recommend [EBseq](http://www.bioconductor.org/packages/release/bioc/html/EBSeq.html).

The network needs to include GPR information and the gene identifiers should to match with those in the data.

# Shout-outs

moominJL uses

https://github.com/JuliaIO/MAT.jl to read .mat files

https://jump.dev/ to build and solve MILPs

https://github.com/jump-dev/CPLEX.jl to interact with CPLEX

The E. coli model used for tests was downloaded from http://systemsbiology.ucsd.edu/Downloads/EcoliCore

# Tips

If you are experiencing long running times, try turning off stoichiometry (`stoichiometry=false`). This should make the MILP easier to solve.

You can export your results in different formats using the `writeOutput` function:

-`format=full` prints a tab-delimited file with the most information.

-`format=json` with `string=false` can be used to visualise results using [Escher](https://escher.github.io/#/).

# Citation

If you use moominJL, please cite

Taneli Pusa, Mariana Galvão Ferrarini, Ricardo Andrade, Arnaud Mary, Alberto Marchetti-Spaccamela, Leen Stougie, Marie-France Sagot,
MOOMIN – Mathematical explOration of ’Omics data on a MetabolIc Network, Bioinformatics, Volume 36, Issue 2, 15 January 2020, Pages 514–523,
https://doi.org/10.1093/bioinformatics/btz584
