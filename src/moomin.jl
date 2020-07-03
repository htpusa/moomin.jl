#=
   MOOMIN is a method for interpreting differential
   expression data, see publication for further details
   https://doi.org/10.1093/bioinformatics/btz584
   Author: T. Pusa
   Date:   06/2020
=#
module moomin

using MAT
using DelimitedFiles
using JuMP
using CPLEX

include("data.jl")
include("weight.jl")
include("genes.jl")
include("reactions.jl")
include("model.jl")
include("MILP.jl")
include("output.jl")

"""
    runMoomin(pathToData, pathToModel;
               delimiter='\\t', geneIDhead="GeneID", PPDEhead="PPDE", logFChead="logFC",
               modelStr="model",
               pThresh=0.9, alpha=3, precision=7,
               stoichiometry=true, enumerate=1,
               printLevel=1, timeLimit=1000)

   Main function to run the full algorithm. Input is a DE data file and metabolic
   network.

#  INPUTS

- `pathToData`:      A String that gives the location of file containing DE results
- `pathToModel`:     A String that gives the location of a .mat file containing
                      the metabolic network to be used

# OPTIONAL INPUTS

- `delimiter`:      A Char that is used as the delimiter in `pathToData` (default: '\\t')
- `geneIDhead`:     A String that gives the header of the column in `pathToData`
                     that contains gene IDs (default: "GeneID")
- `PPDEhead`:       A String that gives the header of the column in `pathToData`
                     that contains PPDEs (default: "PPDE")
- `logFChead`:      A String that gives the header of the column in `pathToData`
                     that contains log fold changes (default: "logFC")
- `modelStr`:       A String that is the name of the model structure in `pathToModel`
                     Default usually works (default: "model")
- `pThresh`:        Double to be used as the threshold for a gene to be DE
                     (default: 0.9)
- `alpha`:          Double, parameter of the weight function, a higher value means less
                    evidence is needed for a change to be inferred (default: 3)
- `precision`:      Int to determine with how many significant digits reaction weights
                     are given
- `stoichiometry`:  Boolean to determine whether to use stoichiometric (true) or
                     simpler, topological constraints (false) (default: true)
- `enumerate`:      Int, how many alternative solutions are sought at most
- `printLevel`:     Int, determine how much progress info to print (0-2, default: 1)
- `timeLimit`:      Int, time limit for the MILP solver, in seconds (default: 1000)

"""
function runMoomin(pathToData, pathToModel;
               delimiter='\t', geneIDhead="GeneID", PPDEhead="PPDE", logFChead="logFC",
               modelStr="model",
            	pThresh=0.9, alpha=3, precision=7,
               stoichiometry=true, enumerate=1,
               printLevel=1, timeLimit=1000)

   model = readModel(pathToData, pathToModel;
   	delimiter=delimiter, geneIDhead=geneIDhead, PPDEhead=PPDEhead, logFChead=logFChead,
   	modelStr=modelStr,
   	pThresh=pThresh, alpha=alpha, precision=precision)

   solve!(model, CPLEX.Optimizer, stoichiometry=stoichiometry, enumerate=enumerate,
            printLevel=printLevel, timeLimit=timeLimit)

   return model
end

export Data, Genes, Reactions, MoominModel, nReactions, readData,
   mapDataToGenes!, weight, colour!,
   getAssGenes, readModel, createMILP, solve!, interpretSolution, runMoomin,
   fillOutputs!, writeOutput

end
