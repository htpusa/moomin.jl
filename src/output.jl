"""
  writeOutput(model::MoominModel, fileName; nSol=1, format="standard", type="output", asString=true)

  Function for exporting outputs of the algorithm into a text file in various formats.

# INPUTS

-`model`:     A MoominModel structure that is an output of `runMoomin`
-`fileName`:  String to determine were to write the results

# OPTIONAL INPUTS

-`nSol`:      Int, number of the solution to be printed (default: 1)
-`format`:    String, one of the three options (default: "standard"):
               "standard": output is a delimited file with two columns
                listing the ID and the selected output for each
                reaction (see later)
               "full": a table that lists various outputs as well as other information
                for each reaction
               "json": a json file with the selected output for each reaction
-`type`:       String, one of the four options (default: "output"):
                "output": output colour from one solution
                "input": a priori colours ie colours based solely on data
                "combined": attempted consensus between all enumerated solutions
                 If colour differs between solutions, 6/yellow (indeterminate)
                 is entered
                "frequency": how often a reaction is coloured in all enumerated
                 solutions
-`asString`:    Boolean, if true, output is given in actual colours (eg "red")
                (default: true)
                
"""


function writeOutput(model::MoominModel, fileName; nSol=1, format="standard", type="output", asString=true)
  if isempty(model.reactions.outputColours)
    @warn "The model contains no solutions."
    return
  elseif nSol > size(model.reactions.outputColours, 2)
    @warn "`nSol` exceeds number of solutions in the model."
    return
  end

  if format=="full"
    toWrite = writeFullOutput(model, nSol, asString)
  else
    if type == "output"
      solution = model.reactions.outputColours[:, nSol]
    elseif type == "input"
      solution = model.reactions.inputColours
    elseif type == "combined"
      solution = model.reactions.combinedOutput
    elseif type == "frequency"
      solution = model.reactions.outputFrequency
    else
      @warn "Cannot recognise type given."
      return
    end
    if asString & (type != "frequency")
      solution = numberToColour.(solution)
    end
    solution = string.(solution)
    if format=="standard"
      toWrite = ["ID"; model.reactions.IDs]
      toWrite = [toWrite [type; solution]]
    elseif format=="json"
      toWrite = "{\n"
      for rInd in 1:nReactions(model)-1
        toWrite = toWrite * "\"" * model.reactions.IDs[rInd] * "\": " *
                    solution[rInd] * ",\n"
      end
      toWrite = toWrite * "\"" * model.reactions.IDs[end] * "\": " *
                  solution[end] * "\n}"
      open(fileName, "w") do f
        write(f, toWrite)
      end
      return
    else
      @warn "Cannot recognise format given."
    end
  end
  writedlm(fileName, toWrite)
end

function numberToColour(n)
  if n == 2
    return "r.red"
  elseif n==1
    return "red"
  elseif n==0
    return "grey"
  elseif n==-1
    return "blue"
  elseif n==-2
    return "r.blue"
  elseif n==6
    return "yellow"
  else
    return "potato"
  end
end

function writeFullOutput(model::MoominModel, nSol, asString)
  toWrite = ["ID"; model.reactions.IDs]
  if asString
    toWrite = [toWrite ["input"; numberToColour.(model.reactions.inputColours)]]
    toWrite = [toWrite ["sol#$nSol"; numberToColour.(model.reactions.outputColours[:, nSol])]]
    toWrite = [toWrite ["combined"; numberToColour.(model.reactions.combinedOutput)]]
  else
    toWrite = [toWrite ["input"; model.reactions.inputColours]]
    toWrite = [toWrite ["sol#$nSol"; model.reactions.outputColours[:, nSol]]]
    toWrite = [toWrite ["combined"; model.reactions.combinedOutput]]
  end
  toWrite = [toWrite ["frequency"; model.reactions.outputFrequency]]
  toWrite = [toWrite ["subsystem"; model.reactions.subsystems]]
  genes = ["NA"; model.genes.IDs]
  PPDEs = ["NA"; model.genes.PPDE]
  FCs = ["NA"; model.genes.logFC]
  adjustedInds = model.reactions.leadingGenes.+1
  toWrite = [toWrite ["leadingGene"; genes[adjustedInds]]]
  toWrite = [toWrite ["PPDE"; PPDEs[adjustedInds]]]
  toWrite = [toWrite ["logFC"; FCs[adjustedInds]]]
end
