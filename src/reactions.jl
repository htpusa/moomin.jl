mutable struct Reactions
  IDs
  names
  S
  reversible
  grRules
  subsystems
  assGenes
  inputColours
  weights
  leadingGenes
  outputColours
  outputFrequency
  combinedOutput

  function Reactions(IDs,
	  				 names,
	  				 S,
					 reversible,
					 grRules,
					 subsystems,
					 assGenes,
					 inputColours,
					 weights,
					 leadingGenes,
					 outputColours,
					 outputFrequency,
					 combinedOutput)
    new(IDs,
		names,
		S,
		reversible,
		grRules,
		subsystems,
		assGenes,
		inputColours,
		weights,
  		leadingGenes,
		outputColours,
		outputFrequency,
		combinedOutput)
  end

  function Reactions(IDs,
	  				 names,
	  				 S,
					 reversible,
					 grRules,
					 subsystems,
					 assGenes)
    Reactions(IDs,
		names,
		S,
		reversible,
		grRules,
		subsystems,
		assGenes,
		[],
		[],
  		[],
		[],
		[],
		[])
  end
end

function nReactions(reactions::Reactions)
	return length(reactions.IDs)
end

function getAssGenes(rules)
	assGenes = repeat([[]], length(rules))
	for i in 1:length(rules)
	  assGenes[i] = unique(parse.(Int64, collect(m.match for m in eachmatch(r"\d+", rules[i]))))
	end

	return assGenes
end

"""
Determine colours and weights for reactions
"""
function colour!(reactions::Reactions, genes::Genes, pThresh, alpha, precision)
	nReac = nReactions(reactions)
	colours = zeros(nReac)
	weights = fill(weight(0, pThresh, alpha, precision), nReac)
	leadingGenes = repeat([zero(Int64)], nReac)

	for reacInd in 1:nReac
		currAssGenes = reactions.assGenes[reacInd]
		if !isempty(currAssGenes)
			currColours = genes.colours[currAssGenes]
			currWeights = genes.weights[currAssGenes]
			if (1 in currColours) & (-1 in currColours)
				weights[reacInd] = weight(0.5, pThresh, alpha, precision)
			else
				colours[reacInd] = sign(sum(currColours))
				(weights[reacInd], leadGeneInd) = findmax(currWeights)
				if weights[reacInd] > 0
					leadingGenes[reacInd] = currAssGenes[leadGeneInd]
				end
			end
		end
	end

	reactions.inputColours = colours
	reactions.weights = weights
	reactions.leadingGenes = leadingGenes
end
