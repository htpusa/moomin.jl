mutable struct MoominModel
	data::Data
	genes::Genes
	reactions::Reactions

	function MoominModel(data::Data, genes::Genes, reactions::Reactions)
		new(data, genes, reactions)
	end
end

function nReactions(model::MoominModel)
	return nReactions(model.reactions)
end

"""
Read metabolic model and DE data from files
"""
function readModel(pathToData, pathToModel;
	delimiter='\t', geneIDhead="GeneID", PPDEhead="PPDE", logFChead="logFC",
	modelStr="model",
	pThresh=0.9, alpha=3, precision=7)

	data = readData(pathToData; delimiter=delimiter, geneIDhead=geneIDhead,
		PPDEhead=PPDEhead, logFChead=logFChead)

	modelFile = matopen(pathToModel)
	vars = matread(pathToModel)
	network = vars[modelStr]

	S = network["S"]
	nReac = size(S, 2)

	genes = Genes(getField(network, "genes", nReac))

	reactionIDs = getField(network, "rxns", nReac)
	reactionNames = getField(network, "rxnNames", nReac)
	reversible = [i<0 for i in network["lb"][:]]
	grRules = getField(network, "grRules", nReac)
	subsystems = getField(network, "subSystems", nReac)
  	assGenes = getAssGenes(getField(network, "rules", nReac))
	reactions = Reactions(reactionIDs, reactionNames, S, reversible, grRules, subsystems,
	 						assGenes)

	mapDataToGenes!(genes, data)
	colour!(genes, pThresh, alpha, precision)
	colour!(reactions, genes, pThresh, alpha, precision)

	close(modelFile)
	return MoominModel(data, genes, reactions)
end

function getField(network, fieldName, nReactions)
	if haskey(network, fieldName)
		if typeof(network[fieldName][1]) == Array{Any,2}
			return convert(Array{String,1}, [row[1] for row in network[fieldName]][:])
		else
			return convert(Array{String,1}, network[fieldName][:])
		end
	else
		@warn "Model has no entry for $fieldName"
		return fill("NA", nReactions)
	end
end
