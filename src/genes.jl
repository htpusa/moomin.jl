mutable struct Genes
  IDs
  PPDE
  logFC
  weights
  colours

  function Genes(IDs, PPDE, logFC, weights, colours)
    new(IDs, PPDE, logFC, weights, colours)
  end

  function Genes(IDs)
    Genes(IDs, [], [], [], [])
  end
end

function nGenes(genes::Genes)
	return length(genes.IDs)
end

"""
Map DE data to model genes
"""
function mapDataToGenes!(genes::Genes, data::Data)
	inData = [findfirst(isequal(name), data.geneIDs) for name in genes.IDs]
	found = [name in data.geneIDs for name in genes.IDs]
	inData = inData[found]
	PPDE = zeros(nGenes(genes))
	logFC = zeros(nGenes(genes))
	PPDE[found] = data.PPDE[inData]
	logFC[found] = data.logFC[inData]

	genes.PPDE = PPDE
	genes.logFC = logFC
end

"""
Determine colours and weights for genes
"""
function colour!(genes::Genes, pThresh, alpha, precision)
	genes.colours = Int.(genes.PPDE .>0.9) .* sign.(genes.logFC)
	genes.weights = weight.(genes.PPDE, pThresh, alpha, precision)
end
