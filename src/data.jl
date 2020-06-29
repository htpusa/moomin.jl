mutable struct Data
	geneIDs
	PPDE
	logFC

	function Data(geneIDs, PPDE, logFC)
		new(geneIDs, PPDE, logFC)
	end
end

"""
Read DE data from a delimited file
"""
function readData(pathToData; delimiter='\t', geneIDhead="GeneID", PPDEhead="PPDE", logFChead="logFC")
	(data, cNames) = readdlm(pathToData, delimiter, header=true)
	geneIDs = data[:, findfirst(isequal(geneIDhead), cNames[:])][:]
	PPDE = data[:, findfirst(isequal(PPDEhead), cNames[:])][:]
	FC = data[:, findfirst(isequal(logFChead), cNames[:])][:]

	return Data(geneIDs, PPDE, FC)
end
