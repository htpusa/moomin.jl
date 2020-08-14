using Test
using moomin
using JuMP
using CPLEX
using DelimitedFiles

exampleData() = Data(["g1"; "g2"; "g3"; "g4"; "g5"],
                        [1; 0.95; 0; 0.5; 0.9],
                        [1; -1; -2; 3; -1])

exampleReactions() = Reactions(["r1"; "r2"; "r3"; "r4"; "r5"],                  #IDs
                               ["reac1"; "reac2"; "reac3"; "reac4"; "reac5"],   #names
                               zeros(5,5),                                      #S
                               fill(false, 5),                                  #reversible
                               ["g2"; "g6"; "g1"; "g3"; "g4"],                  #grRules
                               ["a"; "a"; "b"; "c"; "a"],                       #subSystems
                               [[[1]]; [[2;1]]; [[1;2;3]]; [[4;5]]; [[]]])      #assGenes

exampleGenes() = Genes(["g2"; "g6"; "g1"; "g3"; "g4"])

@testset "Types" begin
    data = exampleData()
    @test data isa Data
    reactions = exampleReactions()
    @test reactions isa Reactions
    genes = exampleGenes()
    @test genes isa Genes
    model = MoominModel(data, genes, reactions)
    @test model isa MoominModel
end

@testset "Weight function" begin
    @test weight(1, 0.9, 1, 7) == 4.605170
    @test weight(0.95, 0.9, 3, 7) == 1.3862940
    @test weight(0.5, 0.9, 1, 3) == -3.220
    @test weight(0, 0.5, 1, 7) == -1.3862940
end

@testset "Reading data from file" begin
    data = readData("data/toyData.txt", delimiter='\t', geneIDhead="GeneID",
        PPDEhead="PPDE", logFChead="logFC")
    @test data isa Data
    @test data.geneIDs == ["g1"; "g2"; "g3"; "g4"; "g5"]
    @test data.PPDE == [1; 0.95; 0; 0.5; 0.9]
    @test data.logFC == [1; -1; -2; 3; -1]
end

@testset "Mapping data" begin
    data = exampleData()
    genes = exampleGenes()
    mapDataToGenes!(genes, data)
    @test genes.PPDE == [0.95; 0; 1; 0; 0.5]
    @test genes.logFC == [-1; 0; 1; -2; 3]

    colour!(genes, 0.9, 1, 7)
    @test genes.colours == [-1; 0; 1; 0; 0]
    @test genes.weights == [1.386294;
                          -4.60517;
                          4.60517;
                          -4.60517;
                          -3.218876]

    reactions = exampleReactions()
    colour!(reactions, genes, 0.9, 1, 7)
    @test reactions.inputColours == [-1; -1; 0; 0; 0]
    @test reactions.weights == [1.386294;
                                1.386294;
                                -3.218876;
                                -3.218876;
                                -4.60517]
    @test reactions.leadingGenes == [1; 1; 0; 0; 0]
end

@testset "Reading models" begin
    rules = ["x(1)";
             "( x(2) | x(3) )";
             "(( x(2) & x(3) ) & ( x(4) & x(5) & x(6) ) & x(75) )";
             ""]
    @test getAssGenes(rules) == [[[1]]; [[2;3]]; [[2;3;4;5;6;75]]; [[]]]

    model = readModel("data/toyData.txt", "data/toyModel.mat",
    	delimiter='\t', geneIDhead="GeneID", PPDEhead="PPDE", logFChead="logFC",
    	modelStr="model",
    	pThresh=0.9, alpha=1, precision=7)
    @test model isa MoominModel
    @test model.reactions.inputColours == [-1; -1; 0; 0; 0; 0; 0]
    @test model.reactions.weights == [1.386294;
                                      1.386294;
                                      -3.218876;
                                      -3.218876;
                                      -4.60517;
                                      -4.60517;
                                      -4.60517]
end

@testset "MILP" begin
    model = readModel("data/toyData.txt", "data/toyModel.mat")
    MILP = createMILP(model, CPLEX.Optimizer, stoichiometry=true)
    @test lower_bound.(MILP[:v]) == fill(-100, 7)
    @test upper_bound.(MILP[:v]) == [100; 0; fill(100, 5)]

    xPlus = [1; 0; 0; 0; 1; 0]
    xMinus = [0; 1; 0; 0; 0; 1]
    inputColours = [0; 0; 0; 0; -1; 0]
    reversible = [false; false; false; false; true; true]
    @test interpretSolution(xPlus, xMinus, inputColours, reversible) ==
        [1; -1; 0; 0; -2; 6]

    outputColours = [1  0  0;
                     0  2  2;
                     -1 -2 0;
                     0  0  0;
                     1  1  1]
    fillOutputs!(model, outputColours)
    @test model.reactions.outputColours == outputColours
    @test model.reactions.outputFrequency == round.([1/3; 2/3; 2/3; 0; 1], digits=3)
    @test model.reactions.combinedOutput == [1; 2; 6; 0; 1]

    model = runMoomin("data/ecoliData.txt", "data/ecoli.mat",
                delimiter='\t', geneIDhead="GeneID", PPDEhead="PPDE", logFChead="logFC",
                modelStr="model",
                pThresh=0.9, alpha=1, precision=7,
                stoichiometry=true, enumerate=10)
    @test size(model.reactions.outputColours, 2) == 1
    sol = readdlm("data/stoich_sol.txt")
    @test sol[:] == (model.reactions.outputColours .!= 0)

    model = runMoomin("data/ecoliData.txt", "data/ecoli.mat",
                stoichiometry=true, enumerate=10)
    @test size(model.reactions.outputColours, 2) == 3
    sol = readdlm("data/stoich_freq.txt")
    @test all(abs.(sol[:] .- model.reactions.outputFrequency) .< 0.01)

    model = runMoomin("data/ecoliData.txt", "data/ecoli.mat",
                stoichiometry=false, enumerate=10)
    @test size(model.reactions.outputColours, 2) == 1
    sol = readdlm("data/topo_sol.txt")
    @test sol[:] == (model.reactions.outputColours .!= 0)

    model = runMoomin("data/ecoliData.txt", "data/ecoli.mat",
                stoichiometry=false, alpha=1, enumerate=10)
    @test size(model.reactions.outputColours, 2) == 4
    sol = readdlm("data/topo_freq.txt")
    @test all(abs.(sol[:] .- model.reactions.outputFrequency) .< 0.01)

    rm("solver.log")
end

@testset "Printing output" begin
    data = exampleData()
    genes = exampleGenes()
    reactions = exampleReactions()
    mapDataToGenes!(genes, data)
    colour!(genes, 0.9, 1, 7)
    colour!(reactions, genes, 0.9, 1, 7)
    model = MoominModel(data, genes, reactions)
    outputColours = [ 1  0;
                     -1  2;
                     -2  0;
                      1  1;
                      6  1]
    fillOutputs!(model, outputColours)

    writeOutput(model, "data/outputExample")
    @test isfile("data/outputExample")
    output = readdlm("data/outputExample", header=true)
    @test output[2] == ["ID" "output"]
    @test output[1][:, 1] == ["r1"; "r2"; "r3"; "r4"; "r5"]
    @test output[1][:, 2] == ["red"; "blue"; "r.blue"; "red"; "yellow"]
    rm("data/outputExample")
    writeOutput(model, "data/outputExample", nSol=2, type="combined", asString=false)
    @test isfile("data/outputExample")
    output = readdlm("data/outputExample", header=true)
    @test output[2] == ["ID" "combined"]
    @test output[1][:, 2] == [1; 6; -2; 1; 6]
    rm("data/outputExample")

    writeOutput(model, "data/outputExample", format="full")
    @test isfile("data/outputExample")
    output = readdlm("data/outputExample", header=true)
    @test output[2] == ["ID" "input" "sol#1" "combined" "frequency" "subsystem" "leadingGene" "PPDE" "logFC"]
    @test output[1][:, 6] == ["a"; "a"; "b"; "c"; "a"]
    @test output[1][:, 7] == ["g2"; "g2"; "NA"; "NA"; "NA"]
    @test output[1][:, 8] == [0.95; 0.95; "NA"; "NA"; "NA"]
    @test output[1][:, 9] == [-1; -1; "NA"; "NA"; "NA"]
    rm("data/outputExample")

    writeOutput(model, "data/outputExample", format="json")
    @test isfile("data/outputExample")
    f = open("data/outputExample")
    output = read(f, String)
    close(f)
    @test output == "{\n\"r1\": red,\n\"r2\": blue,\n\"r3\": r.blue,\n\"r4\": red,\n\"r5\": yellow\n}"
    rm("data/outputExample")

    writeOutput(model, "data/outputExample", type="input")
    @test isfile("data/outputExample")
    rm("data/outputExample")
    writeOutput(model, "data/outputExample", type="frequency")
    @test isfile("data/outputExample")
    rm("data/outputExample")
end
