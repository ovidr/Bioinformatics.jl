include("../src/Bioinformatics.jl")
if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

<<<<<<< HEAD
@test Bioinformatics.verifyDna("ACGT") == true
@test_throws ErrorException Bioinformatics.verifyDna("ACFY")

@test Bioinformatics.countBases("AGTCGACG") == (2, 2, 3, 1)
@test_throws ErrorException Bioinformatics.countBases("AGTGCFF")

@test Bioinformatics.transcribeDnaToRna("GATGGAACTT") == "GAUGGAACUU"
@test Bioinformatics.translateRNA(
        Bioinformatics.transcribeDnaToMRna("AAACTTGAATAAACGT")) == "FELIC"

@test Bioinformatics.reverseComplement("AAAACCCGGT") == "ACCGGGTTTT"

@test round(Bioinformatics.gcContent(
        "CCACCCTCGTGGTATGGCTAGGCATTCAG"), digits=2) == 58.62
@test Bioinformatics.gcContent(
        Dict("DNA_1" => "AGTC", "DNA_2" => "AATG")) ==
        Dict("DNA_1" => 50.0, "DNA_2" => 25.0)

@test Bioinformatics.hammingDist("GAGCCTACTAACGGGAT",
                                 "CATCGTAATGACGGCCT") == 7

@test Bioinformatics.translateRNA("AUGGCCAUGGCGCCCAGAACUGA") == "MAMAPRT"

@test Bioinformatics.findMotif("GATATATGCATATACTT", "ATAT") == [2, 4, 10]

@test Bioinformatics.profileMatrix(
        Dict("Rosalind_1" => "ATCCAGCT", "Rosalind_2" => "GGGCAACT",
             "Rosalind_3" => "ATGGATCT", "Rosalind_4" => "AAGCAACC",
             "Rosalind_5" => "TTGGAACT", "Rosalind_6" => "ATGCCATT",
             "Rosalind_7" => "ATGGCACT")) ==
            [[5.0  1.0  0.0  0.0  5.0  5.0  0.0  0.0]
             [0.0  0.0  1.0  4.0  2.0  0.0  6.0  1.0]
             [1.0  1.0  6.0  3.0  0.0  1.0  0.0  0.0]
             [1.0  5.0  0.0  0.0  0.0  1.0  1.0  6.0]]

@test Bioinformatics.consensusString(
            [[5.0  1.0  0.0  0.0  5.0  5.0  0.0  0.0]
             [0.0  0.0  1.0  4.0  2.0  0.0  6.0  1.0]
             [1.0  1.0  6.0  3.0  0.0  1.0  0.0  0.0]
             [1.0  5.0  0.0  0.0  0.0  1.0  1.0  6.0]]) == "ATGCAACT"

@test round(Bioinformatics.proteinMass("SKADYEK"), digits=2) == 821.39

@test Bioinformatics.longestCommonSubstring(
        Dict("Rosalind_1" => "AATCCACT",
             "Rosalind_2" => "GGCAACT",
             "Rosalind_3" => "AATCT",
             "Rosalind_4" => "AAAACC")) == "AA"
=======
@test verifyDna("ACGT") == true
@test_throws ErrorException verifyDna("ACFY")
@test countBases("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC") == (20, 12, 17, 21)
@test_throws ErrorException countBases("AGTGCFF")
@test transcribeDnaToRna("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"
@test translateRNA(transcribeDnaToMRna("AAACTTGAATAAACGTAACGACGGTTGCTAACGCGTTAGTGCCGCGCCCTCAGCTATAATGACATG")) == "FELICIAANDCAITARESILLY"
@test reverseComplement("AAAACCCGGT") == "ACCGGGTTTT"
@test round(gcContent("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"), 2) == 60.92
@test gcContent(Dict("DNA_1" => "AGTC", "DNA_2" => "AATG")) == Dict("DNA_1" => 50.0, "DNA_2" => 25.0)
@test hammingDist("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7
@test translateRNA("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") == "MAMAPRTEINSTRING"
@test findMotif("GATATATGCATATACTT", "ATAT") == [2, 4, 10]
@test profileMatrix(Dict("Rosalind_1" => "ATCCAGCT", "Rosalind_2" => "GGGCAACT",
                         "Rosalind_3" => "ATGGATCT", "Rosalind_4" => "AAGCAACC",
                         "Rosalind_5" => "TTGGAACT", "Rosalind_6" => "ATGCCATT",
                         "Rosalind_7" => "ATGGCACT")) ==
                                    [[5.0  1.0  0.0  0.0  5.0  5.0  0.0  0.0]
                                     [0.0  0.0  1.0  4.0  2.0  0.0  6.0  1.0]
                                     [1.0  1.0  6.0  3.0  0.0  1.0  0.0  0.0]
                                     [1.0  5.0  0.0  0.0  0.0  1.0  1.0  6.0]]
@test consensusString([[5.0  1.0  0.0  0.0  5.0  5.0  0.0  0.0]
                       [0.0  0.0  1.0  4.0  2.0  0.0  6.0  1.0]
                       [1.0  1.0  6.0  3.0  0.0  1.0  0.0  0.0]
                       [1.0  5.0  0.0  0.0  0.0  1.0  1.0  6.0]]) == "ATGCAACT"
@test round(proteinMass("SKADYEK"), 2) == 821.39
@test longestCommonSubstring(Dict("Rosalind_1" => "AATCCACT",
                                  "Rosalind_2" => "GGCAACT",
                                  "Rosalind_3" => "AATCT",
                                  "Rosalind_4" => "AAAACC")) == "AA"
>>>>>>> e30983c50c7c334acafab3c6066203162ccfa95a
