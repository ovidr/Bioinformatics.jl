include("../src/Bioinformatics.jl")
using Bioinformatics
if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@test verifyDna("ACGT") == true
@test_throws ErrorException verifyDna("ACFY")
@test countBases("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC") == (20, 12, 17, 21)
@test_throws ErrorException countBases("AGTGCFF")
@test transcribeDna("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"
@test reverseComplement("AAAACCCGGT") == "ACCGGGTTTT"
@test round(gcContent("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"), 2) == 60.92
@test gcContent(Dict("DNA_1" => "AGTC", "DNA_2" => "AATG")) == Dict("DNA_1" => 50.0, "DNA_2" => 25.0)
@test hammingDist("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7
