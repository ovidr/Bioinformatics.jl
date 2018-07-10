include("../src/Bioinformatics.jl")
using Bioinformatics
if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@test verifyDna("ACGT") == true
@test_throws ErrorException verifyDna("ACFY") == false
@test countBases("AAGCCCT") == (2, 3, 1, 1)
@test transcribeDna("ACGT") == "ACGU"
@test reverseComplement("AAAACCCGGT") == "ACCGGGTTTT"
