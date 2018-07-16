## Available Functions

The following functions are imported upon `include("Bioinformatics.jl")`

[`verifyDna(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L62)
> Verifies that the string contains only 'A', 'T', 'G', 'C' characters

[`readStringFromFile(filename::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L75)
> Loads a text file that contains a DNA/RNA string

[`countBases(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L87)
> Counts the amount of each base in a DNA chain

[`transcribeDnaToRna(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L113)
> Transcibes a DNA string to RNA

[`transcribeDnaToMRna(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L130)
> Transcibes a DNA string to mRNA

[`reverseComplement(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L144)
> Calculates a reverse complement of a DNA string

[`readFASTA(filename::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L155)
> Reads a FASTA formatted file.

[`gcContent(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L188)
> Calculates G+C ratio

[`hammingDist(s::String, t::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L217)
> Calculates Hamming distance between two strings

[`translateRNA(rna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L233)
> Translates a RNA string to protein string

[`findMotif(s::String, t::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L251)
> Returns all locations of t as a substring of s

[`profileMatrix(dna::Dict)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L266)
> Calculates profile matrix of a collection of DNA strings

[`consensusString(profileMatrix::Array)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L291)
> Calculates consensus string from a profile matrix

[`proteinMass(protein::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L316)
> Calculates mass of given protein string

[`longestCommonSubstring(dna::Dict)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L329)
> Calculates longest common substring of the collection

[`reversePalindrome(dna::String)`](https://github.com/mrtkp9993/Bioinformatics.jl/blob/master/src/Bioinformatics.jl#L361)
> Finds reverse palindrome substrings of a DNA string
