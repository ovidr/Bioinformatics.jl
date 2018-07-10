module Bioinformatics
export verifyDna, readStringFromFile, countBases, transcribeDna, reverseComplement

const dnaAlphabet = ['A', 'C', 'G', 'T']
const rnaAlphabet = ['A', 'C', 'G', 'U']

const dnaComplements = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')

function verifyDna(dna::String)
    if all(base -> base in dnaAlphabet, dna)
        return true
    else
        error("String is not a valid DNA string!")
    end
end

function readStringFromFile(filename::String)
    s = open(filename) do f
        read(f, String)
    end
    return strip(s)
end

function countBases(dna::String)
    countA = 0
    countC = 0
    countG = 0
    countT = 0
    for c in dna
        if c == 'A'
            countA += 1
        elseif c == 'C'
            countC += 1
        elseif c == 'G'
            countG += 1
        elseif c == 'T'
            countT += 1
        else
            error("Char $c is not a valid base!")
        end
    end
    return countA, countC, countG, countT
end

function transcribeDna(dna::String)
    rna = ""
    for c in dna
        if c == 'T'
            rna = string(rna, 'U')
        else
            rna = string(rna, c)
        end
    end
    return rna
end

function reverseComplement(dna::String)
    reversed = reverse(dna)
    reverseComplement = join([dnaComplements[c] for c in reversed])
    return reverseComplement
end

end
