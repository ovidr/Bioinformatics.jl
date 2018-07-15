module Bioinformatics
export verifyDna,
       readStringFromFile,
       countBases,
       transcribeDnaToRna,
       transcribeDnaToMRna,
       reverseComplement,
       readFASTA,
       gcContent,
       hammingDist,
       translateRNA,
       findMotif,
       profileMatrix,
       consensusString,
       proteinMass,
       longestCommonSubstring

const dnaAlphabet = ['A', 'C', 'G', 'T']
const rnaAlphabet = ['A', 'C', 'G', 'U']

const dnaComplements = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
const rnaComplements = Dict('A' => 'U', 'U' => 'A', 'C' => 'G', 'G' => 'C')

const proteinAlphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

const rnaCodonTable = Dict("UUU" => 'F',      "CUU" => 'L',      "AUU" => 'I',
                           "GUU" => 'V',      "UUC" => 'F',      "CUC" => 'L',
                           "AUC" => 'I',      "GUC" => 'V',      "UUA" => 'L',
                           "CUA" => 'L',      "AUA" => 'I',      "GUA" => 'V',
                           "UUG" => 'L',      "CUG" => 'L',      "AUG" => 'M',
                           "GUG" => 'V',      "UCU" => 'S',      "CCU" => 'P',
                           "ACU" => 'T',      "GCU" => 'A',      "UCC" => 'S',
                           "CCC" => 'P',      "ACC" => 'T',      "GCC" => 'A',
                           "UCA" => 'S',      "CCA" => 'P',      "ACA" => 'T',
                           "GCA" => 'A',      "UCG" => 'S',      "CCG" => 'P',
                           "ACG" => 'T',      "GCG" => 'A',      "UAU" => 'Y',
                           "CAU" => 'H',      "AAU" => 'N',      "GAU" => 'D',
                           "UAC" => 'Y',      "CAC" => 'H',      "AAC" => 'N',
                           "GAC" => 'D',      "UAA" => "Stop",   "CAA" => 'Q',
                           "AAA" => 'K',      "GAA" => 'E',      "UAG" => "Stop",
                           "CAG" => 'Q',      "AAG" => 'K',      "GAG" => 'E',
                           "UGU" => 'C',      "CGU" => 'R',      "AGU" => 'S',
                           "GGU" => 'G',      "UGC" => 'C',      "CGC" => 'R',
                           "AGC" => 'S',      "GGC" => 'G',      "UGA" => "Stop",
                           "CGA" => 'R',      "AGA" => 'R',      "GGA" => 'G',
                           "UGG" => 'W',      "CGG" => 'R',      "AGG" => 'R',
                           "GGG" => 'G')

const proteinMassTable = Dict('A' => 71.03711, 'C' => 103.00919,
                              'D' => 115.02694, 'E' => 129.04259,
                              'F' => 147.06841, 'G' => 57.02146,
                              'H' => 137.05891, 'I' => 113.08406,
                              'K' => 128.09496, 'L' => 113.08406,
                              'M' => 131.04049, 'N' => 114.04293,
                              'P' => 97.05276, 'Q' => 128.05858,
                              'R' => 156.10111, 'S' => 87.03203,
                              'T' => 101.04768, 'V' => 99.06841,
                              'W' => 186.07931, 'Y' => 163.06333)

"""
    verifyDna(dna::String)

Verifies that the string contains only 'A', 'T', 'G', 'C' characters.
"""
function verifyDna(dna::String)
    if all(base -> base in dnaAlphabet, dna)
        return true
    else
        error("String is not a valid DNA string!")
    end
end

"""
    readStringFromFile(filename::String)

Loads a text file that contains a DNA/RNA string.
"""
function readStringFromFile(filename::String)
    s = open(filename) do f
        read(f, String)
    end
    return strip(s)
end

"""
    countBases(dna::String)

Counts the amount of each base in a DNA chain.
"""
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

"""
    transcribeDnaToRna(dna::String)

Transcibes a DNA string to RNA.
"""
function transcribeDnaToRna(dna::String)
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

"""
    transcribeDnaToMRna(dna::String)

Transcibes a DNA string to mRNA.
"""
function transcribeDnaToMRna(dna::String)
    rna = transcribeDnaToRna(dna)
    mRNA = ""
    for b in rna
        mRNA = string(mRNA, rnaComplements[b])
    end
    return mRNA
end

"""
    reverseComplement(dna::String)

Calculates a reverse complement of a DNA string.
"""
function reverseComplement(dna::String)
    reversed = reverse(dna)
    reverseComplement = join([dnaComplements[c] for c in reversed])
    return reverseComplement
end

"""
    readFASTA(filename::String)

Reads a FASTA formatted file.
"""
function readFASTA(filename::String)
    data = Dict()
    lines = open(filename) do f
        readlines(f)
    end
    key = ""
    val = ""
    for i = 1:length(lines)
        line = lines[i]
        if startswith(line, '>')
            key = match(r"([^>])+", line).match
        else
            val = string(val, line)
            if i == length(lines)
                data[key] = val
                continue
            end
            nextLine = lines[i+1]
            if startswith(nextLine, '>')
                data[key] = val
                key = ""
                val = ""
            end
        end
    end
    return data
end

"""
    gcContent(dna::String)

Calculates G+C ratio.
"""
function gcContent(dna::String)
    n = length(dna)
    m = 0
    for b in dna
        if b == 'G' || b == 'C'
            m += 1
        end
    end
    return 100 * m / n
end

"""
    gcContent(dna::String)

Calculates G+C ratio.
"""
function gcContent(dna::Dict)
    data = Dict()
    for (k, v) in dna
        data[k] = gcContent(v)
    end
    return data
end

"""
    hammingDist(s::String, t::String)

Calculates Hamming distance between two strings.
"""
function hammingDist(s::String, t::String)
    @assert length(s) == length(t)
    dist = 0
    for i = 1:length(s)
        if s[i] != t[i]
            dist += 1
        end
    end
    return dist
end

"""
    translateRNA(rna::String)

Translates a RNA string to protein string.
"""
function translateRNA(rna::String)
    proteinString = ""
    for i = 1:3:length(rna)-2
        aa = rnaCodonTable[rna[i:i+2]]
        if aa == "Stop"
            continue
        else
            proteinString = string(proteinString, aa)
        end
    end
    return proteinString
end

"""
    findMotif(s::String, t::String)

Returns all locations of t as a substring of s.
"""
function findMotif(s::String, t::String)
  indexes = []
  for i = 1:(length(s) - length(t) + 1)
    if s[i:i+length(t) - 1] == t
      append!(indexes, i)
    end
  end
  return indexes
end

"""
    profileMatrix(dna::Dict)

Calculates profile matrix of a collection of DNA strings.
"""
function profileMatrix(dna::Dict)
  colSize = length(collect(Iterators.take(dna, 1))[1][2])
  profile = zeros(4, colSize)
  for k in keys(dna)
    str = (dna[k])
    for i = 1:length(str)
      if str[i] == 'A'
        profile[1, i] += 1
      elseif str[i] == 'C'
        profile[2, i] += 1
      elseif str[i] == 'G'
        profile[3, i] += 1
      else
        profile[4, i] += 1
      end
    end
  end
  return profile
end

"""
    consensusString(profileMatrix::Array)

Calculates consensus string from a profile matrix.
"""
function consensusString(profileMatrix::Array)
  rowSize, colSize = size(profileMatrix)
  consensus = ""
  for i = 1:colSize
    b = ' '
    ind = argmax(profileMatrix[:, i])
    if ind == 1
      b = 'A'
    elseif ind == 2
      b = 'C'
    elseif ind == 3
      b = 'G'
    else
      b = 'T'
    end
    consensus = string(consensus, b)
  end
  return consensus
end

"""
    proteinMass(protein::String)

Calculates mass of given protein string.
"""
function proteinMass(protein::String)
    mass = 0.0
    for aa in protein
        mass += proteinMassTable[aa]
    end
    return mass
end

"""
    longestCommonSubstring(dna::Dict)

Calculates longest common substring of the collection.
"""
function longestCommonSubstring(dna::Dict)
    v = collect(values(dna))
    sort!(v, by=length)
    motif = ""
    shortestStr = v[1]
    otherStr = v[2:end]
    len = length(shortestStr)
    for i = 1:len
        for j = i:len
            m = shortestStr[i:j]
            found = false
            for dna in otherStr
                if occursin(m, dna)
                    found = true
                else
                    found = false
                    break
                end
            end
            if found & (length(m) > length(motif))
                    motif = m
            end
        end
    end
    return motif
end

end
