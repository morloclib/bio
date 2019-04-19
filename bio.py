from Bio import Entrez
import pandas as pd
import numpy as np

# # example call:
# composition(map_ProteinSeq_to_String(map_GeneID_to_ProteinSeq(map_Integer_to_GeneID(1616443878))))
# # for the morloc composition:
# > composition (1616443878 :: GeneID)

Entrez.email = "zbwrnz@gmail.com"

#  map_Integer_to_GeneID :: map => Integer -> ProteinSeq
def map_Integer_to_GeneID(i):
    geneID = str(i)
    handle = Entrez.esummary(db="protein", id=geneID)
    try:
        Entrez.read(handle)
    except:
        return(None)
    return(str(i))

#  map_GeneID_to_ProteinSeq :: map => GeneID -> ProteinSeq
def map_GeneID_to_ProteinSeq(geneID):
    "Look up the protein sequence for this geneID"
    handle = Entrez.efetch(db="protein", id=geneID, rettype="fasta", retmode="text")
    x = handle.readlines()
    return (geneID, ''.join([line.strip() for line in x[1:]]))

#  map_ProteinSeq_to_String :: map => ProteinSeq -> String
def map_ProteinSeq_to_String(proteinSeq):
    "proteinSeq should be a tuple with two elements, where the second is the sequence"
    return(proteinSeq[1])

#  composition :: String -> Table { Character :: Char, Count :: Int }
def composition(seq):
   return(pd.Series(list(seq)).value_counts().to_frame())
