from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError

def retrieve_sequences_from_ncbi(database, email, term):
  # Set up the Entrez search.
  Entrez.email = email
  handle = Entrez.esearch(db=database, term=term)
  record = Entrez.read(handle)
  handle.close()
  
  # Retrieve the sequences and their metadata.
  sequences = []
  for id in record["IdList"]:
    handle = Entrez.efetch(db=database, id=id, rettype="gb", retmode="text")
    entry = SeqIO.read(handle, "gb")
    handle.close()

    # Extract the sequence and metadata.
    try:
        sequence = str(entry.seq)
    except UndefinedSequenceError:
        sequence = None

    metadata = {
      "id": entry.id,
      "description": entry.description,
      "annotations": entry.annotations,
      "features": entry.features,
    }
    
    sequences.append((metadata, sequence))
  return sequences


import numpy as np

def UPGMA(D, labels):
    """
    Calculates the UPGMA tree from a distance matrix.
    
    Parameters
    ----------
    D : numpy array
        Distance matrix.
    labels : list
        Labels for the objects in the distance matrix.
        
    Returns
    -------
    tree : list
        List representation of the UPGMA tree.
    """
    # Number of objects
    n = len(D)
    
    # Initialize the cluster list with singleton clusters
    clusters = [[i] for i in range(n)]
    
    # Initialize the tree list
    tree = []
    
    # Iterate until there is only one cluster left
    while len(clusters) > 1:
        # Find the two closest clusters
        min_distance = float('inf')
        min_i = 0
        min_j = 0
        for i in range(len(clusters)):
            for j in range(i+1, len(clusters)):
                distance = 0
                for k in clusters[i]:
                    for l in clusters[j]:
                        distance += D[k][l]
                distance /= len(clusters[i]) * len(clusters[j])
                if distance < min_distance:
                    min_distance = distance
                    min_i = i
                    min_j = j
        
        # Merge the two closest clusters
        new_cluster = clusters[min_i] + clusters[min_j]
        clusters.pop(min_j)
        clusters.pop(min_i)
        clusters.append(new_cluster)
        
        # Update the tree list
        tree.append((min_distance/2, min_i, min_j))
    
    # Construct the tree in a more readable format
    tree_str = []
    for i, (distance, left, right) in enumerate(tree):
        left_str = labels[left] if isinstance(left, int) else f'({tree_str[left]})'
        right_str = labels[right] if isinstance(right, int) else f'({tree_str[right]})'
        tree_str.append(f'{left_str}:{distance}:{right_str}')
    
    return tree_str[0]

# Example distance matrix
example_distmat = np.array(
    [[0, 4, 3, 5],
    [ 4, 0, 6, 7],
    [ 3, 6, 0, 8],
    [ 5, 7, 8, 0]])

# Example labels
example_labels = ['A', 'B', 'C', 'D']
