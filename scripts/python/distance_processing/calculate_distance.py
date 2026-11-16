def get_exon_index(pos, exons):
    """
    Helper function. Get the index of the exon that contains the given position.
    """
    for i, (start, end) in enumerate(exons):
        if start <= pos <= end:
            return i
    return None  # Position is not within any exon

def calculate_distance(pos, exons, strand):
    """
    Calculate the distance from the 3'end to the variant's position.
    Parameters:
    pos (int): The position of the variant.
    exons (list of tuples): List of exon positions as (start, end) tuples.
    strand (str): The strand information ('+' or '-').
    
    Returns:
    int: The distance from the 3'end to the variant's position.
    """
    if strand not in ("+", "-"):
        raise ValueError(f"strand must be '+' or '-', got {strand!r}")

    exons = sorted(exons)  # Sort exons by start position
    k = get_exon_index(pos, exons)
    if k == None:
            return None  # Position is not within any exon
    s_k = exons[k][0]  # variant's exon start
    e_k = exons[k][1]  # variant's exon end
    exons_lengths = [end - start + 1 for (start, end) in exons]
    if strand == "+":
        d = (e_k - pos) + sum(exons_lengths[k+1:])
    elif strand == "-":
        d = (pos - s_k) + sum(exons_lengths[:k])
    return d