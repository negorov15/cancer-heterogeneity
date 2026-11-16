import os
from numpy import mean
import pandas as pd
from calculate_distance import calculate_distance


def create_chr_dict(annotation_df):
    """
    Create a dictionary from the annotation DataFrame for faster access.
    Instead of searching through the DataFrame each time, we can use a dictionary to store the exons and strands for each transcript.
    Input:
        annotation_df: DataFrame containing the annotation data with columns "transcript_id", "genomic start", "genomic end", and "strand".
    Output:
        annotation_dict: Dictionary where keys are chromosomes.
    """
    # Group annotation file by chromosomes
    chromosomes = annotation_df.groupby("chr_name")
    annotation_dict = {}
    # create a dictionary where keys are chromosomes and values are corresponding DataFrames
    for chr, values in chromosomes:
        annotation_dict[chr] = values
    return annotation_dict


def get_transcript(var_pos, annotation):
    """
    Get the transcript(s) corresponding to a given variant position and chromosome.
    Input:
        var_pos: Variant position (int).
        chr: Chromosome name (str).
        annotation: DataFrame containing the annotation data for the given chromosome.
    Output:
        transcript_id: 
            Transcript ID(s) corresponding to the variant position. 
            If multiple transcripts are found, return a list of transcript IDs. 
            If no transcript is found, return None.
    """
    if not isinstance(var_pos, int):
        raise TypeError(f"Variant position must be an integer, got {type(var_pos).__name__}")
    variant_transcripts = annotation[(annotation["genomic start"] <= var_pos) & (annotation["genomic end"] >= var_pos)]
    if variant_transcripts.empty:
        return None
    else:
        return variant_transcripts["transcript_id"].unique().tolist()

def get_exons(transcript_id, annotation):
    """
    Get list of exons corresponding to a given transcript ID.
    input:
        transcript_id: Transcript ID (str).
        annotation: DataFrame containing the annotation data for the given chromosome.
    output:
        exons: List of tuples containing the start and end positions of the exons corresponding to the transcript ID. 
        If no exons are found, return None.
    """
    if not isinstance(transcript_id, str):
        raise TypeError(f"transcript_id must be a string, got {type(transcript_id).__name__}")
    exons = annotation[annotation["transcript_id"] == transcript_id]
    if exons.empty:
        return None
    else:
        return list(zip(exons["genomic start"], exons["genomic end"]))

def get_strand(transcript_id, annotation):
    """
    Get the strand corresponding to a given transcript ID.
    Input: 
        transcript_id: Transcript ID (str).
        annotation: DataFrame containing the annotation data for the given chromosome.
    Output:
        strand: Strand corresponding to the transcript ID.
        If no strand is found, return None.
    """
    if not isinstance(transcript_id, str):
        raise TypeError(f"transcript_id must be a string, got {type(transcript_id).__name__}")
    strand = annotation[annotation["transcript_id"] == transcript_id]["genomic strand"]
    if strand.empty:
        return None
    else:
        return strand.iloc[0]

def main(status):
    """
    Calculate distances for expressed variants
    """
    # Load annotation data and create a dictionary for faster access
    annotation_file = "output/gtf/gtf_exons_only.csv"
    annotation_df = pd.read_csv(annotation_file, sep=",")
    annotation_dict = create_chr_dict(annotation_df)
    # Load cell lines
    cell_lines = pd.read_csv("input/cell_lines_modelids.tsv",sep="\t")
    # Loop through each cell line and process its variants
    for line in cell_lines["Cell line"]:
        print(f"Processing cell line: {line}")
        results = []
        # cl_file = f"output/SomaticMutationsFiltered/{line}/{line}_SomaticMutations_expressed.vcf"
        # cl_variants = pd.read_csv(cl_file, sep="\t", comment="#", header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
        if status == "expressed":
            cl_file = f"output/FilteredSomaticMutations/{line}/expressed_variants.csv"
        elif status == "unexpressed":
            cl_file = f"output/FilteredSomaticMutations/{line}/not_expressed_variants.csv"
        cl_variants = pd.read_csv(cl_file, sep="\t")
        # For each variant in the cell line, 
        # find the corresponding transcript(s), strand, and exons, and calculate the distance to the 3' end of the transcript
        for chr, var_pos in zip(cl_variants["Chrom"], cl_variants["Pos"]):
            if chr not in annotation_dict:  # Check if the chromosome is present in the annotation dictionary
                results.append((chr, var_pos, None)) # If not, set the distance to None for that variant
                continue
            annotation = annotation_dict.get(chr) # Get the annotation DataFrame for the chromosome of the variant
            # For each variant, find the corresponding transcripts
            var_transcript_id = get_transcript(var_pos, annotation)
            if not var_transcript_id:  # If no transcript is found for the variant, set the distance to None
                results.append((chr, var_pos, None))
                continue
            transcript_df = annotation[annotation["transcript_id"].isin(var_transcript_id)]
            distances = []
            for transcript_id in var_transcript_id:
                var_strand = get_strand(transcript_id, transcript_df)  # for each variant, find the corresponding strand
                var_exons = get_exons(transcript_id, transcript_df) # for each variant, find the corresponding exons
                distance = calculate_distance(var_pos, var_exons, var_strand) # for each variant, calculate the distance to the 3' end of the transcript
                distances.append(distance)
            final_dist = mean(distances)
            results.append((chr, var_pos, final_dist))

        # Save the results to a CSV file
        results_df = pd.DataFrame(results, columns=["CHROM", "POS", "Distance to 3'end"])
        output_dir = f"output/distances/{line}"
        os.makedirs(output_dir, exist_ok=True)
        results_df.to_csv(f"output/distances/{line}/{line}_SomaticMutations_distance_{status}.csv", sep="\t", index=False)


if __name__ == "__main__":
    main("unexpressed")