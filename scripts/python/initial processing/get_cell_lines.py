import pandas as pd
import sys

def get_cell_lines(input_file):
    """
    Extract unique cell lines from the input CSV file.
    The input CSV is expected to have a column named "Cell line".
    """
    df = pd.read_csv(input_file, sep=";")
    cell_lines = df["Cell line"].unique()
    return cell_lines

def get_flowcells(input_file):
    """
    Extract flowcell IDs for each cell line from the input CSV file.
    The input CSV is expected to have columns "Cell line" and "Cell ID".
    The flowcell ID is extracted from the "Cell ID" column, which has the format:
        <Cell line>_<barcode>.<flowcell ID><Lane>
    """
    df = pd.read_csv(input_file, sep=";")
    
    # Extract flowcell ID (part after the dot from Cell ID column)
    df["Flowcell ID"] = df["Cell ID"].str.split(".").str[1]
    
    # Group by cell line and get unique flowcell IDs
    flowcells_by_line = df.groupby("Cell line")["Flowcell ID"].unique()
    
    return flowcells_by_line


if __name__ == "__main__":
    
    file_path = "input/cell_lines.csv"
    
    # Extract flowcells by cell line
    flowcells_by_line = get_flowcells(file_path)
    
    # Convert to DataFrame and save as CSV
    flowcells_data = []
    for cell_line, flowcells in flowcells_by_line.items():
        flowcells_data.append({
            "Cell line": cell_line,
            "Flowcell IDs": ", ".join(flowcells)
        })
    
    df_flowcells = pd.DataFrame(flowcells_data)
    df_flowcells.to_csv("output/flowcells_by_cellline.csv", index=False)
    print(f"Extracted flowcells by cell line to output/flowcells_by_cellline.csv")