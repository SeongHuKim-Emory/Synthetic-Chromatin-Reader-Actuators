import pandas as pd
import os


def process_genomic_file(input_path, output_path):
    """
    Reads a genomic data file, adds and calculates promoter, distal, and enhancer
    region coordinates, and saves the result to a new file.

    Args:
        input_path (str): The path to the input CSV file.
        output_path (str): The path where the processed CSV file will be saved.
    """
    try:
        # 1. Import the CSV file into a pandas DataFrame
        df = pd.read_csv(input_path)
        print(f"Successfully loaded {os.path.basename(input_path)}. Shape: {df.shape}")

        # 2. Add the new columns, initializing with pandas' NA value
        new_cols = [
            'DSPP_start', 'DSPP_end', 'CP_start', 'CP_end', 'USPP_start', 'USPP_end',
            'Distal_start', 'Distal_end', 'Enhancer_start', 'Enhancer_end'
        ]
        for col in new_cols:
            df[col] = pd.NA

        # Ensure key columns are numeric for calculations
        numeric_cols = ['TSS', 'element_start', 'element_end']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        # 3. Calculate Promoter Regions (DSPP, CP, USPP) based on gene_strand
        # Create boolean masks to identify forward and reverse strands
        forward_mask = df['gene_strand'] == 'forward'
        reverse_mask = df['gene_strand'] == 'reverse'

        # --- Calculations for 'forward' strand ---
        df.loc[forward_mask, 'DSPP_start'] = df.loc[forward_mask, 'TSS'] + 201
        df.loc[forward_mask, 'DSPP_end'] = df.loc[forward_mask, 'TSS'] + 2000
        df.loc[forward_mask, 'CP_start'] = df.loc[forward_mask, 'TSS'] - 200
        df.loc[forward_mask, 'CP_end'] = df.loc[forward_mask, 'TSS'] + 200
        df.loc[forward_mask, 'USPP_start'] = df.loc[forward_mask, 'TSS'] - 2000
        df.loc[forward_mask, 'USPP_end'] = df.loc[forward_mask, 'TSS'] - 201

        # --- Calculations for 'reverse' strand ---
        df.loc[reverse_mask, 'DSPP_start'] = df.loc[reverse_mask, 'TSS'] - 2000
        df.loc[reverse_mask, 'DSPP_end'] = df.loc[reverse_mask, 'TSS'] - 201
        df.loc[reverse_mask, 'CP_start'] = df.loc[reverse_mask, 'TSS'] - 200
        df.loc[reverse_mask, 'CP_end'] = df.loc[reverse_mask, 'TSS'] + 200
        df.loc[reverse_mask, 'USPP_start'] = df.loc[reverse_mask, 'TSS'] + 201
        df.loc[reverse_mask, 'USPP_end'] = df.loc[reverse_mask, 'TSS'] + 2000

        # 4. Calculate Distal Regions
        # This requires creating masks for four distinct conditions.

        # Condition 1: forward strand and element is upstream of TSS
        cond1 = (df['gene_strand'] == 'forward') & (df['element_start'] < df['TSS'])
        df.loc[cond1, 'Distal_start'] = df.loc[cond1, 'element_end'] + 1
        df.loc[cond1, 'Distal_end'] = df.loc[cond1, 'USPP_start'] - 1

        # Condition 2: forward strand and element is downstream of TSS
        cond2 = (df['gene_strand'] == 'forward') & (df['element_start'] > df['TSS'])
        df.loc[cond2, 'Distal_start'] = df.loc[cond2, 'DSPP_end'] + 1
        df.loc[cond2, 'Distal_end'] = df.loc[cond2, 'element_start'] - 1

        # Condition 3: reverse strand and element is genomically upstream of TSS
        cond3 = (df['gene_strand'] == 'reverse') & (df['element_start'] < df['TSS'])
        df.loc[cond3, 'Distal_start'] = df.loc[cond3, 'element_end'] + 1
        df.loc[cond3, 'Distal_end'] = df.loc[cond3, 'DSPP_start'] - 1

        # Condition 4: reverse strand and element is genomically downstream of TSS
        cond4 = (df['gene_strand'] == 'reverse') & (df['element_start'] > df['TSS'])
        df.loc[cond4, 'Distal_start'] = df.loc[cond4, 'USPP_end'] + 1
        df.loc[cond4, 'Distal_end'] = df.loc[cond4, 'element_start'] - 1

        # 5. Define Enhancer Regions
        # This is a direct mapping from the element's coordinates.
        df['Enhancer_start'] = df['element_start']
        df['Enhancer_end'] = df['element_end']

        # Convert all new coordinate columns to nullable integer types
        for col in new_cols:
            df[col] = df[col].astype('Int64')

        # 6. Export the updated DataFrame to a new CSV file
        # Ensure the output directory exists
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        df.to_csv(output_path, index=False)
        print(f"Processing complete. Saved updated data to {os.path.basename(output_path)}")

    except FileNotFoundError:
        print(f"Error: The file was not found at {input_path}")
    except Exception as e:
        print(f"An unexpected error occurred while processing {input_path}: {e}")


if __name__ == '__main__':
    # Define the file paths
    base_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/'

    # --- File Set 1: Non-DEG ---
    non_deg_input = os.path.join(base_path, 'GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool.csv')
    non_deg_output = os.path.join(base_path, 'GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Pro.csv')

    # --- File Set 2: SRA Specific Up-DEG ---
    up_deg_input = os.path.join(base_path, 'GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool.csv')
    up_deg_output = os.path.join(base_path,
                                 'GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Pro.csv')

    # Process both files
    print("--- Starting Genomic Data Processing ---")
    process_genomic_file(non_deg_input, non_deg_output)
    print("-" * 40)
    process_genomic_file(up_deg_input, up_deg_output)
    print("\n--- All files processed successfully. ---")
