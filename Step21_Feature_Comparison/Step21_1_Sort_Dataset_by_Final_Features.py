import polars as pl
import json

# Define file paths
up_deg_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/SRA_Specific_UpDEG_with_GO_CCLE.csv'
non_deg_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/Non_DEG_with_GO_CCLE.csv'
features_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_Final_Model/selected_final_features.json'

# Define output paths
output_up_deg_path = '../../public/Multiomics/Step21_Feature_Comparison/SRA_Specific_UpDEG_with_GO_CCLE_Final_Features.csv'
output_non_deg_path = '../../public/Multiomics/Step21_Feature_Comparison/Non_DEG_with_GO_CCLE_Final_Features.csv'

try:
    # 1. Import files
    # Load the list of selected features from the JSON file
    with open(features_path, 'r') as f:
        selected_features = json.load(f)

    # Load the CSV files into Polars DataFrames
    sra_specific_updeg_df = pl.read_csv(up_deg_path)
    non_deg_df = pl.read_csv(non_deg_path)

    # 2. Make a new table for SRA_Specific_UpDEG_with_GO_CCLE
    sra_specific_updeg_final_features_df = sra_specific_updeg_df.clone()

    # 3. Make a new table for Non_DEG_with_GO_CCLE
    non_deg_final_features_df = non_deg_df.clone()

    # 4. Filter columns for both tables
    # Define the base columns to always keep
    base_columns = [
        "Symbol", "chr", "elementstart", "elementend", "genechr", "genestart", 
        "geneend", "genestrand", "TSS", "025onelement", "ElementTSSDistance", 
        "DSPPstart", "DSPPend", "CPstart", "CPend", "USPPstart", "USPPend", 
        "Distalstart", "Distalend", "Enhancerstart", "Enhancerend"
    ]

    # Get the list of all columns from one of the dataframes (they share the same structure)
    all_columns = sra_specific_updeg_df.columns

    # Identify columns to keep
    columns_to_keep = [col for col in all_columns if col in base_columns or col in selected_features]

    # Select the desired columns
    sra_specific_updeg_final_features_df = sra_specific_updeg_final_features_df.select(columns_to_keep)
    non_deg_final_features_df = non_deg_final_features_df.select(columns_to_keep)

    # 5. Export SRA_Specific_UpDEG_with_GO_CCLE_Final_Features table
    sra_specific_updeg_final_features_df.write_csv(output_up_deg_path)
    print(f"Successfully exported the filtered UpDEG data to: {output_up_deg_path}")

    # 6. Export Non_DEG_with_GO_CCLE_Final_Features table
    non_deg_final_features_df.write_csv(output_non_deg_path)
    print(f"Successfully exported the filtered NonDEG data to: {output_non_deg_path}")

except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the input files are in the correct directory.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")

