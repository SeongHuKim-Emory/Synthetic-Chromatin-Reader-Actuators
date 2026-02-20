import pandas as pd
import os

def calculate_feature_statistics():
    """
    Imports UpDEG and Non-DEG feature sets, calculates descriptive statistics for
    each feature, and exports the results to a new CSV file.
    """
    # 1. Define file paths
    path_updeg = '../../public/Multiomics/Step21_Feature_Comparison/SRA_Specific_UpDEG_with_GO_CCLE_Final_Features.csv'
    path_non_deg = '../../public/Multiomics/Step21_Feature_Comparison/Non_DEG_with_GO_CCLE_Final_Features.csv'
    output_path = '../../public/Multiomics/Step21_Feature_Comparison/Feature_Statistics.csv'

    # Ensure the output directory exists
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 4. Define the list of features to be analyzed
    features = [
        'Mean_DSPP_025', 'Mean_DSPP_CEBPB', 'Max_DSPP_CHD1', 'Mean_DSPP_CHD1', 'Max_DSPP_CLOCK', 'Max_DSPP_CTCF', 
        'Max_DSPP_E4F1', 'Mean_DSPP_EGR1', 'Mean_DSPP_EP300', 'Mean_DSPP_FOS', 'Max_DSPP_GATA3', 'Max_DSPP_H3K4me3', 
        'Max_DSPP_H3K9ac', 'Max_DSPP_H3K9me3', 'Max_DSPP_H3K27ac', 'Mean_DSPP_H3K27ac', 'Max_DSPP_H3K27me3', 
        'Mean_DSPP_H3K27me3', 'Max_DSPP_H3K36me3', 'Mean_DSPP_H3K36me3', 'Max_DSPP_H4K20me1', 'Mean_DSPP_HCFC1', 
        'Mean_DSPP_HDAC2', 'Max_DSPP_HDGF', 'Max_DSPP_HSF1', 'Max_DSPP_JUND', 'Mean_DSPP_LARP7', 'Max_DSPP_MAX', 
        'Max_DSPP_MAZ', 'Mean_DSPP_PAX8', 'Mean_DSPP_PPP1R10', 'Max_DSPP_RAD21', 'Max_DSPP_RAD51', 'Max_DSPP_SIX4', 
        'Mean_DSPP_SREBF1', 'Max_DSPP_SRF', 'Mean_DSPP_TARDBP', 'Max_DSPP_TRIM22', 'Mean_DSPP_ZBTB11', 
        'Mean_DSPP_ZBTB40', 'Mean_DSPP_ZFX', 'Max_DSPP_ZKSCAN1', 'Mean_DSPP_ZKSCAN1', 'Max_DSPP_ZNF217', 
        'Max_DSPP_ZNF512B', 'Max_DSPP_ZNF574', 'Mean_DSPP_ZNF574', 'Max_DSPP_ZNF687', 'Mean_DSPP_ZNF687', 
        'Mean_CP_025', 'Mean_CP_CEBPB', 'Max_CP_COPS2', 'Max_CP_CUX1', 'Max_CP_EGR1', 'Mean_CP_ELK1', 'Max_CP_FOSL2', 
        'Mean_CP_FOXA1', 'Mean_CP_FOXK2', 'Max_CP_GABPA', 'Max_CP_H3K9ac', 'Mean_CP_H3K27ac', 'Max_CP_H3K36me3', 
        'Mean_CP_HDGF', 'Max_CP_HES1', 'Mean_CP_MBD2', 'Max_CP_NFXL1', 'Mean_CP_POLR2A', 'Max_CP_RAD51', 
        'Mean_CP_RAD51', 'Max_CP_RFX1', 'Mean_CP_RFX5', 'Mean_CP_SNIP1', 'Mean_CP_TCF7L2', 'Max_CP_TEAD4', 
        'Mean_CP_ZFX', 'Max_CP_ZHX2', 'Max_CP_ZNF8', 'Mean_CP_ZNF8', 'Max_CP_ZNF24', 'Mean_CP_ZNF512B', 
        'Mean_USPP_ATF7', 'Max_USPP_CREB1', 'Max_USPP_DPF2', 'Max_USPP_EGR1', 'Mean_USPP_ELF1', 'Max_USPP_ELK1', 
        'Mean_USPP_EP300', 'Max_USPP_FOXK2', 'Mean_USPP_FOXM1', 'Max_USPP_H3K4me1', 'Mean_USPP_H3K4me2', 
        'Mean_USPP_H3K27me3', 'Mean_USPP_H3K36me3', 'Max_USPP_HDAC2', 'Mean_USPP_JUN', 'Mean_USPP_JUND', 
        'Max_USPP_MAX', 'Mean_USPP_MAZ', 'Max_USPP_MBD2', 'Mean_USPP_MBD2', 'Max_USPP_MTA3', 'Max_USPP_NCOA3', 
        'Mean_USPP_NFRKB', 'Max_USPP_PPP1R10', 'Max_USPP_RAD51', 'Max_USPP_REST', 'Max_USPP_SMARCA5', 
        'Max_USPP_SNIP1', 'Max_USPP_SREBF1', 'Max_USPP_TEAD4', 'Max_USPP_ZBTB11', 'Mean_USPP_ZBTB33', 'Mean_USPP_ZFX', 
        'Mean_USPP_ZNF8', 'Mean_USPP_ZNF574', 'Mean_USPP_ZNF579', 'Max_USPP_ZNF687', 'Mean_USPP_ZNF687', 
        'Max_Enhancer_ARID3A', 'Max_Enhancer_BMI1', 'Mean_Enhancer_BMI1', 'Mean_Enhancer_COPS2', 'Mean_Enhancer_DPF2', 
        'Mean_Enhancer_E2F8', 'Mean_Enhancer_FOS', 'Max_Enhancer_H3K27me3', 'Mean_Enhancer_H3K27me3', 
        'Max_Enhancer_H3K36me3', 'Max_Enhancer_H4K20me1', 'Mean_Enhancer_H4K20me1', 'Mean_Enhancer_HDGF', 
        'Mean_Enhancer_HSF1', 'Max_Enhancer_LARP7', 'Mean_Enhancer_LARP7', 'Mean_Enhancer_MBD2', 'Max_Enhancer_MTA2', 
        'Max_Enhancer_NEUROD1', 'Mean_Enhancer_NFRKB', 'Max_Enhancer_PAX8', 'Mean_Enhancer_TEAD4', 'Max_Enhancer_ZFX', 
        'Mean_Enhancer_ZNF8', 'Max_Enhancer_ZNF207', 'Max_Enhancer_ZNF592', 'Mean_Enhancer_ZNF592', 'TPM', 
        'GO:0099050', 'GO:0005788', 'GO:2000741', 'MS_Value'
    ]

    try:
        # 1. Import files
        updeg_df = pd.read_csv(path_updeg)
        non_deg_df = pd.read_csv(path_non_deg)

        # 2, 3, 5, & 6: Initialize a list to hold statistics dictionaries
        stats_list = []
        for feature in features:
            # Check if the feature exists in both dataframes
            if feature in non_deg_df.columns and feature in updeg_df.columns:
                # Calculate statistics for the Non-DEG dataframe
                non_deg_stats = non_deg_df[feature].agg(['min', 'median', 'mean', 'max'])
                
                # Calculate statistics for the UpDEG dataframe
                updeg_stats = updeg_df[feature].agg(['min', 'median', 'mean', 'max'])
                
                # Append the combined statistics to the list
                stats_list.append({
                    'Feature': feature,
                    'Non_DEG_Min': non_deg_stats['min'],
                    'Non_DEG_Median': non_deg_stats['median'],
                    'Non_DEG_Mean': non_deg_stats['mean'],
                    'Non_DEG_Max': non_deg_stats['max'],
                    'UpDEG_Min': updeg_stats['min'],
                    'UpDEG_Median': updeg_stats['median'],
                    'UpDEG_Mean': updeg_stats['mean'],
                    'UpDEG_Max': updeg_stats['max']
                })
        
        # Create the Feature_Statistics table from the list of dictionaries
        feature_statistics_df = pd.DataFrame(stats_list)

        # 7. Export the table to a CSV file
        feature_statistics_df.to_csv(output_path, index=False)
        print(f"Successfully generated and saved the statistics file to: {output_path}")

    except FileNotFoundError:
        print(f"Execution failed: One or both input files were not found.")
        print(f"Please ensure the following files are correctly located:")
        print(f"- {path_updeg}")
        print(f"- {path_non_deg}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# To run the script, call the function
if __name__ == '__main__':
    calculate_feature_statistics()
