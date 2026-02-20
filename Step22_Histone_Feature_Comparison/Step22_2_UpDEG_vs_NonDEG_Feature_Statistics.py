import pandas as pd

# 1. Import files
# Define the file paths
updeg_file_path = '../../public/Multiomics/Step22_Histone_Feature_Comparison/SRA_Specific_UpDEG_with_GO_CCLE_Select_Features.csv'
non_deg_file_path = '../../public/Multiomics/Step22_Histone_Feature_Comparison/Non_DEG_with_GO_CCLE_Select_Features.csv'
output_file_path = '../../public/Multiomics/Step22_Histone_Feature_Comparison/Select_Feature_Statistics.csv'

# Read the CSV files into pandas DataFrames
try:
    updeg_data = pd.read_csv(updeg_file_path)
    non_deg_data = pd.read_csv(non_deg_file_path)
except FileNotFoundError as e:
    print(f"Error loading files: {e}")
    print("Please ensure the script is run from a directory where the relative paths to the CSV files are correct.")
    # Exit or handle the error as appropriate
    exit()


# 2. Make a new table named Select_Feature_Statistics
# 4. Transpose the feature list into the 'Feature' column
features = [
    'Max_DSPP_025', 'Mean_DSPP_025', 'Max_DSPP_H3K4me1', 'Mean_DSPP_H3K4me1', 'Max_DSPP_H3K4me2', 'Mean_DSPP_H3K4me2',
    'Max_DSPP_H3K4me3', 'Mean_DSPP_H3K4me3', 'Max_DSPP_H3K9ac', 'Mean_DSPP_H3K9ac', 'Max_DSPP_H3K9me3', 'Mean_DSPP_H3K9me3',
    'Max_DSPP_H3K27ac', 'Mean_DSPP_H3K27ac', 'Max_DSPP_H3K27me3', 'Mean_DSPP_H3K27me3', 'Max_DSPP_H3K36me3', 'Mean_DSPP_H3K36me3',
    'Max_DSPP_H4K20me1', 'Mean_DSPP_H4K20me1', 'Max_CP_025', 'Mean_CP_025', 'Max_CP_H3K4me1', 'Mean_CP_H3K4me1',
    'Max_CP_H3K4me2', 'Mean_CP_H3K4me2', 'Max_CP_H3K4me3', 'Mean_CP_H3K4me3', 'Max_CP_H3K9ac', 'Mean_CP_H3K9ac',
    'Max_CP_H3K9me3', 'Mean_CP_H3K9me3', 'Max_CP_H3K27ac', 'Mean_CP_H3K27ac', 'Max_CP_H3K27me3', 'Mean_CP_H3K27me3',
    'Max_CP_H3K36me3', 'Mean_CP_H3K36me3', 'Max_CP_H4K20me1', 'Mean_CP_H4K20me1', 'Max_USPP_025', 'Mean_USPP_025',
    'Max_USPP_H3K4me1', 'Mean_USPP_H3K4me1', 'Max_USPP_H3K4me2', 'Mean_USPP_H3K4me2', 'Max_USPP_H3K4me3', 'Mean_USPP_H3K4me3',
    'Max_USPP_H3K9ac', 'Mean_USPP_H3K9ac', 'Max_USPP_H3K9me3', 'Mean_USPP_H3K9me3', 'Max_USPP_H3K27ac', 'Mean_USPP_H3K27ac',
    'Max_USPP_H3K27me3', 'Mean_USPP_H3K27me3', 'Max_USPP_H3K36me3', 'Mean_USPP_H3K36me3', 'Max_USPP_H4K20me1', 'Mean_USPP_H4K20me1',
    'Max_Enhancer_025', 'Mean_Enhancer_025', 'Max_Enhancer_H3K4me1', 'Mean_Enhancer_H3K4me1', 'Max_Enhancer_H3K4me2', 'Mean_Enhancer_H3K4me2',
    'Max_Enhancer_H3K4me3', 'Mean_Enhancer_H3K4me3', 'Max_Enhancer_H3K9ac', 'Mean_Enhancer_H3K9ac', 'Max_Enhancer_H3K9me3', 'Mean_Enhancer_H3K9me3',
    'Max_Enhancer_H3K27ac', 'Mean_Enhancer_H3K27ac', 'Max_Enhancer_H3K27me3', 'Mean_Enhancer_H3K27me3', 'Max_Enhancer_H3K36me3', 'Mean_Enhancer_H3K36me3',
    'Max_Enhancer_H4K20me1', 'Mean_Enhancer_H4K20me1'
]

# 3. Add columns to the Select_Feature_Statistics table
# Initialize the statistics DataFrame
select_feature_statistics = pd.DataFrame(features, columns=['Feature'])

# 5. Calculate and record statistics for Non-DEG data
stats_non_deg = []
for feature in features:
    stats_non_deg.append({
        'Non_DEG_Min': non_deg_data[feature].min(),
        'Non_DEG_Median': non_deg_data[feature].median(),
        'Non_DEG_Mean': non_deg_data[feature].mean(),
        'Non_DEG_Max': non_deg_data[feature].max()
    })

# 6. Calculate and record statistics for UpDEG data
stats_updeg = []
for feature in features:
    stats_updeg.append({
        'UpDEG_Min': updeg_data[feature].min(),
        'UpDEG_Median': updeg_data[feature].median(),
        'UpDEG_Mean': updeg_data[feature].mean(),
        'UpDEG_Max': updeg_data[feature].max()
    })

# Combine the statistics into the main DataFrame
stats_non_deg_df = pd.DataFrame(stats_non_deg)
stats_updeg_df = pd.DataFrame(stats_updeg)

select_feature_statistics = pd.concat([select_feature_statistics, stats_non_deg_df, stats_updeg_df], axis=1)

# 7. Export the Select_Feature_Statistics table
select_feature_statistics.to_csv(output_file_path, index=False)

print(f"Successfully generated and saved statistics to {output_file_path}")

