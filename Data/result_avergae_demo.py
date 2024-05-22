import os
import pandas as pd

def average_csv_files(folder_path, prefix, output_file):
    csv_files = [f for f in os.listdir(folder_path) if f.startswith(prefix) and f.endswith('.csv')]

    if not csv_files:
        print(f"No files found with prefix {prefix}")
        return

    df_list = []

    for file in csv_files:
        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path, index_col=0)
        df_list.append(df)

    total_df = sum(df_list)
    average_df = total_df / len(df_list)

    output_path = os.path.join('simu_target_avg', output_file)
    average_df.to_csv(output_path)
    print(f"Average CSV saved to {output_path}")

directory = 'synthetic_result'
Ps = [0.05,0.1,0.3,0.5]
Ks=[1000,3000,5000]
for P in Ps:
    for K in Ks:
        prefix = f"P_{P}K_{K}filter_0.1"
        output_filename = f'{prefix}.csv'
        average_csv_files(directory, prefix, output_filename)
