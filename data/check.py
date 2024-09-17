import pandas as pd
import os

def find_missing_values(directory):
    # Walk through all files and subdirectories in the given directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.csv'):
                # Construct the full file path
                file_path = os.path.join(root, file)
                # Load the CSV file
                df = pd.read_csv(file_path)
                # Check if there are any missing values
                if df.isnull().any().any():
                    print(f"Missing values found in: {file_path}")
                else:
                    print(f"No missing values in: {file_path}")

# Replace 'your_directory_path' with the path to the directory containing your CSV files
find_missing_values('processed')