import pandas as pd
import os
import numpy as np

def get_lambda_values():
    """Get lambda values from user input"""
    print("\nPlease enter the lambda values for each metric:")
    lambda1 = float(input("λ₁ (coefficient for Branch Length) = "))
    lambda2 = float(input("λ₂ (coefficient for Weight) = "))
    lambda3 = float(input("λ₃ (coefficient for Degree) = "))
    lambda4 = float(input("λ₄ (coefficient for Height) = "))
    lambda5 = float(input("λ₅ (coefficient for Hamming Distance) = "))
    return lambda1, lambda2, lambda3, lambda4, lambda5

def calculate_wmfd(row, lambda1, lambda2, lambda3, lambda4, lambda5):
    """Calculate WMFD with debugging information"""
    try:
        # Get penalty value
        penalty = float(row['Penalty'])
        
        # Get common metrics values
        common_bl = float(row['Normalized_Common_BL']) if pd.notna(row['Normalized_Common_BL']) else 0
        common_weight = float(row['Normalized_Common_Weight']) if pd.notna(row['Normalized_Common_Weight']) else 0
        common_degree = float(row['Normalized_Common_Degree']) if pd.notna(row['Normalized_Common_Degree']) else 0
        common_height = float(row['Normalized_Common_Height']) if pd.notna(row['Normalized_Common_Height']) else 0
        
        # Get uncommon metrics values
        uncommon_bl = float(row['Normalized_Uncommon_BL']) if pd.notna(row['Normalized_Uncommon_BL']) else 0
        uncommon_weight = float(row['Normalized_Uncommon_Weight']) if pd.notna(row['Normalized_Uncommon_Weight']) else 0
        uncommon_degree = float(row['Normalized_Uncommon_Degree']) if pd.notna(row['Normalized_Uncommon_Degree']) else 0
        uncommon_height = float(row['Normalized_Uncommon_Height']) if pd.notna(row['Normalized_Uncommon_Height']) else 0
        
        # Get Hamming distance
        hamming_dist = float(row['Normalized_Hamming_Distance']) if pd.notna(row['Normalized_Hamming_Distance']) else 0
        
        # Calculate common nodes part
        common_part = (
            lambda1 * common_bl +
            lambda2 * common_weight +
            lambda3 * common_degree +
            lambda4 * common_height
        )
        
        # Calculate uncommon nodes part
        uncommon_part = (
            lambda1 * uncommon_bl +
            lambda2 * uncommon_weight +
            lambda3 * uncommon_degree +
            lambda4 * uncommon_height
        )
        
        # Calculate final WMFD with penalty only applied to uncommon part
        # Now including lambda5 for Hamming distance
        wmfd = common_part + (penalty * uncommon_part) + (lambda5 * hamming_dist)
        
        return wmfd
        
    except Exception as e:
        print(f"Error calculating WMFD for {row['Tree_Pair']}: {str(e)}")
        return None

def main():
    try:
        # File paths
        input_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/tree_metrics 2.csv')
        output_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/wmfd_results.csv')
        
        # Read input CSV
        print("Reading input file...")
        df = pd.read_csv(input_path)
        
        # Get lambda values from user
        lambda1, lambda2, lambda3, lambda4, lambda5 = get_lambda_values()
        
        print("\nCalculating WMFD values...")
        
        # Calculate WMFD for each row
        results = []
        for idx, row in df.iterrows():
            wmfd = calculate_wmfd(row, lambda1, lambda2, lambda3, lambda4, lambda5)
            results.append({
                'Tree_Pair': row['Tree_Pair'],
                'WMFD': round(wmfd, 4) if wmfd is not None else None
            })
            # Print calculation for verification
            print(f"\nWMFD for {row['Tree_Pair']}: {round(wmfd, 4) if wmfd is not None else None}")
        
        # Create results DataFrame
        results_df = pd.DataFrame(results)
        
        # Save results
        results_df.to_csv(output_path, index=False)
        
        print(f"\nResults have been saved to: {output_path}")
        
        # Print lambda values used
        print("\nLambda values used:")
        print(f"λ₁ (Branch Length) = {lambda1}")
        print(f"λ₂ (Weight) = {lambda2}")
        print(f"λ₃ (Degree) = {lambda3}")
        print(f"λ₄ (Height) = {lambda4}")
        print(f"λ₅ (Hamming Distance) = {lambda5}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        raise

if __name__ == "__main__":
    main()