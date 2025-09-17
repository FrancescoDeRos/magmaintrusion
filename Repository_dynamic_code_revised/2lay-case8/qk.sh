#!/usr/bin/env bash
#!/bin/bash

# Define the input file name
input_file="./input/input_BE.dat"
echo "Starting script..."

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found." >&2 # Output to stderr
    exit 1
fi

# Function to extract numeric values from a specific line identified by a header.
get_numeric_params() {
    local header_pattern="$1"
    local line_content=""

    line_content=$(awk -v header="$header_pattern" '
        found && !/^[[:space:]]*#/ && !/^[[:space:]]*$/ {
            print;
            exit;
        }
        $0 ~ header {
            found=1;
        }
    ' "$input_file")
    
    echo "$line_content"
}

# --- Extraction ---

echo "Extracting MU and NU..."
mu_nu_raw=$(get_numeric_params "MU / NU / rhoR1")
if [ -z "$mu_nu_raw" ]; then
    echo "Error: Could not find MU / NU data line after header. Check file format." >&2
    exit 1
fi
# ADDED sed 's/\r//g' to remove carriage returns!
mu1=$(echo "$mu_nu_raw" | awk '{print $1}' | sed 's/[Dd]/e/g; s/\r//g')
nu1=$(echo "$mu_nu_raw" | awk '{print $2}' | sed 's/[Dd]/e/g; s/\r//g')

echo "Extracting E_prop1 and E_prop2..."
e_prop_raw=$(get_numeric_params "E_prop1 / E_prop2")
if [ -z "$e_prop_raw" ]; then
    echo "Error: Could not find E_prop1 / E_prop2 data line after header. Check file format." >&2
    exit 1
fi
# ADDED sed 's/\r//g' to remove carriage returns!
e_prop1=$(echo "$e_prop_raw" | awk '{print $1}' | sed 's/[Dd]/e/g; s/\r//g')
e_prop2=$(echo "$e_prop_raw" | awk '{print $2}' | sed 's/[Dd]/e/g; s/\r//g')

# --- Validation of extracted values ---
# This regex is correct for numbers without trailing \r.
# With the added sed command, it should now pass.


echo "Extracted values: MU1=$mu1, NU1=$nu1, E_prop1=$e_prop1, E_prop2=$e_prop2"

# --- Calculations ---
echo "Performing calculations..."

# Function to safely calculate a square root using bc
calculate_kc() {
    local e_prop_val="$1"
    local mu_val="$2"
    local nu_val="$3"
    
    local numerator
    local denominator
    local term_inside_sqrt
    local kc_result

    numerator=$(echo "2 * $e_prop_val * $mu_val" | bc -l 2>/dev/null)
    if [ -z "$numerator" ]; then
        echo "Error: bc failed to calculate numerator for E_prop=$e_prop_val, MU=$mu_val." >&2
        return 1
    fi

    denominator=$(echo "1 - $nu_val" | bc -l 2>/dev/null)
    if [ -z "$denominator" ]; then
        echo "Error: bc failed to calculate denominator for NU=$nu_val." >&2
        return 1
    fi

    if (( $(echo "$denominator == 0" | bc -l) )); then
        echo "Error: Denominator is zero (1 - $nu_val = 0). Cannot calculate." >&2
        return 1
    fi
    
    term_inside_sqrt=$(echo "$numerator / $denominator" | bc -l 2>/dev/null)
    if [ -z "$term_inside_sqrt" ]; then
        echo "Error: bc failed to calculate term inside sqrt ($numerator / $denominator)." >&2
        return 1
    fi

    if (( $(echo "$term_inside_sqrt < 0" | bc -l) )); then
        echo "Error: Cannot calculate square root of a negative number ($term_inside_sqrt). Check input values." >&2
        return 1
    fi

    kc_result=$(echo "sqrt($term_inside_sqrt)" | bc -l 2>/dev/null)
    if [ -z "$kc_result" ]; then
        echo "Error: bc failed to calculate square root of $term_inside_sqrt." >&2
        return 1
    fi

    echo "$kc_result"
    return 0
}

# Calculate K_c1
kc1=$(calculate_kc "$e_prop1" "$mu1" "$nu1")
if [ $? -ne 0 ]; then
    echo "Failed to calculate K_c1." >&2
    exit 1
fi

# Calculate K_c2
kc2=$(calculate_kc "$e_prop2" "$mu1" "$nu1")
if [ $? -ne 0 ]; then
    echo "Failed to calculate K_c2." >&2
    exit 1
fi

# Print the results
echo "Calculated K_c1: $(printf "%.4f" "$kc1")"
echo "Calculated K_c2: $(printf "%.4f" "$kc2")"