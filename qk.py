import re
import math
import os
def calculate_kc_values(file_content):
    """
    Reads the file content, extracts MU, NU, E_prop params, computes Kc1 and Kc2 based on the specified header matching.
    input: str: The entire content of the input file as a string.
    output: tuple containing (Kc1, Kc2) or (None, None) if values are not found.
    """
    mu1, nu1 = None, None  # lower layer
    mu2, nu2 = None, None  # upper layer
    e_prop1, e_prop2 = None, None
    fortran_float_pattern = r"[+-]?\d+\.?\d*(?:[Dd][+-]?\d+)?"
    mu_nu_rho_pattern = (
        r"MU\s*/\s*NU\s*/\s*rhoR1.*?\n"  
        r"^\s*(" + fortran_float_pattern + r")"  
        r"\s*(" + fortran_float_pattern + r")"  
        r".*?\n"                               
        r"^\s*(" + fortran_float_pattern + r")"  
        r"\s*(" + fortran_float_pattern + r")"  
    )
    mu_nu_rho_re = re.compile(mu_nu_rho_pattern, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    e_prop_pattern = r"E_prop1\s*/\s*E_prop2.*?\n" + \
                     r"^\s*(" + fortran_float_pattern + r")" + \
                     r"\s*(" + fortran_float_pattern + r")"
    e_prop_re = re.compile(e_prop_pattern, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    print("Attempting to find MU and NU (lower and upper layer data)...")
    match_mu_nu = mu_nu_rho_re.search(file_content)
    if match_mu_nu:
        mu1_str = match_mu_nu.group(1).replace('D', 'e').replace('d', 'e')
        nu1_str = match_mu_nu.group(2).replace('D', 'e').replace('d', 'e')
        mu2_str = match_mu_nu.group(3).replace('D', 'e').replace('d', 'e')
        nu2_str = match_mu_nu.group(4).replace('D', 'e').replace('d', 'e')
        try:
            mu1 = float(mu1_str)
            nu1 = float(nu1_str)
            mu2 = float(mu2_str)
            nu2 = float(nu2_str)
            print(f"Found MU (lower layer): {mu1_str}, NU (lower layer): {nu1_str}")
            print(f"Found MU (upper layer): {mu2_str}, NU (upper layer): {nu2_str}")
        except ValueError:
            print(f"ERROR: MU/NU not convrt to float. Raw strings:")
            print(f"Lower layer: MU='{mu1_str}', NU='{nu1_str}'")
            print(f"Upper layer: MU='{mu2_str}', NU='{nu2_str}'")
            mu1, nu1, mu2, nu2 = None, None, None, None # Invalidate all if any conversion fails
    else:
        print("  FAIL: MU/NU regex did NOT find a match for both data lines after header.")
        header_start_index = file_content.find("MU / NU / rhoR1")
        if header_start_index != -1:
            print(f"Content near MU/NU header: \n{file_content[header_start_index:header_start_index+300]}...")
        else:
            print("MU / NU / rhoR1 header not found in file.")
    print("\nSearching E_prop1 and E_prop2...")
    match_e_prop = e_prop_re.search(file_content)
    if match_e_prop:
        e_prop1_str = match_e_prop.group(1).replace('D', 'e').replace('d', 'e')
        e_prop2_str = match_e_prop.group(2).replace('D', 'e').replace('d', 'e')
        try:
            e_prop1 = float(e_prop1_str)
            e_prop2 = float(e_prop2_str)
            print(f"Found E_prop1: {e_prop1_str}, E_prop2: {e_prop2_str}")
        except ValueError:
            print(f"ERROR: E_prop not convrt to float. Raw strings: E_prop1='{e_prop1_str}', E_prop2='{e_prop2_str}'")
            e_prop1, e_prop2 = None, None
    else:
        print("ERROR: E_prop regex not found.")
        header_start_index = file_content.find("E_prop1 / E_prop2")
        if header_start_index != -1:
            print(f"Content near E_prop header: \n{file_content[header_start_index:header_start_index+100]}...")
        else:
            print("E_prop1 / E_prop2 header not found in file.")
    kc1, kc2 = None, None
    if mu1 is not None and nu1 is not None and e_prop1 is not None and e_prop2 is not None:
        print("\nMU1, NU1, E_prop1, E_prop2 found. Computing...")
        try:
            denominator = (1 - nu1)
            if denominator == 0:
                print("ERROR: D(x) is zero.")
                return None, None
            term_inside_sqrt1 = (2 * e_prop1 * mu1) / denominator
            if term_inside_sqrt1 < 0:
                 print(f"ERROR: Term inside sqrt for K_c1 is negative: {term_inside_sqrt1}.")
                 return None, None
            kc1 = math.sqrt(term_inside_sqrt1)
            term_inside_sqrt2 = (2 * e_prop2 * mu1) / denominator
            if term_inside_sqrt2 < 0:
                 print(f"ERROR: Term inside sqrt for K_c2 is negative: {term_inside_sqrt2}.")
                 return None, None
            kc2 = math.sqrt(term_inside_sqrt2)
        except Exception as ex:
            print(f"ERROR:{ex}")
    else:
        print("\nWarning: Could not find all necessary parameters (MU1, NU1, E_prop1, E_prop2) for calculation. Calculation aborted.")
    return kc1, kc2
#percorso
file_name = "/home/francesco/perceccel/larocca/furst2023code/Repository_dynamic_code_revised/2lay-case7/input/input_BE.dat" 
try:
    with open(file_name, 'r') as file:
        file_content = file.read()
    file_content = file_content.replace('\r\n', '\n').replace('\r', '\n')
    file_content = file_content.replace('\xa0', ' ')
    file_content = file_content.replace('\t', ' ')
    print("File processed.")
    kc1_result, kc2_result = calculate_kc_values(file_content)
    if kc1_result is not None and kc2_result is not None:
        print(f"\nK_c1 computed:: {kc1_result:.4f}")
        print(f"K_c2 computed:: {kc2_result:.4f}")
    else:
        print("\nCalculation failed. : (")
except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found.")
except Exception as e:
    print(f"ERROR: {e}")