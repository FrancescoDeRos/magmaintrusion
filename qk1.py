import re
import math
import os
def calculate_kc_value(file_content):
    """
    Reads the file content, extracts MU, NU, E_prop params, computes Kc based on the specified header matching.
    input: str The entire content of the input file as a string.
    output:float The calculated Kc value, or None if values are not found.
    """
    mu1, nu1 = None, None
    e_prop = None 
    fortran_float_pattern = r"[+-]?\d+\.?\d*(?:[Dd][+-]?\d+)?"
    mu_nu_rho_pattern = r"MU\s*/\s*NU\s*/\s*rhoR1.*?\n" + \
                        r"^\s*(" + fortran_float_pattern + r")" + \
                        r"\s*(" + fortran_float_pattern + r")"
    mu_nu_rho_re = re.compile(mu_nu_rho_pattern, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    e_prop_pattern = r"E_prop\s*\(.*?\)\s*\n" + \
                     r"^\s*(" + fortran_float_pattern + r")"
    e_prop_re = re.compile(e_prop_pattern, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    match_mu_nu = mu_nu_rho_re.search(file_content)
    if match_mu_nu:
        mu1_str = match_mu_nu.group(1).replace('D', 'e').replace('d', 'e')
        nu1_str = match_mu_nu.group(2).replace('D', 'e').replace('d', 'e')
        try:
            mu1 = float(mu1_str)
            nu1 = float(nu1_str)
            print(f"Found MU : {mu1_str}, NU: {nu1_str}")
        except ValueError:
            print(f"ERROR: not converted to float. Raw strings: MU='{mu1_str}', NU='{nu1_str}'")
            mu1, nu1 = None, None
    else:
        print("ERROR: MU/NU regex not found for the first data line after header.")
        header_start_index = file_content.find("MU / NU / rhoR1")
        if header_start_index != -1:
            print(f"Content snippet near MU/NU header: \n{file_content[header_start_index:header_start_index+200]}...")
        else:
            print("ERROR: MU / NU / rhoR1 header not found in the file.")
    print("\nsearching E_prop...")
    match_e_prop = e_prop_re.search(file_content)
    if match_e_prop:
        e_prop_str = match_e_prop.group(1).replace('D', 'e').replace('d', 'e') 
        try:
            e_prop = float(e_prop_str)
            print(f"Found E_prop: {e_prop_str}")
        except ValueError:
            print(f"ERROR: E_prop not converted float.E_prop='{e_prop_str}'")
            e_prop = None
    else:
        print("ERROR: E_prop regex not found.")
        header_start_index = file_content.find("E_prop")
        if header_start_index != -1:
            print(f"Content snippet near E_prop header: \n{file_content[header_start_index:header_start_index+100]}...")
        else:
            print("ERROR:E_prop header not found in the file.")
    kc_result = None
    if mu1 is not None and nu1 is not None and e_prop is not None:
        print("\nAll parameters found. Calculating fracture thoughness..")
        try:
            denominator = (1 - nu1)
            if denominator == 0:
                print("ERROR: D(x) is zero.")
                return None
            term_inside_sqrt = (2 * e_prop * mu1) / denominator
            if term_inside_sqrt < 0:
                print(f"  ERROR: Term inside sqrt for K_c is negative: {term_inside_sqrt}. Cannot calculate.")
                return None
            kc_result = math.sqrt(term_inside_sqrt)
        except Exception as ex:
            print(f"  ERROR: An unexpected error occurred during calculation: {ex}")
    else:
        print("\nERROR: MU, NU, E_prop not found in the file.")
    return kc_result
#percorso
file_name = "/home/francesco/perceccel/larocca/furst2023code/Repository_dynamic_code_revised/1lay-test/input/input_BE.dat"

try:
    with open(file_name, 'r') as file:
        file_content = file.read()
    file_content = file_content.replace('\r\n', '\n').replace('\r', '\n')
    file_content = file_content.replace('\xa0', ' ')
    file_content = file_content.replace('\t', ' ')
    print("File processed.")
    kc_result = calculate_kc_value(file_content)
    if kc_result is not None:
        print(f"\nK_c computed: {kc_result:.4f}")
    else:
        print("\nCalculation failed. : (")
except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found.")
except Exception as e:
    print(f"Error: {e}")