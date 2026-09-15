# Scrapes drug smiles from PubChem and outputs to drugs.xlsx

import requests
import pandas as pd

mappings = {
    # "2169Uganda_rifampicin_resistant": {"Full name": "Rifampicin", "abbrev": "RIF"}, #
    "AMIKACIN": {"Full name": "Amikacin", "abbrev": "AMIKACIN"},
    "AMP": {"Full name": "Ampicillin", "abbrev": "AMP"},
    "AZITHROMYCIN": {"Full name": "Azithromycin", "abbrev": "AZITHROMYCIN"},
    "BDQ": {"Full name": "Bedaquiline", "abbrev": "BDQ"},
    "CAP": {"Full name": "Capreomycin", "abbrev": "CAP"},
    "CEFACLOR": {"Full name": "Cefaclor", "abbrev": "CEF"},
    "CHLORAMPHENICOL": {"Full name": "Chloramphenicol", "abbrev": "CHLORAMPHENICOL"},
    "CIPROFLOXACIN": {"Full name": "Ciprofloxacin", "abbrev": "CIPROFLOXACIN"},
    "CLARYTHROMYCIN": {"Full name": "Clarithromycin", "abbrev": "CLARYTHROMYCIN"},
    "CLOFAZIMINE": {"Full name": "Clofazimine", "abbrev": "CLOFAZIMINE"},
    "CPZ": {"Full name": "Chlorpromazine", "abbrev": "CPZ"},
    "CYCLOSERINED": {"Full name": "Cycloserine", "abbrev": "CYCLOSERINED"},
    "DELx": {"Full name": "Delamanid", "abbrev": "DELx"},
    "DOXYCYCLINE": {"Full name": "Doxycycline", "abbrev": "DOXYCYCLINE"},
    "ECONAZOLE": {"Full name": "Econazole", "abbrev": "ECONAZOLE"},
    "EMB": {"Full name": "Ethambutol", "abbrev": "EMB"},
    "EMBx": {"Full name": "Ethambutol", "abbrev": "EMBx"},
    "ERYTHROMYCIN": {"Full name": "Erythromycin", "abbrev": "ERY"},
    "ETA": {"Full name": "Ethionamide", "abbrev": "ETA"},
    "ETHIDIUMBROMIDE": {"Full name": "Ethidium bromide", "abbrev": "ETHIDIUMBROMIDE"},
    "FUSIDICACID": {"Full name": "Fusidic acid", "abbrev": "FUSIDICACID"},
    # "HN878_rif_resistant": {"Full name": "Rifampicin", "abbrev": "RIF"},
    "INH": {"Full name": "Isoniazid", "abbrev": "INH"},
    "Kanamycin": {"Full name": "Kanamycin", "abbrev": "Kanamycin"},
    "LEVO": {"Full name": "Levofloxacin", "abbrev": "LEVO"},
    "LINEZOLID": {"Full name": "Linezolid", "abbrev": "LINEZOLID"},
    "LZDx": {"Full name": "Linezolid", "abbrev": "LZDx"},
    "MENADIONE": {"Full name": "Menadione", "abbrev": "MENADIONE"},
    "MINOCYCLINE": {"Full name": "Minocycline", "abbrev": "MINOCYCLINE"},
    "MTM": {"Full name": "Mitomycin", "abbrev": "MTM"},
    "Moxifloxacin": {"Full name": "Moxifloxacin", "abbrev": "Moxifloxacin"},
    "NIGERICIN": {"Full name": "Nigericin", "abbrev": "NIGERICIN"},
    "NITROFURANTOIN": {"Full name": "Nitrofurantoin", "abbrev": "NITROFURANTOIN"},
    "NORFLOXACIN": {"Full name": "Norfloxacin", "abbrev": "NORFLOXACIN"},
    "NOVOBIOCIN": {"Full name": "Novobiocin", "abbrev": "NOVOBIOCIN"},
    "OFLOX": {"Full name": "Ofloxacin", "abbrev": "OFLOX"},
    "OFX1": {"Full name": "Ofloxacin", "abbrev": "OFX1"},
    "OXACILLIN": {"Full name": "Oxacillin", "abbrev": "OXA"}, #
    "PA824": {"Full name": "PA-824", "abbrev": "PA824"},
    "PBTZ169x": {"Full name": "Macozinone", "abbrev": "PBTZ169x"},
    "RIF": {"Full name": "Rifampicin", "abbrev": "RIF"}, #
    "ROX": {"Full name": "Roxithromycin", "abbrev": "ROX"},
    "SM": {"Full name": "Streptomycin", "abbrev": "SM"},
    "SPECTINOMYCIN": {"Full name": "Spectinomycin", "abbrev": "SPECTINOMYCIN"},
    "SQ109": {"Full name": "SQ109", "abbrev": "SQ109"},
    "SUTx": {"Full name": "Sulfamethoxazole", "abbrev": "SUTx"},
    "TET": {"Full name": "Tetracycline", "abbrev": "TET"},
    "THZ (1hrMIC)": {"Full name": "Thiazole", "abbrev": "THZ (1hrMIC)"},
    # "TKK_010025_ethionamide_resistant": {"Full name": "Ethionamide", "abbrev": "TKK_010025_ethionamide_resistant"},
    # "TKK_010033_ethionamide_resistant": {"Full name": "Ethionamide", "abbrev": "TKK_010033_ethionamide_resistant"},
    # "TKK_010040_ethionamide_resistant": {"Full name": "Ethionamide", "abbrev": "TKK_010040_ethionamide_resistant"},
    # "TRS_OFX_resistant": {"Full name": "Ofloxacin", "abbrev": "TRS_OFX_resistant"}, #
    "TRZ": {"Full name": "Trimethoprim", "abbrev": "TRZ"},
    "TUNICAMYCIN": {"Full name": "Tunicamycin", "abbrev": "TUNICAMYCIN"}, # MF: C39H64N4O16  MW: 844.9 g/mol  
    "VANCOMYCIN": {"Full name": "Vancomycin", "abbrev": "VANCOMYCIN"},
    "VERAPAMIL": {"Full name": "Verapamil", "abbrev": "VERAPAMIL"},
    "VERx": {"Full name": "Verapamil", "abbrev": "VERx"},
}

# Function to get SMILES from PubChem using compound name
def get_smiles_by_name(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/SMILES/TXT"
    response = requests.get(url)

    if response.status_code == 200:
        return response.text.strip()
    else:
        print(f"Failed to retrieve SMILES for {drug_name}: {response.status_code}")
        return None

drug_smiles_dict = {}

for original_name, info in mappings.items():
    full_name = info["Full name"]
    abbrev = info["abbrev"]
    smiles = get_smiles_by_name(full_name)

    drug_smiles_dict[original_name] = {"Full name": full_name, "abbrev": abbrev, "SMILES": smiles}

df = pd.DataFrame.from_dict(drug_smiles_dict, orient="index").reset_index(drop=True)[["Full name", "abbrev", "SMILES"]]
df.to_excel("drugs_mtb.xlsx", index=False)

print("Data saved to drugs_mtb.xlsx")