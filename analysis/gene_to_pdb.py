'''
This file is used to search the PDB for experimental structures for our genes of interest.

Note: written in part using gen AI.
'''

import requests
import pandas as pd
import time
from analysis import PROJ_DIR

INPUT_CSV = PROJ_DIR / "data/uniprot.csv"
OUTPUT_CSV = PROJ_DIR / "data/gene_to_pdb.csv"

# Load UniProt accessions
df_input = pd.read_csv(INPUT_CSV)
df_input = df_input.dropna(subset=["UniProt"])  # drop rows with missing UniProt IDs

# query RCSB for a UniProt ID
def get_pdb_ids(uniprot_acc):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
  "query": {
    "type": "group",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
          "operator": "in",
          "negation": False,
          "value": [
            uniprot_acc
          ]
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
          "operator": "exact_match",
          "value": "UniProt",
          "negation": False
        }
      }
    ],
    "logical_operator": "and",
    "label": "nested-attribute"
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {
      "start": 0,
      "rows": 25
    },
    "results_content_type": [
      "experimental"
    ],
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ],
    "scoring_strategy": "combined"
  }
}

    response = requests.post(url, json=query)
    if response.status_code == 200:
        data = response.json()
        return [entry["identifier"] for entry in data.get("result_set", [])]
    else:
        print(f"Error {response.status_code} for UniProt {uniprot_acc}")
        return []

results = []

for _, row in df_input.iterrows():
    gene = row["Gene"]
    uniprot_id = row["UniProt"]
    pdb_ids = get_pdb_ids(uniprot_id)

    if pdb_ids:
        for pdb in pdb_ids:
            results.append((gene, uniprot_id, pdb))
    else:
        results.append((gene, uniprot_id, None))
    
    time.sleep(0.1) 

df_out = pd.DataFrame(results, columns=["Gene", "UniProt", "PDB_ID"])
df_out.to_csv(OUTPUT_CSV, index=False)
print(f"Results saved to {OUTPUT_CSV}")

