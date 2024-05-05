def mlc_search_entrez(config, query):
    """
    example:
        config = dict(db="nuccore", retmax=20, mindate="2021/01/01", maxdate="2021/12/31")
        query = "Influenza+A+Virus[Organism]+H5N1"
        result = mlc_search_entrez(config, query)
    """
    # The import is placed inside the function to avoid loading the module when
    # this 
    import requests
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": config["db"],
        "term": query,
        "retmode": "json",
        "retmax": str(config["retmax"]),
        "datetype": "pdat",
        "mindate": config["mindate"],
        "maxdate": config["maxdate"],
        "idtype": "acc",
    }

    result = dict()

    req = requests.get(base, params=params)
    result = req.json()["esearchresult"]

    return result["idlist"]

def mlc_nucleotide_accession_to_json(config, gb_ids):
    """
    Lookup json metadata for a list of ids in entrez.
    """
    from Bio import Entrez 
    import time
    import sys

    Entrez.email = config["email"]

    attempt = 0
    while attempt < 5:
        try:
            h = Entrez.efetch(db="nucleotide", id=gb_ids, retmode="xml")
            x = Entrez.read(h)
            h.close()
            return(x)
        except Exception as err:
            attempt += 1
            print(f"Received error from server {err}", file=sys.stderr)
            print(f"Attempt {str(attempt)} of 5 attempts", file=sys.stderr)
            print(str(config), file=sys.stderr)
            print(str(gb_ids), file=sys.stderr)
            time.sleep(15)
    raise ValueError("Failed to retrieve ids")
