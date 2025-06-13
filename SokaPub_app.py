import pandas as pd
from Bio import Entrez
import urllib.error

# Helper function to parse publication date
def parse_pub_date(pub_date):
    year = pub_date.get('Year', 'Unknown')
    month = pub_date.get('Month', '01')
    day = pub_date.get('Day', '01')
    return f"{year}-{month}-{day}"

# Set email for NCBI Entrez
Entrez.email = 'adegbesamson@gmail.com'

def search_pubmed(query, db="pubmed"):
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=50)
        record = Entrez.read(handle)
        id_list = record.get('IdList', [])

        if not id_list:
            return pd.DataFrame()

        rows = []
        for pmid in id_list:
            try:
                handle = Entrez.efetch(db=db, id=pmid, retmode='xml')
                records = Entrez.read(handle)
                for record in records['PubmedArticle']:
                    article = record['MedlineCitation']['Article']
                    title = article.get('ArticleTitle', 'Title Not Available')
                    abstract = ' '.join(article['Abstract']['AbstractText']) if 'Abstract' in article else 'Abstract Not Available'
                    authors_list = ', '.join(
                        f"{a.get('ForeName', '')} {a.get('LastName', '')}" for a in article.get('AuthorList', [])
                    ) or 'Authors Not Available'
                    journal = article['Journal'].get('Title', 'Journal Not Available')
                    pub_date = parse_pub_date(article['Journal']['JournalIssue']['PubDate'])
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    rows.append({
                        'PMID': pmid,
                        'Title': title,
                        'Abstract': abstract,
                        'Authors': authors_list,
                        'Journal': journal,
                        'URL': url,
                        'Publication Date': pub_date
                    })
            except Exception:
                continue
        return pd.DataFrame(rows)
    except Exception:
        return pd.DataFrame()

queries = [
    '"Olono"[All Fields]',
    '"genomic capacity"[All Fields] AND "precision health"[All Fields]',
    '"Africa"[All Fields] AND "genomic"[All Fields]',
]

results = pd.DataFrame()
for query in queries:
    df_pubmed = search_pubmed(query, db="pubmed")
    results = pd.concat([results, df_pubmed], ignore_index=True)
    if df_pubmed.empty:
        df_pmc = search_pubmed(query, db="pmc")
        results = pd.concat([results, df_pmc], ignore_index=True)

results.drop_duplicates(subset='PMID', inplace=True)

if not results.empty:
    results.to_excel("Broad_Search_Results.xlsx", index=False)
    print("Results saved to Broad_Search_Results.xlsx")
else:
    print("No results found.")
