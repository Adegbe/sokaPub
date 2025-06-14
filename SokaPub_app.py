import pandas as pd
import streamlit as st
from Bio import Entrez
import urllib.error

# Set your email for Entrez access
Entrez.email = 'adegbesamson@gmail.com'

# Helper function to parse publication date
def parse_pub_date(pub_date):
    year = pub_date.get('Year', 'Unknown')
    month = pub_date.get('Month', '01')
    day = pub_date.get('Day', '01')
    return f"{year}-{month}-{day}"

# Function to search PubMed or PMC
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

# Streamlit UI
st.title("PubMed/PMC Literature Search")
queries = [
    '"Olono"[All Fields]',
    '"genomic capacity"[All Fields] AND "precision health"[All Fields]',
    '"Africa"[All Fields] AND "genomic"[All Fields]',
]

selected_query = st.selectbox("Select a query to run:", queries)

if st.button("Search"):
    with st.spinner("Searching PubMed..."):
        df = search_pubmed(selected_query)
        if df.empty:
            st.warning("No results found.")
        else:
            st.success(f"Found {len(df)} articles.")
            st.dataframe(df[['Title', 'Authors', 'Journal', 'Publication Date', 'URL']])
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("Download results as CSV", csv, "search_results.csv", "text/csv")
