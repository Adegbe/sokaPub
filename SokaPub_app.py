import pandas as pd
from Bio import Entrez
import streamlit as st

# Set your email for NCBI Entrez (required by NCBI)
Entrez.email = 'adegbesamson@gmail.com'

# Helper function to parse publication date
def parse_pub_date(pub_date):
    year = pub_date.get('Year', 'Unknown')
    month = pub_date.get('Month', '01')
    day = pub_date.get('Day', '01')
    return f"{year}-{month}-{day}"

# Function to search PubMed
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
                        'Publication Date': pub_date,
                        'URL': url
                    })
            except Exception:
                continue
        return pd.DataFrame(rows)
    except Exception:
        return pd.DataFrame()

# Streamlit UI
st.title("PubMed Literature Search")

# Text input for dynamic PubMed query
query = st.text_input("Enter your PubMed search query", value='"Africa"[All Fields] AND "genomic"[All Fields]')

if st.button("Search"):
    if query.strip():
        st.info(f"Running query: {query}")
        df = search_pubmed(query)
        if not df.empty:
            st.success(f"Found {len(df)} results")
            st.dataframe(df)
        else:
            st.warning("No results found.")
    else:
        st.warning("Please enter a search query.")
