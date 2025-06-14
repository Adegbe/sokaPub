import streamlit as st
import pandas as pd
from Bio import Entrez

# Set Entrez email
Entrez.email = 'adegbesamson@gmail.com'

# Helper: parse publication date
def parse_pub_date(pub_date):
    year = pub_date.get('Year', 'Unknown')
    month = pub_date.get('Month', '01')
    day = pub_date.get('Day', '01')
    return f"{year}-{month}-{day}"

# Function: PubMed search
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

# Streamlit Interface
st.title("PubMed Boolean Search Tool")

# Search input sections
st.subheader("Build Your Query")

field_map = {
    "All Fields": "All Fields",
    "Title": "Title",
    "Author": "Author",
    "Abstract": "Abstract",
    "Journal": "Journal"
}

# Query parts
keyword1 = st.text_input("Keyword 1", "Africa")
field1 = st.selectbox("Field for Keyword 1", list(field_map.keys()))

bool_operator = st.selectbox("Boolean Operator", ["AND", "OR", "NOT"])

keyword2 = st.text_input("Keyword 2", "genomic")
field2 = st.selectbox("Field for Keyword 2", list(field_map.keys()))

# Execute search
if st.button("Search"):
    query = f'"{keyword1}"[{field_map[field1]}] {bool_operator} "{keyword2}"[{field_map[field2]}]'
    st.write(f"**Formatted Query:** `{query}`")

    results = search_pubmed(query)
    if not results.empty:
        st.success(f"✅ Found {len(results)} results")
        st.dataframe(results)
        results.to_excel("Broad_Search_Results.xlsx", index=False)
        st.download_button("Download Results as Excel", data=results.to_csv(index=False), file_name="pubmed_results.csv")
    else:
        st.warning("No results found.")
