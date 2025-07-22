import streamlit as st
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import py3Dmol
import streamlit.components.v1 as components
import io
import zipfile

# --- Page Configuration & Theming ---
st.set_page_config(
    page_title="Protein Structure Finder",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# -- Inject Custom CSS for Pro-Level UI --
st.markdown(
    """
    <style>
    .reportview-container, .main {
      background: linear-gradient(180deg, #ffffff 0%, #f9f9f9 100%);
    }
    .sidebar .sidebar-content {
      background-color: #2C3E50;
      color: #ecf0f1;
    }
    .css-18e3th9 h1 {
      color: #1abc9c;
      text-align: center;
      font-weight: bold;
    }
    .stButton>button {
      background-color: #1abc9c;
      color: white;
      border-radius: 8px;
      height: 3em;
    }
    [role="tab"] {
      color: #555;
      background-color: #ecf0f1;
      border-radius: 8px 8px 0 0;
      padding: 0.5em;
      margin-right: 0.5em;
    }
    [role="tab"][aria-selected="true"] {
      background-color: #1abc9c;
      color: white;
    }
    .css-1n76uvr th {
      background-color: #1abc9c !important;
      color: white !important;
    }
    .stDownloadButton>button {
      background-color: #3498db;
      color: white;
      border-radius: 8px;
      height: 3em;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# --- Hybrid Search Function ---
@st.cache_data(show_spinner=False)
def search_entries(query):
    def precise_query():
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        payload = {
            "query": {
                "type": "group",
                "logical_operator": "or",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "struct.title",
                            "operator": "contains_phrase",
                            "value": query
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity.pdbx_description",
                            "operator": "contains_words",
                            "value": query
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.gene.rcsb_gene_name.value",
                            "operator": "contains_words",
                            "value": query.upper()
                        }
                    }
                ]
            },
            "return_type": "entry",
            "request_options": {"return_all_hits": True}
        }
        r = requests.post(url, json=payload)
        return [rec["identifier"] for rec in r.json().get("result_set", [])] if r.ok else []

    def fallback_full_text():
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        payload = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": query}
            },
            "return_type": "entry",
            "request_options": {"return_all_hits": True}
        }
        r = requests.post(url, json=payload)
        return [rec["identifier"] for rec in r.json().get("result_set", [])] if r.ok else []

    results = precise_query()
    return results if results else fallback_full_text()

# --- Metadata Fetcher ---
@st.cache_data(show_spinner=False)
def fetch_all_metadata(ids, max_threads):
    def fetch_one(eid):
        try:
            e = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{eid}", timeout=10)
            ed = e.json()
            method = ed.get("exptl", [{}])[0].get("method", "Unknown")
            resolution = ed.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
            title = ed.get("struct", {}).get("title", "No Title")
            orgs, chains = set(), []
            for ent in ed.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
                pe = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{eid}/{ent}", timeout=10)
                pdj = pe.json()
                src = pdj.get("rcsb_entity_source_organism", [])
                if src:
                    orgs.add(src[0].get("scientific_name", "Unknown"))
                chains += pdj.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", [])
            return {
                "PDB ID": eid,
                "Title": title,
                "Method": method,
                "Resolution (Å)": resolution,
                "Organism": ", ".join(sorted(orgs)) or "Unknown",
                "Chain Count": len(set(chains))
            }
        except:
            return None
    rows = []
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [executor.submit(fetch_one, eid) for eid in ids]
        for future in as_completed(futures):
            md = future.result()
            if md:
                rows.append(md)
    return pd.DataFrame(rows)

# --- 3D Viewer ---
def view_structure_3d(pdb_id):
    pdb_txt = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb").text
    view = py3Dmol.view(width=900, height=600)
    view.addModel(pdb_txt, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    return view

# --- Sidebar Controls ---
with st.sidebar:
    st.header("🧬 Structure Finder Filters")
    query = st.text_input("🔍 Protein/Gene Name", help="e.g. EGFR, BRCA1, TP53")
    only_human = st.checkbox("✅ Only Homo sapiens", value=True)
    monomer_only = st.checkbox("🔗 Monomeric Only (1 chain)", value=True)
    max_res = st.slider("📉 Max Resolution (Å)", 1.0, 5.0, 3.0, 0.1)
    method_filter = st.selectbox("🧪 Method", ["Any", "X-RAY", "ELECTRON MICROSCOPY", "NMR"])
    threads = st.slider("⚙️ Fetch Threads", 1, 20, 10)
    search_btn = st.button("🚀 Search PDB Structures")

# --- Main Logic ---
if search_btn:
    if not query:
        st.warning("⚠️ Please enter a protein or gene name.")
    else:
        with st.spinner("🔎 Searching RCSB PDB..."):
            entries = search_entries(query)
        if not entries:
            st.error("❌ No entries found for your query.")
        else:
            st.success(f"✅ Found {len(entries)} entries. Fetching metadata...")
            df_all = fetch_all_metadata(entries, threads)
            df = df_all.copy()
            if only_human:
                df = df[df["Organism"].str.contains("Homo sapiens", na=False)]
            if monomer_only:
                df = df[df["Chain Count"] == 1]
            if method_filter != "Any":
                df = df[df["Method"].str.contains(method_filter, case=False, na=False)]
            df = df[df["Resolution (Å)"].notnull() & (df["Resolution (Å)"] <= max_res)]
            if df.empty:
                st.warning("⚠️ No structures matched your filters.")
            else:
                st.session_state.df = df
                st.session_state.raw_count = len(entries)
                st.session_state.filtered_count = len(df)

# --- Display Results & Viewer ---
if "df" in st.session_state:
    df = st.session_state.df
    raw_n = st.session_state.raw_count
    filt_n = st.session_state.filtered_count

    col1, col2 = st.columns([1, 1])
    col1.metric("Total Entries", raw_n)
    col2.metric("Filtered Entries", filt_n)

    tab1, tab2 = st.tabs(["📋 Results Table", "🧬 3D Molecular Viewer"])

    with tab1:
        st.subheader("📄 Matching Structures")
        st.dataframe(df.sort_values("Resolution (Å)"), use_container_width=True)

        top3 = (
            df[df["Method"].str.contains("X-RAY", case=False)]
            .sort_values("Resolution (Å)")
            .head(3)
        )
        if not top3.empty:
            st.markdown("### 🏆 Top 3 X-ray Structures")
            st.table(top3)

            buf_all = io.BytesIO()
            with zipfile.ZipFile(buf_all, "w") as zipf:
                for pdb_id in df["PDB ID"]:
                    pdb_txt = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb").text
                    zipf.writestr(f"{pdb_id}.pdb", pdb_txt)
            buf_all.seek(0)
            st.download_button("📥 Download All Filtered PDBs (ZIP)", buf_all, "filtered_structures.zip", "application/zip")

            buf_top = io.BytesIO()
            with zipfile.ZipFile(buf_top, "w") as zipf:
                for pdb_id in top3["PDB ID"]:
                    pdb_txt = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb").text
                    zipf.writestr(f"{pdb_id}.pdb", pdb_txt)
            buf_top.seek(0)
            st.download_button("📥 Download Top 3 PDBs (ZIP)", buf_top, "top3_structures.zip", "application/zip")
        else:
            st.warning("⚡ No X-ray structures available in filtered results.")

    with tab2:
        st.subheader("🧬 Interactive 3D Viewer")
        choice = st.selectbox("🔬 Select a PDB ID:", df["PDB ID"].tolist())
        if choice:
            with st.spinner(f"Loading 3D model for {choice}..."):
                view = view_structure_3d(choice)
                components.html(view._make_html(), height=700)
