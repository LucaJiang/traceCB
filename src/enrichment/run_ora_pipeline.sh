#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/py312/bin/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_ROOT="${RESULT_ROOT:-/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea}"
GMT_DIR="${GMT_DIR:-/home/wjiang49/group/wjiang49/data/gsea_gmt}"
WORKERS="${WORKERS:-24}"
SKIP_GMT_PREP="${SKIP_GMT_PREP:-0}"
PREPARE_GMT_ONLY="${PREPARE_GMT_ONLY:-0}"
FORCE_GMT="${FORCE_GMT:-0}"
MSIGDB_RELEASE="${MSIGDB_RELEASE:-2026.1.Hs}"
MSIGDB_BASE_URL="${MSIGDB_BASE_URL:-https://data.broadinstitute.org/gsea-msigdb/msigdb/release/${MSIGDB_RELEASE}}"
PROTEIN_ATLAS_URL="${PROTEIN_ATLAS_URL:-https://www.proteinatlas.org/download/proteinatlas.tsv.zip}"

download_file() {
  local url="$1"
  local out="$2"
  local tmp="${out}.tmp.$$"
  if [[ "${FORCE_GMT}" != "1" && -s "${out}" ]]; then
    echo "[GMT] exists: ${out}"
    return
  fi
  echo "[GMT] download: ${url}"
  mkdir -p "$(dirname "${out}")"
  curl -L --fail --retry 3 --connect-timeout 20 -o "${tmp}" "${url}"
  mv "${tmp}" "${out}"
}

download_msigdb_gmt() {
  local file="$1"
  download_file "${MSIGDB_BASE_URL}/${file}" "${GMT_DIR}/${file}"
}

download_enrichr_gmt() {
  local library="$1"
  local out="${GMT_DIR}/${library}.gmt"
  if [[ "${FORCE_GMT}" != "1" && -s "${out}" ]]; then
    echo "[GMT] exists: ${out}"
    return
  fi
  echo "[GMT] gseapy download: ${library}"
  "${PYTHON_BIN}" - "${library}" "${out}.tmp.$$" <<'PY'
import sys
import gseapy as gp

library, out = sys.argv[1], sys.argv[2]
gp.get_library(name=library, organism="Human", save=out)
PY
  mv "${out}.tmp.$$" "${out}"
}

build_hpa_blood_immune_gmt() {
  local zip_path="${GMT_DIR}/proteinatlas.tsv.zip"
  local out="${GMT_DIR}/hpa_blood_immune_2025.gmt"
  if [[ "${FORCE_GMT}" != "1" && -s "${out}" ]]; then
    echo "[GMT] exists: ${out}"
    return
  fi
  download_file "${PROTEIN_ATLAS_URL}" "${zip_path}"
  echo "[GMT] build: ${out}"
  "${PYTHON_BIN}" - "${zip_path}" "${out}.tmp.$$" <<'PY'
import re
import sys
import zipfile
from collections import defaultdict

import pandas as pd

zip_path, out_path = sys.argv[1], sys.argv[2]
columns = [
    "Gene",
    "RNA blood cell specific nTPM",
    "RNA blood lineage specific nTPM",
]

def clean_term(prefix: str, name: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_").upper()
    return f"{prefix}_{text}"

def add_terms(term_to_genes, prefix: str, gene: str, value: object) -> None:
    if not isinstance(value, str) or not value.strip():
        return
    for item in value.split(";"):
        if ":" not in item:
            continue
        label, amount = item.rsplit(":", 1)
        try:
            numeric = float(amount.strip())
        except ValueError:
            continue
        if numeric <= 0:
            continue
        term_to_genes[clean_term(prefix, label)].add(gene)

with zipfile.ZipFile(zip_path) as archive:
    with archive.open("proteinatlas.tsv") as handle:
        df = pd.read_csv(handle, sep="\t", usecols=columns)

term_to_genes = defaultdict(set)
for row in df.itertuples(index=False):
    gene = str(row.Gene).strip()
    if not gene or gene == "nan":
        continue
    add_terms(term_to_genes, "HPA_BLOOD_CELL", gene, getattr(row, "_1"))
    add_terms(term_to_genes, "HPA_BLOOD_LINEAGE", gene, getattr(row, "_2"))

with open(out_path, "w") as handle:
    for term in sorted(term_to_genes):
        genes = sorted(term_to_genes[term])
        if genes:
            handle.write(f"{term}\t\t" + "\t".join(genes) + "\n")
PY
  mv "${out}.tmp.$$" "${out}"
}

prepare_gmt_files() {
  if [[ "${SKIP_GMT_PREP}" == "1" ]]; then
    echo "[GMT] skip GMT preparation"
    return
  fi
  mkdir -p "${GMT_DIR}"

  # MSigDB release GMTs used by the significant-pathway and two-tier ORA scripts.
  download_msigdb_gmt "h.all.v2026.1.Hs.symbols.gmt"
  download_msigdb_gmt "c5.go.bp.v2026.1.Hs.symbols.gmt"
  download_msigdb_gmt "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
  download_msigdb_gmt "c7.all.v2026.1.Hs.symbols.gmt"
  download_msigdb_gmt "c8.all.v2026.1.Hs.symbols.gmt"

  # Enrichr libraries are downloaded through GSEApy to local GMT format.
  download_enrichr_gmt "GTEx_Tissues_V8_2023"
  download_enrichr_gmt "GWAS_Catalog_2025"

  # HPA blood/immune GMT is derived from the Protein Atlas TSV.
  build_hpa_blood_immune_gmt

  echo "[GMT] required GMT files:"
  for file in \
    h.all.v2026.1.Hs.symbols.gmt \
    c5.go.bp.v2026.1.Hs.symbols.gmt \
    c2.cp.reactome.v2026.1.Hs.symbols.gmt \
    c7.all.v2026.1.Hs.symbols.gmt \
    c8.all.v2026.1.Hs.symbols.gmt \
    GTEx_Tissues_V8_2023.gmt \
    GWAS_Catalog_2025.gmt \
    hpa_blood_immune_2025.gmt; do
    test -s "${GMT_DIR}/${file}"
    printf '[GMT] %s\t%s bytes\n' "${file}" "$(stat -c%s "${GMT_DIR}/${file}")"
  done
}

echo "[ORA] result root: ${RESULT_ROOT}"
echo "[ORA] workers: ${WORKERS}"
echo "[ORA] GMT dir: ${GMT_DIR}"

prepare_gmt_files

if [[ "${PREPARE_GMT_ONLY}" == "1" ]]; then
  echo "[GMT] preparation complete; PREPARE_GMT_ONLY=1 so ORA is not run"
  exit 0
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/06_ora_significant_pathways.py" \
  --out-dir "${RESULT_ROOT}/ora/significant_pathways" \
  --workers "${WORKERS}" \
  --top-terms-per-library 4

echo "[ORA] done"
