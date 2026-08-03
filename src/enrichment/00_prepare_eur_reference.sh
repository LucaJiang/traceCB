#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RESULT_DIR="${RESULT_DIR:-${REPO_ROOT}/output/sldsc_gsea_eur_release_matched}"
REFERENCE_DIR="${REFERENCE_DIR:-${RESULT_DIR}/reference}"
BASELINE_ARCHIVE="${REFERENCE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores.tgz"
BASELINE_DIR="${REFERENCE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores"
BASELINE_URL="${BASELINE_URL:-https://zenodo.org/records/10515792/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz?download=1}"
BASELINE_MD5="b261e0caf06a003e7522938e01b3d349"
WEIGHTS_ARCHIVE="${REFERENCE_DIR}/1000G_Phase3_weights_hm3_no_MHC.tgz"
WEIGHTS_URL="${WEIGHTS_URL:-https://zenodo.org/records/10515792/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1}"
WEIGHTS_MD5="a98ac0f089ee285177544a3e6e721ca3"
HM3_SOURCE="${REFERENCE_DIR}/hm3_no_MHC.source.list.txt"
HM3_URL="${HM3_URL:-https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt?download=1}"
HM3_MD5="65a34c68833eb4a764d0707b5505b508"
BFILE_PREFIX="${BFILE_PREFIX:-${REFERENCE_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.}"
PLINK_ARCHIVE="${REFERENCE_DIR}/1000G_Phase3_plinkfiles.tgz"
PLINK_URL="${PLINK_URL:-https://zenodo.org/records/10515792/files/1000G_Phase3_plinkfiles.tgz?download=1}"
PLINK_MD5="a7773ab485827b533cb300c76356d76b"
FRQ_ARCHIVE="${REFERENCE_DIR}/1000G_Phase3_frq.tgz"
FRQ_URL="${FRQ_URL:-https://zenodo.org/records/10515792/files/1000G_Phase3_frq.tgz?download=1}"
FRQ_MD5="ac29686ffd5b6378789857a522ebca77"
PYTHON_BIN="${PYTHON_BIN:-python}"

verify_md5() {
  local path="$1"
  local expected="$2"
  local observed
  observed="$(md5sum "${path}" | awk '{print $1}')"
  [[ "${observed}" == "${expected}" ]]
}

ensure_download() {
  local label="$1"
  local target="$2"
  local url="$3"
  local expected_md5="$4"
  local partial="${target}.download"

  if [[ -s "${target}" ]] && verify_md5 "${target}" "${expected_md5}"; then
    return 0
  fi
  echo "[download] ${label}"
  curl --location --fail --retry 3 --continue-at - \
    --output "${partial}" "${url}"
  if ! verify_md5 "${partial}" "${expected_md5}"; then
    echo "[error] Published MD5 mismatch for ${partial}" >&2
    return 1
  fi
  mv -f "${partial}" "${target}"
}

mkdir -p "${REFERENCE_DIR}"

ensure_download \
  "EUR baseline-LD v2.2" "${BASELINE_ARCHIVE}" "${BASELINE_URL}" "${BASELINE_MD5}"
tar -tzf "${BASELINE_ARCHIVE}" >/dev/null

ensure_download \
  "release-matched EUR PLINK reference" "${PLINK_ARCHIVE}" "${PLINK_URL}" "${PLINK_MD5}"
tar -tzf "${PLINK_ARCHIVE}" >/dev/null
if [[ ! -f "${BFILE_PREFIX}1.bim" ]]; then
  echo "[extract] release-matched EUR PLINK reference"
  tar -xzf "${PLINK_ARCHIVE}" -C "${REFERENCE_DIR}"
fi

ensure_download \
  "release-matched EUR allele frequencies" "${FRQ_ARCHIVE}" "${FRQ_URL}" "${FRQ_MD5}"
tar -tzf "${FRQ_ARCHIVE}" >/dev/null
if [[ ! -f "${REFERENCE_DIR}/1000G_Phase3_frq/1000G.EUR.QC.1.frq" ]]; then
  echo "[extract] release-matched EUR allele frequencies"
  tar -xzf "${FRQ_ARCHIVE}" -C "${REFERENCE_DIR}"
fi

if [[ ! -f "${BASELINE_DIR}/baselineLD.1.l2.ldscore.gz" ]]; then
  echo "[extract] EUR baseline-LD v2.2"
  mkdir -p "${BASELINE_DIR}"
  tar -xzf "${BASELINE_ARCHIVE}" -C "${BASELINE_DIR}"
fi

ensure_download \
  "EUR HapMap3 non-MHC weights" "${WEIGHTS_ARCHIVE}" "${WEIGHTS_URL}" "${WEIGHTS_MD5}"
tar -tzf "${WEIGHTS_ARCHIVE}" >/dev/null
if [[ ! -f "${REFERENCE_DIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.1.l2.ldscore.gz" ]]; then
  echo "[extract] EUR HapMap3 non-MHC regression weights"
  tar -xzf "${WEIGHTS_ARCHIVE}" -C "${REFERENCE_DIR}"
fi

ensure_download \
  "HapMap3 non-MHC SNP list" "${HM3_SOURCE}" "${HM3_URL}" "${HM3_MD5}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/00_build_compatible_hm3_list.py" \
  --source-list "${HM3_SOURCE}" \
  --bfile-prefix "${BFILE_PREFIX}" \
  --baseline-prefix "${BASELINE_DIR}/baselineLD." \
  --weights-prefix "${REFERENCE_DIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." \
  --output "${REFERENCE_DIR}/hm3_no_MHC.list.txt"

{
  printf 'File\tZenodoMD5\tObservedSHA256\tSourceURL\n'
  for record in \
    "${BASELINE_ARCHIVE}|${BASELINE_MD5}|${BASELINE_URL}" \
    "${PLINK_ARCHIVE}|${PLINK_MD5}|${PLINK_URL}" \
    "${FRQ_ARCHIVE}|${FRQ_MD5}|${FRQ_URL}" \
    "${WEIGHTS_ARCHIVE}|${WEIGHTS_MD5}|${WEIGHTS_URL}" \
    "${HM3_SOURCE}|${HM3_MD5}|${HM3_URL}"
  do
    IFS='|' read -r path published_md5 source_url <<< "${record}"
    observed_sha256="$(sha256sum "${path}" | awk '{print $1}')"
    printf '%s\t%s\t%s\t%s\n' \
      "$(basename "${path}")" "${published_md5}" "${observed_sha256}" "${source_url}"
  done
  compatible_sha256="$(sha256sum "${REFERENCE_DIR}/hm3_no_MHC.list.txt" | awk '{print $1}')"
  printf '%s\t%s\t%s\t%s\n' \
    "hm3_no_MHC.list.txt" "derived" "${compatible_sha256}" \
    "intersection of the five verified reference components"
} > "${REFERENCE_DIR}/reference_manifest.tsv"

cat > "${REFERENCE_DIR}/SOURCE.txt" <<EOF
Reference population: EUR
Genome build: GRCh37/hg19
Baseline model: 1000 Genomes Phase 3 EUR baseline-LD v2.2 (97 annotations)
Reference archive: S-LDSC reference files, version 4
Reference DOI: 10.5281/zenodo.10515792
Retrieval/verification date (UTC): $(date -u +%F)
Baseline source: ${BASELINE_URL}
EUR PLINK source: ${PLINK_URL}
EUR allele-frequency source: ${FRQ_URL}
Regression weights source: ${WEIGHTS_URL}
Regression SNP-list source: ${HM3_URL}
Regression SNP list: hm3_no_MHC.list.txt (derived exact-stack intersection; MHC excluded)
PLINK prefix: ${BFILE_PREFIX}
Allele-frequency prefix: ${REFERENCE_DIR}/1000G_Phase3_frq/1000G.EUR.QC.
EOF

echo "[done] EUR reference resources prepared under ${REFERENCE_DIR}"
