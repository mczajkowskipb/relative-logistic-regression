#!/usr/bin/env python3
"""GSE72094 GEO source-header inventory. NO models / no expression retrieval.

Reads the public NCBI GEO Series Matrix header only. It cannot establish
full compressed-file checksum, patient independence or exact KRAS labels.
Run from repo root: python3 scripts/rlmt_dc_geo_header_intake.py
"""
from __future__ import annotations

import collections
import csv
import datetime
import gzip
import hashlib
import io
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/rlmt_dc_gse72094_header_intake_v1"
URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/"
    "GSE72094/matrix/GSE72094_series_matrix.txt.gz"
)
EXPECTED_ACCESSION = "GSE72094"


def fail(message: str) -> None:
    raise SystemExit("INTAKE_STOP: " + message)


def parse_geo_line(line: str) -> list[str]:
    try:
        return next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
    except (csv.Error, StopIteration) as exc:
        fail(f"invalid GEO sample metadata TSV: {exc}")


def main() -> None:
    if Path.cwd().resolve() != ROOT:
        fail("run from repository root")
    if OUT.exists():
        fail(f"refusing to overwrite existing report: {OUT}")
    if not (ROOT / "docs/RLMT_DC_GSE72094_SOURCE_INTAKE.md").is_file():
        fail("source intake protocol has not been recorded")
    if not re.fullmatch(
        r"https://ftp\.ncbi\.nlm\.nih\.gov/geo/series/"
        r"GSE72nnn/GSE72094/matrix/GSE72094_series_matrix\.txt\.gz",
        URL,
    ):
        fail("unexpected source URL")

    sample_fields: dict[str, list[list[str]]] = collections.defaultdict(list)
    non_sample_fields: dict[str, list[str]] = {}
    seen_table = False
    header_lines = 0
    header_bytes = 0
    header_digest = hashlib.sha256()
    table_columns: list[str] | None = None
    request = urllib.request.Request(
        URL,
        headers={
            "User-Agent": "RLMT-DC-research-provenance-header-audit/0.1 "
                          "(NCBI-GEO-public-data)"
        },
    )

    print("=== GEO HEADER ONLY: GSE72094 ===", flush=True)
    print("SOURCE_URL=" + URL, flush=True)
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            final_url = response.geturl()
            if final_url != URL:
                fail(f"unexpected redirect; source requires approval: {final_url}")
            with gzip.GzipFile(fileobj=response, mode="rb") as decompressed:
                with io.TextIOWrapper(
                    decompressed, encoding="utf-8", errors="strict"
                ) as stream:
                    for line in stream:
                        if header_lines > 3000 or header_bytes > 40_000_000:
                            fail("metadata header exceeds bounded safety cap")
                        if line.startswith("!series_matrix_table_begin"):
                            seen_table = True
                            table_header = stream.readline()
                            table_columns = parse_geo_line(table_header)
                            break
                        header_lines += 1
                        encoded = line.encode("utf-8")
                        header_bytes += len(encoded)
                        header_digest.update(encoded)
                        if line.startswith("!Sample_"):
                            fields = parse_geo_line(line)
                            sample_fields[fields[0]].append(fields[1:])
                        elif line.startswith("!Series_") or line.startswith(
                            "!Platform_"
                        ):
                            fields = parse_geo_line(line)
                            if fields[0] in (
                                "!Series_geo_accession", "!Series_platform_id",
                                "!Series_sample_id", "!Series_title",
                                "!Platform_geo_accession",
                            ):
                                non_sample_fields.setdefault(fields[0], []).extend(
                                    fields[1:]
                                )
    except (urllib.error.URLError, OSError, UnicodeError, EOFError) as exc:
        fail(f"official GEO matrix could not be streamed: {exc}")

    if not seen_table or not table_columns:
        fail("missing expression table boundary/header: input incomplete")
    if not table_columns or table_columns[0].strip('"') != "ID_REF":
        fail("unexpected matrix expression-table column schema")
    ids = table_columns[1:]
    if not ids or any(not re.fullmatch(r"GSM[0-9]+", x) for x in ids):
        fail("matrix sample columns are not all GSM accession IDs")
    if len(ids) != len(set(ids)):
        fail("duplicate GSM IDs in matrix")
    if len(ids) != 442:
        fail(f"unexpected GEO matrix count {len(ids)}; expected 442")

    gsms = sample_fields.get("!Sample_geo_accession")
    if not gsms or len(gsms) != 1:
        fail("GEO sample accession row missing or duplicated")
    if ids != gsms[0]:
        fail("expression columns are not aligned with metadata GSM IDs")
    for name, values in sample_fields.items():
        for record in values:
            if len(record) != len(ids):
                fail(f"uneven sample metadata width for {name}")

    series_geo = non_sample_fields.get("!Series_geo_accession", [])
    if EXPECTED_ACCESSION not in series_geo:
        fail("series accession missing in header")
    platforms = non_sample_fields.get("!Series_platform_id", [])
    if "GPL15048" not in platforms:
        fail("unexpected GEO platform; expected GPL15048")

    columns = {}
    for field, values in sample_fields.items():
        for i, record in enumerate(values):
            key = field.removeprefix("!Sample_")
            if len(values) > 1:
                key += f"_{i+1}"
            columns[key] = record
    values_by_sample = [
        {key: value[row_idx] for key, value in columns.items()}
        for row_idx in range(len(ids))
    ]

    characteristic_counts: dict[str, collections.Counter] = collections.defaultdict(
        collections.Counter
    )
    mutation_fields: set[str] = set()
    for row in values_by_sample:
        for key, entry in row.items():
            if not key.startswith("characteristics_ch"):
                continue
            if ":" not in entry:
                continue
            k, v = entry.split(":", 1)
            k, v = k.strip().lower(), v.strip()
            characteristic_counts[k][v] += 1
            if "kras" in k.lower() or "mutation" in k.lower():
                mutation_fields.add(k)

    overview = {
        "accession": EXPECTED_ACCESSION,
        "source_url": URL,
        "retrieved_utc": datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat(),
        "status": "SOURCE_HEADER_ONLY_NOT_MODEL_ELIGIBLE",
        "samples_in_expression_header": len(ids),
        "sample_accession_unique": True,
        "platforms": platforms,
        "header_line_count": header_lines,
        "header_decompressed_bytes": header_bytes,
        "metadata_header_sha256_NOT_FULL_FILE": header_digest.hexdigest(),
        "full_source_file_sha256": None,
        "sample_metadata_columns": sorted(columns),
        "characteristics_name_counts": {
            k: sum(v.values()) for k, v in sorted(characteristic_counts.items())
        },
        "possible_mutation_characteristic_names_NOT_VERIFIED_TARGETS": sorted(
            mutation_fields
        ),
        "mutation_status_source_verified": False,
        "unique_patient_ids_verified": False,
        "full_expression_matrix_validated": False,
        "licensing_and_redistribution_review_complete": False,
        "cohort_overlap_review_complete": False,
        "model_training_authorized": False,
    }

    OUT.mkdir(parents=True)
    with (OUT / "sample_metadata_header.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["geo_accession", *sorted(columns)]
        )
        writer.writeheader()
        for idx, sid in enumerate(ids):
            writer.writerow({"geo_accession": sid, **values_by_sample[idx]})
    (OUT / "source_header_audit.json").write_text(
        json.dumps(overview, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (OUT / "README.txt").write_text(
        "Header-only source audit. NO expression matrix was loaded, "
        "NO patient-level independence confirmed, NO KRAS class labels "
        "created, NO predictive fitting performed.\n"
        "The SHA256 in source_header_audit.json is for a decoded metadata "
        "header and is NOT a checksum of the full GEO series matrix file.\n"
        "Review candidate phenotype keys and obtain a verified primary "
        "mutation table before training.\n",
        encoding="utf-8",
    )
    print("GEO_SAMPLE_COUNT=" + str(len(ids)))
    print("GEO_PLATFORM=GPL15048")
    print(
        "POSSIBLE_MUTATION_FIELDS="
        + (";".join(sorted(mutation_fields)) if mutation_fields else "NONE")
    )
    print("PATIENT_ID_VERIFIED=NO")
    print("TARGET_LABELS_VERIFIED=NO")
    print("MODEL_TRAINING_AUTHORIZED=NO")
    print("RLMT_DC_GEO_HEADER_INTAKE=PASS")
    print("REPORT=" + str(OUT))


if __name__ == "__main__":
    main()
