#!/usr/bin/env python3
"""
Phase A — AlgPred 2.0 Batch Submitter (Allergenicity)
=====================================================
Submits short peptides to the AlgPred 2.0 batch server for
allergenicity prediction. Uses AAC-RF model which is more
appropriate for short peptides than Hybrid.

Usage:
    python phaseA_run_algpred.py

Output:
    data/phaseA_inputs/results_algpred.csv
    Columns: peptide_id, sequence, ml_score, prediction
"""

import csv
import os
import re
import ssl
import time
import urllib.request

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(BASE_DIR, "..", "..", "data", "phaseA_inputs"))

FASTA_MASTER = os.path.join(DATA_DIR, "moso_short_2to6.fasta")
OUTPUT_CSV = os.path.join(DATA_DIR, "results_algpred.csv")

POST_URL = "https://webs.iiitd.edu.in/raghava/algpred2/batch_action.php"
CHUNK_SIZE = 50
SLEEP_BETWEEN = 3
MAX_RETRIES = 2
TIMEOUT = 120

ssl_ctx = ssl.create_default_context()
ssl_ctx.check_hostname = False
ssl_ctx.verify_mode = ssl.CERT_NONE


def build_multipart(fields: list) -> tuple[bytes, str]:
    boundary = f"----PhaseABoundary{int(time.time())}"
    body_parts = []
    for name, value in fields:
        body_parts.append(f"--{boundary}".encode())
        body_parts.append(f'Content-Disposition: form-data; name="{name}"'.encode())
        body_parts.append(b"")
        body_parts.append(value.encode() if isinstance(value, str) else value)
    body_parts.append(f"--{boundary}--".encode())
    data = b"\r\n".join(body_parts)
    return data, boundary


def parse_results(html_text: str) -> list[dict]:
    """Parse AlgPred results table. Look for the prediction results table."""
    rows = []
    # Find the table with prediction results (typically the second table)
    # Look for <tr> rows containing prediction info
    for tr_match in re.finditer(
        r"<tr[^>]*>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>",
        html_text, re.DOTALL | re.IGNORECASE,
    ):
        cells = []
        for g in tr_match.groups():
            text = re.sub(r"<[^>]+>", "", g).strip()
            cells.append(text)

        # Validate: first cell starts with pep_ or is a short ID
        if cells[0] and (cells[0].startswith("pep_") or cells[0].isdigit()):
            rows.append({
                "peptide_id": cells[0],
                "sequence": "",  # AlgPred doesn't return sequence
                "ml_score": cells[1] if len(cells) > 1 else "",
                "prediction": cells[2] if len(cells) > 2 else "",
            })
    return rows


def submit_chunk(fasta_chunk: str, chunk_num: int, total: int) -> list[dict]:
    label = f"chunk {chunk_num}/{total}"
    fields = [
        ("name", f"batch_{chunk_num}"),
        ("seq", fasta_chunk),
        ("terminus", "0"),      # AAC-RF: more appropriate for short peptides
        ("svm_th", "0.0"),
        ("submit", "Submit"),
    ]

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            data, boundary = build_multipart(fields)
            req = urllib.request.Request(POST_URL, data=data, method="POST")
            req.add_header("Content-Type", f"multipart/form-data; boundary={boundary}")
            resp = urllib.request.urlopen(req, timeout=TIMEOUT, context=ssl_ctx)
            html = resp.read().decode("utf-8", errors="replace")

            # Look for results in the body - sometimes there's JavaScript, check content
            rows = parse_results(html)
            if rows:
                print(f"  ✓  {label}: {len(rows)} results")
                return rows

            # Maybe the result is at a different part - check if there's a redirect/refresh
            redirect = re.search(
                r"<meta\s+http-equiv=[\"']refresh[\"']\s+content=[\"']\d+;url=([^\"']+)[\"']",
                html, re.IGNORECASE
            )
            if redirect:
                time.sleep(2)
                url = redirect.group(1)
                if not url.startswith("http"):
                    url = "https://webs.iiitd.edu.in/raghava/algpred2/" + url.lstrip("/")
                resp2 = urllib.request.urlopen(url, timeout=TIMEOUT, context=ssl_ctx)
                html = resp2.read().decode("utf-8", errors="replace")
                rows = parse_results(html)
                if rows:
                    print(f"  ✓  {label}: {len(rows)} results (after redirect)")
                    return rows

            print(f"  ⚠  {label}: no results (attempt {attempt}), HTML size={len(html)}")
            if attempt < MAX_RETRIES:
                time.sleep(5)

        except Exception as e:
            print(f"  ⚠  {label} (attempt {attempt}): {e}")
            time.sleep(5)

    print(f"  ✗  {label}: FAILED after {MAX_RETRIES} attempts")
    return []


def main():
    print("=" * 60)
    print("Phase A — AlgPred 2.0 Batch Submitter (Allergenicity)")
    print("=" * 60)

    if not os.path.exists(FASTA_MASTER):
        print(f"ERROR: {FASTA_MASTER} not found")
        return

    with open(FASTA_MASTER) as f:
        fasta_lines = f.read().strip().split("\n")

    # Split into chunks
    chunks = []
    current = []
    count = 0
    for line in fasta_lines:
        current.append(line)
        if line.startswith(">"):
            count += 1
            if count % CHUNK_SIZE == 0:
                chunks.append("\n".join(current))
                current = []
    if current:
        chunks.append("\n".join(current))

    total_chunks = len(chunks)
    print(f"Master: {len(fasta_lines)//2} peptides, {total_chunks} chunks of {CHUNK_SIZE}\n")

    all_results = []
    for i, chunk in enumerate(chunks, 1):
        print(f"[chunk {i}/{total_chunks}] Submitting...")
        rows = submit_chunk(chunk, i, total_chunks)
        all_results.extend(rows)
        if i < total_chunks:
            time.sleep(SLEEP_BETWEEN)

    print(f"\nTotal results: {len(all_results)}")

    # Dedup
    seen = {}
    for r in all_results:
        key = r["peptide_id"]
        seen[key] = r
    deduped = list(seen.values())
    print(f"After dedup: {len(deduped)} unique")

    if deduped:
        os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
        with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["peptide_id", "sequence", "ml_score", "prediction"])
            writer.writeheader()
            for r in deduped:
                writer.writerow(r)
        print(f"Saved → {OUTPUT_CSV}")
    else:
        print("No results obtained. Check server availability.")

    print("AlgPred phase complete.")


if __name__ == "__main__":
    main()
