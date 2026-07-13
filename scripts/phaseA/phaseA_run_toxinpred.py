#!/usr/bin/env python3
"""
Phase A — ToxinPred3 Batch Submitter (v2)
==========================================
Submits short peptides (2-6aa) to the ToxinPred batch web server.
Server limit: ~50 peptides per submission.
Splits the master input into chunks of 50, submits sequentially,
saves aggregated results as results_toxinpred.csv.

Usage:
    python phaseA_run_toxinpred.py

Output:
    data/phaseA_inputs/results_toxinpred.csv
    Columns: peptide_id, sequence, svm_score, prediction
"""

import csv
import os
import re
import ssl
import time
import urllib.parse
import urllib.request

# ── Config ────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(BASE_DIR, "..", "..", "data", "phaseA_inputs"))

FASTA_MASTER = os.path.join(DATA_DIR, "moso_short_2to6.fasta")
OUTPUT_CSV = os.path.join(DATA_DIR, "results_toxinpred.csv")

POST_URL = "http://www.raghavagps.net/raghava/toxinpred/multiple_test.php"
CHUNK_SIZE = 50            # peptides per submission — 50 works reliably
SLEEP_BETWEEN = 3          # seconds between chunks
MAX_RETRIES = 2
TIMEOUT = 120

# SSL — ToxinPred uses HTTP, but we may follow redirects
ssl_ctx = ssl.create_default_context()
ssl_ctx.check_hostname = False
ssl_ctx.verify_mode = ssl.CERT_NONE


def build_form_data(fasta_chunk: str) -> bytes:
    """Build URL-encoded POST body."""
    params = {
        "seq": fasta_chunk,
        "method": "1",
        "thval": "0.0",
        "eval": "0.0001",
        "field[]": ["4", "7", "9"],
    }
    return urllib.parse.urlencode(params, doseq=True).encode()


def extract_redirect(html_text: str) -> str | None:
    m = re.search(
        r'<meta\s+http-equiv=[\'"]refresh[\'"]\s+content=[\'"]\d+;url=([^\'"]+)[\'"]',
        html_text, re.IGNORECASE,
    )
    if m:
        url = m.group(1)
        if url.startswith("http"):
            return url
        # Relative path — construct absolute URL preserving query string
        parsed = urllib.parse.urlparse(POST_URL)
        base_path = "/".join(parsed.path.split("/")[:-1]) + "/"
        return f"{parsed.scheme}://{parsed.netloc}{base_path}{url}"
    return None


def parse_results(html_text: str) -> list[dict]:
    """
    Parse ToxinPred result table.
    Columns: Peptide ID, Peptide Sequence, SVM Score, Prediction
    Some rows have the sequence embedded in an <a> tag.
    """
    rows = []
    # Find the tbody
    tbody_match = re.search(
        r"<tbody[^>]*>(.*?)</tbody>", html_text, re.DOTALL | re.IGNORECASE
    )
    if not tbody_match:
        return []

    tbody = tbody_match.group(1)

    # Find all <tr> inside tbody
    for tr_match in re.finditer(
        r"<tr[^>]*>(.*?)</tr>", tbody, re.DOTALL | re.IGNORECASE
    ):
        tds = re.findall(
            r"<td[^>]*>(.*?)</td>", tr_match.group(1), re.DOTALL | re.IGNORECASE
        )
        if len(tds) >= 4:
            # Strip all HTML tags from each field
            fields = []
            for td in tds[:4]:
                text = re.sub(r"<[^>]+>", "", td).strip()
                # Decode HTML entities
                text = text.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
                fields.append(text)

            pid, seq, score, pred = fields
            if pid or seq:
                rows.append({
                    "peptide_id": pid,
                    "sequence": seq.upper(),
                    "svm_score": score,
                    "prediction": pred,
                })
    return rows


def submit_chunk(fasta_chunk: str, label: str) -> list[dict]:
    """Submit one chunk to ToxinPred, return parsed results."""
    data = build_form_data(fasta_chunk)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            req = urllib.request.Request(POST_URL, data=data)
            resp = urllib.request.urlopen(req, timeout=TIMEOUT)
            html = resp.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"  ⚠  {label} (attempt {attempt}/{MAX_RETRIES}): POST failed — {e}")
            time.sleep(5)
            continue

        # Follow meta-refresh redirect
        redirect_url = extract_redirect(html)
        if redirect_url:
            time.sleep(2)
            try:
                resp2 = urllib.request.urlopen(redirect_url, timeout=TIMEOUT)
                html = resp2.read().decode("utf-8", errors="replace")
            except Exception as e:
                print(f"  ⚠  {label} (attempt {attempt}/{MAX_RETRIES}): redirect failed — {e}")
                time.sleep(5)
                continue

        rows = parse_results(html)
        if rows:
            print(f"  ✓  {label}: {len(rows)} results")
            return rows
        else:
            print(f"  ⚠  {label}: no results table (attempt {attempt})")
            time.sleep(5)

    print(f"  ✗  {label}: FAILED after {MAX_RETRIES} attempts")
    return []


def main():
    print("=" * 60)
    print("Phase A — ToxinPred Batch Submitter v2")
    print("=" * 60)

    # Read master FASTA
    if not os.path.exists(FASTA_MASTER):
        print(f"ERROR: {FASTA_MASTER} not found")
        return

    with open(FASTA_MASTER) as f:
        fasta_lines = f.read().strip().split("\n")

    # Split into chunks of CHUNK_SIZE (each peptide = 2 lines: header + seq)
    # FASTA format: >header\nSEQUENCE\n
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

    print(f"Master: {len(fasta_lines)//2} peptides, split into {len(chunks)} chunks of {CHUNK_SIZE}\n")

    all_results = []
    for i, chunk in enumerate(chunks, 1):
        chunk_label = f"chunk {i}/{len(chunks)} ({len(chunk)} chars)"
        print(f"[{chunk_label}] Submitting...")
        rows = submit_chunk(chunk, chunk_label)
        all_results.extend(rows)
        if i < len(chunks):
            time.sleep(SLEEP_BETWEEN)

    print(f"\nTotal results: {len(all_results)} / {len(fasta_lines)//2}")

    # Deduplicate (keep last entry per sequence)
    seen = {}
    for r in all_results:
        seen[r["sequence"]] = r
    deduped = list(seen.values())
    print(f"After dedup: {len(deduped)} unique peptides")

    # Write CSV
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["peptide_id", "sequence", "svm_score", "prediction"])
        writer.writeheader()
        for r in deduped:
            writer.writerow(r)
    print(f"Saved → {OUTPUT_CSV}")
    print("ToxinPred phase complete.")


if __name__ == "__main__":
    main()
