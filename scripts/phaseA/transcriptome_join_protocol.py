# -*- coding: utf-8 -*-
"""
Idea 3 transcriptome expression layer — quantitative-join protocol
(for environments with a BLAST / MMseqs2 host)
====================================================================
This script is **not run in the current sandbox** (no BLAST, no large-file
download capability). It packages the per-candidate shoot TPM quantitative join
as a one-click reproducible workflow, to be completed in a resource-rich env.

Two quantitative-join paths:
  PATH A (de novo shoot transcriptome, e.g. Peng 2013 Suppl / SRA):
      use our 253 proteins (query) to sequence-align against the shoot
      transcriptome unigenes (library); take the best-hit identity/coverage
      as evidence that "this source protein is detectable in the shoot transcriptome";
      if the transcriptome provides TPM/FPKM, carry the expression level through.
  PATH B (genome gene-model-keyed shoot expression matrix, e.g. MDPI 2020 Suppl):
      requires a UniProt accession -> PH0100... gene-model mapping table
      (currently missing, see below); once available, join per-source-protein TPM by gene model.

Note: UniProt TrEMBL (tax 38705) entries currently do **not** carry PH0100...
gene-model cross-references, so PATH B first needs an accession<->PH0100
mapping (derivable from the genome GFF proteinId<->geneId relation, or
reverse-derived from PATH A best-hits).

Dependencies: mmseqs2 (recommended) or blast+; both must be on PATH.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA_253 = os.path.join(HERE, "..", "..", "data", "moso_253.fasta")
# the following paths must be supplied by the runing env:
SHOOT_TRANSCRIPTOME_FASTA = os.environ.get("SHOOT_FASTA", "SHOOT_TRANSCRIPTOME.fasta")
SHOOT_TPM_MATRIX = os.environ.get("SHOOT_TPM", "SHOOT_TPM.tsv")  # optional, PATH B
OUT = os.path.join(HERE, "..", "..", "data", "phaseA", "transcriptome_join_result.tsv")


def run(cmd):
    print("RUN:", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr)
        sys.exit(r.returncode)
    return r.stdout


def path_a_mmseqs():
    """Use mmseqs2 for protein-level alignment against the de novo shoot transcriptome."""
    # 1) build databases
    run(["mmseqs", "createdb", SHOOT_TRANSCRIPTOME_FASTA, "shoot_db"])
    run(["mmseqs", "createdb", FASTA_253, "query_db"])
    # 2) search (protein sequence alignment, sensitive mode)
    run(["mmseqs", "search", "query_db", "shoot_db", "result_db", "tmp",
         "--sensitivity", "5.7", "-e", "1e-5"])
    # 3) convert to TSV
    run(["mmseqs", "convertalis", "query_db", "shoot_db", "result_db", "align.tsv",
         "--format-output", "query,target,pident,alnlen,qlen,tlen,bits"])
    print("-> align.tsv generated; per query (= our accession) take the best-hit as the "
          "shoot-detection evidence for that source protein.")


def path_b_genome_model():
    """Join a genome-gene-model-keyed shoot expression matrix (needs accession<->PH0100 mapping)."""
    if not os.path.exists(SHOOT_TPM_MATRIX):
        print("Missing SHOOT_TPM matrix (env var SHOOT_TPM). PATH B unavailable.")
        return
    # assume acc2model.tsv exists: accession \t PH0100...G....
    acc2model = os.environ.get("ACC2MODEL", "acc2model.tsv")
    if not os.path.exists(acc2model):
        print("Missing accession<->PH0100 mapping (ACC2MODEL). Must first derive from genome GFF.")
        return
    # pseudo-code: read acc2model -> {acc: model}; read SHOOT_TPM (model-keyed)
    # -> write {acc, model, shoot_tpm_mean}
    print("Execute PATH B join (reference implementation in script comments).")


if __name__ == "__main__":
    print("This script is a protocol template, not run in the sandbox.")
    print("Usage: on an mmseqs2 host -> set SHOOT_FASTA=shoot_transcriptome.fasta then python transcriptome_join_protocol.py")
    if os.path.exists(SHOOT_TRANSCRIPTOME_FASTA) and SHOOT_TRANSCRIPTOME_FASTA != "SHOOT_TRANSCRIPTOME.fasta":
        path_a_mmseqs()
    else:
        print("SHOOT_FASTA not detected -> protocol print only. PATH A command example:")
        print("  mmseqs createdb shoot_transcriptome.fasta shoot_db && mmseqs createdb",
              FASTA_253, "query_db && mmseqs search query_db shoot_db result_db tmp -e 1e-5")
