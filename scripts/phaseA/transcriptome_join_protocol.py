# -*- coding: utf-8 -*-
"""
思路③ 转录组表达层 — 定量连接协议（留给有 BLAST / MMseqs2 主机的环境）
====================================================================
本脚本【不在当前沙箱运行】(无 BLAST、无大文件下载能力)。它把"逐候选笋 TPM"
的定量连接封装为一键可复现流程, 供有计算资源的环境补齐。

两种定量连接路径:
  PATH A (de novo 笋转录组, 如 Peng 2013 Suppl / SRA):
      用我们的 253 条蛋白(查询) 对 笋转录组 unigene(库) 做序列比对,
      取 best-hit 的 identity/coverage 作为"该源蛋白在笋转录组中可被检出"的证据;
      若转录组提供 TPM/FPKM, 则直接带出表达量。
  PATH B (基因组基因模型为键的笋表达矩阵, 如 MDPI 2020 Suppl):
      需要 UniProt accession -> PH0100... 基因模型的映射表(当前缺失, 见下);
      一旦有该映射, 按基因模型 join 出逐源蛋白 TPM。

注意: UniProt TrEMBL (tax 38705) 条目当前【不含】PH0100... 基因模型交叉引用,
因此 PATH B 还需先生成 accession<->PH0100 映射(可由基因组 GFF 的
proteinId<->geneId 关系派生, 或 PATH A 的 best-hit 反向得到)。

依赖: mmseqs2 (推荐) 或 blast+; 均需在 PATH 中。
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA_253 = os.path.join(HERE, "..", "..", "data", "moso_253.fasta")
# 以下路径需在使用环境提供:
SHOOT_TRANSCRIPTOME_FASTA = os.environ.get("SHOOT_FASTA", "SHOOT_TRANSCRIPTOME.fasta")
SHOOT_TPM_MATRIX = os.environ.get("SHOOT_TPM", "SHOOT_TPM.tsv")  # 可选, PATH B
OUT = os.path.join(HERE, "..", "..", "data", "phaseA", "transcriptome_join_result.tsv")


def run(cmd):
    print("RUN:", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr)
        sys.exit(r.returncode)
    return r.stdout


def path_a_mmseqs():
    """用 mmseqs2 对 de novo 笋转录组做蛋白级比对。"""
    # 1) 建库
    run(["mmseqs", "createdb", SHOOT_TRANSCRIPTOME_FASTA, "shoot_db"])
    run(["mmseqs", "createdb", FASTA_253, "query_db"])
    # 2) 搜索 (蛋白序列比对, 敏感模式)
    run(["mmseqs", "search", "query_db", "shoot_db", "result_db", "tmp",
         "--sensitivity", "5.7", "-e", "1e-5"])
    # 3) 转 TSV
    run(["mmseqs", "convertalis", "query_db", "shoot_db", "result_db", "align.tsv",
         "--format-output", "query,target,pident,alnlen,qlen,tlen,bits"])
    print("-> align.tsv 生成; 用其按 query(=我们的 accession) 取 best-hit 即源蛋白笋检出证据。")


def path_b_genome_model():
    """基因组基因模型为键的笋表达矩阵 join (需 accession<->PH0100 映射)。"""
    if not os.path.exists(SHOOT_TPM_MATRIX):
        print("缺少 SHOOT_TPM 矩阵 (环境变量 SHOOT_TPM)。PATH B 不可用。")
        return
    # 假设存在 acc2model.tsv: accession \t PH0100...G....
    acc2model = os.environ.get("ACC2MODEL", "acc2model.tsv")
    if not os.path.exists(acc2model):
        print("缺少 accession<->PH0100 映射 (ACC2MODEL)。需先由基因组 GFF 派生。")
        return
    # 简易 join 伪代码: 读 acc2model -> {acc: model}; 读 SHOOT_TPM (model keyed)
    # -> 写出 {acc, model, shoot_tpm_mean}
    print("执行 PATH B join (参考实现见脚本注释)。")


if __name__ == "__main__":
    print("本脚本为协议模板, 不在沙箱运行。")
    print("用法: 在有 mmseqs2 的主机 -> 设 SHOOT_FASTA=笋转录组.fasta 后 python transcriptome_join_protocol.py")
    if os.path.exists(SHOOT_TRANSCRIPTOME_FASTA) and SHOOT_TRANSCRIPTOME_FASTA != "SHOOT_TRANSCRIPTOME.fasta":
        path_a_mmseqs()
    else:
        print("未检测到 SHOOT_FASTA -> 仅打印协议。PATH A 命令示例:")
        print("  mmseqs createdb 笋转录组.fasta shoot_db && mmseqs createdb",
              FASTA_253, "query_db && mmseqs search query_db shoot_db result_db tmp -e 1e-5")
