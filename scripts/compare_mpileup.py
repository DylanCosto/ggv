"""Compare `ggv query pileup --format jsonl` against `samtools mpileup` position by position.

usage: compare_mpileup.py GGV_EXE BAM_WINDOWS_PATH LOCUS [REF_WINDOWS_PATH]
samtools runs in WSL with the settings GGV's defaults are meant to match.
"""
import json
import subprocess
import sys
from collections import Counter


def wsl_path(p):
    p = p.replace("\\", "/")
    return f"/mnt/{p[0].lower()}{p[2:]}"


def parse_bases(s):
    c = Counter()
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "^":
            i += 2
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            j = i + 1
            while s[j].isdigit():
                j += 1
            n = int(s[i + 1 : j])
            seq = s[j : j + n]
            if ch == "+":
                c["ins_" + ("rev" if seq.islower() else "fwd")] += 1
            i = j + n
            continue
        if ch in "><":
            c["skip"] += 1
        elif ch == "*":
            c["DEL_fwd"] += 1
        elif ch == "#":
            c["DEL_rev"] += 1
        elif ch == ".":
            c["REF_fwd"] += 1
        elif ch == ",":
            c["REF_rev"] += 1
        else:
            b = ch.upper()
            if b not in "ACGT":
                b = "N"
            c[b + ("_rev" if ch.islower() else "_fwd")] += 1
        i += 1
    return c


def main():
    ggv, bam, locus = sys.argv[1:4]
    ref = sys.argv[4] if len(sys.argv) > 4 else None

    cmd = ["samtools", "mpileup", "-A", "-B", "-Q", "0", "-q", "0", "-d", "0", "-x",
           "--reverse-del", "--ff", "UNMAP,SECONDARY,QCFAIL,DUP", "-r", locus]
    if ref:
        cmd += ["-f", wsl_path(ref)]
    cmd.append(wsl_path(bam))
    sam = subprocess.run(["wsl.exe", "-e", *cmd], capture_output=True, text=True, check=True).stdout

    ggv_cmd = [ggv, "query", "pileup", bam, locus, "--format", "jsonl"]
    if ref:
        ggv_cmd += ["-f", ref]
    ours = {}
    for line in subprocess.run(ggv_cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        s = json.loads(line)
        ours[s["pos"]] = s

    checked = mismatches = 0
    seen = set()
    for line in sam.splitlines():
        chrom, pos, refb, depth, bases = line.split("\t")[:5]
        pos = int(pos)
        seen.add(pos)
        c = parse_bases(bases)
        s = ours.get(pos)
        problems = []
        if s is None:
            problems.append("missing in GGV")
        else:
            # samtools prints '.'/',' for reference matches; attribute them to the ref base.
            rb = refb.upper()
            if c["REF_fwd"] or c["REF_rev"]:
                key = rb if rb in "ACGT" else "N"
                c[key + "_fwd"] += c.pop("REF_fwd", 0)
                c[key + "_rev"] += c.pop("REF_rev", 0)
            if int(depth) != s["depth"] + s["ref_skips"]:
                problems.append(f"depth sam={depth} ggv={s['depth']}+{s['ref_skips']}")
            if c["skip"] != s["ref_skips"]:
                problems.append(f"ref_skips sam={c['skip']} ggv={s['ref_skips']}")
            for b in "ACGTN":
                for strand in ("fwd", "rev"):
                    if c[f"{b}_{strand}"] != s["bases"][b][strand]:
                        problems.append(f"{b}_{strand} sam={c[f'{b}_{strand}']} ggv={s['bases'][b][strand]}")
            for strand in ("fwd", "rev"):
                if c[f"DEL_{strand}"] != s["deletions"][strand]:
                    problems.append(f"DEL_{strand} sam={c[f'DEL_{strand}']} ggv={s['deletions'][strand]}")
                if c[f"ins_{strand}"] != s["insertions_after"][strand]:
                    problems.append(f"INS_{strand} sam={c[f'ins_{strand}']} ggv={s['insertions_after'][strand]}")
            if ref and (s["ref"] or "N").upper() != rb:
                problems.append(f"ref sam={rb} ggv={s['ref']}")
        checked += 1
        if problems:
            mismatches += 1
            if mismatches <= 10:
                print(f"{chrom}:{pos}: " + "; ".join(problems))

    extra = sorted(set(ours) - seen)
    print(f"{bam.split('/')[-1]} {locus} ref={'yes' if ref else 'no'}: "
          f"{checked} positions checked, {mismatches} mismatched, {len(extra)} only in GGV")
    if extra[:5]:
        print("  e.g. only in ggv:", extra[:5])
    sys.exit(1 if mismatches or extra else 0)


main()
