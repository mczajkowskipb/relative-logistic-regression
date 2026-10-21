#!/usr/bin/env python3
"""Strict, no-overwrite launcher for independent 95-fold TSPDT-inspired run.

Preflight is read-only; execution launches detached R and uses a separate
TSPDT worktree/result root. Never alters the 285 RLR/RLMT fit artifacts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BRANCH = "feature/tspdt-baseline-evidence"
MAIN = "4079a5ac400f1d2cfded523c7000f655470925ce"
CFG = ROOT / "configs/tspdt/tspdt_reconstruction_19x5_v1.json"
OUT = ROOT / "results/tspdt_recon_19x5_v1"
ORIGINAL_SHA = "aa4dae241ba583394501f2f11efd757c8d74ad46"
REFERENCE_OUT = "results/rlmt_dc_published_data_5fold_v1"
AUDITS = {
    "published_5fold_audit_v1/independent_fold_metrics.csv": 285,
    "published_5fold_firstfold_parity_v1/firstfold_prediction_parity.csv": 57,
    "published_5fold_summary_v1/per_fold.csv": 95,
    "published_5fold_summary_v1/per_dataset.csv": 19,
}

def stop(message: str) -> None:
    raise SystemExit("TSPDT_LAUNCH_STOP: " + message)

def git(*args: str, cwd: Path = ROOT) -> str:
    return subprocess.run(["git", *args], cwd=cwd, check=True,
                          capture_output=True, text=True).stdout.strip()

def hashfile(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def csvrows(path: Path) -> list[dict[str,str]]:
    if not path.is_file():
        stop("missing CSV: "+str(path))
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-sha",required=True)
    parser.add_argument("--reference-repo",required=True,
                        help="Full absolute path of completed RLR worktree")
    parser.add_argument("--preflight-only",action="store_true")
    args=parser.parse_args()
    ref=Path(args.reference_repo).resolve(strict=True)
    if Path.cwd().resolve()!=ROOT:
        stop("run from separate TSPDT Git worktree root")
    expected=args.expected_sha.lower()
    if not re.fullmatch(r"[0-9a-f]{40}",expected):
        stop("invalid expected Git SHA")
    if git("rev-parse","HEAD")!=expected:
        stop("TSPDT HEAD differs from registered execution SHA")
    if git("rev-parse","origin/"+BRANCH)!=expected:
        stop("remote TSPDT branch differs from pinned SHA")
    if git("rev-parse","origin/main")!=MAIN:
        stop("published main moved")
    if git("status","--porcelain","--untracked-files=no"):
        stop("working tree contains tracked changes")
    if ref==ROOT or not (ref/".git").exists():
        stop("reference repository missing or is not separate")
    if git("rev-parse","origin/main",cwd=ref)!=MAIN:
        stop("reference repo publication main mismatch")
    # No simultaneous writes to original RLR/RLMT experiment.
    refout=ref/REFERENCE_OUT
    status=refout/"status.txt"
    if not status.is_file() or status.read_text().strip()!="RUNNER_EXIT=0":
        stop("285-fit RLR/RLMT experiment still incomplete")
    if not (refout/"CODE_SHA.txt").is_file():
        stop("original RLR 285-fit CODE_SHA missing")
    code_sha=(refout/"CODE_SHA.txt").read_text().strip()
    if not re.fullmatch(r"[0-9a-f]{40}",code_sha):
        stop("original RLR code SHA malformed")
    if (refout/"PUBLICATION_MAIN_SHA.txt").read_text().strip()!=MAIN:
        stop("reference results not linked to publication main")
    if not (refout/"published_5fold_accepted_replay_v1/accepted_models.csv").is_file():
        stop("serialized RLMT accepted-tree replay missing")
    for rel,count in AUDITS.items():
        rows=csvrows(refout/rel)
        if len(rows)!=count:
            stop(f"reference audit {rel} has {len(rows)} rows, expected {count}")
    desc=refout/"published_5fold_summary_v1/descriptive_aggregates.json"
    if not desc.exists():
        stop("dataset-level RLR/RLMT summary missing")
    summary=json.loads(desc.read_text())
    if any(summary.get(key)!=val for key,val in {
        "n_datasets":19,"n_fold_triplets":95,"n_model_fits":285,
        "n_fits_with_error":0}.items()):
        stop("RLR/RLMT summary not fully complete")
    if OUT.exists():
        stop("refusing to overwrite existing independent TSPDT result root")
    cfg=json.loads(CFG.read_text())
    if (cfg.get("stage")!="TSPDT_2011_INSPIRED_95FOLD_EXPLORATORY_V1"
        or cfg.get("model")!="TSPDT2011_INSPIRED_RECONSTRUCTION_V1"
        or cfg.get("publication_main_sha")!=MAIN
        or cfg.get("reference_result_root")!=REFERENCE_OUT
        or cfg.get("output_root")!="results/tspdt_recon_19x5_v1"
        or cfg.get("dataset_count")!=19
        or cfg.get("outer_folds_per_dataset")!=5
        or cfg.get("total_fits")!=95
        or cfg.get("fold_seed")!=42
        or cfg.get("top_n")!=100
        or cfg.get("min_leaf")!=5
        or cfg.get("max_depth")!=3
        or cfg.get("pruning")!="wilson"
        or cfg.get("confidence")!=0.25
        or cfg.get("laplace_alpha")!=1):
        stop("registered TSPDT model hyperparameters differ")
    registered=[(shard["world"],dataset)
                for shard in cfg["datasets"] for dataset in shard["names"]]
    if len(registered)!=19 or len(set(registered))!=19:
        stop("not 19 independently registered datasets")
    original=csvrows(refout/"input_manifest.csv")
    if (len(original)!=19 or
        [(r["world"],r["dataset"]) for r in original]!=registered):
        stop("reference source inventory differs from registered 19 datasets")
    all_models=csvrows(refout/"run_manifest.csv")
    if len(all_models)!=285 or any(r["status"]!="OK" for r in all_models):
        stop("original reference 285 model attempts not all successful")
    canonical=git("show",MAIN+":checksums/DATA_SHA256_EXPECTED.txt",
                  cwd=ref)
    hashes={}
    for rawline in canonical.splitlines():
        match=re.fullmatch(r"([a-fA-F0-9]{64})\s+(\S+)",rawline.strip())
        if not match: stop("invalid publication input hash manifest")
        hashes[match.group(2)]=match.group(1).lower()
    locked=[]
    for r in original:
        name=r["dataset"]
        collection="comparison" if r["world"]=="legacy8" else "CompCancer"
        key=f"data/{collection}/{name}.arff"
        expected_file_hash=hashes.get(key)
        if expected_file_hash!=r["sha256"]:
            stop("reference input SHA differs from publication: "+name)
        source=Path(r["source_path"]).resolve(strict=True)
        if hashfile(source)!=expected_file_hash:
            stop("source ARFF changed since published run: "+name)
        locked.append({
            "world":r["world"],"dataset":name,
            "source_path":str(source),"sha256":expected_file_hash,
            "paper_manifest_path":key
        })
    # Link the synthetic test evidence to the exact source code bytes.
    gate=ROOT/"results/tspdt_reconstruction_v1_synthetic_gate"
    if (not (gate/"status.txt").is_file() or
        (gate/"status.txt").read_text().strip()!="TSPDT_DEV_EXIT=0" or
        not (gate/"test.log").is_file() or
        "TSPDT_RECON_DEV_GATE=PASS" in (gate/"test.log").read_text()):
        # The shell writes the final gate PASS to stdout; test.log
        # contains TSPDT_RECON_SYNTHETIC=PASS, not DEV_GATE.
        pass
    if (not (gate/"status.txt").is_file() or
        (gate/"status.txt").read_text().strip()!="TSPDT_DEV_EXIT=0" or
        not (gate/"test.log").is_file() or
        "TSPDT_RECON_SYNTHETIC=PASS" not in
        (gate/"test.log").read_text()):
        stop("independent synthetic TSPDT DEV gate has not passed")
    old=git("show",ORIGINAL_SHA+":R/baselines/tspdt_reconstruction_v1.R")
    curr=(ROOT/"R/baselines/tspdt_reconstruction_v1.R").read_text()
    if old.rstrip("\n")!=curr.rstrip("\n"):
        stop("TSPDT algorithm changed since successful synthetic DEV")
    subprocess.run([
        "Rscript","--vanilla","-e",
        'parse(file="R/baselines/tspdt_reconstruction_v1.R");'
        'parse(file="R/baselines/50_tspdt_published_19x5_run.R");'
        'parse(file="R/baselines/51_tspdt_published_19x5_audit.R")'
    ],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    print("TSPDT_PUBLISHED_R_PARSE=PASS",flush=True)
    print("TSPDT_PUBLISHED_PREFLIGHT=PASS",flush=True)
    print("SOURCE_DATASETS=19",flush=True)
    print("MATCHED_OUTER_FOLDS=95",flush=True)
    print("TSPDT_MODEL_ATTEMPTS=95",flush=True)
    print("PAST_RLR_RLMT_FITS_REUSED_READ_ONLY=285",flush=True)
    print("BIGTSP=NOT_INCLUDED_NOT_INSTALLED",flush=True)
    if args.preflight_only:
        print("TSPDT_OUTPUT_DIRECTORY_NOT_CREATED=TRUE",flush=True)
        return

    OUT.mkdir(parents=True)
    with (OUT/"input_manifest.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(locked[0]))
        writer.writeheader()
        writer.writerows(locked)
    (OUT/"CODE_SHA.txt").write_text(expected+"\n")
    (OUT/"REFERENCE_MODEL_SHA.txt").write_text(code_sha+"\n")
    (OUT/"PUBLICATION_MAIN_SHA.txt").write_text(MAIN+"\n")
    (OUT/"CONFIG_SHA256.txt").write_text(hashfile(CFG)+"\n")
    (OUT/"REFERENCE_REPO.txt").write_text(str(ref)+"\n")
    (OUT/"README.txt").write_text(
        "2011 TSPDT-INSPIRED RECONSTRUCTION, not historical original. "
        "19 published RLR datasets x five paired outer folds. "
        "Same sample IDs as completed 285-fit comparison. "
        "Train-only MAD f100, depth3, minleaf5, Wilson cf0.25, "
        "Laplace alpha1 scores. NEVER independent external validation.\n")
    runner=OUT/"run.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        f"cd '{ROOT}' || exit 2\n"
        f"if [ \"$(git rev-parse HEAD)\" != \"{expected}\" ] || "
        f"[ -n \"$(git status --porcelain --untracked-files=no)\" ]; then\n"
        f"  echo RUNNER_EXIT=90 > '{OUT}/status.txt'\n"
        "  exit 90\n"
        "fi\n"
        "Rscript --vanilla R/baselines/50_tspdt_published_19x5_run.R "
        f"'{CFG}' '{ref}' > '{OUT}/run.log' 2>&1\n"
        "rc=$?\n"
        f"printf 'RUNNER_EXIT=%s\\n' \"$rc\" > '{OUT}/status.txt'\n"
        "exit \"$rc\"\n")
    with (OUT/"launcher.log").open("w") as handle:
        proc=subprocess.Popen(["bash",str(runner)],cwd=ROOT,
                              stdin=subprocess.DEVNULL,stdout=handle,
                              stderr=subprocess.STDOUT,start_new_session=True)
    (OUT/"runner.pid").write_text(str(proc.pid)+"\n")
    print("TSPDT_PUBLISHED_STARTED_PID="+str(proc.pid),flush=True)
    print("TSPDT_RUN_LOG="+str(OUT/"run.log"),flush=True)
    print("TSPDT_MODEL_MANIFEST="+str(OUT/"run_manifest.csv"),flush=True)

if __name__=="__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        stop("preflight subprocess failed: "+str(e))
