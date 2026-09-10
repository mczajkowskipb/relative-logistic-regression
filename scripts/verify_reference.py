#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

TARGETS=['compcancer_published_state','ptb_fixed_published_state','ptb_loco_published_state','legacy8_published_state_x8']
FILES=['paper_tables/table1_world_performance.csv','paper_tables/table2_world_structure.csv','paper_tables/table3_dataset_level_selected_models.csv','paper_stats/stats_tests_all.csv','paper_stats/stats_tests_paper_focus.csv']

def load(p):
    with p.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
def norm(r):
    o={}
    for k,v in r.items():
        s='' if v is None else str(v).strip()
        if k in {'p_value','p_adjust_holm','holm_p'}:
            try:s=f'{float(s):.12g}'
            except:pass
        o[k]=s
    return tuple(sorted(o.items()))

ap=argparse.ArgumentParser(description='Compare selected publication CSVs between two reference-output roots.')
ap.add_argument('--reference-root', required=True, type=Path)
ap.add_argument('--candidate-root', required=True, type=Path)
a=ap.parse_args()
all_ok=True
for state in TARGETS:
    print('\n'+state)
    for rel in FILES:
        rp=a.reference_root/state/rel; cp=a.candidate_root/state/rel
        if not rp.exists() or not cp.exists():
            print(rel, 'MISSING', rp.exists(), cp.exists()); all_ok=False; continue
        r,c=load(rp),load(cp); rs={norm(x) for x in r}; cs={norm(x) for x in c}
        ok=rs==cs
        print(rel, 'OK' if ok else 'DIFF', f'rows={len(r)}/{len(c)} only_ref={len(rs-cs)} only_candidate={len(cs-rs)}')
        all_ok &= ok
raise SystemExit(0 if all_ok else 1)
