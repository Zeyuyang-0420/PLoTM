#!/usr/bin/env python3
"""export_network_data.py — slim the explorer's demo_data.json down to just what the React network
app needs, and write it to src/networks.json (bundled at build time so the app is self-contained).

    python3 scripts/export_network_data.py [--disease RA]
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, "..", "..", "PLoTM_explorer", "data", "demo_data.json"))
OUT = os.path.abspath(os.path.join(HERE, "..", "src", "networks.json"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disease", default=None, help="disease id; defaults to the first in demo_data.json")
    ap.add_argument("--src", default=SRC)
    args = ap.parse_args()

    with open(args.src) as fh:
        demo = json.load(fh)
    did = args.disease if (args.disease and args.disease in demo["diseases"]) else next(iter(demo["diseases"]))
    d = demo["diseases"][did]
    nets = d["networks"]
    groups = [{"id": g, "label": nets["rest"][g]["label"]} for g in nets["rest"]]
    out = {
        "disease": d["id"],
        "name": d["name"],
        "cell_type": d["cell_type"],
        "conditions": d["conditions"],
        "groups": groups,
        "networks": nets,   # { cond: { group: { label, nodes, edges } } }
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    ng = sum(len(nets[c][g]["nodes"]) for c in nets for g in nets[c])
    ne = sum(len(nets[c][g]["edges"]) for c in nets for g in nets[c])
    print(f"[export] {did}: {len(groups)} groups x {len(d['conditions'])} conditions, "
          f"{ng} nodes / {ne} edges -> {OUT} ({os.path.getsize(OUT)//1024} KB)")


if __name__ == "__main__":
    main()
