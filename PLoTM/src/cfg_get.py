#!/usr/bin/env python
"""cfg_get.py <config.yaml> <dotted.key> — emit a config value for shell consumption.

Special keys:
  perturbseq.conditions   -> one condition per line
  lof.traits.names        -> one trait name per line
Otherwise a dotted path into the YAML (scalars printed as-is; lists one-per-line).
Kept dependency-light (pyyaml only) so the shell orchestrator can read the single source of truth.
"""
import sys

def _load(path):
    """Prefer pyyaml; fall back to a minimal indentation parser (no external deps) so the shell
    orchestrator can read the config in any environment."""
    try:
        import yaml
        with open(path) as fh:
            return yaml.safe_load(fh)
    except ImportError:
        return _mini_yaml(path)

def _mini_yaml(path):
    import re
    root = {}
    stack = [(-1, root)]
    with open(path) as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].rstrip()
            if not line.strip():
                continue
            indent = len(line) - len(line.lstrip())
            s = line.strip()
            while stack and indent <= stack[-1][0]:
                stack.pop()
            parent = stack[-1][1]
            if s.startswith("- "):                       # list item (scalar or inline dict)
                parent.setdefault("__list__", [])
                item = s[2:].strip()
                if item.startswith("{"):                  # {k: v, ...}
                    d = {}
                    for kv in item.strip("{}").split(","):
                        if ":" in kv:
                            k, v = kv.split(":", 1); d[k.strip()] = _scalar(v.strip())
                    parent["__list__"].append(d)
                else:
                    parent["__list__"].append(_scalar(item))
            elif ":" in s:
                k, v = s.split(":", 1); k = k.strip(); v = v.strip()
                if v == "":                               # nested mapping/list follows
                    child = {}
                    parent[k] = child
                    stack.append((indent, child))
                elif v.startswith("["):                   # inline list
                    parent[k] = [_scalar(x.strip()) for x in v.strip("[]").split(",") if x.strip()]
                else:
                    parent[k] = _scalar(v)
    return _delist(root)

def _scalar(v):
    v = v.strip().strip('"').strip("'")
    if re.fullmatch(r"-?\d+", v): return int(v)
    if re.fullmatch(r"-?\d*\.\d+", v): return float(v)
    if v in ("true", "True"): return True
    if v in ("false", "False"): return False
    return v

def _delist(node):
    if isinstance(node, dict):
        if set(node.keys()) == {"__list__"}:
            return [_delist(x) for x in node["__list__"]]
        return {k: _delist(v) for k, v in node.items()}
    return node

import re

def main():
    cfg_path, key = sys.argv[1], sys.argv[2]
    cfg = _load(cfg_path)
    if key == "perturbseq.conditions":
        print("\n".join(map(str, cfg["perturbseq"]["conditions"]))); return
    if key == "lof.traits.names":
        print("\n".join(t["name"] for t in cfg["lof"]["traits"])); return
    node = cfg
    for part in key.split("."):
        node = node[part]
    if isinstance(node, list):
        print("\n".join(map(str, node)))
    else:
        print(node)

if __name__ == "__main__":
    main()
