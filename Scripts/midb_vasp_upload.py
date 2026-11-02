#!/usr/bin/env python3
"""Upload VASP calculations to MIDb (midb.cloud) through the public API.

    python midb_vasp_upload.py check   CALC_DIR [CALC_DIR ...]
    python midb_vasp_upload.py package --dataset NAME [--info info.json] -o pacote.zip CALC_DIR [...]
    python midb_vasp_upload.py upload  --url https://www.midb.cloud --metadata metadata.json pacote.zip
    python midb_vasp_upload.py status  --url https://www.midb.cloud --dataset NAME

Each CALC_DIR is one VASP run and becomes one MIDb entry. MIDb recognizes it by
vasprun.xml. Token: env MIDB_TOKEN (MIDb dashboard > API token). Requires `requests`.
"""
import argparse
import gzip
import json
import os
import re
import shutil
import sys
import zipfile

# Files MIDb keeps for a VASP entry (exact names; anything else is ignored by MIDb).
KEEP = ["vasprun.xml", "INCAR", "KPOINTS", "POSCAR", "CONTCAR", "OUTCAR", "OSZICAR",
        "DOSCAR", "EIGENVAL", "PROCAR", "IBZKPT", "XDATCAR"]
# Large files: only with --include-large. POTCAR is never sent (VASP license).
LARGE = ["CHGCAR", "CHG", "LOCPOT", "ELFCAR", "PARCHG"]
NEVER = ["POTCAR", "WAVECAR", "WAVEDER", "vaspwave.h5", "TMPCAR"]


def _source(calc, name):
    """Path of NAME in CALC, or of NAME.gz (decompressed while packaging)."""
    for cand in (name, name + ".gz"):
        p = os.path.join(calc, cand)
        if os.path.isfile(p):
            return p
    return None


def _read_text(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        return f.read()


def check_calc(calc, include_large=False):
    """Return (problems, kept_files, notes) for one calculation folder."""
    problems, notes = [], []
    if not os.path.isdir(calc):
        return [f"{calc}: not a directory"], [], notes
    vr = _source(calc, "vasprun.xml")
    if vr is None:
        problems.append("vasprun.xml missing (MIDb only recognizes VASP folders by it)")
    else:
        text = _read_text(vr)
        if "</modeling>" not in text[-2000:]:
            problems.append("vasprun.xml incomplete (run did not finish: no </modeling>)")
    names = KEEP + (LARGE if include_large else [])
    kept = [n for n in names if _source(calc, n)]
    kp = _source(calc, "KPOINTS")
    kp_lines = _read_text(kp).splitlines() if kp else []
    if len(kp_lines) > 2 and kp_lines[2].strip()[:1] in ("l", "L"):
        notes.append("line-mode KPOINTS: band structure + DOS plot on MIDb")
    elif vr and kp:
        notes.append("not line-mode: DOS shown, band labels guessed by MIDb")
    elif vr:
        notes.append("no KPOINTS: MIDb shows no band/DOS plot (structure only)")
    for n in NEVER:
        if os.path.exists(os.path.join(calc, n)):
            notes.append(f"{n} will NOT be sent")
    return problems, kept, notes


def _parse_calc_spec(calc_spec):
    """Split 'CALC_DIR:ENTRY_NAME' or default to basename of CALC_DIR."""
    if ":" in calc_spec:
        p, n = calc_spec.split(":", 1)
        if os.path.exists(p):
            return p, re.sub(r"[^A-Za-z0-9._-]+", "_", n)
    return calc_spec, re.sub(r"[^A-Za-z0-9._-]+", "_", os.path.basename(os.path.normpath(calc_spec)))


def cmd_check(a):
    bad = 0
    for calc_spec in a.calcs:
        calc, name = _parse_calc_spec(calc_spec)
        problems, kept, notes = check_calc(calc, a.include_large)
        status = "ERRO" if problems else "ok"
        label = f"{name} ({calc})" if name != os.path.basename(os.path.normpath(calc)) else calc
        print(f"[{status}] {label}\n    envia: {', '.join(kept) or '-'}")
        for p in problems:
            print(f"    problema: {p}")
        for n in notes:
            print(f"    nota: {n}")
        bad += bool(problems)
    return 1 if bad else 0


def cmd_package(a):
    folders = {}
    for calc_spec in a.calcs:
        calc, name = _parse_calc_spec(calc_spec)
        if name in folders:
            sys.exit(f"duplicate folder name {name!r} ({calc} and {folders[name]}); rename one")
        folders[name] = calc
    info = json.load(open(a.info)) if a.info else {}
    top = re.sub(r"[^A-Za-z0-9._-]+", "_", a.dataset)
    bad = 0
    with zipfile.ZipFile(a.output, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name, calc in folders.items():
            problems, kept, _ = check_calc(calc, a.include_large)
            if problems:
                print(f"skip {calc}: {'; '.join(problems)}")
                bad += 1
                continue
            for fname in kept:
                src = _source(calc, fname)
                arc = f"{top}/{name}/{fname}"
                if src.endswith(".gz"):
                    with gzip.open(src, "rb") as fin, z.open(arc, "w") as fout:
                        shutil.copyfileobj(fin, fout)
                else:
                    z.write(src, arc)
            sup = {**info.get("*", {}), **info.get(name, {})}
            sup.setdefault("calculation", name)
            z.writestr(f"{top}/{name}/entry.json",
                       json.dumps({"supplemental_information": sup, "images": [], "files": []}, indent=1))
            print(f"ok   {calc} -> {top}/{name} ({len(kept)} files)")
    size = os.path.getsize(a.output) / 2**20
    print(f"{a.output}: {size:.1f} MiB, {len(folders) - bad} entries")
    return 1 if bad else 0


def cmd_upload(a):
    import requests
    token = os.environ.get("MIDB_TOKEN") or sys.exit("set MIDB_TOKEN")
    meta = json.load(open(a.metadata))
    for key in ("authors", "contact_email", "keywords", "dataset", "description"):
        if not meta.get(key):
            sys.exit(f"metadata.json: missing {key}")
    with open(a.archive, "rb") as f:
        # MIDb reads the metadata as a FILE part named 'form' and decodes it as ASCII
        files = {"form": ("form.json", json.dumps(meta), "application/json"),
                 "file": (os.path.basename(a.archive), f, "application/zip")}
        r = requests.post(a.url.rstrip("/") + "/public_api/upload/dataset",
                          headers={"Authorization": f"Bearer {token}"}, files=files, timeout=3600)
    print(r.status_code, r.text.strip())
    return 0 if r.ok else 1


def cmd_status(a):
    import requests
    token = os.environ.get("MIDB_TOKEN") or sys.exit("set MIDB_TOKEN")
    q = {"filter": {"entry_type": "VASP", "dataset": a.dataset},
         "fields": ["formula", "space_group", "is_public", "upload_date"]}
    r = requests.post(a.url.rstrip("/") + "/public_api/query",
                      headers={"Authorization": f"Bearer {token}"}, json=q, timeout=120)
    if not r.ok:
        print(r.status_code, r.text)
        return 1
    rows = r.json()
    for e in rows:
        print(f"{e['_id']}  {e.get('formula', '?'):12s} {e.get('space_group', ''):10s} public={e.get('is_public')}")
    print(f"{len(rows)} VASP entries in dataset {a.dataset!r}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check"); c.add_argument("calcs", nargs="+"); c.add_argument("--include-large", action="store_true")
    p = sub.add_parser("package"); p.add_argument("calcs", nargs="+"); p.add_argument("--dataset", required=True)
    p.add_argument("--info", help="JSON {folder_name or '*': {supplemental key: value}}")
    p.add_argument("-o", "--output", required=True); p.add_argument("--include-large", action="store_true")
    u = sub.add_parser("upload"); u.add_argument("archive"); u.add_argument("--url", required=True)
    u.add_argument("--metadata", required=True)
    s = sub.add_parser("status"); s.add_argument("--url", required=True); s.add_argument("--dataset", required=True)
    a = ap.parse_args()
    sys.exit({"check": cmd_check, "package": cmd_package, "upload": cmd_upload, "status": cmd_status}[a.cmd](a))


if __name__ == "__main__":
    main()
