#!/usr/bin/env python3
"""Regenerate tests/fixtures/vl_reference.json and vl_ladder.json from the
original unplage/humanization Kabat light-chain engine.

Run from the repository root:

    python3 tools/gen_vl_fixtures.py

Reads docs/reference/humanization_numbering.py (MIT) and stubs the optional
ANARCI/AbRSA adapters, then emits the golden references used by
tests/kabat-light.test.mjs.
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "docs", "reference", "humanization_numbering.py")
FIX = os.path.join(ROOT, "tests", "fixtures")

VL = {
    "vl1": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
    "vl2": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK",
    "vl_insert": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQWYSSPPTFGQGTKVEIK",
    "lambda": "QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL",
}


def load_number_light():
    src = open(REF).read()
    src = src.replace(
        "from .abrsa import is_abrsa_available, number_with_abrsa",
        "is_abrsa_available = lambda: False; number_with_abrsa = lambda *a, **k: None",
    )
    src = src.replace(
        "from .anarci_adapter import is_anarci_available, number_with_anarci_imgt",
        "is_anarci_available = lambda: False; number_with_anarci_imgt = lambda *a, **k: None",
    )
    ns = {}
    exec(compile(src, "numbering.py", "exec"), ns)
    return ns["number_light"]


def main():
    number_light = load_number_light()
    out = {
        name: [[r.pos, r.aa, r.region] for r in number_light(seq).residues]
        for name, seq in VL.items()
    }
    json.dump(out, open(os.path.join(FIX, "vl_reference.json"), "w"))

    base = VL["vl1"]
    ladder = {}
    for extra in range(0, 6):
        cdr3 = "QQHYTTP" + "A" * extra + "PT"
        seq = base.replace("QQHYTTPPT", cdr3)
        ladder["cdr3_%d" % extra] = [[r.pos, r.aa] for r in number_light(seq).residues]
    json.dump(ladder, open(os.path.join(FIX, "vl_ladder.json"), "w"))
    print("wrote vl_reference.json and vl_ladder.json")


if __name__ == "__main__":
    main()
