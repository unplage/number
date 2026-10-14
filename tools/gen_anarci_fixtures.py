#!/usr/bin/env python3
"""Regenerate tests/fixtures/anarci_all_schemes.json from real ANARCI.

Requires the ANARCI Python package and the HMMER binaries (hmmscan/hmmbuild):

    pip install anarci          # plus hmmer on PATH

Run from the repository root:

    python3 tools/gen_anarci_fixtures.py

The output is the golden reference used by tests/schemes.test.mjs. It numbers a
set of VH and VL sequences with every scheme ANARCI supports.
"""
import json
import os
from anarci import anarci

VH = {
    "trast_vh": "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
    "tm_vh": "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS",
    "long_vh": "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARRRRRRRRRRRSSSSSSSSSYYYYYYYYYWGQGTLVTVSS",
    "vhh": "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAAGDYWGQGTLVTVSS",
    "vh2": "QVQLQQSGAELARPGASVKLSCKASGYTFTSYWMQWVKQRPGQGLEWIGAIYPGDGDTRYTQKFKGKATLTADKSSSTAYMQLSSLASEDSAVYYCARGGTGYFDYWGQGTTLTVSS",
    "vh_short": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYGISWVRQAPGQGLEWMGWISAYNGNTNYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARDLGYSYYYYGMDVWGQGTTVTVSS",
}
VL = {
    "trast_vl": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
    "vk2": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK",
    "vl_lambda": "QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL",
    "vl_long_cdr3": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQWYSSPPPTFGQGTKVEIK",
}
SCHEMES = ["kabat", "chothia", "martin", "aho", "wolfguy"]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests", "fixtures", "anarci_all_schemes.json")


def main():
    out = {}
    for seqs in (VH, VL):
        for name, seq in seqs.items():
            for scheme in SCHEMES:
                res = anarci([(name, seq)], scheme=scheme, output=False)
                if not res[0] or not res[0][0]:
                    out[f"{name}::{scheme}"] = {"error": "not numbered"}
                    continue
                dom = res[0][0][0][0]
                details = res[1][0][0]
                if isinstance(details, list):
                    details = details[0]
                out[f"{name}::{scheme}"] = {
                    "chain": details["chain_type"],
                    "numbering": [[str(n) + ic, aa] for (n, ic), aa in dom],
                }
    with open(OUT, "w") as f:
        json.dump(out, f)
    print(f"wrote {OUT} ({len(out)} entries)")


if __name__ == "__main__":
    main()
