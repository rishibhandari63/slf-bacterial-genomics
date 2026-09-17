#!/usr/bin/env python3
"""
One-time patch for a known pyani 0.2.x rough edge: ANIm's parse_delta()
divides by the total aligned base count with no zero-check. When two
genomes are too distantly related for NUCmer to find ANY alignable
sequence (common when your genome set spans multiple phyla, like this
project's), that sum is legitimately zero, and pyani crashes instead of
just reporting "no meaningful ANI for this pair" and continuing.

This only touches YOUR conda env's private copy of pyani - not a shared
system file, not the pyani package on PyPI. Safe to run repeatedly:
it's a no-op if already patched, and aborts without changing anything
if the expected source line isn't found verbatim (e.g. a different
pyani version than the one this was checked against).

Usage:
    python3 10_patch_pyani_zerodiv.py /project/YOUR_PROJECT_HERE/envs/pyani_env
"""
import sys
from pathlib import Path


def find_anim_py(env_prefix: Path) -> Path:
    matches = list(env_prefix.glob("lib/python*/site-packages/pyani/anim.py"))
    if not matches:
        sys.exit(f"Couldn't find pyani/anim.py under {env_prefix} - is the path right?")
    return matches[0]


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 10_patch_pyani_zerodiv.py <path to pyani_env>")
    env_prefix = Path(sys.argv[1])
    target = find_anim_py(env_prefix)

    old = "avrg_ID = sum(weighted_identical_bases) / sum(aligned_bases)"
    new = (
        "avrg_ID = (\n"
        "        sum(weighted_identical_bases) / sum(aligned_bases)\n"
        "        if sum(aligned_bases)\n"
        "        else 0.0\n"
        "    )"
    )

    content = target.read_text()
    if new.split("\n")[0] in content and "if sum(aligned_bases)" in content:
        print(f"Already patched: {target}")
        return
    if old not in content:
        sys.exit(
            f"Expected line not found verbatim in {target} - this pyani version's "
            "source may differ from what this patch was checked against. Not "
            "changing anything; open the file and look for the "
            "'avrg_ID = sum(...) / sum(...)' line in parse_delta() by hand."
        )

    content = content.replace(old, new, 1)
    target.write_text(content)
    print(f"Patched: {target}")
    print("Pairs with zero alignable bases will now report ANI = 0.0 instead of crashing the whole run.")


if __name__ == "__main__":
    main()
