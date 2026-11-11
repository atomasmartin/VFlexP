#!/usr/bin/env python3
"""Git textconv driver for Simulink .slx archives.

Configure once per repository clone:
    git config diff.slx.textconv "python tools/slx_textconv.py"

Together with `.gitattributes` (`*.slx diff=slx`), `git diff` on any `.slx`
prints structural block, parameter, mask, and netlist changes while ignoring
purely visual layout noise (Position, ZOrder, ZoomFactor, Points, SIDs).
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

VISUAL_PARAMS = {
    "Position",
    "ZOrder",
    "Points",
    "Location",
    "Open",
    "ZoomFactor",
    "ScrollbarOffset",
    "ScreenColor",
    "PaperPosition",
    "PaperSize",
    "TiledPageScale",
    "ModelBrowserVisibility",
    "ModelBrowserWidth",
}


def dump_system_xml(sys_name: str, xml_bytes: bytes, out: list[str]) -> None:
    root = ET.fromstring(xml_bytes)
    blocks = []
    for blk in root.findall(".//Block"):
        btype = blk.attrib.get("BlockType", "")
        bname = " ".join(blk.attrib.get("Name", "").split())
        params = []
        for p in blk.findall("P"):
            pname = p.attrib.get("Name", "")
            if pname in VISUAL_PARAMS or pname.startswith("Font"):
                continue
            val = " ".join((p.text or "").split())
            if val:
                params.append(f"    .{pname} = {val}")
        mask = blk.find("Mask")
        if mask is not None:
            for mp in mask.findall(".//MaskParameter"):
                mname = mp.attrib.get("Name", "")
                mtype = mp.attrib.get("Type", "")
                val_el = mp.find("Value")
                mval = " ".join((val_el.text if val_el is not None and val_el.text else "").split())
                params.append(f"    .mask.{mname} ({mtype}) = {mval}")
            init_el = mask.find("Initialization")
            if init_el is not None and init_el.text and init_el.text.strip():
                init_txt = "; ".join(line.strip() for line in init_el.text.strip().splitlines() if line.strip())
                params.append(f"    .MaskInitialization = {init_txt}")
        params.sort()
        blocks.append((bname, btype, params))

    blocks.sort(key=lambda x: (x[0], x[1]))
    out.append(f"=== {sys_name} ({len(blocks)} blocks) ===")
    for bname, btype, params in blocks:
        out.append(f"  [{bname}] ({btype})")
        out.extend(params)


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python tools/slx_textconv.py <file.slx>", file=sys.stderr)
        return 1
    slx_path = Path(sys.argv[1])
    out: list[str] = [f"# SLX Textconv: {slx_path.name}"]
    with zipfile.ZipFile(slx_path, "r") as zf:
        sys_files = sorted(
            n for n in zf.namelist()
            if n.startswith("simulink/systems/system_") and n.endswith(".xml")
        )
        for name in sys_files:
            dump_system_xml(name, zf.read(name), out)
    sys.stdout.buffer.write(("\n".join(out) + "\n").encode("utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
