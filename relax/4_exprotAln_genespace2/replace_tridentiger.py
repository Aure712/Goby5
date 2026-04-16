#!/usr/bin/env python3
import sys
from pathlib import Path

OLD = "(Tridentiger_barbatus{T},Tridentiger_radiatus{T})"
NEW = "(Tridentiger_barbatus{R},Tridentiger_radiatus{R})"

def main(dir_path="Relax_T_only"):
    base = Path(dir_path)
    if not base.is_dir():
        print(f"Error: directory not found: {base}")
        return 1

    files_changed = 0
    total_replacements = 0

    for f in sorted(base.glob("Group*.txt")):
        try:
            text = f.read_text(encoding="utf-8")
        except Exception as e:
            print(f"Skip {f}: cannot read ({e})")
            continue

        count = text.count(OLD)
        if count > 0:
            new_text = text.replace(OLD, NEW)
            try:
                f.write_text(new_text, encoding="utf-8")
            except Exception as e:
                print(f"Error writing {f}: {e}")
                continue
            files_changed += 1
            total_replacements += count
            print(f"{f}: {count} replacement(s)")

    print(f"Done. Files changed: {files_changed}, total replacements: {total_replacements}")
    return 0

if __name__ == "__main__":
    # Usage: python3 replace_tridentiger.py [Relax_T_only]
    path = sys.argv[1] if len(sys.argv) > 1 else "Relax_T_only"
    raise SystemExit(main(path))
