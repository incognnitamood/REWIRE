"""Download and unpack STRING gene name mapping file.

Run with:
    python download_gene_map.py
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import requests
from tqdm import tqdm


URLS = [
    "https://stringdb-static.org/download/9606.protein.info.v12.0.txt.gz",
    "https://stringdb-static.org/download/9606.protein.info.v11.5.txt.gz",
    "https://stringdb-static.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
    "https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz",
]


def _download_with_progress(url: str, dest_path: Path) -> bool:
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            with open(dest_path, "wb") as f, tqdm(
                total=total,
                unit="B",
                unit_scale=True,
                desc="Downloading",
            ) as bar:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if not chunk:
                        continue
                    f.write(chunk)
                    bar.update(len(chunk))
        return True
    except requests.exceptions.SSLError:
        print("WARNING: SSL verification failed. Retrying with verify=False.")
        with requests.get(url, stream=True, verify=False) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            with open(dest_path, "wb") as f, tqdm(
                total=total,
                unit="B",
                unit_scale=True,
                desc="Downloading (no SSL)",
            ) as bar:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if not chunk:
                        continue
                    f.write(chunk)
                    bar.update(len(chunk))
        return True
    except Exception as exc:
        print(f"ERROR: Download failed: {exc}")
        return False


def _unzip_gz(gz_path: Path, out_path: Path) -> bool:
    try:
        with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
        return True
    except Exception as exc:
        print(f"ERROR: Failed to unzip {gz_path}: {exc}")
        return False


def _convert_if_needed(src_path: Path, out_path: Path) -> bool:
    try:
        with open(src_path, "r", encoding="utf-8") as src:
            header = src.readline().rstrip("\n")
            cols = header.split("\t")

        if "preferred_name" in cols and "#string_protein_id" in cols:
            # Already in desired format
            os.replace(src_path, out_path)
            return True

        if "preferred_name" in cols and "string_protein_id" in cols:
            id_idx = cols.index("string_protein_id")
            name_idx = cols.index("preferred_name")
        elif "preferred_name" in cols and "protein_external_id" in cols:
            id_idx = cols.index("protein_external_id")
            name_idx = cols.index("preferred_name")
        else:
            print("ERROR: Unrecognized header format:")
            print(header)
            return False

        with open(src_path, "r", encoding="utf-8") as src, open(out_path, "w", encoding="utf-8") as dst:
            src.readline()
            dst.write("#string_protein_id\tpreferred_name\n")
            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(id_idx, name_idx):
                    continue
                dst.write(f"{parts[id_idx]}\t{parts[name_idx]}\n")
        return True
    except Exception as exc:
        print(f"ERROR: Failed to convert mapping file: {exc}")
        return False


def _print_preview_and_count(path: Path) -> None:
    line_count = 0
    preview = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line_count += 1
            if len(preview) < 5:
                preview.append(line.rstrip("\n"))

    print("\nFirst 5 lines:")
    for line in preview:
        print(line)

    print(f"\nTotal rows: {line_count}")


def main() -> int:
    repo_root = Path(__file__).parent
    out_dir = repo_root / "REWIRE" / "data" / "processed"
    out_dir.mkdir(parents=True, exist_ok=True)

    gz_path = out_dir / "human.name_2_string.tsv.gz"
    raw_path = out_dir / "human.name_2_string.raw.tsv"
    out_path = out_dir / "human.name_2_string.tsv"

    downloaded = False
    for url in URLS:
        print(f"Trying: {url}")
        if _download_with_progress(url, gz_path):
            downloaded = True
            break

    if not downloaded:
        return 1

    if not _unzip_gz(gz_path, raw_path):
        return 1

    if not _convert_if_needed(raw_path, out_path):
        return 1

    _print_preview_and_count(out_path)

    print(f"\nDOWNLOAD COMPLETE — saved to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
