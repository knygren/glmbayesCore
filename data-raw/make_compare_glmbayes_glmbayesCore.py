#!/usr/bin/env python3
"""Compare glmbayes vs glmbayesCore (R/, src/, inst/cl/) and write assessment README.

Run from package root or via the R wrapper:
  source("data-raw/make_compare_glmbayes_glmbayesCore.R")
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from collections import OrderedDict
from datetime import date
from pathlib import Path

CORE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GB = CORE_ROOT.parent / "glmbayes"


def git_sha(repo: Path) -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip()
    except Exception:
        return "unknown"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def normalize_text(text: str) -> str:
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = [ln.rstrip() for ln in text.split("\n")]
    return "\n".join(lines)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def file_hashes(path: Path) -> tuple[str, str]:
    raw = path.read_bytes()
    byte_h = sha256_bytes(raw)
    norm_h = sha256_bytes(normalize_text(raw.decode("utf-8", errors="replace")).encode("utf-8"))
    return byte_h, norm_h


def list_r_files(root: Path) -> list[str]:
    d = root / "R"
    return sorted(p.name for p in d.glob("*.R"))


def list_src_toplevel(root: Path) -> list[str]:
    d = root / "src"
    out = []
    for p in d.iterdir():
        if not p.is_file():
            continue
        if p.suffix.lower() not in {".cpp", ".h", ".c"}:
            continue
        out.append(p.name)
    return sorted(out)


def list_nmath(root: Path) -> list[str]:
    d = root / "src" / "nmath"
    if not d.is_dir():
        return []
    rels = []
    for p in d.rglob("*"):
        if p.is_file() and p.suffix.lower() in {".cpp", ".h", ".c"}:
            rels.append(p.relative_to(root / "src").as_posix())
    return sorted(rels)


def list_cl_relevant(root: Path) -> list[str]:
    d = root / "inst" / "cl"
    if not d.is_dir():
        return []
    rels = []
    for p in d.rglob("*"):
        if not p.is_file():
            continue
        suf = p.suffix.lower()
        name = p.name.lower()
        if suf == ".cl" or suf in {".tsv", ".rds"} or name == "readme.md":
            rels.append(p.relative_to(d).as_posix())
    return sorted(rels)


# ---------------------------------------------------------------------------
# Function parsers
# ---------------------------------------------------------------------------

def extract_r_functions(text: str) -> dict[str, str]:
    """Top-level name <- function(...) { ... } bodies (normalized)."""
    text = normalize_text(text)
    lines = text.split("\n")
    # Match start lines
    start_re = re.compile(
        r"^([A-Za-z.][A-Za-z0-9._]*)\s*(<-|=)\s*function\s*\("
    )
    funcs: dict[str, str] = {}
    i = 0
    while i < len(lines):
        m = start_re.match(lines[i])
        if not m:
            i += 1
            continue
        name = m.group(1)
        # Collect from 'function' to matching close brace after args
        # Find start of function token
        chunk_lines = [lines[i]]
        # Track paren depth for args, then brace depth for body
        rest = lines[i]
        fn_pos = rest.find("function")
        after = rest[fn_pos + len("function") :]
        paren = 0
        brace = 0
        in_body = False
        started_paren = False
        j = i
        # scan character by character across lines
        buf = []
        done = False
        while j < len(lines) and not done:
            line = lines[j] if j > i else after
            # For first line, we already sliced after "function"
            if j == i:
                scan = after
                prefix_for_buf = rest  # full first line goes into body text later
            else:
                scan = line
            k = 0
            while k < len(scan):
                ch = scan[k]
                if not in_body:
                    if ch == "(":
                        paren += 1
                        started_paren = True
                    elif ch == ")":
                        paren -= 1
                    elif ch == "{" and started_paren and paren == 0:
                        in_body = True
                        brace = 1
                else:
                    if ch == "{":
                        brace += 1
                    elif ch == "}":
                        brace -= 1
                        if brace == 0:
                            done = True
                            k += 1
                            break
                k += 1
            if j == i:
                buf.append(rest)
            else:
                buf.append(line)
            j += 1
            if done:
                break
        body = normalize_text("\n".join(buf))
        # Strip package-rename noise for "identical after rename" optional note
        funcs[name] = body
        i = j
    return funcs


def _strip_cpp_noise(s: str) -> str:
    s = normalize_text(s)
    # Collapse package rename differences for comparison helper
    s2 = s
    s2 = s2.replace("glmbayesCore", "PKG")
    s2 = s2.replace("glmbayes", "PKG")
    s2 = s2.replace("_glmbayesCore_", "_PKG_")
    s2 = s2.replace("_glmbayes_", "_PKG_")
    s2 = s2.replace('"glmbayesCore"', '"PKG"')
    s2 = s2.replace('"glmbayes"', '"PKG"')
    return s2


def extract_cpp_functions(text: str, is_header: bool) -> dict[str, str]:
    """Heuristic extraction of C++ function definitions / declarations."""
    text = normalize_text(text)
    funcs: dict[str, str] = {}

    # Rcpp::export blocks: capture following function signature name
    export_re = re.compile(
        r"//\s*\[\[Rcpp::export\]\]\s*\n"
        r"(?:[^\n]*\n)*?"
        r"(?:[\w:<>,\s\*&]+?)\s+(\w+)\s*\([^;]*?\)\s*\{",
        re.M,
    )
    for m in export_re.finditer(text):
        name = m.group(1)
        start = m.start()
        # find body end from opening brace in match
        brace_start = text.find("{", m.end() - 1)
        if brace_start < 0:
            continue
        body = _extract_brace_block(text, brace_start)
        if body is not None:
            funcs[name] = normalize_text(text[start : brace_start + len(body)])

    # Namespace-qualified definitions: ReturnType name( ... ) {
    # Also bare definitions at column 0-ish
    def_re = re.compile(
        r"(?m)^(?:template\s*<[^>]+>\s*)?"
        r"(?:inline\s+|static\s+|virtual\s+|constexpr\s+)*"
        r"(?:[\w:<>\*&\s]+?)\s+"
        r"([A-Za-z_]\w*)\s*\(([^;]*)\)\s*(?:const\s*)?\{",
    )
    for m in def_re.finditer(text):
        name = m.group(1)
        # skip control-like
        if name in {"if", "for", "while", "switch", "catch", "else"}:
            continue
        brace_start = text.find("{", m.end() - 1)
        if brace_start < 0:
            continue
        body = _extract_brace_block(text, brace_start)
        if body is None:
            continue
        key = name
        # disambiguate overloads lightly
        if key in funcs and funcs[key] != normalize_text(text[m.start() : brace_start + len(body)]):
            args = re.sub(r"\s+", " ", m.group(2).strip())[:60]
            key = f"{name}({args})"
        funcs[key] = normalize_text(text[m.start() : brace_start + len(body)])

    if is_header:
        # Declarations ending with );
        decl_re = re.compile(
            r"(?m)^(?:\s*)(?:inline\s+|static\s+)*"
            r"(?:[\w:<>\*&\s]+?)\s+"
            r"([A-Za-z_]\w*)\s*\(([^;]*)\)\s*;"
        )
        for m in decl_re.finditer(text):
            name = m.group(1)
            if name in {"if", "for", "while", "switch", "return"}:
                continue
            key = f"decl:{name}"
            if key not in funcs:
                funcs[key] = normalize_text(m.group(0))

    return funcs


def _extract_brace_block(text: str, brace_start: int) -> str | None:
    if brace_start >= len(text) or text[brace_start] != "{":
        return None
    depth = 0
    i = brace_start
    while i < len(text):
        ch = text[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return text[brace_start : i + 1]
        i += 1
    return None


def extract_cl_symbols(text: str) -> dict[str, str]:
    text = normalize_text(text)
    funcs: dict[str, str] = {}
    # @provides tags
    for m in re.finditer(r"@provides\s+(\S+)", text):
        funcs[f"provides:{m.group(1)}"] = m.group(0)
    # __kernel void name(
    for m in re.finditer(r"__kernel\s+void\s+(\w+)\s*\(", text):
        name = m.group(1)
        brace_start = text.find("{", m.end())
        if brace_start < 0:
            funcs[name] = m.group(0)
            continue
        body = _extract_brace_block(text, brace_start)
        if body:
            funcs[name] = normalize_text(text[m.start() : brace_start + len(body)])
        else:
            funcs[name] = m.group(0)
    return funcs


def compare_func_maps(
    a: dict[str, str], b: dict[str, str]
) -> dict[str, list[str]]:
    sa, sb = set(a), set(b)
    both = sorted(sa & sb)
    identical = []
    different = []
    rename_only = []
    for name in both:
        if a[name] == b[name]:
            identical.append(name)
        elif _strip_cpp_noise(a[name]) == _strip_cpp_noise(b[name]):
            rename_only.append(name)
        else:
            different.append(name)
    return {
        "identical": identical,
        "rename_only": rename_only,
        "different": different,
        "core_only": sorted(sa - sb),
        "glmbayes_only": sorted(sb - sa),
    }


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_shared(
    names: list[str],
    path_core: Path,
    path_gb: Path,
    rel_fn,
) -> dict:
    identical_byte = []
    identical_norm_only = []
    differing = []
    details = {}
    for name in names:
        pc = path_core / rel_fn(name)
        pg = path_gb / rel_fn(name)
        if not pc.exists() or not pg.exists():
            continue
        hb_c, hn_c = file_hashes(pc)
        hb_g, hn_g = file_hashes(pg)
        details[name] = {
            "core_byte": hb_c,
            "gb_byte": hb_g,
            "core_norm": hn_c,
            "gb_norm": hn_g,
        }
        if hb_c == hb_g:
            identical_byte.append(name)
        elif hn_c == hn_g:
            identical_norm_only.append(name)
        else:
            differing.append(name)
    return {
        "identical_byte": identical_byte,
        "identical_norm_only": identical_norm_only,
        "differing": differing,
        "details": details,
    }


def md_list(items: list[str], empty: str = "*(none)*") -> str:
    if not items:
        return empty
    return "\n".join(f"- `{x}`" for x in items)


def md_table_counts(rows: list[tuple[str, int]]) -> str:
    lines = ["| Category | Count |", "|----------|------:|"]
    for k, v in rows:
        lines.append(f"| {k} | {v} |")
    return "\n".join(lines)


def main() -> None:
    import os

    gb_root = Path(os.environ.get("GLMBAYES_COMPARE_ROOT", str(DEFAULT_GB)))
    if not gb_root.is_dir():
        raise SystemExit(f"glmbayes root not found: {gb_root}")

    core_sha = git_sha(CORE_ROOT)
    gb_sha = git_sha(gb_root)
    today = date.today().isoformat()

    # --- inventories ---
    r_core = list_r_files(CORE_ROOT)
    r_gb = list_r_files(gb_root)
    r_both = sorted(set(r_core) & set(r_gb))
    r_core_only = sorted(set(r_core) - set(r_gb))
    r_gb_only = sorted(set(r_gb) - set(r_core))
    r_cls = classify_shared(
        r_both, CORE_ROOT, gb_root, lambda n: Path("R") / n
    )

    src_core = [n for n in list_src_toplevel(CORE_ROOT) if n != "RcppExports.cpp"]
    src_gb = [n for n in list_src_toplevel(gb_root) if n != "RcppExports.cpp"]
    src_both = sorted(set(src_core) & set(src_gb))
    src_core_only = sorted(set(src_core) - set(src_gb))
    src_gb_only = sorted(set(src_gb) - set(src_core))
    src_cls = classify_shared(
        src_both, CORE_ROOT, gb_root, lambda n: Path("src") / n
    )

    nmath_core = list_nmath(CORE_ROOT)
    nmath_gb = list_nmath(gb_root)

    cl_core = list_cl_relevant(CORE_ROOT)
    cl_gb = list_cl_relevant(gb_root)
    cl_both = sorted(set(cl_core) & set(cl_gb))
    cl_core_only = sorted(set(cl_core) - set(cl_gb))
    cl_gb_only = sorted(set(cl_gb) - set(cl_core))
    cl_cls = classify_shared(
        cl_both, CORE_ROOT, gb_root, lambda n: Path("inst") / "cl" / n
    )

    # --- function-level for differing files ---
    fn_reports: OrderedDict[str, dict] = OrderedDict()

    for name in r_cls["differing"]:
        if name in ("RcppExports.R",):
            fn_reports[f"R/{name}"] = {
                "kind": "R",
                "note": "Generated by Rcpp::compileAttributes(); skip detailed body compare.",
                "result": None,
            }
            continue
        a = extract_r_functions(read_text(CORE_ROOT / "R" / name))
        b = extract_r_functions(read_text(gb_root / "R" / name))
        fn_reports[f"R/{name}"] = {
            "kind": "R",
            "result": compare_func_maps(a, b),
            "counts": {"core": len(a), "glmbayes": len(b)},
        }

    for name in src_cls["differing"]:
        if name == "RcppExports.cpp":
            continue
        pc = CORE_ROOT / "src" / name
        pg = gb_root / "src" / name
        is_h = name.endswith(".h")
        a = extract_cpp_functions(read_text(pc), is_header=is_h)
        b = extract_cpp_functions(read_text(pg), is_header=is_h)
        fn_reports[f"src/{name}"] = {
            "kind": "cpp",
            "result": compare_func_maps(a, b),
            "counts": {"core": len(a), "glmbayes": len(b)},
        }

    for name in cl_cls["differing"]:
        if not name.endswith(".cl"):
            fn_reports[f"inst/cl/{name}"] = {
                "kind": "cl",
                "note": "Non-.cl differing file; hash-only.",
                "result": None,
            }
            continue
        a = extract_cl_symbols(read_text(CORE_ROOT / "inst" / "cl" / name))
        b = extract_cl_symbols(read_text(gb_root / "inst" / "cl" / name))
        fn_reports[f"inst/cl/{name}"] = {
            "kind": "cl",
            "result": compare_func_maps(a, b),
            "counts": {"core": len(a), "glmbayes": len(b)},
        }

    # Kernel entry names from identical f2_f3 files (reference)
    kernel_entries = []
    for rel in sorted(x for x in cl_both if x.startswith("src/f2_f3_") and x.endswith(".cl")):
        syms = extract_cl_symbols(read_text(CORE_ROOT / "inst" / "cl" / rel))
        kernels = [k for k in syms if not k.startswith("provides:")]
        kernel_entries.append((rel, kernels))

    # DESCRIPTION OpenCL bits
    def desc_opencl(root: Path) -> dict:
        t = read_text(root / "DESCRIPTION")
        linking = "opencltools" in re.search(
            r"LinkingTo:\s*((?:.|\n)*?)(?:\n\S|\Z)", t + "\n"
        ).group(1) if re.search(r"LinkingTo:", t) else False
        # simpler
        m = re.search(r"LinkingTo:\s*((?:[^\n]|\n(?![A-Za-z]))*)", t)
        lt = m.group(1) if m else ""
        m2 = re.search(r"Imports:\s*((?:[^\n]|\n(?![A-Za-z]))*)", t)
        im = m2.group(1) if m2 else ""
        return {
            "LinkingTo_opencltools": "opencltools" in lt,
            "Imports_opencltools": "opencltools" in im,
            "Imports_nmathopencl": "nmathopencl" in im,
        }

    ocl_core = desc_opencl(CORE_ROOT)
    ocl_gb = desc_opencl(gb_root)

    # Build markdown
    lines: list[str] = []
    L = lines.append

    L("# glmbayes vs glmbayesCore — comparison assessment")
    L("")
    L(f"**Generated:** {today}  ")
    L(f"**glmbayesCore:** `{CORE_ROOT}` @ `{core_sha}`  ")
    L(f"**glmbayes:** `{gb_root}` @ `{gb_sha}`  ")
    L("")
    L("Regenerate: `source(\"data-raw/make_compare_glmbayes_glmbayesCore.R\")`")
    L("")
    L("---")
    L("")
    L("## 1. Purpose and method")
    L("")
    L("This document inventories and compares the **iid GLM/LM + OpenCL** surfaces of")
    L("**glmbayes** and **glmbayesCore** before Batch 7b (OpenCL loader alignment).")
    L("Mixed-model code temporarily staged in **lmebayesCore** is out of scope.")
    L("")
    L("### Identity rules")
    L("")
    L("- **Byte-identical:** SHA256 of raw file bytes.")
    L("- **Normalized-identical:** SHA256 after UTF-8 decode, CRLF→LF, strip trailing")
    L("  whitespace per line (catches Windows line-ending-only drift).")
    L("- **Differing:** neither byte nor normalized hash matches.")
    L("- **Package-rename noise:** when comparing function bodies, a secondary check")
    L("  replaces `glmbayes`/`glmbayesCore` and `_glmbayes_`/`_glmbayesCore_` with a")
    L("  placeholder; matches after that are labeled **rename-only**.")
    L("")
    L("### Scope")
    L("")
    L("| Tree | Included |")
    L("|------|----------|")
    L("| `R/` | all `*.R` |")
    L("| `src/` | top-level `*.cpp` / `*.h` / `*.c` (excl. `RcppExports.cpp` from engine tables) |")
    L("| `src/nmath/` | Core vendored Mathlib (absent in glmbayes) |")
    L("| `inst/cl/` | `*.cl`, dependency manifests (`.tsv`/`.rds`), `README.md` |")
    L("")
    L("`RcppExports.R` / `RcppExports.cpp` are listed where present but treated as")
    L("generated noise for engine-parity narrative.")
    L("")
    L("---")
    L("")
    L("## 2. File inventory by type")
    L("")

    # R
    L("### 2.1 `R/`")
    L("")
    L(
        md_table_counts(
            [
                ("glmbayesCore `R/*.R`", len(r_core)),
                ("glmbayes `R/*.R`", len(r_gb)),
                ("Shared basenames", len(r_both)),
                ("Byte-identical", len(r_cls["identical_byte"])),
                ("Normalized-identical only", len(r_cls["identical_norm_only"])),
                ("Differing", len(r_cls["differing"])),
                ("glmbayesCore-only", len(r_core_only)),
                ("glmbayes-only", len(r_gb_only)),
            ]
        )
    )
    L("")
    L("#### Identical (byte)")
    L("")
    L(md_list(r_cls["identical_byte"]))
    L("")
    L("#### Identical after normalize only")
    L("")
    L(md_list(r_cls["identical_norm_only"]))
    L("")
    L("#### Differing (shared basename)")
    L("")
    L(md_list(r_cls["differing"]))
    L("")
    L("#### glmbayesCore-only")
    L("")
    L(md_list(r_core_only))
    L("")
    L("#### glmbayes-only")
    L("")
    L(md_list(r_gb_only))
    L("")

    # src
    L("### 2.2 `src/` (top-level, excl. `RcppExports.cpp`)")
    L("")
    L(
        md_table_counts(
            [
                ("glmbayesCore top-level", len(src_core)),
                ("glmbayes top-level", len(src_gb)),
                ("Shared basenames", len(src_both)),
                ("Byte-identical", len(src_cls["identical_byte"])),
                ("Normalized-identical only", len(src_cls["identical_norm_only"])),
                ("Differing", len(src_cls["differing"])),
                ("glmbayesCore-only", len(src_core_only)),
                ("glmbayes-only", len(src_gb_only)),
            ]
        )
    )
    L("")
    L("#### Identical (byte)")
    L("")
    L(md_list(src_cls["identical_byte"]))
    L("")
    L("#### Identical after normalize only")
    L("")
    L(md_list(src_cls["identical_norm_only"]))
    L("")
    L("#### Differing (shared basename)")
    L("")
    L(md_list(src_cls["differing"]))
    L("")
    L("#### glmbayesCore-only")
    L("")
    L(md_list(src_core_only))
    L("")
    L("#### glmbayes-only")
    L("")
    L(md_list(src_gb_only))
    L("")

    # nmath
    L("### 2.3 `src/nmath/` (Core-only vendored Mathlib)")
    L("")
    L(
        md_table_counts(
            [
                ("glmbayesCore `src/nmath` files", len(nmath_core)),
                ("glmbayes `src/nmath` files", len(nmath_gb)),
            ]
        )
    )
    L("")
    if nmath_core and not nmath_gb:
        L("**glmbayes** does not vendor `src/nmath/`; it relies on **nmathopencl**")
        L("(Imports) for OpenCL Mathlib pieces and does not ship the R Mathlib `.c`")
        L("tree under `src/`.")
        L("")
        L("<details><summary>glmbayesCore src/nmath file list</summary>")
        L("")
        L(md_list(nmath_core))
        L("")
        L("</details>")
        L("")
    else:
        L(md_list(nmath_core) if nmath_core else "*(none)*")
        L("")

    # cl
    L("### 2.4 `inst/cl/`")
    L("")
    L(
        md_table_counts(
            [
                ("glmbayesCore relevant paths", len(cl_core)),
                ("glmbayes relevant paths", len(cl_gb)),
                ("Shared relative paths", len(cl_both)),
                ("Byte-identical", len(cl_cls["identical_byte"])),
                ("Normalized-identical only", len(cl_cls["identical_norm_only"])),
                ("Differing", len(cl_cls["differing"])),
                ("glmbayesCore-only", len(cl_core_only)),
                ("glmbayes-only", len(cl_gb_only)),
            ]
        )
    )
    L("")
    L("#### Identical (byte)")
    L("")
    cl_only = [x for x in cl_cls["identical_byte"] if x.endswith(".cl")]
    L(
        f"All **{len(cl_only)}** shared `*.cl` files are byte-identical. "
        f"Total byte-identical shared paths (incl. manifests): "
        f"**{len(cl_cls['identical_byte'])}**."
    )
    L("")
    L("<details><summary>Byte-identical path list</summary>")
    L("")
    L(md_list(cl_cls["identical_byte"]))
    L("")
    L("</details>")
    L("")
    L("#### Identical after normalize only")
    L("")
    L(md_list(cl_cls["identical_norm_only"]))
    L("")
    L("#### Differing")
    L("")
    L(md_list(cl_cls["differing"]))
    L("")
    L("#### glmbayesCore-only")
    L("")
    L(md_list(cl_core_only))
    L("")
    L("#### glmbayes-only")
    L("")
    L(md_list(cl_gb_only))
    L("")
    L("#### Shared f2/f3 kernel entry points (reference)")
    L("")
    for rel, kernels in kernel_entries:
        L(f"- `{rel}`: {', '.join(f'`{k}`' for k in kernels) if kernels else '*(no __kernel parsed)*'}")
    L("")

    # DESCRIPTION OpenCL
    L("### 2.5 OpenCL-related DESCRIPTION fields")
    L("")
    L("| Field | glmbayesCore | glmbayes |")
    L("|-------|:------------:|:--------:|")
    L(f"| Imports `opencltools` | {ocl_core['Imports_opencltools']} | {ocl_gb['Imports_opencltools']} |")
    L(f"| Imports `nmathopencl` | {ocl_core['Imports_nmathopencl']} | {ocl_gb['Imports_nmathopencl']} |")
    L(f"| LinkingTo `opencltools` | {ocl_core['LinkingTo_opencltools']} | {ocl_gb['LinkingTo_opencltools']} |")
    L("")

    L("---")
    L("")
    L("## 3. Function-level report (differing files)")
    L("")
    L("For each shared file that is not byte- or normalize-identical, top-level")
    L("functions / kernels are parsed and compared.")
    L("")

    if not cl_cls["differing"]:
        L("### 3.0 `inst/cl`")
        L("")
        L("N/A — all shared `inst/cl` paths are byte-identical; no per-file kernel")
        L("body diffs. See kernel entry list in §2.4.")
        L("")

    def emit_fn_section(title: str, prefix: str):
        L(f"### {title}")
        L("")
        keys = [k for k in fn_reports if k.startswith(prefix)]
        if not keys:
            L("*(no differing files)*")
            L("")
            return
        for key in keys:
            rep = fn_reports[key]
            L(f"#### `{key}`")
            L("")
            if rep.get("note"):
                L(rep["note"])
                L("")
                continue
            res = rep["result"]
            if res is None:
                L("*(skipped)*")
                L("")
                continue
            L(
                f"Parsed symbols — Core: {rep['counts']['core']}, "
                f"glmbayes: {rep['counts']['glmbayes']}."
            )
            L("")
            L(
                md_table_counts(
                    [
                        ("Identical body", len(res["identical"])),
                        ("Rename-only (package string)", len(res["rename_only"])),
                        ("Different body", len(res["different"])),
                        ("Core-only", len(res["core_only"])),
                        ("glmbayes-only", len(res["glmbayes_only"])),
                    ]
                )
            )
            L("")
            if res["identical"]:
                L("**Identical:** " + ", ".join(f"`{x}`" for x in res["identical"]))
                L("")
            if res["rename_only"]:
                L(
                    "**Rename-only:** "
                    + ", ".join(f"`{x}`" for x in res["rename_only"])
                )
                L("")
            if res["different"]:
                L(
                    "**Different:** "
                    + ", ".join(f"`{x}`" for x in res["different"])
                )
                L("")
            if res["core_only"]:
                L(
                    "**Core-only:** "
                    + ", ".join(f"`{x}`" for x in res["core_only"])
                )
                L("")
            if res["glmbayes_only"]:
                L(
                    "**glmbayes-only:** "
                    + ", ".join(f"`{x}`" for x in res["glmbayes_only"])
                )
                L("")

    emit_fn_section("3.1 Differing `R/` files", "R/")
    emit_fn_section("3.2 Differing `src/` files", "src/")
    if cl_cls["differing"]:
        emit_fn_section("3.3 Differing `inst/cl/` files", "inst/cl/")

    # Aggregate rename-only vs different across src/R for summary
    def _bucket_totals(prefix: str) -> tuple[int, int, int]:
        ident = ren = diff = 0
        for k, v in fn_reports.items():
            if not k.startswith(prefix):
                continue
            res = v.get("result")
            if not res:
                continue
            ident += len(res["identical"])
            ren += len(res["rename_only"])
            diff += len(res["different"])
        return ident, ren, diff

    r_id, r_ren, r_diff = _bucket_totals("R/")
    s_id, s_ren, s_diff = _bucket_totals("src/")
    kl = fn_reports.get("src/kernel_loader.cpp", {}).get("result") or {}

    cl_cl_files = [x for x in cl_both if x.endswith(".cl")]
    cl_cl_ident = [
        x for x in cl_cls["identical_byte"] if x.endswith(".cl")
    ]

    # Summary — polished narrative (kept in script so regenerates)
    L("---")
    L("")
    L("## 4. Summary overview")
    L("")
    L("### Snapshot")
    L("")
    L("| Layer | Shared | Byte-identical | Norm-only | Differing |")
    L("|-------|-------:|---------------:|----------:|----------:|")
    L(
        f"| `R/` | {len(r_both)} | {len(r_cls['identical_byte'])} | "
        f"{len(r_cls['identical_norm_only'])} | {len(r_cls['differing'])} |"
    )
    L(
        f"| `src/` (top-level) | {len(src_both)} | {len(src_cls['identical_byte'])} | "
        f"{len(src_cls['identical_norm_only'])} | {len(src_cls['differing'])} |"
    )
    L(
        f"| `inst/cl/` | {len(cl_both)} | {len(cl_cls['identical_byte'])} | "
        f"{len(cl_cls['identical_norm_only'])} | {len(cl_cls['differing'])} |"
    )
    L("")
    L("### What is already the same")
    L("")
    L(f"- **OpenCL program sources:** all **{len(cl_cl_ident)}** shared `*.cl` files")
    L("  are **byte-identical** (entry kernels under `src/f2_f3_*.cl` plus")
    L("  nmath/shim prelude). Manifests (`.tsv`/`.rds`) also match.")
    L("  The only shared `inst/cl` text that differs is **`README.md`**")
    L("  (documentation wording, not kernels).")
    L(f"- **C++ already in lockstep (byte):** {len(src_cls['identical_byte'])} files,")
    L("  notably **`kernel_runners.cpp`**, `configure_OpenCL.cpp`,")
    L("  `opencl_detect.cpp`, `rNormalReg.cpp`, `rnorm_ct.cpp`, `invgamma_ct.cpp`,")
    L("  `EnvelopeEval.cpp`, `EnvelopeSort.cpp`, `Set_Grid.cpp`, `Set_LogP.cpp`,")
    L("  `famfuncs.h`, `famfuncs_poisson.cpp`, `cuda_probe.cpp`, `rng_utils.h`.")
    L(f"- **C++ normalize-only (line endings):** {', '.join(f'`{x}`' for x in src_cls['identical_norm_only']) or '*(none)*'}.")
    L(f"- **R already in lockstep (byte):** {', '.join(f'`{x}`' for x in r_cls['identical_byte'])}.")
    L(f"- **R normalize-only:** {', '.join(f'`{x}`' for x in r_cls['identical_norm_only']) or '*(none)*'}.")
    L("")
    L("### Structural / packaging differences (not logic forks)")
    L("")
    L("| Topic | glmbayesCore | glmbayes |")
    L("|-------|--------------|----------|")
    L("| Package identity | `package_ns.h`, Core symbol prefixes | No `package_ns.h`; glmbayes prefixes |")
    L("| Formula / S3 UX | Matrix API focus (`rglmb`/`rlmb`) | Full `glmb()`/`lmb()` + many S3 methods |")
    L("| Multi-response helpers | `multi_*` R files present | Not in this tree |")
    L("| Vendored CPU nmath under `src/` | Yes (%d files) | No |" % len(nmath_core))
    L("| OpenCL Mathlib at runtime | Imports `nmathopencl`; still vendors `inst/cl/nmath` | Same Imports; thin loader via opencltools C API |")
    L(f"| `LinkingTo: opencltools` | **{ocl_core['LinkingTo_opencltools']}** | **{ocl_gb['LinkingTo_opencltools']}** |")
    L("| `kernel_loader.cpp` | Fat in-tree `system.file` / dependency assembly | Thin opencltools C-API wrapper |")
    L("| Blocks leftover | Stripped | Still has `rNormalGLMBlocks.cpp` |")
    L("")
    L("### Function-level deltas on differing files (§3 aggregates)")
    L("")
    L("| Layer | Identical bodies | Rename-only | Different bodies |")
    L("|-------|-----------------:|------------:|-----------------:|")
    L(f"| Differing `R/` | {r_id} | {r_ren} | {r_diff} |")
    L(f"| Differing `src/` | {s_id} | {s_ren} | {s_diff} |")
    L("")
    L("Rename-only means the function text matches after substituting package")
    L("name / dynlib prefixes. **Different** is the residual that needs a human")
    L("merge decision when Core becomes the shared backend.")
    L("")
    L("### Substantive engine diffs (high level)")
    L("")
    if kl:
        kl_line = (
            f"- **OpenCL loader (Batch 7b target):** `kernel_loader.cpp` is a real "
            f"fork. Parsed symbols — identical {len(kl.get('identical', []))}, "
            f"different {len(kl.get('different', []))}, "
            f"Core-only {len(kl.get('core_only', []))}, "
            f"glmbayes-only {len(kl.get('glmbayes_only', []))}."
        )
    else:
        kl_line = (
            "- **OpenCL loader (Batch 7b target):** `kernel_loader.cpp` is a real "
            "fork (see §3.2)."
        )
    L(kl_line)
    L("  Meanwhile **`kernel_runners.cpp` is already byte-identical**, so the")
    L("  enqueue/runtime path does not need a port — only the program-assembly")
    L("  loader and `LinkingTo: opencltools`.")
    L("- **Envelope / famfuncs / samplers:** several large `src/` files still")
    L("  differ beyond rename (see §3.2 `different` lists for")
    L("  `EnvelopeBuild*.cpp`, `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp`,")
    L("  `export_wrappers.cpp`, …). Treat these as a separate sync track from")
    L("  OpenCL loader alignment.")
    L("- **Core-only R capabilities:** `rindepNormalGamma_reg_with_envelope`,")
    L("  multi-response samplers (`multi_*`), `ing_prior_guard` helpers.")
    L("- **glmbayes-only R layer:** formula modelling + diagnostics S3")
    L("  (`glmb.R`, `lmb.R`, `summary.glmb.R`, `predict.glmb.R`, …) — stays in")
    L("  glmbayes when Core is the engine.")
    L("")
    L("### Recommended Batch 7b actions (code not done here)")
    L("")
    L("1. Port glmbayes’ thin `kernel_loader.cpp` into Core; archive the fat")
    L("   loader under `src/backup/` (same pattern as glmbayes).")
    L("2. Add `LinkingTo: opencltools` to Core `DESCRIPTION`.")
    L("3. Keep shared `*.cl` trees as-is (already identical); optionally sync")
    L("   `inst/cl/README.md` wording. Later: decide whether to stop shipping")
    L("   vendored nmath/shims and rely solely on **nmathopencl**.")
    L("4. Do **not** port `rNormalGLMBlocks.cpp` back into Core as part of OpenCL")
    L("   alignment (mixed-model path belongs with gradual lmebayesCore merge).")
    L("")

    L("---")
    L("")
    L("## 5. Appendix")
    L("")
    L("### Regenerate")
    L("")
    L("```r")
    L("# From glmbayesCore package root")
    L('source("data-raw/make_compare_glmbayes_glmbayesCore.R")')
    L("# Optional override:")
    L('# Sys.setenv(GLMBAYES_COMPARE_ROOT = "C:/path/to/glmbayes")')
    L("```")
    L("")
    L("### Machine-readable sidecar")
    L("")
    L("`data-raw/compare_glmbayes_glmbayesCore.json` — file hashes and function")
    L("bucket lists for the differing files.")
    L("")

    out_md = CORE_ROOT / "data-raw" / "COMPARE_glmbayes_vs_glmbayesCore.md"
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    sidecar = {
        "generated": today,
        "glmbayesCore": {"path": str(CORE_ROOT), "sha": core_sha},
        "glmbayes": {"path": str(gb_root), "sha": gb_sha},
        "R": {
            "identical_byte": r_cls["identical_byte"],
            "identical_norm_only": r_cls["identical_norm_only"],
            "differing": r_cls["differing"],
            "core_only": r_core_only,
            "glmbayes_only": r_gb_only,
        },
        "src": {
            "identical_byte": src_cls["identical_byte"],
            "identical_norm_only": src_cls["identical_norm_only"],
            "differing": src_cls["differing"],
            "core_only": src_core_only,
            "glmbayes_only": src_gb_only,
            "nmath_core_count": len(nmath_core),
        },
        "cl": {
            "identical_byte": cl_cls["identical_byte"],
            "identical_norm_only": cl_cls["identical_norm_only"],
            "differing": cl_cls["differing"],
            "core_only": cl_core_only,
            "glmbayes_only": cl_gb_only,
        },
        "opencl_description": {"glmbayesCore": ocl_core, "glmbayes": ocl_gb},
        "functions": {
            k: (v["result"] if v.get("result") is not None else {"note": v.get("note")})
            for k, v in fn_reports.items()
        },
    }
    out_json = CORE_ROOT / "data-raw" / "compare_glmbayes_glmbayesCore.json"
    out_json.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {out_md}")
    print(f"Wrote {out_json}")
    print(
        f"R: ident={len(r_cls['identical_byte'])} diff={len(r_cls['differing'])}; "
        f"src: ident={len(src_cls['identical_byte'])} diff={len(src_cls['differing'])}; "
        f"cl: ident={len(cl_cls['identical_byte'])} diff={len(cl_cls['differing'])}"
    )


if __name__ == "__main__":
    main()
