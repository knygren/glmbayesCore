"""Regenerate glmbayesCore/R/lmb.R from glmbayes/R/lmb.R (Phase 1 staging)."""
from pathlib import Path

src = Path(r"C:/Rpackages/glmbayes/R/lmb.R").read_text(encoding="utf-8")
lines = src.splitlines(keepends=True)
parts = [
    "".join(lines[154:279]),
    "".join(lines[564:823]),
    "".join(lines[856:980]),
    "".join(lines[1086:1133]),
]
out = "".join(parts)
out = out.replace("{glmbayes}", "{glmbayesCore}")
out = out.replace("@example inst/examples/Ex_lmb.R", "@examples ## See tests/testthat/test-lmb.R")
header = (
    "#' Temporary Phase 1 staging in glmbayesCore (same API as glmbayes \\code{lmb()}).\n"
    "#' S3 for \\code{lmb} remains in glmbayes until Phase 3.\n\n"
)
Path(r"C:/Rpackages/glmbayesCore/R/lmb.R").write_text(header + out, encoding="utf-8")
print("wrote R/lmb.R", len(out), "bytes")
