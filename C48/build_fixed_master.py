"""Build Magma's fixed integer input table, requiring a complete independent source export."""
import hashlib
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent / "data"
source = HERE / "source_partition_inputs_complete.py"
a = runpy.run_path(str(source))
assert a.get("complete") is True and a.get("complete_by_newton_count") is True
assert set(a["representative_pairs"]) == set(range(3,9))
assert [a["noncuspidal_counts"][d] for d in range(3,9)] == [16,64,656,2712,10784,48992]
assert all(sum(a["partition_orbit_sizes"][d]) == a["noncuspidal_counts"][d] for d in range(3,9))
assert all(len(a["representative_pairs"][d]) == len(a["partition_orbit_sizes"][d]) for d in range(3,9))
sha = hashlib.sha256(source.read_bytes()).hexdigest()
out = HERE / "fixed_source_inputs.m"
lines = ["/* Fixed actual order48 moduli; degree8 completeness uses Newton cardinality. */",
         f'C48MasterSHA := "{sha}";', "C48Degrees := [3,4,5,6,7,8];"]
for key, name in (("field_polynomials", "C48FieldPolynomials"),
                  ("representative_pairs", "C48RepresentativePairs"),
                  ("partition_orbit_sizes", "C48OrbitSizes")):
    lines.append(f"{name} := {repr([a[key][d] for d in range(3,9)])};")
lines.append("C48NoncuspidalCounts := [16,64,656,2712,10784,48992];")
text = "\n".join(lines) + "\n"
if out.exists():  # never change the committed fixed table; a rerun must reproduce it
    assert out.read_text() == text, "fixed input table differs"
else:
    out.write_text(text)
print("C48_FIXED_SOURCE_MASTER_BUILT", sha,
      "representative_counts", [len(a["representative_pairs"][d]) for d in range(3,9)])
