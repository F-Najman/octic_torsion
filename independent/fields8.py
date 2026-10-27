# fields8.py -- LMFDB identification of the 12 realization fields and of a rational curve with j = 1331/8.
# Reads the FIELD lines of logs/examples8.log (Magma OptimizedRepresentation polynomials), computes the
# PARI polredabs polynomial (LMFDB normal form), and looks it up in nf_fields on the LMFDB mirror.
import re, sys, json
import psycopg2
from sage.all import pari, QQ, PolynomialRing

R = PolynomialRing(QQ, "t")
log = open("independent/logs/examples8.log").read()
fields = re.findall(r"^FIELD (\S+): optimized polynomial (.*?) ; disc (\S+) = .*?; signature (\d+) ; Galois group (\S+)", log, re.M)
conn = psycopg2.connect(host="devmirror.lmfdb.xyz", dbname="lmfdb", user="lmfdb", password="lmfdb", connect_timeout=60)
cur = conn.cursor()
out = []
for name, poly, disc, r1, gal in fields:
    f = R(poly.replace("^", "**"))
    pr = pari(f).polredabs()
    coeffs = [int(c) for c in pr.Vec()][::-1]          # constant term first, as in LMFDB
    cur.execute("select label, disc_sign*disc_abs, r2, galois_label, class_number from nf_fields where coeffs = %s::numeric[]", (coeffs,))
    rows = cur.fetchall()
    lab = rows[0] if rows else None
    D = int(pari(pr).nfdisc())
    out.append(dict(entry=name, polredabs=str(pr), disc=D, r1=int(r1), magma_gal=gal, lmfdb=None if lab is None else
                    dict(label=lab[0], disc=int(lab[1]), r2=int(lab[2]), galois=lab[3], h=int(lab[4]) if lab[4] is not None else None)))
    print("FIELDLMFDB %s: polredabs %s ; disc %d (Magma %s) ; signature (%d,%d) ; Galois %s ; LMFDB %s" % (
        name, pr, D, disc, int(r1), (8 - int(r1)) // 2, gal, "not found" if lab is None else
        "%s (disc %s, r2 %s, galois %s, h %s)" % lab))
    assert D == int(disc)
# rational elliptic curves with j = 1331/8
cur.execute("select lmfdb_label, conductor, ainvs, torsion_structure, cm from ec_curvedata where jinv = '{1331,8}'::numeric[] order by conductor limit 6")
rows = cur.fetchall()
print("ECJ 1331/8: LMFDB curves with this j (smallest conductors):", rows)
json.dump(out, open("independent/data/fields8.json", "w"), indent=1)
print("FIELDS8_DONE")
