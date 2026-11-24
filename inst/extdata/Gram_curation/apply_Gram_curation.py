#!/usr/bin/env python3
"""
apply_Gram_curation.py

Re-derives the Gram column of the JAMS taxonomy table (data(JAMStaxtable)) from an explicit,
cited rule set, and writes a per-taxid provenance file so that every value can be traced.

Usage (Python >= 3.8, standard library only):

    python3 -I apply_Gram_curation.py \
        --taxtable  JAMStaxtable_<date>.tsv \
        --rules     Gram_evidence_rules.tsv \
        --out_table JAMStaxtable_<date>_Gram_curated.tsv \
        --out_provenance JAMStaxtable_Gram_provenance.tsv.gz \
        --out_review Gram_curation_review.tsv \
        --out_summary Gram_curation_summary.txt

The input is the table exported from R with
    write.table(JAMStaxtable, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
Taxid must be the first column (it is also the row name in R) and Gram the last.
The output table is byte-identical to the input except for the Gram field.

See README_Gram_curation.md for the definitions, the order of precedence and its rationale.
"""

import argparse
import collections
import csv
import gzip
import re
import sys

#Precedence ("tiers") - see README_Gram_curation.md
#  tier 1 : clade-wide N_A rules (non-bacterial domains, Mollicutes). Beat everything.
#  tier 3 : species-level evidence from described isolates (rules with match_column Species or Taxid).
#           Beats tier 4, by decision of J. McCulloch (published isolate evidence beats the MAG rule).
#  tier 4 : no isolate: "uncultured" in the name, Candidatus in the name (switchable), or
#           Proportion_of_MAG_in_taxid == 1 -> undetermined.
#  tier 5 : clade rules (Phylum ... Genus), the deepest matching rank wins.
#  tier 6 : nothing applies -> undetermined.
#  tier 7 : bacterial rank nodes (Genus and above) take the consensus of their terminal
#           descendants: one value -> that value; several -> variable; none known -> undetermined.
#           N_A counts as a known value here (so Bacteria as a whole is variable). Rows decided at
#           tier 1 or tier 3 are never overridden. Intermediate NCBI ranks (subfamily, clade, species
#           group, suborder...) are not rank nodes: their descendants cannot be resolved from the
#           lineage columns, so they keep the value from tiers 1 to 6.
#  tier 8 : rows that share an LKT with the taxid named in that LKT take that taxid's value.

LINEAGE_DEEPEST_FIRST = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]
RANK_NODE_COLUMNS = ["Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]
#NCBI ranks at or above genus level, or otherwise not a terminal organism; never counted as descendants in tier 7
SUPRASPECIFIC_RANKS = {"superkingdom", "domain", "kingdom", "subkingdom", "phylum", "subphylum", "superclass", "class",
                       "subclass", "infraclass", "cohort", "subcohort", "superorder", "order", "suborder", "infraorder",
                       "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus",
                       "section", "subsection", "series", "species group", "species subgroup", "clade", "root"}
ALLOWED_GRAM = {"positive", "negative", "variable", "N_A", "undetermined"}
UNINFORMATIVE_SPECIES = {"s__Unclassified", "s__Missing"}

TAXID_IN_NAME = re.compile(r"^(?:(?:LKT|CSB)__)?[a-z0-9]+__([0-9]+)_")


def taxid_in_name(name):
    m = TAXID_IN_NAME.match(name or "")
    return m.group(1) if m else None


def load_rules(path):
    rules = {}
    notes = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = {"rule_id", "tier", "match_column", "match_value", "assigned_Gram", "evidence_type", "basis", "ref_ids", "verification", "review_flag", "notes"} - set(rd.fieldnames)
        if missing:
            sys.exit("Rules file is missing column(s): " + ", ".join(sorted(missing)))
        for r in rd:
            r = {k: (v if v is not None else "") for k, v in r.items()}
            tier = r["tier"].strip()
            key = (r["match_column"].strip(), r["match_value"].strip())
            g = r["assigned_Gram"].strip()
            if tier == "note":
                notes[key] = r
                continue
            if g not in ALLOWED_GRAM and g != "pass":
                sys.exit("Rule " + r["rule_id"] + " has an invalid assigned_Gram: " + g)
            if (tier, key) in rules:
                sys.exit("Duplicate rule for " + str(key) + " in tier " + tier)
            rules[(tier, key)] = r
    return rules, notes


def find_rule(rules, tier, row, columns):
    """Deepest matching rule of the given tier, or None."""
    for col in columns:
        r = rules.get((tier, (col, row[col])))
        if r is not None:
            return r
    return None


def is_mag_only(row):
    v = row["Proportion_of_MAG_in_taxid"]
    try:
        return float(v) >= 0.99999
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxtable", required=True)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--out_table", required=True)
    ap.add_argument("--out_provenance", required=True)
    ap.add_argument("--out_review", required=True)
    ap.add_argument("--out_summary", required=True)
    ap.add_argument("--candidatus_is_undetermined", default="yes", choices=["yes", "no"],
                    help="Treat Candidatus taxa (no pure culture by definition) as undetermined. Default yes.")
    args = ap.parse_args()

    rules, notes = load_rules(args.rules)

    #Read raw lines so that every field other than Gram is written back byte for byte
    with open(args.taxtable, "r", newline="") as fh:
        raw = fh.read().split("\n")
    trailing_newline = raw[-1] == ""
    if trailing_newline:
        raw = raw[:-1]
    header = raw[0].split("\t")
    if header[0] != "Taxid" or header[-1] != "Gram":
        sys.exit("Expected Taxid as the first column and Gram as the last.")
    ncol = len(header)
    rows = []
    for i, line in enumerate(raw[1:], start=2):
        f = line.split("\t")
        if len(f) != ncol:
            sys.exit("Line " + str(i) + " has " + str(len(f)) + " fields, expected " + str(ncol))
        rows.append(dict(zip(header, f)))
    taxids = [r["Taxid"] for r in rows]
    if len(set(taxids)) != len(taxids):
        sys.exit("Duplicate Taxid values in the input table.")

    new_gram = [None] * len(rows)
    prov = [None] * len(rows)

    #---- tiers 1 to 6, row by row ----
    for i, row in enumerate(rows):
        decided = None
        r1 = find_rule(rules, "1", row, LINEAGE_DEEPEST_FIRST)
        if r1 is not None and r1["assigned_Gram"] != "pass":
            decided = (r1["assigned_Gram"], r1["rule_id"], r1)
        if decided is None:
            r3 = rules.get(("3", ("Taxid", row["Taxid"]))) or rules.get(("3", ("Species", row["Species"])))
            if r3 is not None:
                decided = (r3["assigned_Gram"], r3["rule_id"], r3)
        if decided is None:
            lkt_low = row["LKT"].lower()
            cand = args.candidatus_is_undetermined == "yes" and any("Candidatus" in row[c] for c in ("LKT", "Species", "Genus"))
            #Name-based reasons first, as they are the more informative label; all three give undetermined
            if "uncultured" in lkt_low:
                decided = ("undetermined", "G4_uncultured", None)
            elif cand:
                decided = ("undetermined", "G4_Candidatus", None)
            elif is_mag_only(row):
                decided = ("undetermined", "G4_MAG_only", None)
        if decided is None:
            r5 = find_rule(rules, "5", row, LINEAGE_DEEPEST_FIRST)
            if r5 is not None:
                decided = (r5["assigned_Gram"], r5["rule_id"], r5)
        if decided is None:
            decided = ("undetermined", "G6_no_rule", None)
        new_gram[i] = decided[0]
        prov[i] = {"decided_by": decided[1], "rule": decided[2], "descendants": "",
                   "tier": decided[2]["tier"] if decided[2] is not None else ""}

    #---- tier 7: bacterial rank nodes take the consensus of their terminal descendants ----
    def is_terminal(row):
        rk = row["NCBI_taxonomic_rank"]
        if rk in SUPRASPECIFIC_RANKS:
            return False
        if rk == "no rank" and row["Species"] in UNINFORMATIVE_SPECIES:
            return False
        return True

    node_of = {}
    for i, row in enumerate(rows):
        if row["Domain"] != "d__2_Bacteria":
            continue
        for col in RANK_NODE_COLUMNS:
            if taxid_in_name(row[col]) == row["Taxid"]:
                node_of[i] = (col, row[col])
                break
    wanted = set(node_of.values())
    desc = collections.defaultdict(collections.Counter)
    for i, row in enumerate(rows):
        if row["Domain"] != "d__2_Bacteria" or not is_terminal(row):
            continue
        for col in RANK_NODE_COLUMNS:
            key = (col, row[col])
            if key in wanted:
                desc[key][new_gram[i]] += 1
    for i, key in node_of.items():
        #Tier 1 (N_A clades) and tier 3 (isolate evidence for this very taxid or species) are never overridden
        if prov[i]["tier"] in ("1", "3"):
            continue
        c = desc.get(key)
        if not c:
            #No terminal descendants in the table: keep the value from tiers 1 to 6
            prov[i]["descendants"] = "no terminal descendants in table"
            continue
        known = {g for g in c if g != "undetermined"}
        if not known:
            g = "undetermined"
        elif len(known) == 1:
            g = next(iter(known))
        else:
            g = "variable"
        new_gram[i] = g
        prov[i]["decided_by"] = "G7_consensus"
        prov[i]["rule"] = None
        prov[i]["descendants"] = "; ".join(k + "=" + str(v) for k, v in sorted(c.items()))

    #---- tier 8: one LKT, one Gram ----
    #Several taxids can carry the same LKT (e.g. an NCBI "no rank" or "clade" node placed directly
    #under a phylum is labelled LKT__p__<phylum>). JAMS reads Gram onto features by the taxid embedded
    #in the LKT name, so the taxid named in the LKT is the canonical one, and every other row with that
    #LKT takes its value. This guarantees that a feature called, say, LKT__c__909932_Negativicutes has the
    #Gram of taxid 909932 whichever row of the table one looks at.
    row_of_taxid = {r["Taxid"]: i for i, r in enumerate(rows)}
    for i, row in enumerate(rows):
        tid = taxid_in_name(row["LKT"])
        if (tid is None) or (tid == row["Taxid"]) or (tid not in row_of_taxid):
            continue
        j = row_of_taxid[tid]
        if new_gram[i] != new_gram[j]:
            prov[i]["descendants"] = "value of taxid " + tid + " (" + rows[j]["LKT"] + "); own rules gave " + new_gram[i]
            new_gram[i] = new_gram[j]
            prov[i]["decided_by"] = "G8_same_LKT"
            prov[i]["rule"] = None

    #---- write the curated table ----
    gidx = ncol - 1
    out_lines = [raw[0]]
    for i, line in enumerate(raw[1:]):
        f = line.split("\t")
        f[gidx] = new_gram[i]
        out_lines.append("\t".join(f))
    with open(args.out_table, "w", newline="") as fh:
        fh.write("\n".join(out_lines) + ("\n" if trailing_newline else ""))

    #---- provenance, one line per taxid ----
    pcols = ["Taxid", "LKT", "NCBI_taxonomic_rank", "Domain", "Phylum", "Proportion_of_MAG_in_taxid", "Num_isolate",
             "Gram_previous", "Gram_new", "changed", "decided_by", "evidence_type", "basis", "ref_ids", "verification",
             "review_flag", "descendant_values", "notes"]
    basis_g = {"G4_MAG_only": "No isolate: every genome assigned to this taxid is a MAG (Proportion_of_MAG_in_taxid = 1)",
               "G4_uncultured": "No isolate: 'uncultured' in the taxon name",
               "G4_Candidatus": "No pure culture: Candidatus taxon",
               "G6_no_rule": "No evidence-backed rule covers this lineage",
               "G7_consensus": "Rank node: consensus of its terminal descendants in this table",
               "G8_same_LKT": "Shares its LKT with the taxid named in that LKT, and takes that taxid's value"}
    opener = gzip.open if args.out_provenance.endswith(".gz") else open
    with opener(args.out_provenance, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(pcols)
        for i, row in enumerate(rows):
            p = prov[i]
            r = p["rule"]
            notes_txt = []
            if r is not None and r.get("notes"):
                notes_txt.append(r["notes"])
            #Note rules attach to the taxid itself and to everything below the taxon they name
            for key in [("Taxid", row["Taxid"])] + [(c, row[c]) for c in LINEAGE_DEEPEST_FIRST]:
                nt = notes.get(key)
                if nt is not None:
                    notes_txt.append("[" + nt["rule_id"] + "] " + nt["basis"] + " Refs: " + nt["ref_ids"])
            w.writerow([row["Taxid"], row["LKT"], row["NCBI_taxonomic_rank"], row["Domain"], row["Phylum"],
                        row["Proportion_of_MAG_in_taxid"], row["Num_isolate"], row["Gram"], new_gram[i],
                        "yes" if row["Gram"] != new_gram[i] else "no", p["decided_by"],
                        r["evidence_type"] if r is not None else ("consensus" if p["decided_by"] == "G7_consensus" else "same_LKT" if p["decided_by"] == "G8_same_LKT" else "no_isolate" if p["decided_by"].startswith("G4") else "none"),
                        r["basis"] if r is not None else basis_g.get(p["decided_by"], ""),
                        r["ref_ids"] if r is not None else "", r["verification"] if r is not None else "",
                        r["review_flag"] if r is not None else "", p["descendants"], " | ".join(notes_txt)])

    #---- summary ----
    lines = []
    def tab(title, counter):
        lines.append(title)
        for k, v in sorted(counter.items(), key=lambda kv: (-kv[1], str(kv[0]))):
            lines.append("    " + str(k) + "\t" + str(v))
        lines.append("")
    tab("Gram, previous table:", collections.Counter(r["Gram"] for r in rows))
    tab("Gram, curated table:", collections.Counter(new_gram))
    tab("Curated Gram by domain:", collections.Counter((r["Domain"], new_gram[i]) for i, r in enumerate(rows)))
    tab("Decided by:", collections.Counter(p["decided_by"] for p in prov))
    tab("Transitions previous -> curated (all rows):", collections.Counter((r["Gram"], new_gram[i]) for i, r in enumerate(rows) if r["Gram"] != new_gram[i]))
    bact = [i for i, r in enumerate(rows) if r["Domain"] == "d__2_Bacteria"]
    tab("Bacteria, curated Gram by phylum (top 40 phyla):", collections.Counter())
    byp = collections.defaultdict(collections.Counter)
    for i in bact:
        byp[rows[i]["Phylum"]][new_gram[i]] += 1
    for ph, c in sorted(byp.items(), key=lambda kv: -sum(kv[1].values()))[:40]:
        lines.append("    " + ph + "\t" + ", ".join(k + "=" + str(v) for k, v in sorted(c.items())))
    lines.append("")
    with open(args.out_summary, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    #---- review list ----
    with open(args.out_review, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "item", "n_taxids", "n_with_isolates", "detail", "examples"])
        #(a) rules flagged for review, with how many taxids each decides
        used = collections.Counter(p["decided_by"] for p in prov)
        for (tier, key), r in sorted(rules.items(), key=lambda kv: kv[1]["rule_id"]):
            if r["review_flag"] == "yes":
                w.writerow(["A_rule_flagged_for_review", r["rule_id"] + " " + key[0] + "=" + key[1] + " -> " + r["assigned_Gram"],
                            used.get(r["rule_id"], 0), "", r["verification"] + ". " + r["notes"], ""])
        for key, r in notes.items():
            w.writerow(["A_note_rule", r["rule_id"] + " " + key[0] + "=" + key[1], "", "", r["basis"], ""])
        #(b) direct positive <-> negative flips, by deciding rule
        flips = collections.defaultdict(list)
        for i, r in enumerate(rows):
            a, b = r["Gram"], new_gram[i]
            if {a, b} == {"positive", "negative"}:
                flips[(a, b, prov[i]["decided_by"])].append(i)
        for (a, b, rid), idx in sorted(flips.items(), key=lambda kv: -len(kv[1])):
            iso = sum(1 for i in idx if rows[i]["Num_isolate"] not in ("NA", "0", ""))
            ex = "; ".join(rows[i]["LKT"] for i in idx[:8])
            w.writerow(["B_flip_" + a + "_to_" + b, rid, len(idx), iso, "", ex])
        #(c) bacterial phyla with no rule: taxids left undetermined that have isolates
        noru = collections.defaultdict(list)
        for i in bact:
            if prov[i]["decided_by"] == "G6_no_rule":
                noru[rows[i]["Phylum"]].append(i)
        for ph, idx in sorted(noru.items(), key=lambda kv: -len(kv[1])):
            iso = [i for i in idx if rows[i]["Num_isolate"] not in ("NA", "0", "")]
            w.writerow(["C_phylum_without_rule", ph, len(idx), len(iso), "Left undetermined; needs an evidence-backed rule", "; ".join(rows[i]["LKT"] for i in iso[:6])])
        #(d) previously positive/negative, now undetermined, by reason
        lost = collections.defaultdict(list)
        for i, r in enumerate(rows):
            if r["Gram"] in ("positive", "negative") and new_gram[i] == "undetermined":
                lost[(r["Gram"], prov[i]["decided_by"])].append(i)
        for (a, rid), idx in sorted(lost.items(), key=lambda kv: -len(kv[1])):
            iso = sum(1 for i in idx if rows[i]["Num_isolate"] not in ("NA", "0", ""))
            w.writerow(["D_" + a + "_to_undetermined", rid, len(idx), iso, "", "; ".join(rows[i]["LKT"] for i in idx[:8])])
        #(e) MAG-only taxids that may nonetheless rest on a described isolate (decision needed)
        sp_row = {r["Species"]: i for i, r in enumerate(rows) if r["NCBI_taxonomic_rank"] == "species"}
        e_species, e_infra = [], []
        for i, r in enumerate(rows):
            if prov[i]["decided_by"] != "G4_MAG_only":
                continue
            if r["NCBI_taxonomic_rank"] == "species":
                if not re.search(r"Unclassified|bacterium|symbiont|metagenome|Candidatus|uncultured", r["Species"], re.I):
                    e_species.append(i)
            elif r["Species"] in sp_row and not is_mag_only(rows[sp_row[r["Species"]]]):
                e_infra.append(i)
        w.writerow(["E_MAG_only_named_species", "species-rank binomials with Proportion_of_MAG_in_taxid = 1", len(e_species), "",
                    "A validly published binomial normally implies a deposited type strain, so these may have described isolates whose genomes are simply not in NCBI. A strain-collection source (e.g. BacDive) would settle them.",
                    "; ".join(rows[i]["LKT"] for i in e_species[:12])])
        w.writerow(["E_MAG_only_strain_of_cultured_species", "strain-level taxids that are MAG-only, whose species has isolate genomes", len(e_infra), "",
                    "Kept undetermined by the MAG rule as written; they could instead follow their species.",
                    "; ".join(rows[i]["LKT"] for i in e_infra[:12])])

    print("Done. Rows: " + str(len(rows)) + "; changed: " + str(sum(1 for i, r in enumerate(rows) if r["Gram"] != new_gram[i])))


if __name__ == "__main__":
    main()
