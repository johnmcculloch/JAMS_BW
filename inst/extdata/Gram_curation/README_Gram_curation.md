# Gram curation of the JAMS taxonomy table (JAMStaxtable)

These rules were applied to `data/JAMStaxtable.rda` on 2026-11-24.

This folder documents how the `Gram` column of `data(JAMStaxtable)` was derived. Nothing here is loaded into an R session by JAMS. Read it if you want to know why a taxon carries the Gram value it does.

| File | What it is |
|---|---|
| `README_Gram_curation.md` | This document: definitions, rules, order of precedence, limitations. |
| `Gram_evidence_rules.tsv` | Every rule used, one row each, with the evidence and the references behind it. |
| `Gram_references.tsv` | The references, and what was actually read of each. |
| `apply_Gram_curation.py` | The script that applies the rules to the table. Python 3, standard library only, so anyone can rerun it. |
| `JAMStaxtable_Gram_provenance.tsv.gz` (to be added once approved) | One line per taxid: previous value, new value, the rule that decided it, the type of evidence (envelope, stain, definition, consensus...), its basis, references, and any notes. |

## Disclaimer

Gram is not a taxonomic property. It is not even meaningful for many organisms. Assigning it from taxonomy is an inference that holds for many taxa, but not all of them. The values in JAMStaxtable are intended for annotating plots and for broad, exploratory grouping of microbiome features. They are not a substitute for staining, or for a genome or lipid analysis of a particular strain. Where the evidence is weak or absent, the table says "undetermined" instead of guessing. Any JAMS user relying on Gram for a conclusion should check the provenance of the taxa involved.

## What the values mean

The column describes cell-envelope architecture, with the human immune system in mind, **except where a Gram stain has been published for the taxon, in which case the published stain wins**. The reason is practical: a reader who looks up the description of an organism and finds "Gram-stain-negative" must not find "positive" in JAMS. Where published stains conflict, or are reported as Gram-variable, the value is `variable`.

| Value | Meaning |
|---|---|
| `negative` | A published Gram-negative stain; otherwise, diderm: an outer membrane potentially bearing lipopolysaccharide (LPS) or lipooligosaccharide (LOS). |
| `positive` | A published Gram-positive stain; otherwise, monoderm: no outer membrane. |
| `variable` | (a) An outer membrane without LPS, for example Borrelia, Treponema, Sphingomonas (glycosphingolipids), Anaplasmataceae, Deinococcota, Thermotogota, Synergistota and the mycolic-acid layer of Mycobacteriales; (b) published stains that conflict or are reported as Gram-variable, for example Faecalibacterium, Gemmiger, Eisenbergiella, Acutalibacter, Intestinimonas; or (c) a rank above species whose descendants include more than one category. |
| `N_A` | No Gram-type cell wall expected: viruses, Eukaryota (fungi, protists), all Archaea, the wall-less Mycoplasmatota (Mollicutes, Mycoplasmoidales, and the Izimaplasma/Izemoplasmatales clade), and Dehalococcoides (no peptidoglycan). |
| `undetermined` | No defensible call. This covers taxa known only from metagenome-assembled genomes, "uncultured" taxa, Candidatus taxa, and lineages for which no evidence-backed rule exists. |

Decisions taken by J. McCulloch on 2026-11-23/24, in consultation with his supervisor:

1. Envelope architecture with LPS potential defines positive and negative, but a published Gram stain wins over it (revised 2026-11-24; the first draft used envelope alone).
2. A taxid whose genomes are all MAGs (`Proportion_of_MAG_in_taxid` = 1) is undetermined. So is a higher clade made only of such taxids.
3. Taxa named "uncultured" are undetermined.
4. N_A means "not applicable": no Gram-type envelope is expected (see the table above).
5. An outer membrane without LPS is `variable`.
6. Consequently, monoderm Clostridia with a published Gram-negative stain (Pseudoruminococcus, Oscillibacter, Subdoligranulum, Ruthenibacterium, Dysosmobacter, Vescimonas and others) are `negative`. An outer membrane without LPS stays `variable` even where a stain is published (e.g. Mycobacteriales, Treponema), because `variable` does not contradict either stain.
7. Published evidence for a described isolate beats the MAG-only rule.

Two further choices were made during implementation and are open to revision:

- **Candidatus taxa** are treated like "uncultured" ones, because Candidatus names are given to organisms without a pure culture. Switch this off with `--candidatus_is_undetermined no`. This choice decides 2,566 taxids. Its boundary is nomenclatural rather than biological: 552 of those taxids have non-MAG genomes (for example Candidatus Pelagibacter, which has since been cultured), and unculturable endosymbionts with valid names such as Buchnera receive their clade value while Candidatus endosymbionts such as Blochmanniella do not.
- **Taxids 0 and 1** (root, unclassified) keep `N_A`, as before.

## Order of precedence

Each taxid is evaluated in this order. The first step that applies decides.

1. **Tier 1, N_A clades.** Non-bacterial domains (Archaea, Eukaryota, Viruses, unclassified root) and the whole phylum Mycoplasmatota are N_A. For the Izimaplasma/Izemoplasmatales clade the wall-less envelope is shown by cryo-EM and by the isolate Hujiaoplasma nucleasis.
2. **Tier 3, species-level evidence from described isolates.** These are rules matched on the Species column or the Taxid, including published stains of a species (e.g. Pseudoruminococcus massiliensis). Strains inherit them through their Species column.
3. **Tier 4, no isolate.** "uncultured" in the name, Candidatus in the name (switchable), or `Proportion_of_MAG_in_taxid` = 1: undetermined. The name-based reasons are tested first only because they make the more informative label in the provenance file; all three give the same value.
4. **Tier 5, clade rules.** Phylum, Class, Order, Family or Genus; the deepest matching rank wins, so exceptions override the phylum default. For example, the Negativicutes rule (negative) overrides the Bacillota rule (positive), and the Oscillibacter stain rule (negative) overrides the Bacillota rule too. Each rule's `evidence_type` says whether it rests on envelope evidence, a published stain, or both. A genus-level stain rule comes from the description of the type species (and, usually, the genus), and so applies to the species of that genus that have isolates; species known only from MAGs, uncultured or Candidatus species remain undetermined at tier 4.
5. **Tier 6.** If nothing applies, the taxid is undetermined.
6. **Tier 7, bacterial rank nodes.** A genus, family, order, class, phylum, kingdom or domain row takes the consensus of its terminal descendants in the table (species and below): one value gives that value, several give `variable`, and none known gives undetermined. `N_A` counts as a known value here, so Bacteria as a whole is `variable`. Rows decided at tier 1 or tier 3 are never overridden. Intermediate NCBI ranks (subfamily, clade, species group, suborder, subclass) are not rank nodes, because their descendants cannot be resolved from the lineage columns; they keep the value from tiers 1 to 6. This is the same logic `agglomerate_features()` uses when it agglomerates features. Large mixed clades therefore come out `variable` (for example Bacteria, Pseudomonadota because of Anaplasmataceae and Sphingomonas, Bacillota and Clostridia because of the diderm lineages and the Gram-negative-staining genera, and Actinomycetota because of Mycobacteriales).
7. **Tier 8, one LKT, one Gram.** Several taxids can carry the same LKT; for example an NCBI "no rank" node placed directly under a phylum is labelled `LKT__p__<phylum>`. JAMS reads Gram onto features through the taxid embedded in the LKT name, so every row sharing an LKT takes the value of the taxid named in it. This changed 95 rows, and guarantees that, say, `LKT__c__909932_Negativicutes` has the value of taxid 909932 on every row of the table.

Unnamed "sp." isolates (names containing `_Unclassified_`) and named species without a genome in NCBI (`Proportion_of_MAG_in_taxid` = NA) receive their clade's value at tier 5. That value is an inference from the clade, and the provenance file says so.

## The rules in brief

The full detail, evidence and references are in `Gram_evidence_rules.tsv`.

**Negative (diderm with LPS/LOS):**

- Pseudomonadota
- Bacteroidota
- Campylobacterota
- Thermodesulfobacteriota (Desulfovibrio and relatives)
- Myxococcota
- Bdellovibrionota
- Fusobacteriota
- Chlamydiota
- Verrucomicrobiota
- Planctomycetota
- Cyanobacteriota
- Acidobacteriota
- Aquificota
- Chlorobiota
- Fibrobacterota
- The Spirochaetota families Leptospiraceae and Brachyspiraceae
- Within Bacillota: the classes Negativicutes and Limnochordia, and the order Halanaerobiales (diderm with LPS, and they also stain Gram-negative)
- Deltaproteobacteria that NCBI leaves without a phylum (class rule; lineage inference)

**Negative by published stain (monoderm, but described as Gram-stain-negative):**

- Bacillota: Pseudoruminococcus massiliensis, Oscillibacter, Pseudoflavonifractor, Subdoligranulum, Ruthenibacterium, Dysosmobacter, Vescimonas, Hydrogenoanaerobacterium, Anaeromassilibacillus, Negativibacillus
- Chloroflexota: the classes Anaerolineae, Caldilineae, Chloroflexia and Ardenticatenia, and the genus Dehalogenimonas
- Thermomicrobiota (secondary source; flagged)

**Variable by published stain (conflicting or Gram-variable reports):** Faecalibacterium (flagged for review), Gemmiger, Eisenbergiella, Acutalibacter, Intestinimonas, Sphaerobacter.

**Positive (monoderm):**

- Bacillota, except the three diderm lineages above
- Actinomycetota, except Mycobacteriales
- Chloroflexota by default, and Ktedonobacteria, which are also described as Gram-positive
- The six Corynebacterium species without mycolic acids: C. amycolatum, C. caspium, C. ciconiae, C. kroppenstedtii, C. lactis and C. otitidis
- The genus Hoyosella, which lacks mycolic acids (it includes the former Amycolicicoccus)

**Variable (outer membrane without LPS):**

- Borreliaceae
- Treponema
- Anaplasmataceae
- Sphingomonas
- Mycobacteriales (Lawsonella keeps this value: it produces mycolic acids)
- Deinococcota
- Thermotogota
- Synergistota

**Undetermined by rule:** Dictyoglomota, where published reports on lipid A genes conflict.

**N_A by rule:** Dehalococcoides, which has no peptidoglycan and no outer membrane.

**No rule, so undetermined:** about twenty smaller phyla (25 entries in the review file, counting p__Missing and p__Unclassified), including Nitrospirota, Gemmatimonadota, Lentisphaerota, Kiritimatiellota, Elusimicrobiota, Deferribacterota (which includes the gut genus Mucispirillum), Balneolota, Rhodothermota, Ignavibacteriota and Armatimonadota; the other Spirochaetota; and all Candidatus phyla. No evidence about their envelopes was read for this curation; Victivallis (Lentisphaerota) and Elusimicrobium are described as staining Gram-negative (R26), which under the envelope definition is not enough. They are listed in the review file with examples, so that rules can be added once evidence is in hand.

## How confident each rule is

Each rule records what was actually read (`verification` column):

- **direct:** the abstract, or an excerpt, of a primary paper stating the property.
- **secondary:** an encyclopaedic source summarising primary work.
- **lineage inference:** the clade was formerly part of a clade with documented properties.

Rules that are secondary, lineage-only or conflicting carry `review_flag = yes`. Full texts were not available while this was prepared, so supplementary tables were not read; for example, the per-phylum LPS gene survey in Sutcliffe 2010, Table S1. Checking those flagged rules against the full texts is the most useful next step.

## Limits

- **Clade rules are generalisations.** An undocumented exception inside a clade inherits the clade's value. The rule set is designed to grow: adding a row to `Gram_evidence_rules.tsv` and rerunning the script is all it takes.
- **`Proportion_of_MAG_in_taxid` is taken as given.** It reflects NCBI assembly metadata. Some strain-level taxids of well-cultured species are flagged MAG-only (for example Pseudomonas aeruginosa C3719). Some validly named species also have only MAG genomes in NCBI although a type strain exists. Both groups are listed in the review file.
- **Stain coverage is partial.** Published stains were looked up for the taxa where the previous hand curation and the envelope rules disagreed, and for a few others found along the way (about 30 genera or classes). Any other taxon whose description reports a stain contrary to its clade's envelope keeps the clade value until a rule is added. Clostridial cluster IV in particular contains Gram-stain-positive, -negative and -variable species, and a 2026 study (R49) reports monoderm lineages of Bacillaceae that stain Gram-negative; neither is yet covered comprehensively. A systematic pass with a strain-collection source such as BacDive, which records the stain of tens of thousands of type strains, is the most valuable next step, and the rule format is designed to take its results as additional rows.
- **Stain checks that failed.** No published stain could be found for Waltera, Massiliimalia, Thermoflexus, Tepidiforma, Thermobaculum and Phototrophicus, all of which the previous hand curation called negative. They keep the envelope value (positive) and are listed in the review file.

## Reproducing the curated table

```bash
python3 -I apply_Gram_curation.py \
    --taxtable JAMStaxtable_<date>.tsv \
    --rules Gram_evidence_rules.tsv \
    --out_table JAMStaxtable_<date>_Gram_curated.tsv \
    --out_provenance JAMStaxtable_Gram_provenance.tsv.gz \
    --out_review Gram_curation_review.tsv \
    --out_summary Gram_curation_summary.txt
```

The input TSV is exported from R with:

```r
write.table(JAMStaxtable, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
```

The output is identical to the input except for the Gram field; the script checks the layout and refuses malformed input.

To put the result back into the package without disturbing the class of any other column, replace only the Gram column of the existing object:

```r
load("data/JAMStaxtable.rda")
cur <- read.delim("JAMStaxtable_<date>_Gram_curated.tsv", colClasses = "character", quote = "", comment.char = "", na.strings = character(0))
stopifnot(nrow(cur) == nrow(JAMStaxtable), identical(as.character(JAMStaxtable$Taxid), cur$Taxid))
old_gram <- JAMStaxtable$Gram
JAMStaxtable$Gram <- if (is.factor(old_gram)) factor(cur$Gram) else cur$Gram
stopifnot(all(as.character(JAMStaxtable$Gram) %in% c("positive", "negative", "variable", "N_A", "undetermined")))
print(table(previous = old_gram, curated = JAMStaxtable$Gram, useNA = "ifany"))
save(JAMStaxtable, file = "data/JAMStaxtable.rda", compress = "xz")
```

Check the compression of the original first with `tools::checkRdaFiles("data")`, and use the same.

## When JAMStaxtable is rebuilt

`libexec/JAMSbuildk2db` builds JAMStaxtable from the NCBI taxonomy dump and creates no Gram column. After a rebuild, export the new table with any Gram column (an empty one will do), run `apply_Gram_curation.py` on it, and merge the result back as above. New taxids are then covered by the same rules automatically, and taxa with no evidence come out as undetermined.

## Effect on the 2026-11-23 table

On the 2026-11-23 table, 43,903 of 212,123 Gram values change:

| Value | Previous | Curated |
|---|---|---|
| N_A | 147,406 | 113,131 |
| negative | 37,218 | 48,928 |
| positive | 25,827 | 33,661 |
| undetermined | 1,672 | 10,355 |
| variable | 0 | 6,048 |

Most of the change is 25,492 bacterial taxa that had been N_A for lack of curation and now receive a positive or negative value (24,601 from their clade or stain rule, e.g. Agathobacter rectalis: positive; 891 rank nodes from the consensus of their descendants). A further 1,327 went from N_A to variable. Only 61 values go directly from negative to positive and 60 from positive to negative; each is listed in the review file with its rule. See `Gram_curation_summary.txt` and `Gram_curation_review.tsv` for the complete accounting.

Rule identifiers are not contiguous: C16 (a phylum-wide Spirochaetota rule) was withdrawn in favour of family and genus rules, M03 was folded into M02, and N01 (the Pseudoruminococcus note) and C30 (Thermomicrobiota, positive by lineage) were replaced by the stain rules S01 and S20. Rule prefixes: D = domain, M = Mycoplasmatota, T = species-level envelope evidence, C = clade envelope evidence, S = published stain, E = other envelope evidence, N = notes.
