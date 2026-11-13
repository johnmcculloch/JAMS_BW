---
title: "A Practical Guide to Downstream Analysis with JAMS2"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 2
    theme: cosmo
    highlight: textmate
---

# A Practical Guide to Downstream Analysis with JAMS2

*Working interactively in an R session produced by JAMSbeta*

This guide walks through how to explore, filter, and plot your data **after** JAMSbeta has run. It assumes you have a JAMSbeta output folder containing an `.RData` workspace (and/or the saved `..._Expvec_list_of_SummarizedExperiment_objects.rds` file) and that you can open that workspace on a laptop with the `JAMS` package installed. Everything below is meant to be run in a live R session.

> **Scope.** This is about analysis *downstream* of JAMSbeta. It does not cover running JAMSalpha or JAMSbeta themselves. Where a concept from the upstream pipeline is needed to understand a downstream choice, it is summarised briefly.

---

## Table of contents

1. [Prerequisites and opening the workspace](#1-prerequisites-and-opening-the-workspace)
2. [A worked example dataset to practise on](#2-a-worked-example-dataset-to-practise-on)
3. [Getting oriented: the key objects](#3-getting-oriented-the-key-objects)
4. [Metadata: `pheno`, `ctable`, and replacing metadata safely](#4-metadata-pheno-ctable-and-replacing-metadata-safely)
5. [Analysis spaces: taxonomic and functional](#5-analysis-spaces-taxonomic-and-functional)
6. [Core concepts: LKT, CGB, and CSB](#6-core-concepts-lkt-cgb-and-csb)
7. [The shared filtering vocabulary of JAMS functions](#7-the-shared-filtering-vocabulary-of-jams-functions)
8. [Choosing taxonomic filters with `titrate_feature_counts`](#8-choosing-taxonomic-filters-with-titrate_feature_counts)
9. [Arguments common to all the plotting functions](#9-arguments-common-to-all-the-plotting-functions)
10. [Ordination plots](#10-ordination-plots)
11. [Relative abundance heatmaps](#11-relative-abundance-heatmaps)
12. [Alpha diversity](#12-alpha-diversity)
13. [Putting it together with `compare_within_subsets`](#13-putting-it-together-with-compare_within_subsets)
14. [Worked examples on the MR50 dataset](#14-worked-examples-on-the-mr50-dataset)
15. [Troubleshooting](#15-troubleshooting)
16. [Quick reference](#16-quick-reference)

---

## 1. Prerequisites and opening the workspace

### What you need

- **R** (a recent version) on a laptop is fine for all downstream analysis; the heavy lifting was done upstream by JAMSalpha/JAMSbeta on the cluster.
- The **`JAMS` package** installed, along with its dependencies. Several downstream functions rely on Bioconductor packages (notably `SummarizedExperiment` and `ComplexHeatmap`). If they are missing, JAMS provides a helper:
  ```r
  library(JAMS)
  install_biocon_deps()   # installs dada2, ComplexHeatmap, genefilter, SummarizedExperiment, ...
  ```
- The **JAMSbeta output** for your project. This is either the `.RData` workspace image written into the JAMSbeta output folder, or the saved `..._Expvec_list_of_SummarizedExperiment_objects.rds` object, or both.

### Loading the package

```r
library(JAMS)
```

Loading `JAMS` also makes the tidyverse, ggplot2, SummarizedExperiment and the other plotting dependencies available, so you generally do not need to attach them yourself.

### Option A — open the JAMSbeta workspace image

If you have the `.RData` workspace that JAMSbeta saved, restoring it brings back `expvec`, `pheno`, `ctable`, `variable_list` and everything else exactly as JAMSbeta left them. The safe wrapper (which uses multiple cores via `fastSave` if available) is:

```r
# From within the JAMSbeta output folder:
IO_jams_workspace_image(workspaceimage = ".RData", operation = "load", verbose = TRUE)
```

You can also just use base R's `load(".RData")`; the wrapper is simply more informative and faster on large images.

### Option B — load only the experiment list

If all you kept is the saved experiment list, load it directly and reconstruct `pheno` from any object's `colData`:

```r
expvec <- readRDS("MyProject_Expvec_list_of_SummarizedExperiment_objects.rds")
pheno  <- as.data.frame(colData(expvec[[1]]), optional = TRUE)   # optional = TRUE keeps column names exactly as they are
```

### Sanity checks before you start

```r
names(expvec)                                   # which analysis spaces are present?
dim(expvec$ConsolidatedGenomeBin)               # features x samples
head(as.data.frame(colData(expvec$ConsolidatedGenomeBin), optional = TRUE))   # the metadata
metadata(expvec$ConsolidatedGenomeBin)$version  # which JAMS version built this?
```

Working directory matters, because every plot is written to a PDF in the current folder. Set it deliberately:

```r
setwd("~/path/to/where/I/want/my/figures")
```

---

## 2. A worked example dataset to practise on

Before touching your own data, it is well worth getting fluent on the small example dataset shipped with the package. It is a real mouse experiment with **stark, easy-to-see effects**, so it is ideal for learning what each argument does.

```r
library(JAMS)
data(expvec_MR50)

names(expvec_MR50)
# "Contig_LKT" "ConsolidatedGenomeBin" "ECNumber" "Product" "resfinder"

# Pull the metadata out of any object
pheno_MR50 <- as.data.frame(colData(expvec_MR50$ConsolidatedGenomeBin), optional = TRUE)
head(pheno_MR50)
```

### The MR50 experiment in one paragraph

Three groups of mice were followed over time under different antibiotic regimens: **G1** received no antibiotic, **G2** received vancomycin, and **G3** received vancomycin + neomycin. Faecal/caecal samples were collected at four timepoints — **Day 0** (baseline, pre-treatment), **Day 5**, **Day 15**, and **Day 56** (endpoint). Five mice per group, individually tracked, so you can look at group differences, time series within a mouse, and dose–response (0 vs 1 vs 2 antibiotics) all in one dataset.

The metadata columns you will use most:

| Column | Meaning | Good for |
|---|---|---|
| `Group` | G1 / G2 / G3 | the primary comparison |
| `Treatment` | `none` / `vancomycin` / `vancomycin+neomycin` | human-readable version of Group |
| `Timepoint` | `Day_0` / `Day_05` / `Day_15` / `Day_56` | subsetting; time series |
| `Timepoint_days` | numeric 0/5/15/56 | continuous time (correlations, ordering) |
| `Mouse_number` | individual mouse ID | pairing/connecting a subject over time |
| `Vancomycin_used`, `Neomycin_used` | Used / Not_Used | simple binary contrasts |

Because `expvec_MR50` already carries this metadata in its `colData`, you can use it immediately — no metadata replacement needed. (When you *do* need to replace metadata on your own project, see [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely).)

> **Note on the example object.** `expvec_MR50` still contains the deprecated `Contig_LKT` space so that older material keeps working. As with any project, ignore it and use `ConsolidatedGenomeBin` (see [§5](#5-analysis-spaces-taxonomic-and-functional)).

Fully worked, copy-pasteable analyses on this dataset are collected in [§14](#14-worked-examples-on-the-mr50-dataset). If you are new to JAMS, read [§3](#3-getting-oriented-the-key-objects)–[§13](#13-putting-it-together-with-compare_within_subsets) first, then run §14 to see it all in action.

---

## 3. Getting oriented: the key objects

When you load a JAMSbeta workspace, a handful of objects do almost all the work. You will spend your whole analysis manipulating these.

### `expvec`

`expvec` is a **named list of `SummarizedExperiment` objects**, one per analysis space. This is the heart of everything. Each element is a self-contained experiment: a counts matrix, a feature (row) table, the sample (column) metadata, and some extra bits tucked into the object's metadata.

```r
# What analysis spaces do I have?
names(expvec)
# e.g. "ConsolidatedGenomeBin" "Product" "ECNumber" "Pfam" "Interpro" "GO" ...

# Pull out a single experiment to work with
SEobj <- expvec[["ConsolidatedGenomeBin"]]
SEobj
```

Anatomy of one `SummarizedExperiment` (SEobj):

- **`assays(SEobj)`** — the numeric matrices, features (rows) × samples (columns). The primary one is always `BaseCounts` (sequencing depth in bases; see [§4](#4-core-concepts-lkt-cgb-and-csb)). Taxonomic objects also carry `GenomeCompleteness` and `GenomeContamination`; functional objects carry `GeneCounts` and `GeneLengths`.
- **`rowData(SEobj)`** — the feature table. For taxonomic spaces this holds the taxonomic lineage columns (Domain … Species, LKT, Gram, Taxid). For functional spaces it holds the accession and its description.
- **`colData(SEobj)`** — the sample metadata (this is your `pheno`, stored inside the object).
- **`metadata(SEobj)`** — everything else: `analysis` (the space name), `version`, `TotalBasesSequenced`, the colour table `ctable`, and, for taxon-stratified functional objects, the sparse feature-by-taxon matrices.

### `pheno`

`pheno` is the **metadata table** (a data frame): one row per sample, one column per variable, with row names equal to the sample names. It exists both as a standalone object in the workspace and as a copy inside every SEobj (`colData`). Keeping those two in sync is important — see [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely). Column names in `pheno` follow JAMS style (letters, digits and underscores, e.g. `IL_6_pg_mL`), because JAMSbeta cleans the metadata with `enforce_JAMS_style_metadata()` when it builds the objects; clean any curated replacement the same way (also §4).

### `ctable`

`ctable` is the **colour dictionary**: a table (columns `Name`, `Colour`, `Hex`) mapping metadata class names to fixed hex colours, so that (say) a given treatment arm is the same colour in every plot across your whole project. It lives in each SEobj's metadata and is consulted automatically by the plotting functions.

You build a `ctable` with `make_colour_dictionary()` and, if you wish, hand-tune it interactively with `tweak_ctable()` — both are covered in [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely). Two things worth knowing up front:

- A `ctable` holds **one colour per class name**, globally. A class such as `Pittsburgh` that appears in more than one variable (say `Cohort` and `Study`) therefore has a single colour shared across both.
- `make_colour_dictionary()` now generates a **distinct, harmonious palette within each variable** and actively **prevents two classes of the same variable from sharing a colour**, so what you get is plot-ready by default.

---

## 4. Metadata: `pheno`, `ctable`, and replacing metadata safely

### Why you often need to replace metadata

A common pattern: JAMSbeta is run with a **placeholder/dummy metadata** file (just the sample column and one variable) so the experiment objects can be built, and the **real, curated metadata** is swapped in afterwards. You may also want to add derived variables (e.g. dichotomising a continuous variable, adding a grouping column) after the fact. Either way, you must get the new metadata **into every SEobj in `expvec` without corrupting the object**.

> If you are just practising on `expvec_MR50`, its metadata is already curated and loaded, so you can skip straight to plotting — this section matters when you move to your own project.

### Step zero: clean the metadata with `enforce_JAMS_style_metadata()`

Badly named metadata columns are the single most common cause of baffling downstream errors. A space, hyphen, bracket, `%`, `/` or `+` in a column name looks harmless in a spreadsheet. But R quietly renames such columns in some places (e.g. `IL-6 (pg/mL)` becomes `IL.6..pg.mL.`) and not in others, so a variable you pass to `compareby` or `subsetby` can go missing halfway through a function. You then get errors like `undefined columns selected`, or a PERMANOVA formula that will not parse. Stray whitespace, blank cells and duplicated sample names cause similar grief.

JAMS has **one function that conforms a metadata table to JAMS style**: `enforce_JAMS_style_metadata()`. Crucially, **JAMSbeta and JAMS16 run this same function on your metadata file when they build the objects**, via `load_metadata_from_file()`. So a table you clean downstream ends up with exactly the same column names as the ones inside `expvec`. Always clean curated metadata this way before pushing it into the objects.

#### The easiest way: read the file the way JAMSbeta does

If your curated metadata is a file, let `load_metadata_from_file()` read *and* clean it. This is exactly the code path JAMSbeta uses, so the result is guaranteed to match:

```r
md    <- load_metadata_from_file(xlsxFile = "MyProject_curated_metadata.xlsx")
# or:  md <- load_metadata_from_file(phenotable_tsv = "MyProject_curated_metadata.tsv")
pheno <- md$phenotable
attr(pheno, "JAMS_column_name_changes")   # table of every column that was renamed
```

`md` also contains `phenolabels`, and `ctable` if the spreadsheet has a colour sheet.

#### Cleaning a data frame you already have in R

```r
pheno <- enforce_JAMS_style_metadata(pheno = pheno)
attr(pheno, "JAMS_column_name_changes")   # Original -> New, one row per renamed column
```

> **Watch how you read files yourself.** `read.table()`/`read.delim()` (default `check.names = TRUE`) and `openxlsx::read.xlsx()` (default `sep.names = "."`) both turn spaces into dots *before* `enforce_JAMS_style_metadata()` ever sees the names. A column read that way stays dotted (`IL.6..pg.mL.`) instead of becoming `IL_6_pg_mL`, and will not match the objects. Either use `load_metadata_from_file()` as above, or read TSVs with `read.delim(file, check.names = FALSE, colClasses = "character")`.

#### What it does

In order, `enforce_JAMS_style_metadata()`:

1. **Checks the input.** It must be a data frame; tibbles and data.tables are converted to a plain data frame.
2. **Makes column names safe:**
   - `%` becomes `pct` and `+` becomes `plus`.
   - Any other character that is not a letter, digit, `_` or `.` becomes `_`.
   - Repeated underscores are collapsed, and leading or trailing ones removed.
   - Names starting with a digit get an `X` in front, and R reserved words (`TRUE`, `NA`, `if`, …) get `_var` appended.
   - Names that end up identical get `_1`, `_2` added, with a warning.
3. **Finds the `Sample` column.** A column called `sample` or `SAMPLE` is renamed to `Sample`. If there is no such column it stops, and suggests likely candidates.
4. **Trims whitespace** from every cell.
5. **Enforces the JAMS convention:** every cell is character, and blank or `NA` cells become `N_A`.
6. **Vets the sample names:**
   - Completely empty rows are dropped quietly.
   - Rows without a sample name are dropped with a warning.
   - Duplicated sample names make it **stop**, listing them, because there is no safe way to choose between them.
   - Sample names themselves are **never changed**, since they must match the JAMS file prefixes; odd characters in them only trigger a warning.
7. **Flags the `§` character** in any cell. JAMS uses `§` internally to separate subset levels.
8. **Sets the row names to the sample names** (`rownames(pheno) == pheno$Sample`), as required for building and updating the SummarizedExperiment objects.

Some typical renames:

| Original column name | JAMS-style name |
|---|---|
| `Olink_IL_1 alpha_Serum` | `Olink_IL_1_alpha_Serum` |
| `IL-6 (pg/mL)` | `IL_6_pg_mL` |
| `Fat %` | `Fat_pct` |
| `CD4+ count` | `CD4plus_count` |
| `2nd_visit` | `X2nd_visit` |
| `sample` | `Sample` |

Every change is reported through the JAMS logger (one line per renamed column), so nothing happens silently.

Useful arguments:

- **`sample_column`** — name of the column holding the sample names. The default is `"Sample"`.
- **`class_to_ignore`** — the missing-data sentinel. The default is `"N_A"`.
- **`fix_column_names = FALSE`** — check-only mode: unsafe names are reported, with suggested fixes, but not changed.
- **`on_duplicate_samples = "drop"`** — drop every row of a duplicated sample instead of stopping. This is what JAMSbeta uses.
- **`phenolabels`** — if you still use a phenolabels table, pass it as well and its `Var_label` entries are renamed in step. The function then returns a list with `phenotable` and `phenolabels`.

#### What JAMSbeta and JAMS16 record when they clean your metadata

When JAMSbeta or JAMS16 read the metadata file, everything above is written to the project's **`<project>_JAMS.log`**:

- a summary of how many columns were renamed;
- one line per rename, e.g. `metadata column renamed: "IL-6 (pg/mL)" -> "IL_6_pg_mL"`;
- any dropped empty rows, unnamed rows and duplicated samples.

When anything was renamed, the full table of original and new names is also written to **`<project>_metadata_column_name_changes.tsv`** in the JAMSbeta output folder, and kept in the workspace as `opt$metadata_column_name_changes`. JAMSbeta always calls the sample column `Sample`, whatever it was called in the file.

Existing phenolabels files keep working. A label still matches its column whether it was written with the original name or with the dotted name produced by earlier JAMS versions.

> **"My column names changed!"** Projects built with earlier JAMS versions may have dotted names such as `IL.6..pg.mL.`, while newer runs produce `IL_6_pg_mL`. This is intentional: one set of rules for every way metadata enters JAMS. Look in the `_JAMS.log` or the `_metadata_column_name_changes.tsv` file for the exact mapping, and use the new names in `compareby`, `subsetby`, `colourby` and so on.

### The golden rule

Never assign to `colData()` piecemeal in a way that lets sample order drift. The `SummarizedExperiment` class will not let you casually swap in a reordered or differently-sized metadata table, and silently mismatching rows to columns is the worst kind of bug because nothing errors. JAMS provides a safe helper.

### `update_ExpObj_metadata()`

This function rebuilds the SEobj from scratch with the new phenotable, after checking that every sample in the new table exists in the object. Use it per analysis space.

```r
# 1. Read and clean the curated metadata the same way JAMSbeta does (see "Step zero" above)
pheno <- load_metadata_from_file(xlsxFile = "MyProject_curated_metadata.xlsx")$phenotable
#    ...or, for a data frame already in R:  pheno <- enforce_JAMS_style_metadata(pheno = pheno)

# 2. Nothing to do for row names: the cleaned pheno already has rownames == pheno$Sample,
#    which is how the objects match samples.

# 3. Rebuild a fresh colour dictionary from the curated metadata.
#    You do NOT need to pre-subset pheno by hand. Name the columns you want
#    coloured with columnsToUse (or the ones to skip with columnsToIgnore).
#    "Sample" is always protected: it is never coloured, whatever you pass.
#    Use the cleaned (JAMS-style) column names here.
ctable <- make_colour_dictionary(
    pheno         = pheno,
    columnsToUse  = c("Group", "Treatment", "Timepoint", "Response", "Sex"),
    legacy_ctable = "ctable_legacy.tsv",   # optional; see "Legacy colours" below
    object_to_return = "ctable"
)

# 4. Push the curated pheno (and ctable) into every experiment, safely
expvec_backup <- expvec          # always keep a backup
expvec <- list()
for (anal in names(expvec_backup)){
    # Subset the object's columns to the samples in pheno, in pheno's order,
    # THEN replace the metadata. Order is enforced by the [ , rownames(pheno)] step.
    expvec[[anal]] <- update_ExpObj_metadata(
        ExpObj    = expvec_backup[[anal]][ , rownames(pheno)],
        phenotable = pheno
    )
    # Re-attach the colour dictionary
    metadata(expvec[[anal]])$ctable <- ctable
}
```

Key points:

- **Clean first.** `enforce_JAMS_style_metadata()` (directly, or through `load_metadata_from_file()`) replaces the old pair of calls to `trim_whitespace_from_df()` and `replace_NAs_with_character()`. It does both, plus the column-name and sample-name checks.
- **Subset columns before replacing** (`expvec_backup[[anal]][ , rownames(pheno)]`). This guarantees the object's columns and your `pheno` rows are in the same order and the same set.
- `update_ExpObj_metadata()` **stops with an error** if `pheno` contains samples not present in the object — a feature, not a nuisance. It protects you from typos and stale sample names.
- Rebuild `ctable` **from the curated metadata**, then re-attach it in the loop, so colours are consistent across every space.
- The `N_A` token matters: it is the default `class_to_ignore` value throughout JAMS, so samples marked `N_A` in a comparison variable are automatically dropped from *that* comparison (see [§7](#7-the-shared-filtering-vocabulary-of-jams-functions)).
- **Getting `pheno` back out of an object?** Use `as.data.frame(colData(SEobj), optional = TRUE)`. Without `optional = TRUE`, `as.data.frame()` renames any column whose name is not a valid R name.

### What `make_colour_dictionary()` now does for you

The colour builder has been substantially improved. You mostly get good results with the call above, but it is worth knowing what it is doing, because the defaults are deliberate.

**Choose columns without pre-subsetting.** Use `columnsToUse` to name the columns to colour, or `columnsToIgnore` to name the ones to skip; you can pass the full `pheno`. `"Sample"` is force-protected — never coloured, and never counted — no matter what you pass to either argument. This replaces the old idiom of subsetting `pheno` down to the wanted columns.

**Per-variable, divergent palettes.** With `within_variable_palette = TRUE` (the default), each variable gets its *own* harmonious palette:

- A **two-class** variable gets a clean pair of contrasting colours (`binary_positive_colour` / `binary_negative_colour`).
- A **many-class** variable gets colours spread across a diverging ramp (`diverging_palette`, default `"Spectral"`).

It does not matter that, say, `Study` and `Cohort` end up with similar-looking palettes — they are independent variables and are never plotted on the same axis.

**A best-effort guess at binary polarity.** For two-class variables, `make_colour_dictionary()` tries to place the *affirmative/positive* class (Responder, Good, Yes, Y, CR, PR, Sensitive, …) on `binary_positive_colour` and the *negative* class (Non_Responder, Bad, No, N, PD, Resistant, Refractory, …) on `binary_negative_colour`. The defaults are:

- `binary_positive_colour = "#0004e9"` (blue) — Responder / Good / Yes
- `binary_negative_colour = "#E41A1C"` (red)  — Non_Responder / Bad / No

Class labels are normalised before matching (lower-cased, punctuation/spaces/underscores stripped), so `Non_Responder`, `Non-Responder` and `Non Responder` are treated identically. You can extend the vocabulary with `extra_positive_tokens` / `extra_negative_tokens`, or turn the guess off entirely with `guess_binary_polarity = FALSE` (which falls back to a stable alphabetical assignment).

> **Heads-up: some binary labels are valence-neutral.** `Yes`/`No` and `Y`/`N` mean "affirmative/negative", not "good/bad". For a column like `Death_event` or `Prog_event`, `Y` is the *bad* outcome, yet the guesser will still paint `Y` with the affirmative colour. This is an educated guess that is deliberately wrong-but-consistent in such cases. If it bothers you for a particular column, either drop that column from `columnsToUse`, or fix it afterwards in `tweak_ctable()` (below). Genuinely context-dependent pairs like `High`/`Low` are **not** guessed at all and fall through to alphabetical.

**Legacy colours (fixed project palette).** To keep a class the same colour across successive studies (e.g. `Responder` always `#0004e9`), pass a `legacy_ctable` — either a data frame or a path to a tsv. The file is column-name agnostic: it looks for a class column among `Name` / `Class_label` and a colour column among `Hex` / `Colour` / `Class_colour`, so your existing `ctable_legacy.tsv` works untouched. Wherever a class name matches, the legacy colour **overrides** the generated one. If a legacy file is supplied but nothing matches, you get a warning rather than silence.

> The old `colour_table` argument still works as a deprecated alias for `legacy_ctable`.

**Automatic collision avoidance.** With `avoid_collisions = TRUE` (the default), the builder guarantees that **within each variable, every class ends up with a distinct colour** in the final flat `ctable`. This fixes a subtle artefact of the flat, one-colour-per-name model: two different variables can independently land a class on the same colour (e.g. two 5-class variables both taking the `Spectral` midpoint `#FFFFBF`), and the global de-duplication could then leave two classes of one variable looking identical on a plot. The pass nudges the more "local" class (the one appearing in fewer variables) to a nearby free colour, and it **never moves a legacy-locked or ignore-class colour**. Intended cross-variable reuse — e.g. every binary-positive class sharing one blue — is preserved.

**Returning a per-variable dictionary instead.** `object_to_return = "cdict"` returns a named list with one colour table per variable (rather than a single collapsed `ctable`). Most workflows want the default `"ctable"`.

### Interactively fine-tuning colours with `tweak_ctable()`

Automatic generation gets you 95% of the way, but sometimes you want to hand-pick a shade. `tweak_ctable()` opens a small **Shiny** colour editor, lets you adjust any class with a colour wheel, live-reports any colour clashes, and returns the edited `ctable` straight back to your R session.

```r
# Edit the ctable, then keep the edited result
ctable <- tweak_ctable(ctable)

# Group the pickers by metadata variable, so you can see what plots together
ctable <- tweak_ctable(ctable, pheno = pheno)

# Group AND also save the result to disk on "Done"
ctable <- tweak_ctable(ctable, pheno = pheno, outfile = "ctable_curated.tsv")
```

After editing, re-attach the returned `ctable` to every object exactly as in step 4 above (the replacement loop), or simply:

```r
for (anal in names(expvec)) metadata(expvec[[anal]])$ctable <- ctable
```

What it does:

- **Live duplicate report.** The panel at the top flags any classes that currently share a hex. When you pass `pheno`, the check becomes **per-variable** (the clash that actually matters — two classes of one variable that would be indistinguishable on the same plot); without `pheno`, it reports global duplicates.
- **Grouping by variable (`pheno` supplied).** Each metadata variable becomes a titled box containing its classes, so you can see, at a glance, the set of colours that will appear together on a plot. Columns that are not colour classes (Sample, IDs, dates, continuous variables) are skipped automatically. Any `ctable` class not found in the metadata is collected into a trailing box so nothing is hidden.
- **Shared classes stay in sync.** Because a `ctable` stores one colour per class name, a class that appears in several variables (e.g. `Pittsburgh` in `Cohort` and `Study`) shows up in each of its boxes and is kept locked together — editing it in one box updates it everywhere. That is an honest preview of how the flat `ctable` will behave on your figures.
- **Optional file output.** Set `outfile` to also write the edited table to a tsv when you click Done.
- **Cancel is safe.** Clicking Cancel returns the original `ctable` unchanged, so it is fine to open speculatively.

> **Dependencies.** `tweak_ctable()` needs the `shiny` and `colourpicker` packages. These are **not** JAMS dependencies (the tool is used rarely), so they are not installed with JAMS. If either is missing, the function prints a short install hint and returns your `ctable` untouched — nothing breaks. Install them once with:
> ```r
> install.packages(c("shiny", "colourpicker"))
> ```

A good end-to-end recipe is therefore: **generate** with `make_colour_dictionary()` (per-variable palettes, legacy overrides, collision avoidance all on), then **hand-tune only the rare leftover** with `tweak_ctable(ctable, pheno = pheno)`, then **re-attach** the result across `expvec`.

### Adding a derived variable

If you add or edit a column, do it on `pheno` and then re-run the replacement loop above so the change propagates into `expvec`. Do not edit `colData` inside one object and expect the others to follow. If the new column introduces new classes, rebuild `ctable` (and, if you like, re-open `tweak_ctable()`) before re-attaching. Give new columns JAMS-style names (letters, digits and underscores), or simply run `enforce_JAMS_style_metadata()` again after adding them.

A common derived variable is a **dichotomised continuous variable**, e.g. a cytokine level split into High and Low. `dichotomize_variable_in_md()` does this in a JAMS-safe way:

```r
# Split at the median (the default)
pheno$IL_6_group <- dichotomize_variable_in_md(md = pheno, column_to_dichotomize = "IL_6_pg_mL")

# ...or at a fixed threshold, or with any function returning one number
pheno$IL_6_group <- dichotomize_variable_in_md(md = pheno, column_to_dichotomize = "IL_6_pg_mL", threshold = 5)
pheno$IL_6_group <- dichotomize_variable_in_md(md = pheno, column_to_dichotomize = "IL_6_pg_mL",
                                               func_for_threshold = function(x) quantile(x, 0.75))
table(pheno$IL_6_group)
```

How it handles the values:

- Values **at or above** the threshold become `"High"`, values below it `"Low"`. Ties therefore go to High, which matters for medians of variables with many repeated values.
- Numbers stored as text (the JAMS norm) are handled, and so are factors (read through their labels).
- `N_A` cells, and anything else in `class_to_ignore`, are left out of the threshold calculation and come back as `N_A`.
- Any other value that cannot be read as a number also comes back as `N_A`, with a warning naming it.
- `high_label` / `low_label` change the labels.

If you only need the split for an ordination, `plot_Ordination()` can do it on the fly instead; see [§10](#10-ordination-plots).

---

## 5. Analysis spaces: taxonomic and functional

Each element of `expvec` is one analysis "space". They fall into two families.

### Taxonomic spaces

These describe **who is there**. In JAMS2 the taxonomic space you should use is:

- **`ConsolidatedGenomeBin`** — the consolidated, quality-vetted taxonomic entities (see [§4](#4-core-concepts-lkt-cgb-and-csb)). **This is the taxonomic space to use.**

> **Deprecation notice — do not use `Contig_LKT`.**
> `expvec$Contig_LKT` is deprecated in JAMS2 and should not be used for analysis. If it is present in an older workspace, drop it:
> ```r
> expvec$Contig_LKT <- NULL
> ```
> Older objects may also contain an `MB2bin` space; it too has been phased out. Standardise on `ConsolidatedGenomeBin`.

A convenient feature of the taxonomic space is that you can **agglomerate up the ranks on the fly** with the `glomby` argument (e.g. `glomby = "Species"`, `glomby = "Genus"`), so you do not need a separate object per rank.

### Functional spaces

These describe **what the community can do**. Which ones exist depends on what annotation was available, but common spaces include:

- `Product` — gene product annotation (free-text descriptions)
- `ECNumber` — Enzyme Commission numbers
- `Pfam`, `Interpro`, `NCBIfam`, `CDD`, `PRINTS` — protein family / domain databases
- `GO` — Gene Ontology terms
- BLAST-based special spaces such as `resfinder` (antimicrobial resistance) and `vfdb` (virulence factors)

For functional spaces, the counts represent sequencing depth attributed to each functional feature, and (if the object was built with taxon stratification) you can ask *which taxa carry a given function* — see [§11](#11-relative-abundance-heatmaps) and the stratification note in [§10](#10-ordination-plots).

---

## 6. Core concepts: LKT, CGB, and CSB

These three acronyms recur everywhere in taxonomic analysis. Understanding them is the single biggest step toward using JAMS well.

### Abundance is sequencing depth, not read counts

Before the acronyms: throughout JAMS, **abundance means DNA sequencing depth** (the per-base coverage of assembled contigs, stored as `NumBases`), typically expressed as **PPM** (parts per million of the sample's total). It is *not* a raw read count. This is why relative abundances are directly comparable across samples with different sequencing effort.

### LKT — Last Known Taxon

An **LKT** is the *most resolved taxonomic label* that can confidently be attributed to something. It may sit at any rank — a well-resolved entity might be labelled at species or infraspecies level, while a poorly-resolved one might only reach family or phylum. LKT is the "leaf" label of a taxonomic feature. Tags follow the convention `d__` (domain), `p__` (phylum), `c__` (class), `o__` (order), `f__` (family), `g__` (genus), `s__` (species), `is1__` (infraspecies), and the terminal `LKT__`.

### CGB — Consolidated Genome Bin

A **CGB** is a per-sample taxonomic entity produced by JAMSalpha's consolidation step. In brief: high-quality metagenome-assembled genomes (MAGs) are kept as-is, lower-quality bins are enhanced where possible, and leftover contigs are grouped by taxonomy. The result is that as much sequencing depth as possible is attributed to defined taxonomic entities. Each row of the `ConsolidatedGenomeBin` object is (after cross-sample assembly) one such consolidated taxonomic feature. CGBs come with genome **completeness** and **contamination** estimates, which is what makes the `GenomeCompletenessCutoff` filter meaningful.

### CSB — Consolidated/Contextualized Species Bin

When JAMSbeta builds the taxonomic experiment, high-quality CGBs that share the same species are compared by their **functional gene content** and clustered. Each resulting functionally-coherent cluster is relabelled as a **CSB** (its feature name begins with `CSB_` / `CSB__` and carries an `FC` functional-cluster suffix). A CSB is, roughly, a strain-level grouping: it distinguishes genomically/functionally distinct populations of the same species.

The practical upshot:

- A **CSB is a feature that came from an actual MetaBAT2 bin** (a genuinely binned genome), as opposed to a feature assembled only from leftover, taxonomy-labelled contigs.
- This is exactly what the `only_allow_CSBs` argument keys on (see [§7](#7-the-shared-filtering-vocabulary-of-jams-functions)): setting it `TRUE` restricts a `ConsolidatedGenomeBin` analysis to *bin-derived* features and discards the taxonomy-only remainder.

### When to restrict to CSBs

- Use **`only_allow_CSBs = FALSE`** (the default, all taxa) when you want the most complete picture of the community, including sequence that was classified but not binned.
- Use **`only_allow_CSBs = TRUE`** when you want higher-confidence, genuinely-binned genomes only — e.g. when strain-level resolution or genome quality matters more than completeness of the census.

Running an analysis both ways and comparing is often informative.

---

## 7. The shared filtering vocabulary of JAMS functions

Almost every JAMS plotting function accepts the same filtering arguments, and they mean the same thing everywhere. Learn them once.

### `featcutoff` — abundance/prevalence filter

A **numeric vector of length 2**: `c(minPPM, percentOfSamples)`.

> Keep a feature only if it reaches at least `minPPM` parts-per-million **in at least** `percentOfSamples` % of the samples.

```r
featcutoff = c(250, 10)   # >= 250 PPM (0.025%) in >= 10% of samples
featcutoff = c(0, 0)      # no abundance filtering (the default)
```

This is a **prevalence-aware** abundance filter: it removes features that are either universally rare or present in only a handful of samples. The percentage is evaluated *within the sample set being plotted* (so within a subset, if you are subsetting).

### `GenomeCompletenessCutoff` — genome-quality filter (taxonomic only)

Also a **length-2 numeric vector**: `c(minCompletenessPercent, percentOfSamples)`.

> Keep a taxonomic feature only if its estimated genome completeness is at least `minCompletenessPercent` % **in at least** `percentOfSamples` % of samples.

```r
GenomeCompletenessCutoff = c(30, 10)   # >= 30% complete in >= 10% of samples
```

This only makes sense for **taxonomic** spaces (it relies on the `GenomeCompleteness` assay). If you pass it to a functional space it is ignored with a warning. It is a powerful denoiser: a taxon detected only as a sliver of a genome in every sample is usually noise.

### `applyfilters` — shorthand presets

Instead of setting the two cutoffs by hand you can pass `applyfilters = "light"`, `"moderate"`, or `"stringent"`. Each expands to sensible `featcutoff`, `GenomeCompletenessCutoff`, and `minl2fc` values, and the presets differ between taxonomic and functional objects. For taxonomic objects, roughly:

| preset | featcutoff | GenomeCompletenessCutoff | minl2fc |
|---|---|---|---|
| light | c(50, 5) | c(5, 5) | 1 |
| moderate | c(250, 15) | c(10, 5) | 1 |
| stringent | c(2000, 15) | c(30, 10) | 2 |

Any value you set explicitly **overrides** the preset for that one argument, so you can say `applyfilters = "moderate", GenomeCompletenessCutoff = c(20, 5)` and only that cutoff changes.

### `minl2fc` / `maxl2fc` — effect-size filters

Minimum (and maximum) absolute **log2 fold change** between the two groups of a binary comparison. Used to show only features whose difference is large enough to care about. Only relevant in comparative contexts where a fold change is computed.

### `PPM_normalize_to_bases_sequenced`

Controls the denominator when converting to PPM.

- `FALSE` (default in most plotting functions): for functional objects, normalise to the bases **in that analysis** (i.e. relative to coding regions annotated in that space).
- `TRUE`: normalise to the **total** bases sequenced in the sample (including non-coding / unannotated content).

For **taxonomic** objects the two are effectively the same, because the taxonomic census already spans the whole sample. For **functional** objects the choice changes interpretation: `TRUE` gives "fraction of the whole metagenome", `FALSE` gives "fraction of the annotated coding space". Pick one and be consistent within a project.

### `only_allow_CSBs`

As explained in [§6](#6-core-concepts-lkt-cgb-and-csb): restricts a `ConsolidatedGenomeBin` object to bin-derived (CSB) features only. Ignored (with a warning) for non-taxonomic spaces.

### `class_to_ignore`

Default `"N_A"`. Any sample whose value in the **comparison variable** equals this token is dropped from that comparison. This is how you exclude "not applicable / unknown" samples on a per-variable basis without deleting them from the object.

---

## 8. Choosing taxonomic filters with `titrate_feature_counts`

How do you pick good `featcutoff` and `GenomeCompletenessCutoff` values rather than guessing? **Titrate them.** `titrate_feature_counts` sweeps a grid of PPM thresholds × genome-completeness thresholds (at a fixed prevalence) and, for each cell, reports both how many features survive **and how much sequencing depth those survivors represent.** It is fast — no correlations, no permutations, just repeated feature counting — so it runs in well under a second on a laptop.

The depth accounting is the clever part. Cutting 90% of your features might cost only 5% of your sequencing depth (you removed low-abundance noise — good) or 50% of it (you removed real community — bad), and which one it is depends entirely on each sample's evenness. The function makes that visible directly.

```r
pdf("Taxonomic_feature_titration.pdf", paper = "a4r")

# Titrate at the Species level, all taxa
tit_species <- titrate_feature_counts(
    ExpObj                = expvec$ConsolidatedGenomeBin,
    glomby                = "Species",
    only_allow_CSBs       = FALSE,
    n_ppm                 = 10,     # number of PPM thresholds to sweep
    n_comp                = 10,     # number of completeness thresholds to sweep
    prevalence_pct        = 5,      # the "% of samples" leg of the cutoff, held fixed
    PPM_normalize_to_bases_sequenced = FALSE,
    ignoreunclassified    = TRUE,
    class_to_ignore       = "N_A",
    return_plots          = TRUE,
    verbose               = TRUE
)
print(tit_species$plots)      # prints all four titration panels

# Titrate restricted to CSBs (bin-derived genomes only)
tit_csb <- titrate_feature_counts(
    ExpObj                = expvec$ConsolidatedGenomeBin,
    glomby                = NULL,
    only_allow_CSBs       = TRUE,
    n_ppm                 = 10,
    n_comp                = 10,
    prevalence_pct        = 5,
    return_plots          = TRUE,
    verbose               = TRUE
)
print(tit_csb$plots)

dev.off()
```

### What it returns

With `return_plots = TRUE` you get a list with two elements:

- **`$grid`** — a data frame with one row per (PPM, completeness) cell, carrying `n_features`, `pct_retained`, `pct_depth_median`, `pct_depth_min`, and a local `decline_rate`.
- **`$plots`** — a list of four ggplot tile-plots (so these *can* be captured and re-arranged, unlike the heatmap):
  1. **`feature_count`** — surviving feature count (and % of starting features) per cell.
  2. **`decline_rate`** — the local rate of feature loss per grid step; pale = a stable plateau, dark = a steep face. Turning points (e.g. a recurring ~10% completeness floor) show up here as a drop-off in the rate.
  3. **`depth_retained`** — the % of sequencing depth the survivors represent (median across samples, with the worst-hit sample's value in parentheses).
  4. **`count_vs_depth_divergence`** — depth% minus feature%: **blue** means you axed mostly noise (kept far more depth than features — good), near-**white** means features and depth fell together, **red** means you lost more depth than features (axed real community — bad).

### Useful arguments

- **`glomby`** — titrate at the level you will actually analyse (`"Species"`, `"Genus"`, or `NULL` for the native CGB level).
- **`only_allow_CSBs`** — titrate the CSB-restricted view separately; its curve usually looks quite different.
- **`n_ppm` / `n_comp`** — grid resolution; the PPM ladder is auto-built with roughly log spacing so the informative low-PPM region is well sampled. You can supply explicit `PPM_thresholds` / `completeness_thresholds` vectors instead.
- **`prevalence_pct`** — the fixed "% of samples" leg, held constant across the whole grid. For small cohorts consider **`prevalence_min_samples`** instead: an absolute sample count that overrides the percentage (setting it to `1` recovers "present in at least one sample", handy for spotting single-sample artefacts while still censoring via the completeness axis).

### How to choose

Look for where tightening the threshold stops buying you much: on the **decline-rate** panel, the point past which the surface goes pale; on the **divergence** panel, the most aggressive cell that is still comfortably blue (and whose depth **min** is not alarmingly low). That combination — kills noise, keeps depth, keeps the worst sample intact — is your denoising cutoff. Record your decisions so every downstream plot uses the same values:

```r
filtration_thresholds <- list(
    denoise  = list(featcutoff = c(25, 5),  GenomeCompletenessCutoff = c(10, 5)),
    moderate = list(featcutoff = c(25, 5),  GenomeCompletenessCutoff = c(20, 5))
)
```

Remember these cutoffs apply across the **whole** sample set of the object (or within each subset when you subset). Titrate on the data you will actually analyse.

---

## 9. Arguments common to all the plotting functions

The three main plotting functions (`plot_Ordination`, `plot_relabund_heatmap`, `plot_alpha_diversity`) and the comparator share a large common vocabulary. Here is the core set.

### Selecting samples and features

- **`samplesToKeep`** — vector of sample names to include. `NULL` = all samples.
- **`featuresToKeep`** — vector of feature names to include. `NULL` = all features. When you also agglomerate (`glomby`), the names here must be **post-agglomeration** names (e.g. family names like `f__Enterobacteriaceae` if `glomby = "Family"`). *Names that are mostly absent from the object are treated as a likely namespace mistake: a high miss rate warns, a very high one errors, so a wrong or stale feature list is never plotted silently.*

### Grouping and comparing

- **`compareby`** — the metadata variable that defines the groups you are testing (e.g. treatment vs control). This drives p-values, fold changes, PERMANOVA, etc.
- **`subsetby`** — a metadata variable (or up to three, for tiered subsetting) that splits the data into separate plots, one per class. Analysis then runs **within** each subset. For example `subsetby = "Timepoint"` produces one plot per timepoint; `subsetby = c("Arm", "Timepoint")` produces one plot per Arm×Timepoint combination.

The difference is important: `compareby` compares groups **within** a plot; `subsetby` makes **separate** plots.

### Agglomeration

- **`glomby`** — for taxonomic objects, the rank to roll features up to (`"Genus"`, `"Species"`, …). `NULL` = no agglomeration (use features as-is). Terminal levels like `LKT` need no agglomeration.

### Cosmetic / annotation mapping

- **`colourby`** (US spelling `colorby` accepted in ordination) — metadata variable mapped to colour.
- **`shapeby`** — metadata variable mapped to point shape.
- **`cdict`** — a colour dictionary to enforce specific colours. If `NULL`, the `ctable` stored in the object is used automatically.
- **`addtit`** — extra text appended to the plot title.

### Normalisation and filtering

All the filtering arguments from [§7](#7-the-shared-filtering-vocabulary-of-jams-functions) (`applyfilters`, `featcutoff`, `GenomeCompletenessCutoff`, `only_allow_CSBs`, `PPM_normalize_to_bases_sequenced`, `class_to_ignore`) appear here too, with identical meaning.

---

## 10. Ordination plots

`plot_Ordination` answers **"are the communities different overall?"** It reduces the full sample-by-feature matrix to two dimensions and (optionally) tests group separation with PERMANOVA.

### Minimal call

```r
ord <- plot_Ordination(
    ExpObj    = expvec$ConsolidatedGenomeBin,
    glomby    = "Species",
    algorithm = "tUMAP",
    compareby = "Group",
    colourby  = "Group",
    applyfilters = "light"
)
# plot_Ordination returns a LIST of ggplot objects
print(ord[[1]])
```

Unlike the heatmap, **ordination returns ggplot objects in a list**, so you can capture, further edit, or arrange them.

### Most-used arguments

- **`algorithm`** — `"tUMAP"` (default), `"PCoA"`, `"PCA"`, or `"tSNE"`. tUMAP and tSNE emphasise local structure; PCoA/PCA give interpretable axes with % variance and pair naturally with a dissimilarity index.
- **`distmethod`** — dissimilarity index for PCoA/PCA/plot distances (default `"bray"` for Bray–Curtis).
- **`compareby`** — group variable for PERMANOVA and centroids.
- **`colourby` / `shapeby`** — map groups to colour and shape. Very often you set `colourby = compareby`.
- **`permanova`** (default `TRUE`) and **`permanova_permutations`** — whether to compute PERMANOVA and how many permutations.
- **`plotcentroids`** — draw a centroid per group with spokes to each sample; **`show_centroid_distances`** adds a between-centroid distance table.
- **`ellipseby`** — draw normal-data ellipses around a grouping.
- **`textby`**, **`connectby`**, **`sizeby`** — label points, connect points sharing a class (e.g. a subject across timepoints), or size by a variable.

### Highlighting a subset of samples

Two related mechanisms let you show *where* a group of samples sits relative to the whole cohort:

- **`samplesToHighlight`** — a vector of sample names to draw at full opacity while all others fade into the background. The ordination is still computed on the full set (or on `samplesToKeep`), so positions are preserved; only opacity changes. `unhighlighted_transparency_percentage` controls how faded the rest are.

```r
STH <- rownames(subset(as.data.frame(colData(expvec$ConsolidatedGenomeBin), optional = TRUE),
                       Region == "North"))
plot_Ordination(ExpObj = expvec$ConsolidatedGenomeBin, glomby = "Species",
                samplesToHighlight = STH, shapeby = "Region",
                permanova = FALSE, applyfilters = "light")
```

- **`highlight_subset_in_context = TRUE`** (used together with `subsetby`) — for each subset class, the **entire cohort** is plotted as a faded backdrop and only that subset's samples are highlighted, while **PERMANOVA and centroids are computed on the highlighted samples only**. This keeps every subset panel in the same shared coordinate space, so panels are directly comparable, and preserves each subset's position relative to the whole. When this mode is on, an explicitly passed `samplesToHighlight` is ignored (with a warning), because the highlight channel is being used by the subset mechanism.

### Continuous `compareby`: dichotomising on the fly

PERMANOVA, centroids and group colours need **categorical** groups. If you pass a continuous variable to `compareby` (all of its values are numbers, more than two distinct ones), `plot_Ordination` will not guess what you want. By default it issues a warning explaining the options and **returns an empty list without plotting**. To compare High against Low groups instead, set `dichotomize_continuous_compareby_threshold`:

```r
ord <- plot_Ordination(
    ExpObj    = expvec$ConsolidatedGenomeBin,
    glomby    = "Species",
    compareby = "IL_6_pg_mL",
    colourby  = "IL_6_pg_mL",
    dichotomize_continuous_compareby_threshold = "median",   # or a number, e.g. 5, or function(x) quantile(x, 0.75)
    applyfilters = "light"
)
```

How the split works:

- **What you can pass:** a number, the name of a function (`"median"`, `"mean"`), or your own function returning one number.
- **One threshold for the whole cohort.** It is computed once, across all samples being plotted, before any `subsetby` split, so High and Low mean the same thing in every subset panel. The threshold and the number of samples in each group are printed in the plot title.
- **A new column for the split.** The High/Low labels go into a new variable called `<compareby>_dichotomized`, which is used for PERMANOVA and centroids. If `colourby`, `shapeby` or `ellipseby` were the same variable as `compareby`, they switch to the new column too; `textby`, `sizeby` and `connectby` keep the original numbers. Your objects are not modified.
- **Colours.** High and Low take their colours from `cdict` or the object's `ctable` when **both** labels are defined there. Otherwise they always get the same fixed pair: red for High, blue for Low.
- **Numerically coded groups.** If the variable really is categorical but coded with numbers (e.g. timepoints 0, 7, 14, 28), pass `dichotomize_continuous_compareby_threshold = FALSE` to treat each value as its own group.
- **Inherited `compareby`.** If you did not pass `compareby` and it was taken over from `colourby`, `shapeby` and so on, a continuous variable is plotted as before, without refusing.

For a split you want to keep and reuse across functions, add it to `pheno` with `dichotomize_variable_in_md()` instead ([§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely)).

### Loops that do not crash

When something about the call cannot work, `plot_Ordination` issues a warning naming the problem and **returns an empty list** instead of stopping. This is what you want when looping over many variables. It happens when:

- a variable named in the call (`compareby`, `subsetby`, `colourby`, `shapeby`, …) does not exist in the object;
- the `compareby` variable is entirely `N_A` for the samples being plotted;
- fewer than 3 samples remain;
- a continuous `compareby` was given without a dichotomisation threshold (see above).

Check `length(ord) == 0` to detect a skipped variable.

### A note on interpretation

Points close together have similar community composition (under the chosen dissimilarity). A significant PERMANOVA p-value means the groups differ in composition; the F-model statistic and the visual separation tell you how strongly. Remember tUMAP/tSNE distances between *clusters* are not quantitatively meaningful the way PCoA axis distances are — use PCoA/PCA when you need axes with interpretable variance.

---

## 11. Relative abundance heatmaps

`plot_relabund_heatmap` answers **"which features differ, and by how much?"** It shows features (rows) across samples (columns) as a colour-coded relative-abundance matrix, optionally with per-feature statistics.

### The single most important caveat

> **`plot_relabund_heatmap` is NOT a ggplot2 function.** It draws directly to the active graphics device (via the ComplexHeatmap package) as a side effect. **You cannot capture its plot into an object** the way you can with ordination or alpha diversity. Therefore you must wrap calls in your own `pdf(...)` / `dev.off()`, and the return value (when `returnstats = TRUE`) is a **stats list**, not a plot.

```r
pdf("Heatmaps.pdf", paper = "a4r")
plot_relabund_heatmap(
    ExpObj    = expvec$ConsolidatedGenomeBin,
    glomby    = "Species",
    hmtype    = "comparative",
    compareby = "Group",
    applyfilters = "moderate",
    showonlypbelow = 0.05
)
dev.off()
```

### `hmtype`: exploratory vs comparative

- **`"exploratory"`** — features ranked by **variance** across all samples. Good for "what varies most?" with no group hypothesis. Cap the number shown with `ntop`.
- **`"comparative"`** — requires `compareby`. Computes per-feature statistics between groups: Mann–Whitney–Wilcoxon (2 classes) or ANOVA (3+), plus log2 fold change for binary comparisons. Good for "what differs between my groups?"

### Most-used arguments

- **`colcategories`** — which metadata variables to display as the coloured **annotation bars along the top** of the heatmap (the sample legend strips). Choose the variables that help you read the plot (e.g. group, timepoint, sex). If `NULL`, JAMS auto-includes variables with 2–10 classes. Setting this deliberately keeps the header readable.
- **`showonlypbelow`** — show only features with a p-value below this threshold (comparative heatmaps). The workhorse for keeping a comparative heatmap to a digestible size. Pair with `adj_pval_for_threshold = TRUE` to threshold on the **FDR-adjusted** p-value instead of the raw one.
- **`splitcolsby`** — split the columns into blocks by a metadata variable, so samples are grouped visually (e.g. `splitcolsby = "Group"` puts all controls together, all treated together). Often you set this to your `compareby`.
- **`scaled`** — if `TRUE`, plot per-feature **z-scores** (row-scaled) instead of absolute relative abundance. Use when you care about *pattern* across samples rather than absolute level; a low-abundance feature and a high-abundance one become visually comparable.
- **`minl2fc` / `maxl2fc`** — restrict to features with a large (or bounded) log2 fold change between the two groups.
- **`ntop`** and **`max_rows_in_heatmap`** — cap total features and features-per-heatmap page; **`maxnumheatmaps`** caps how many pages are produced.
- **`show_GenomeCompleteness`** (taxonomic) — add a companion genome-completeness heatmap so you can judge how trustworthy each taxon's abundance is.
- **`showGram`** (taxonomic) — annotate rows with Phylum and predicted Gram status.
- **`ignoreunclassified`** (default `TRUE`) — hide the fully-unclassified/unannotated feature from display (statistics still account for it, so p-values do not change).
- **`no_underscores`** — strip underscores from taxonomic names for readability.

### How to read the heatmap

- **Rows** are features (taxa or functions); **columns** are samples.
- **Cell colour** encodes relative abundance. For taxonomic objects the scale is a fixed percentage spectrum (so colours mean the same thing across heatmaps); the legend gives the breakpoints. When `scaled = TRUE`, colour instead encodes a z-score centred at zero (blue below the feature's mean, red above).
- **Row-label colour** in comparative binary heatmaps encodes direction: by default one colour marks features increased in one group, another colour the other group. The title states the directionality (e.g. "positive l2fc means increased in *X*"), which is affected by `invertbinaryorder`.
- **Right-hand annotations** show the p-value (and FDR-adjusted p-value) and the log2 fold change per feature, when applicable.
- **Top annotation bars** are your `colcategories` — the sample metadata, colour-keyed.
- If columns are split (`splitcolsby`), each block is one class of that variable.
- The title also reports the filter settings and the number of features assessed, so a heatmap is self-documenting.

If a comparative heatmap comes back empty, it usually means nothing survived the combined p-value / fold-change / abundance filters — loosen `showonlypbelow`, `minl2fc`, or the feature cutoffs.

---

## 12. Alpha diversity

`plot_alpha_diversity` answers **"is within-sample richness/evenness different between groups?"** It computes standard diversity measures per sample and compares them across groups with box plots and tests.

```r
pdf("AlphaDiversity.pdf", paper = "a4r")
adiv <- plot_alpha_diversity(
    ExpObj    = expvec$ConsolidatedGenomeBin,
    measures  = c("Observed", "Shannon", "InvSimpson"),
    compareby = "Group",
    colourby  = "Group",
    fillby    = "Group",
    applyfilters = "light"          # NOTE: do not leave this unfiltered — see below
)
dev.off()
```

### The critical point: never run alpha diversity unfiltered

> **Alpha diversity is extremely sensitive to noise.** An unfiltered object contains a long tail of spurious, ultra-low-abundance, sliver-of-a-genome features. Richness measures (Observed, Chao1) count those phantom features and **massively over-inflate diversity**, and even evenness-based measures are distorted. **Always apply a sensible filter** before estimating alpha diversity.

For this reason `plot_alpha_diversity` **defaults `applyfilters` to `"light"`** rather than to no filtering — a deliberate guard rail. Use at least `"light"`, and prefer cutoffs you have justified by titration ([§8](#8-choosing-taxonomic-filters-with-titrate_feature_counts)). Do **not** pass `featcutoff = c(0,0)` / `GenomeCompletenessCutoff = c(0,0)` for a real diversity analysis.

### Most-used arguments

- **`measures`** — which indices to compute. Options include `"Observed"`, `"Chao1"`, `"Shannon"`, `"Simpson"`, `"InvSimpson"`, and `"GeneCount"` (number of distinct genes; functional flavour of richness).
- **`compareby`** — group variable; drives the Wilcoxon (2 groups) or the pairwise comparisons drawn on the boxes.
- **`stratify_by_domains`** (default `TRUE`, taxonomic) — produce a separate diversity panel per Domain (Bacteria, Archaea, …) as well, provided each domain has enough features. Useful when a shift is driven by one domain.
- **`colourby` / `shapeby` / `fillby`** — cosmetic mappings; commonly all set to `compareby`.
- **`pairby`** — pair samples for a **paired** Wilcoxon test (e.g. pre/post on the same subject).
- **`connectby`** — draw a line linking a subject's samples across the boxes (e.g. one line per participant across timepoints); jitter is disabled so lines meet their points.
- **`applyfilters` / `featcutoff` / `GenomeCompletenessCutoff`** — as always; keep them meaningful (see the warning above).
- **`returnstats`** — if `TRUE`, the returned list also contains the per-sample measure tables, handy for your own downstream stats or export.

### Reading it

Each box summarises the distribution of a diversity measure within a group; the printed p-value tests whether groups differ. Higher Observed/Chao1 = more distinct taxa; higher Shannon/InvSimpson = more even, diverse communities. Because the y-axis meaning depends on the filter, always state your filter in the caption (the plot title records it for you).

---

## 13. Putting it together with `compare_within_subsets`

Once the individual functions make sense, `compare_within_subsets` lets you run the **whole standard trio — ordination, heatmap, alpha diversity — for one analysis space, organised subset by subset, in a single call.** For each subset tier it emits, in order:

1. **Ordination** — *is the community different overall?*
2. **Heatmap** — *which features drive/differ?*
3. **Alpha diversity** — *is richness different?*

keeping one subset's complete story together before moving to the next.

### It does not open a device

Because the heatmap draws directly to the active device (see [§11](#11-relative-abundance-heatmaps)), `compare_within_subsets` **does not open a graphics device itself**. You must wrap it in your own `pdf(...)` / `dev.off()`. Ordination and alpha ggplots are printed on the fly so everything lands in your open PDF in subset order.

### A fixed, global feature set across all subsets

This is the single most important thing to understand about how the comparator filters.

The abundance and completeness cutoffs (`featcutoff`, `GenomeCompletenessCutoff`) are **prevalence-aware**: they keep a feature only if it clears a threshold *in at least X % of samples*. That "% of samples" leg depends on how many samples are in the set being filtered. So if each subset were filtered independently — fewer samples per subset, different composition — the set of surviving features would **drift from one subset panel to the next**, and the panels would no longer be comparable.

To prevent this, `compare_within_subsets` establishes a **single global feature set once, up front**, across all retained samples, using your filtration arguments. It then holds that set fixed for every subset and every plot type, handing it to each plotting function with the plotting functions' own filtration pinned to non-filtering values. Every panel, in every subset, therefore considers exactly the same features. (This mirrors the global-feature machinery already built into `plot_relabund_heatmap` and `plot_Ordination` when they subset internally — the comparator subsets *outside* those functions, so it reproduces the same guarantee itself.)

The global feature set is computed in the same vetted order the plotting functions use — any up-front sample/feature cull first, then CSB restriction, then agglomeration — so the feature names line up with what each plot actually sees.

### Up-front sample and feature culling

Two arguments let you remove things you never want considered, **before** the global feature set is computed or any subsetting happens:

- **`samplesToKeep`** — sample names to retain up front. Use this to drop known outliers or samples irrelevant to the question. The global feature set and all subsetting are scoped to what remains.
- **`featuresToKeep`** — feature names to retain up front. Use this to drop noise, or to focus on a shortlist of features of interest, before global filtration runs.

> **`featuresToKeep` is in the POST-agglomeration namespace — identical to the standalone plotting functions.** If you also pass `glomby`, give `featuresToKeep` the **agglomerated-level** names, i.e. the names you would see in `rownames()` *after* agglomerating to that level. For example, with `glomby = "Species"` pass species names (`s__165179_Segatella_copri`); with `glomby = "Family"` pass family names (`f__Enterobacteriaceae`). **Do not** pass the underlying LKT/CGB names that compose those groups. (This matches the `featuresToKeep` contract in [§9](#9-arguments-common-to-all-the-plotting-functions); the two are now the same, so there is nothing special to remember for the comparator.)

> **Wrong namespaces and stale lists are caught, not silently honoured.** Because `featuresToKeep` is matched *after* agglomeration, a vector given in the wrong namespace (e.g. LKT names with `glomby = "Species"`) will mostly fail to match. JAMS now judges this by how many of your requested names are actually present, scaled by how many you asked for:
> - **0 present** → stops with an error (nothing to plot; clearest sign of a bad vector).
> - **a small request** (≈10 or fewer names) → never errors on absence; it is perfectly fine to ask for 4 species and have 2 be absent — you get an informative note and a plot of the rest.
> - **a large request with a high miss rate** → a strong warning above ~20% absent (it still plots), and an **error** above ~50% absent (almost always a wrong namespace or a stale list).
>
> The thresholds are deliberate internal defaults (not user arguments): the aim is simply to guarantee that a wrong or degenerate feature set never lands *silently* in a PDF. If you are knowingly applying a large curated list to a dataset that shares few of its features, expect the error and either trim the list to this dataset or run that case through `plot_relabund_heatmap` directly.

### Anatomy of a call

```r
# Per-function argument bundles (override the comparator's shared defaults)
hm_args  <- list(showonlypbelow = 0.05, adj_pval_for_threshold = FALSE,
                 minl2fc = 1, maxnumheatmaps = 2, max_rows_in_heatmap = 80,
                 show_GenomeCompleteness = TRUE, no_underscores = TRUE,
                 colcategories = c("Group", "Timepoint"))

ord_args <- list(algorithm = "tUMAP", plotcentroids = TRUE,
                 highlight_centroids = TRUE, grid = FALSE, forceaspectratio = 1)

alp_args <- list(measures = c("Observed", "Shannon", "InvSimpson"),
                 stratify_by_domains = TRUE)

pdf("Comparisons_Species.pdf", paper = "a4r")
compare_within_subsets(
    ExpObj    = expvec$ConsolidatedGenomeBin,
    compareby = "Group",
    subsetby  = c("Arm", "Timepoint"),     # tiered subsets: one story per Arm x Timepoint
    samplesToKeep = NULL,                   # optional up-front sample cull (e.g. drop outliers)
    featuresToKeep = NULL,                  # optional up-front feature cull (POST-glom names!)
    glomby    = "Species",
    only_allow_CSBs = FALSE,
    do_ordination = TRUE,
    ordination_highlight_subset_in_context = TRUE,
    do_heatmap = TRUE,
    do_alpha   = TRUE,
    ordination_args = ord_args,
    heatmap_args    = hm_args,
    alpha_args      = alp_args,
    featcutoff = c(25, 5),
    GenomeCompletenessCutoff = c(10, 5),
    PPM_normalize_to_bases_sequenced = TRUE,
    class_to_ignore = "N_A",
    return_stats = FALSE
)
dev.off()
```

Want to restrict to a handful of species up front? Because `featuresToKeep` is post-glom, you name the species directly — no need to look up their constituent LKTs:

```r
pdf("Comparisons_fewspecies.pdf", paper = "a4r")
compare_within_subsets(
    ExpObj    = expvec$ConsolidatedGenomeBin,
    compareby = "Group",
    subsetby  = "Timepoint",
    glomby    = "Species",
    featuresToKeep = c("s__853_Faecalibacterium_prausnitzii",
                       "s__165179_Segatella_copri",
                       "s__562_Escherichia_coli"),   # species names, matched post-glom
    featcutoff = c(25, 5), GenomeCompletenessCutoff = c(10, 5),
    PPM_normalize_to_bases_sequenced = TRUE,
    return_stats = FALSE
)
dev.off()
```

### How its arguments work

- **Shared filtration arguments** (`applyfilters`, `featcutoff`, `GenomeCompletenessCutoff`) are used **once**, up front, to derive the global feature set (see above). They are **not** re-applied inside each subset; the plotting functions receive the fixed feature set and non-filtering values instead. This is what guarantees the same features across every panel.
- **Other shared arguments** (`PPM_normalize_to_bases_sequenced`, `glomby`, `only_allow_CSBs`, `cdict`, `class_to_ignore`) are forwarded to **all three** plotting functions, so normalisation and agglomeration are consistent across the trio.
- **`samplesToKeep` / `featuresToKeep`** are applied **first**, before global filtration and subsetting, as described above (`featuresToKeep` in the post-glom namespace).
- **`ordination_args` / `heatmap_args` / `alpha_args`** are named lists that **override or extend** the shared defaults for one function only. This is how you set, say, a heatmap-specific `showonlypbelow` or an ordination-specific `algorithm` without affecting the others.
- **`do_ordination` / `do_heatmap` / `do_alpha`** switch each plot type on or off.
- **`ordination_highlight_subset_in_context`** turns on the in-context highlighting described in [§10](#10-ordination-plots): the ordination is computed on all (retained) samples and each subset is highlighted against the faded cohort, with PERMANOVA/centroids on the highlighted samples only.
- **`split_heatmap_columns_by_compareby`** (default `TRUE`) automatically sets the heatmap's `splitcolsby` to your `compareby`.
- **`return_stats`** — if `TRUE`, the function invisibly returns a collected list of the heatmap and alpha stats per subset (plots still go to the device).

### Typical workflow

1. Clean the curated metadata with `enforce_JAMS_style_metadata()`, swap it in and rebuild `ctable` ([§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely)).
2. Drop deprecated spaces (`expvec$Contig_LKT <- NULL`) ([§5](#5-analysis-spaces-taxonomic-and-functional)).
3. Titrate to pick filters ([§8](#8-choosing-taxonomic-filters-with-titrate_feature_counts)); store them in a list.
4. Explore with individual functions to sanity-check arguments ([§10](#10-ordination-plots)–[§12](#12-alpha-diversity)).
5. Run `compare_within_subsets` to produce the full comparison PDF(s), once with `only_allow_CSBs = FALSE` and once `TRUE` if strain-level confidence matters. Optionally drop known outlier samples or noise features up front with `samplesToKeep` / `featuresToKeep` (post-glom names).
6. Save the workspace (`IO_jams_workspace_image(...)`) so the state is reproducible.

---

## 14. Worked examples on the MR50 dataset

Everything below runs on the bundled example data with no metadata setup — load it and go. Each block is self-contained (opens its own PDF), so you can run them one at a time, inspect the figure, and keep the ones you want for the tutorial. The MR50 effects are large, so these make clear teaching figures.

```r
library(JAMS)
data(expvec_MR50)

# Use the ConsolidatedGenomeBin space; ignore the deprecated Contig_LKT.
SE <- expvec_MR50$ConsolidatedGenomeBin

# Peek at the design
pheno_MR50 <- as.data.frame(colData(SE), optional = TRUE)
table(pheno_MR50$Group, pheno_MR50$Timepoint)
```

Throughout, we denoise with a titration-justified taxonomic filter. Feel free to titrate first to see where these came from:

```r
pdf("MR50_00_titration.pdf", paper = "a4r")
tit <- titrate_feature_counts(ExpObj = SE, glomby = "Species",
                              only_allow_CSBs = FALSE, prevalence_pct = 5,
                              return_plots = TRUE, verbose = TRUE)
print(tit$plots)
dev.off()

# Filter used in the examples below
fc  <- c(25, 5)     # >= 25 PPM in >= 5% of samples
gcc <- c(10, 5)     # >= 10% genome completeness in >= 5% of samples
```

### Example 1 — Ordination: do the treatment groups separate? (endpoint)

The clearest single figure: at the Day 56 endpoint, colour by treatment group and test with PERMANOVA.

```r
pdf("MR50_01_ordination_endpoint.pdf", paper = "a4r")
ord <- plot_Ordination(
    ExpObj    = SE,
    glomby    = "Species",
    samplesToKeep = rownames(subset(pheno_MR50, Timepoint == "Day_56")),
    algorithm = "PCoA",              # interpretable axes + Bray-Curtis
    distmethod = "bray",
    compareby = "Group",
    colourby  = "Group",
    shapeby   = "Group",
    ellipseby = "Group",
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    permanova = TRUE, plotcentroids = TRUE, highlight_centroids = TRUE,
    grid = FALSE
)
print(ord[[1]])
dev.off()
```

### Example 2 — Ordination as a time series: watch the community move and recover

Colour by timepoint, connect each mouse across time, and see the antibiotic groups swing away from baseline and (partially) return. One panel per group via `subsetby`.

```r
pdf("MR50_02_ordination_timeseries_by_group.pdf", paper = "a4r")
ord_ts <- plot_Ordination(
    ExpObj    = SE,
    glomby    = "Species",
    subsetby  = "Group",             # one plot per group
    algorithm = "PCoA", distmethod = "bray",
    compareby = "Timepoint",
    colourby  = "Timepoint",
    connectby = "Mouse_number",      # trace each mouse
    connection_orderby = "Timepoint_days",   # in true time order
    textby    = "Timepoint",
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    permanova = TRUE, grid = FALSE
)
for (g in ord_ts) print(g)
dev.off()
```

### Example 3 — Comparative heatmap: which species drive the vancomycin effect?

At Day 5 (peak perturbation), contrast the untreated group against the vancomycin group and show only significant, sizeable differences. Remember the heatmap draws straight to the device.

```r
pdf("MR50_03_heatmap_G1_vs_G2_day5.pdf", paper = "a4r")
plot_relabund_heatmap(
    ExpObj    = SE,
    glomby    = "Species",
    hmtype    = "comparative",
    samplesToKeep = rownames(subset(pheno_MR50,
                    Timepoint == "Day_05" & Group %in% c("G1", "G2"))),
    compareby = "Group",
    splitcolsby = "Group",               # group the columns visually
    colcategories = c("Group", "Timepoint"),
    showonlypbelow = 0.05,
    minl2fc   = 1,
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    show_GenomeCompleteness = TRUE,       # trust check per taxon
    showGram  = TRUE,                     # vancomycin hits Gram-positives
    no_underscores = TRUE
)
dev.off()
```

The `showGram = TRUE` annotation is pedagogically nice here: vancomycin targets Gram-positive organisms, and the direction of change (row-label colour + Log2FC column) should line up with Gram status.

### Example 4 — Exploratory heatmap across the whole experiment

No hypothesis — just show the most variable species across all samples, split by group, annotated with timepoint. Good for a "big picture" opening figure.

```r
pdf("MR50_04_heatmap_exploratory_all.pdf", paper = "a4r")
plot_relabund_heatmap(
    ExpObj    = SE,
    glomby    = "Species",
    hmtype    = "exploratory",
    splitcolsby = "Group",
    colcategories = c("Group", "Timepoint"),
    ntop      = 40,
    scaled    = TRUE,                     # z-scores: compare patterns, not levels
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    no_underscores = TRUE
)
dev.off()
```

### Example 5 — Alpha diversity: collapse and recovery over time

The signature antibiotic result: richness crashes with treatment. Compare groups within each timepoint. Note we filter (never run this unfiltered).

```r
pdf("MR50_05_alpha_by_group_within_timepoint.pdf", paper = "a4r")
adiv <- plot_alpha_diversity(
    ExpObj    = SE,
    glomby    = "Species",
    measures  = c("Observed", "Shannon", "InvSimpson"),
    subsetby  = "Timepoint",             # one set of panels per timepoint
    compareby = "Group",
    colourby  = "Group", fillby = "Group",
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    stratify_by_domains = FALSE
)
for (g in adiv) print(g)
dev.off()
```

### Example 6 — Alpha diversity as a time course within one mouse line

Treat timepoint as the x-axis and connect each mouse, within the double-antibiotic group, to see individual trajectories of recovery.

```r
pdf("MR50_06_alpha_timecourse_G3.pdf", paper = "a4r")
adiv_ts <- plot_alpha_diversity(
    ExpObj    = SE,
    glomby    = "Species",
    measures  = c("Observed", "Shannon"),
    samplesToKeep = rownames(subset(pheno_MR50, Group == "G3")),
    compareby = "Timepoint",
    connectby = "Mouse_number",          # one line per mouse
    colourby  = "Timepoint", fillby = "Timepoint",
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    stratify_by_domains = FALSE
)
for (g in adiv_ts) print(g)
dev.off()
```

### Example 7 — A functional question: antimicrobial resistance genes

Switch to a functional space. Does the resistance-gene profile differ by treatment at the endpoint? (Functional spaces take no genome-completeness filter.)

```r
pdf("MR50_07_resfinder_endpoint.pdf", paper = "a4r")
plot_relabund_heatmap(
    ExpObj    = expvec_MR50$resfinder,
    hmtype    = "comparative",
    samplesToKeep = rownames(subset(pheno_MR50, Timepoint == "Day_56")),
    compareby = "Group",
    splitcolsby = "Group",
    colcategories = c("Group", "Timepoint"),
    showonlypbelow = 0.1,
    featcutoff = c(0, 0),               # functional: keep it light, no GC filter
    PPM_normalize_to_bases_sequenced = TRUE
)
dev.off()
```

You can do the same with `expvec_MR50$ECNumber` or `expvec_MR50$Product` to look at enzyme or gene-product shifts.

### Example 8 — The whole story in one call with `compare_within_subsets`

Finally, let the comparator emit ordination + heatmap + alpha for every timepoint, comparing groups within each. This is the workflow you would actually use day-to-day.

```r
hm_args  <- list(showonlypbelow = 0.05, minl2fc = 1,
                 max_rows_in_heatmap = 60, maxnumheatmaps = 2,
                 show_GenomeCompleteness = TRUE, showGram = TRUE,
                 no_underscores = TRUE,
                 colcategories = c("Group", "Timepoint"))

ord_args <- list(algorithm = "PCoA", distmethod = "bray",
                 plotcentroids = TRUE, highlight_centroids = TRUE, grid = FALSE)

alp_args <- list(measures = c("Observed", "Shannon", "InvSimpson"),
                 stratify_by_domains = FALSE)

pdf("MR50_08_compare_within_timepoints.pdf", paper = "a4r")
compare_within_subsets(
    ExpObj    = SE,
    compareby = "Group",
    subsetby  = "Timepoint",
    glomby    = "Species",
    only_allow_CSBs = FALSE,
    do_ordination = TRUE, do_heatmap = TRUE, do_alpha = TRUE,
    ordination_highlight_subset_in_context = TRUE,
    ordination_args = ord_args, heatmap_args = hm_args, alpha_args = alp_args,
    featcutoff = fc, GenomeCompletenessCutoff = gcc,
    PPM_normalize_to_bases_sequenced = TRUE,
    class_to_ignore = "N_A", return_stats = FALSE
)
dev.off()
```

Two nice variations for the tutorial:

- Swap `subsetby = "Timepoint"` for `subsetby = "Group"` with `compareby = "Timepoint"` to tell the story the other way round (how each group changes over time).
- Re-run Example 8 with `only_allow_CSBs = TRUE` to show how restricting to genuinely-binned genomes changes (or confirms) the picture.

> **Tip for building the tutorial figures.** Run these top to bottom in a fresh session; each writes one PDF into your working directory. Because MR50 effects are so pronounced, the PERMANOVA p-values, heatmap fold-changes, and alpha-diversity drops should all be unambiguous — ideal for annotating in the wiki.

---

## 15. Troubleshooting

**"It says my object only has `Contig_LKT` / I get odd taxonomic results."**
You are probably on a deprecated space. Use `ConsolidatedGenomeBin` and drop the rest: `expvec$Contig_LKT <- NULL`. See [§5](#5-analysis-spaces-taxonomic-and-functional).

**A comparative heatmap comes back with "No features fulfilling…".**
Nothing survived the combined filters. Loosen, in this order: raise `showonlypbelow` (e.g. to 0.1), lower `minl2fc` (or drop it), then relax `featcutoff` / `GenomeCompletenessCutoff`. Confirm your `compareby` variable actually has ≥2 classes with ≥2 samples each in the sample set being plotted.

**My heatmap didn't appear / I got an empty or corrupt PDF.**
`plot_relabund_heatmap` (and `compare_within_subsets`) draw directly to the device and return no plot object. You must open a device first and close it after: `pdf("x.pdf", paper="a4r"); plot_relabund_heatmap(...); dev.off()`. Trying to do `p <- plot_relabund_heatmap(...)` and `print(p)` will not work. See [§11](#11-relative-abundance-heatmaps).

**`update_ExpObj_metadata` errors that samples are missing.**
Your `pheno` contains sample names not present in the object (typo, stale name, or wrong sample column). Reconcile them: `setdiff(rownames(pheno), colnames(expvec$ConsolidatedGenomeBin))`. This error is protective — fix the metadata rather than bypassing it. See [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely).

**`undefined columns selected`, or a variable "not found", with a column you can plainly see.**
Almost always a column name with a space, hyphen, bracket or other special character, which R quietly renames in some places and not others. Clean the metadata with `enforce_JAMS_style_metadata()` (or read it with `load_metadata_from_file()`), push it into the objects with the replacement loop, and use the new names. See [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely).

**My metadata column names changed between JAMS runs.**
Since metadata cleaning was unified, JAMSbeta and JAMS16 make every column name JAMS-style (e.g. `IL.6..pg.mL.` or `IL-6 (pg/mL)` becomes `IL_6_pg_mL`), and the sample column is always called `Sample`. The exact mapping is in the project's `_JAMS.log` (one line per rename) and in `<project>_metadata_column_name_changes.tsv` in the JAMSbeta output folder. See [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely).

**`plot_Ordination` returned an empty list.**
It refused to plot and said why in a warning. The usual causes are a continuous `compareby` with no `dichotomize_continuous_compareby_threshold`, a variable missing from the object, a `compareby` that is entirely `N_A`, or fewer than 3 samples. See [§10](#10-ordination-plots).

**PERMANOVA is skipped with a warning.**
A subset had too few samples in the smallest `compareby` class (default minimum is 3 per class), or only one class survived after `class_to_ignore` removals. The ordination is still plotted; only the p-value is withheld. Increase samples, merge classes, or lower `min_samples_per_class_for_permanova`.

**Alpha diversity looks implausibly high / groups all look the same.**
You almost certainly ran it too loosely. Never use `featcutoff = c(0,0)` for diversity; unfiltered noise inflates richness. Use at least `applyfilters = "light"` or a titration-based cutoff. See [§12](#12-alpha-diversity).

**Colours are inconsistent between plots, or random.**
The object has no `ctable`, or you passed a `cdict` in some calls but not others. Rebuild a colour dictionary with `make_colour_dictionary()` and attach it to every object's metadata (see [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely)); once `ctable` is in the object's metadata the plotting functions use it automatically. To keep specific classes fixed across studies, pass a `legacy_ctable`. To keep colours consistent across successive re-runs, note that `make_colour_dictionary()` is deterministic (no random palette) unless you set `shuffle = TRUE`.

**Two categories of one variable look identical on a plot.**
Two classes of the same variable have been assigned the same colour. Regenerate with `make_colour_dictionary(..., avoid_collisions = TRUE)` (the default), which guarantees within-variable distinctness, or fix the offending pair directly in `tweak_ctable(ctable, pheno = pheno)` — its per-variable duplicate report points straight at the clash. See [§4](#4-metadata-pheno-ctable-and-replacing-metadata-safely).

**`tweak_ctable()` says it needs shiny / colourpicker.**
Those packages are intentionally not JAMS dependencies. Install them once with `install.packages(c("shiny", "colourpicker"))` and try again; until then the function returns your `ctable` unchanged.

**`GenomeCompletenessCutoff` "isn't doing anything" / warns it was ignored.**
You applied it to a functional space (it only works on taxonomic objects with a `GenomeCompleteness` assay), or the object predates that assay. Use it on `ConsolidatedGenomeBin`.

**`only_allow_CSBs = TRUE` warns and is ignored.**
That restriction only applies to a `ConsolidatedGenomeBin` object. It is meaningless for functional spaces and for other taxonomic objects. See [§6](#6-core-concepts-lkt-cgb-and-csb).

**`featuresToKeep` doesn't match after agglomeration.**
When you set `glomby`, feature names must be **post-agglomeration** names (e.g. `f__Enterobacteriaceae` for `glomby = "Family"`), not the original CGB/LKT names. See [§9](#9-arguments-common-to-all-the-plotting-functions).

**Relative abundances don't match between two functional plots.**
Check `PPM_normalize_to_bases_sequenced`. `TRUE` normalises to total sample depth; `FALSE` to in-analysis (coding) depth. For functional objects these differ. Be consistent across a project. See [§7](#7-the-shared-filtering-vocabulary-of-jams-functions).

**Everything is very slow or runs out of memory on big projects.**
Taxon-stratified functional objects hold large sparse matrices. Restrict to the analyses you need (fewer spaces), subset samples, agglomerate, and prefer running big comparative sweeps on a machine with more RAM. The example dataset and typical taxonomic work are laptop-friendly.

---

## 16. Quick reference

### Objects

| Object | What it is |
|---|---|
| `expvec` | Named list of `SummarizedExperiment` objects, one per analysis space |
| `pheno` | Sample metadata data frame (rownames = sample names) |
| `ctable` | Colour dictionary; lives in each SEobj's metadata |
| `expvec[[X]]` | One experiment: `assays`, `rowData`, `colData`, `metadata` |

### Taxonomic space

- **Use `ConsolidatedGenomeBin`.** Agglomerate with `glomby` (`"Genus"`, `"Species"`, …).
- **Do not use `Contig_LKT`** (deprecated) or `MB2bin` (phased out).
- **LKT** = last known taxon (leaf label). **CGB** = consolidated per-sample entity with completeness/contamination. **CSB** = functionally-clustered, genuinely-binned genome; `only_allow_CSBs = TRUE` keeps these only.

### Filtering vocabulary

| Argument | Meaning |
|---|---|
| `featcutoff = c(PPM, pct)` | keep feature if ≥ PPM in ≥ pct % of samples |
| `GenomeCompletenessCutoff = c(comp, pct)` | (taxonomic) keep if ≥ comp% complete in ≥ pct % of samples |
| `applyfilters` | `"light"`/`"moderate"`/`"stringent"` presets; explicit values override |
| `minl2fc` / `maxl2fc` | min/max absolute log2 fold change |
| `only_allow_CSBs` | restrict `ConsolidatedGenomeBin` to bin-derived genomes |
| `PPM_normalize_to_bases_sequenced` | total-bases (TRUE) vs in-analysis (FALSE) denominator |
| `class_to_ignore` | drop samples with this token in the comparison variable (default `"N_A"`) |

### The three plots

| Function | Question | Returns | Device? |
|---|---|---|---|
| `plot_Ordination` | different overall? | list of ggplots | you may print/edit |
| `plot_relabund_heatmap` | which features differ? | stats list (not a plot) | **must** wrap in pdf/dev.off |
| `plot_alpha_diversity` | richness different? | list of ggplots (+stats) | you may print/edit |
| `compare_within_subsets` | all three, per subset; one fixed global feature set across subsets | invisible stats (optional) | **must** wrap in pdf/dev.off |


### Don't-forget list

- Build `ctable` with `make_colour_dictionary()` from curated metadata (use `columnsToUse` / `columnsToIgnore`; `Sample` is always protected), and re-attach it after `update_ExpObj_metadata`.
- Reuse a `legacy_ctable` to keep fixed project colours across studies; hand-tune the rest with `tweak_ctable(ctable, pheno = pheno)`.
- Collision avoidance is on by default — every class within a variable gets a distinct colour.
- Clean curated metadata with `enforce_JAMS_style_metadata()` (or read it with `load_metadata_from_file()`, the same code path JAMSbeta uses) before pushing it into the objects; it also sets `rownames(pheno)` to the sample names for you.
- Keep `pheno` and every `colData` in sync — always propagate via the replacement loop.
- To pull metadata out of an object, use `as.data.frame(colData(SEobj), optional = TRUE)` so column names are not mangled.
- A continuous `compareby` in `plot_Ordination` needs `dichotomize_continuous_compareby_threshold` (a number, `"median"`, or a function); for a reusable split, add a column with `dichotomize_variable_in_md()`.
- Titrate before trusting a filter.
- **Never** run alpha diversity unfiltered.
- Wrap heatmaps and `compare_within_subsets` in your own `pdf(...)` / `dev.off()`.
- `compare_within_subsets` computes **one global feature set up front** so panels stay comparable across subsets; use `samplesToKeep` / `featuresToKeep` to drop outliers/noise before that step. `featuresToKeep` uses **post-agglomeration** names (same as the standalone plotting functions) — e.g. species names when `glomby = "Species"`. A wrong namespace or stale list is caught (warns at a high miss rate, errors at a very high one), never plotted silently.
