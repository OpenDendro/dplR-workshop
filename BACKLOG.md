# Backlog

Working notes for the rewrite. Not part of the book.

## Start here (left open 2026-10-06)

State of the repo: the Quarto conversion, the restructure into parts, the stubs and the dataset work are on the `quarto-rewrite` branch. `main` is still the 2024 bookdown version. The folder moved from `workshopMaterials/dplR-workshop` to `openDendro/learningToLoveDplR`; the remote is still `OpenDendro/dplR-workshop`. The writing spec is a local file and is not in the repo.

### Structure suggestions, undecided

- Split "Standardization and Chronologies" in two. Detrending, Chronology Statistics and Chronology Building are the core every reader needs. Regional Curves, Signal-Free and Basal Area are specialist chapters, and basal area is neither standardization nor a chronology.
- Give the ecological material one home. Basal Area, `treeMean` and Events serve readers asking about trees and disturbance, and they sit in three places. Grouping Basal Area and Events in one part, with tree-level series there too, would also take Events out of Time Series.
- Trim the time series chapters (14, 15) against the time series book. Keep what is specific to tree rings (the spline and its frequency response, the residual chronology, `redfit`, wavelets) and point to that book for the rest.
- Crossdating now comes after standardization, so the book's order differs from the order of a project. The Preface and A First Pass should say that dating comes first in practice, and Checking a Collection should point forward to the crossdating chapters.

### Open decisions

- dplR version in CI. `DESCRIPTION` installs dplR from CRAN. Add `OpenDendro/dplR` under `Remotes:` to build on dev.
- Deployment. The workflow in `.github/workflows/publish.yml` deploys with GitHub Actions and `docs/` is gitignored but still tracked. Before pushing: switch Pages to "workflow" and run `git rm -r --cached docs`. Until both are done the old bookdown site stays live.
- Cover page, favicon, `LICENSE.md`, `CITATION.cff`, `.zenodo.json`: not started. `index.qmd` is a placeholder.
- The 54 chunks that were unnamed in the 2022 chapters are still unlabeled.
- Whether to split the chapter skeleton and code conventions out of the writing spec into a `CRC-POLISH.md`, as in the other books.

### Suggested next step

Write one chapter all the way through to the skeleton and the spec before adjusting the plan further. Chronology Statistics (06) is the candidate: the running rbar and EPS section is critical, it does not depend on the signal-free paper, and it will settle the dataset and plotting conventions.

### Outside this repo

- dplR: `strip.rwl()` stops with a `data.frame` row-count error when it removes a series (1.8.1 dev; CRAN 1.8.0 not checked). Likely cause: the `[` methods for `rwl` and `rwi` trim years, so the reinsert step binds objects of different lengths. Not fixed.
- dplR: `ssf(return.info = TRUE)` returns a list with class `crn`.
- dplR `TODO` has a new item: revamp the plots.

## Standing notes

- Plots: include new ggplot figures where they are useful. Most of dplR's own plot methods date from before 2010 and are due a revamp (see dplR's TODO), so don't lean on them where a ggplot shows the data better.

## Running datasets (decided 2026-10-06)

- `co021` (Mesa Verde Douglas-fir) and `wa082` (Hurricane Ridge Pacific silver fir) are carried through the book. The reasons, with the evidence, are written up in A2DataSources.
- `co021` is the strong-signal collection: 14th of 6,764 ITRDB ring-width files for interseries correlation, no flagged segments until an error is planted, 716 absent rings.
- `wa082` is the ordinary one: 44th percentile, EPS near 0.85, and the sample-depth and SSS cutoffs disagree. Its raw file has a `-999` gap in 712011 at 1900 that the on-board copy holds as zero.
- `gp.rwl` for Regional Curves and Basal Area (pith offsets and diameters; `zof.rwl` has them too). `ca533` for Signal-Free only. `cana157` for Events.
- `data/` holds `co021.rwl` and `wa082.rwl` as served by the ITRDB. Chapters 02, 04, 05 and 06 still use `ca533` and need switching when they are rewritten.
- Still needed: a deliberately broken file for the Tucson aside; restore the "zip `data/` into the site" step in `publish.yml`; an entry for the two files in a permissions log if the book goes to a publisher.
- Rule adopted with this decision: don't describe what a dataset is good for without checking it against the data. Unchecked claims are listed in the callout at the top of A2DataSources.

## Coverage decisions (2026-10-06)

- `net` is left out of the book.
- `glk` and `sgc` get an aside (11xAsideSignAgreement).
- `sea` and `pointer` get a chapter (16Events).
- `treeMean` is framed the way it is used: ecological work on tree-to-tree differences.
- The smaller functions (`series.rwl.plot`, `skel.plot`, `xskel.plot`, `common.interval`, `read.crn`, `bakker`, `gini.coef`) go in only where they fit. Not if clunky.
- Chronology Statistics now comes before Chronology Building, because the SSS cutoff, EPS stripping and variance stabilization depend on rbar, EPS and SSS.

## Planned revisions to carried-over chapters

- 07Chronologies, `chron.ci`: say what the interval means. It is uncertainty in the mean of these series, not in a reconstruction.
- 07Chronologies, `chron.stabilized`: show the problem first (variance tracking sample depth and rbar), then the fix. dplPy treats stabilization as a step applied to any finished chronology; dplR's function builds one from indices.
- 11Crossdating: the chapter never mentions `xdate.report` or `write.xdate.report`.
- 14SmoothingFiltering: lead with `pass.filt` and fold the manual `signal` padding code. The splines section still names `ffcsaps`, which is defunct.

## Known problems in carried-over chapters

- 07Chronologies: `strip.rwl()` stops with an error in dplR 1.8.x (dplR bug). The chunk has `error: true`.
- 09SignalFree: written against the 2024 output of `ssf(return.info = TRUE)`; nine chunks fail. The chapter has `error: true`.
- A3Acknowledgements fetches the author list from crandb at build time.
