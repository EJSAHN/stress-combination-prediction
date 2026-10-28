# Plant stress-combination prediction

Reproducible analysis of combined-stress transcriptome prediction and gene-wise prediction gains in *Arabidopsis thaliana* and citrus.

The pipeline uses author-processed expression matrices from three GEO studies. It compares predictions of combined-stress responses, evaluates conditional noise sensitivity, and examines how gene-wise prediction gains repeat within and between genotypes. Models are fitted separately within each genotype; no fitted model is transferred between species or genotypes.

## Run

Use an existing 64-bit Python with NumPy, pandas, matplotlib, Pillow and openpyxl. The launchers do not create environments or install packages. SciPy is needed only for separately rerunning historical QR audits.

This GitHub distribution stores one frozen reference in three byte chunks and excludes the four original GEO matrices. First restore the reference locally, then download and verify the original matrices:

```sh
python scripts/restore_reference.py
python scripts/fetch_inputs.py
python scripts/fetch_inputs.py --check-only
python scripts/run_pipeline.py --output runs/confirmation_01
```

The restoration step uses only the standard library, makes no downloads, verifies every chunk and the assembled file, and preserves an existing verified reference. It refuses to overwrite an existing mismatched file. Run `python scripts/restore_reference.py --check-only` to verify restoration.

Choose a fresh output directory. Downloads must match recorded sizes and SHA256 hashes. The complete local ZIP includes all required inputs and runs offline.

On Windows, run the same restoration and input commands with the existing Python environment, then use the launcher:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run.ps1
```

The launcher never downloads inputs automatically. See [Windows instructions](docs/WINDOWS_RUN.txt) for selecting an existing interpreter and setting the project directory.

## Outputs

| Path below the output directory | Contents |
|---|---|
| `analysis/tables/` | Recomputed scientific results |
| `analysis/validation/` | Comparisons with pinned reference results |
| `figures/` | Four main and four supplementary figures, vector PDF, 300 dpi PNG and numerical sources |
| `tables/Results.xlsx` | Results workbook with 19 sheets |
| `logs/` | Execution logs |
| `validation.json` | Overall status, environment and completed stages |

## Reproduction status

Release **1.0.1** completed the full pipeline on Linux and native Windows. The Windows PowerShell run on 28 October 2026 passed analysis, figure generation and Excel export in approximately 6 minutes 23 seconds using an existing Anaconda environment.

- 61 reference tables and 194 intermediate-array checks passed.
- All 36,000 conditional simulation draws and 1,024 complementary-pool refits were recomputed.
- All eight PNG figures were byte-identical between the reviewed Windows and Linux runs.
- All 19 workbook sheets agreed numerically. Relative source-path separators differ by operating system.

The [validation record](docs/validation/native_windows_confirmation.json) distinguishes the tested runtime from subsequent documentation and packaging changes. Reproduction confirms computation on these inputs; it is not a new biological validation experiment.

The release corrects auxiliary noise diagnostics at an exact-zero projection boundary and compares two declared JSON singular-value vectors numerically. Prediction and simulation engines, source inputs, seeds and numerical tolerances are unchanged. Raw diagnostics and explicit correction records are retained. See [the correction note](docs/NUMERICAL_PORTABILITY_CORRECTION.txt) and [change history](CHANGELOG.md).

## Data and provenance

Source accessions: **GSE147962, GSE330879 and GSE203331**. Exact file URLs, hashes, sizes and expression units are recorded in [data sources](config/data_sources.json). The pipeline starts from deposited processed matrices, not raw sequencing reads.

`scripts/` contains the entry points, `src/stress_prediction/` coordinates analysis, and `config/` records inputs and reference mappings. `provenance/historical/` preserves audited engines and fixed records. Original reference tables remain under `provenance/reference_tables/`; two explicitly amended diagnostic references are stored separately under `provenance/portability_reference_tables/`.

## Interpretation

Treatment cells contain three RNA pools. Genes, plants within pools and overlapping refit partitions are not additional independent RNA replicates. The 1,024 fits use complementary pools within each partition, but different partitions reuse observations.

This release reproduces known outcomes. Repeating fit-before-evaluation does not make previously inspected targets unseen again. The discovery selection failure and unfavorable conditional noise result remain included. An observed projection residual is not identified as a pure higher-order biological interaction. Conditional simulation exceedance fractions are not biological P-values, and refit quantiles are not confidence intervals.

See [interpretation limits](docs/CLAIMS_AND_LIMITS.txt) and [figure legends](docs/figure_legends.txt).

## License and attribution

Original code is distributed under the [MIT license](LICENSE). This license does not relicense deposited expression data, third-party annotations or cited publications. Cite the original studies and accessions when reusing data. See [data and software attribution](docs/DATA_AND_SOFTWARE_LICENSES.txt).
