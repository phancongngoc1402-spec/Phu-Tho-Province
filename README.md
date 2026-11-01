# Temporal transferability of mass-calibrated concentration thresholds

Self-contained reproducibility package for **Temporal Transferability of Mass-Calibrated Concentration Thresholds for Pollutant-Mass Prioritization in High-Frequency Industrial Effluent Monitoring in Phu Tho Province, Vietnam**.

## Clean-environment quick start
Python 3.11 is the tested interpreter. Dataset.xlsx is bundled. No API key, external repository checkout, hard-coded machine path, intermediate pickle or network reference lookup is required to execute the analysis after installing dependencies.

macOS / Linux:
```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
python validate_final.py
```
Windows PowerShell:
```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python run_all.py
python validate_final.py
```
If activation is unavailable, invoke `.venv/Scripts/python.exe` directly; do not alter system policy. All paths are resolved from __file__, not the working directory. `--output /another/path` is supported. Dependencies are pinned in requirements.txt. Full calculation uses 1,000 primary and 500 secondary replicates and can take several minutes.

`outputs/final/` contains ONLY final manuscript tables and figures. OUTPUT_MAP.csv is the file-to-table map. `legacy/uploaded_outputs/` is labelled NOT FOR MANUSCRIPT, retained for provenance and never read by the final pipeline. `reference_final/` is the numerical checkpoint for validate_final.py; `provenance/` stores execution and validation logs. The internal `outputs/base/` and `outputs/work/` are pipeline staging areas, not manuscript results. `--base-only` is an internal resumable stage, not a complete run; a normal run recomputes the base and all remaining steps. Optional --resume reuses a base stage only if its dataset, code and dependency signature matches.

## Data and calibration
Data sheet Data_Tong_hop has 148,906 unique five-minute timestamps from 2025-01-01 00:00 through 2026-07-31 23:55 (UTC+7).
SHA-256: 253fd20f2e26f9e5b5cae57a4ed1140433f01ef141805e20b7b1f11b226fa283.
Accepted statuses are Hoạt động tổt and Vượt qui chuẩn. Primary pairs require finite positive chemistry/flow. Missing slots are not imputed or converted to zero flow. The legacy NH4+ field is treated as NH4-N per study documentation, without chemical conversion. Load rate=C×Q/1000 (kg/h); interval mass=load rate×5/60. No retrospective legal compliance is adjudicated.

analysis/calibration.py is the ONE canonical Algorithm 1 implementation for observed data, weighted burden bootstrap and nested re-estimation. Stable target/concentration ranks locate the required target hit; that fraction is rounded upward on the 0.5-percentage-point grid. All concentration ties are included AFTER selection. This is not an alternative inclusive-tie scan of earlier grid points. Weighted integer multiplicities are tested against explicit stable row-copy expansion. Moving-block evaluation uses stable mass sorting.

## Final streams and script responsibilities
- Primary joint/frozen-target Table 3/S20: seed20261026, shared COD→TSS→NH4-N draws.
- Fixed historical Table5/S7b/c: seed20260828, shared draws; both observed 2025 boundaries fixed. Target-copy and unique-support summaries use identical draws.
- Monthly S4a: 20261026+20×pollutant_index+month.
- Table2/S19 burden: 20260831+10×pollutant_index+year_offset (2025=0,2026=100).
- Nested S8:725911+pollutant_index; moving blocks:725911+100×block_length+pollutant_index.
- JS S21:20261026+200+pollutant_index; hourly/event seeds in analysis/config.py.
Indices COD=0,TSS=1,NH4-N=2. Full software/seed metadata is outputs/final/run_metadata.json.

revision_statistics.py contains the sole FINAL pairwise-Pearson ACF implementation; exact calendar gaps remain gaps, vectors are separately centred/scaled at every lag. The obsolete globally centred ACF approximation is absent from the audit script. Legacy conditional.csv and unified_fixed_ci.csv must NOT be used as S20/S7b. Internally generated alternative outputs stay in work/, and final mappings select the specified streams only.

Monthly/quarterly targets are ranked within each window and are not pooled-target decompositions. Pooled_target_by_month.csv provides a separate pooled diagnostic. May NH4-N contains no pooled Jan–Jul target but has38 month-specific targets. S6a removes13July while holding full-2025 concentration thresholds fixed and reranking remaining target mass; coverage-restricted calibration is separate in S6b.

COD fixed-target100–100% is a degenerate conditional interval: all26 observed targets were detected;6/1000 target-free replicates have undefined recall and are excluded only from recall/mass percentiles. No perfect future detection is guaranteed. ACF mean-ESS is descriptive, not tail-recall information; missingness mechanisms and measurement error are not modeled.

## Validation and provenance
`python validate_final.py` checks final table checkpoints, primary count identities, exclusion partitions, ACF bounds, fixed-target support and independent calibration-oracle tests. `python validate_key_results.py --input Dataset.xlsx --full --revision-bootstrap` checks155 deterministic/historical endpoints. Passing checks does not prove unmeasured assumptions.
Original base provenance: author repository commit088cdb9ad58c6c0f7b98a264a3d604acff5460b9. No remote repository update or DOI deposit is claimed. CHECKSUMS.sha256 records release files. Data/code are provided for verification; no additional licence is inferred. The author should choose an explicit licence before broader reuse/deposit.
