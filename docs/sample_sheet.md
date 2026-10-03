# Sample sheet (`config/samples.tsv`)

Tab-separated, one row per **sequenced sample** (a specimen sequenced on two platforms has two rows sharing
`specimen_id`). Validated at start-up by `mbstats.load_samples`; errors name the offending rows/columns.

| column | required | description |
|---|---|---|
| `sample_id` | yes | unique; letters, digits, `_ . -` only |
| `specimen_id` | yes | biological specimen; links the 16S and shotgun sample of the same stool/biopsy |
| `subject_id` | yes* | participant (random effect, CV grouping). *`NA` allowed for negative controls |
| `data_type` | yes | `16S` or `MGX` |
| `batch` | yes | sequencing run/batch (DADA2 error models are learned per batch; also a model covariate) |
| `sample_type` | yes | `sample` or `negative_control` (16S controls feed `decontam`) |
| `week` | yes* | study week/time of collection (numeric) |
| `diagnosis` | yes* | one of `design.diagnosis_levels` (default `nonIBD`, `CD`, `UC`); one per subject |
| `fastq_1`, `fastq_2` | yes | local path **or** http(s)/ftp URL (downloaded with curl); paired-end only |
| `age`, `bmi` | no | numeric covariates |
| `sex`, `antibiotics`, `immunosuppressants` | no | categorical covariates (`Yes/No` normalised) |
| `hbi`, `sccai` | no | activity scores (CD: Harvey-Bradshaw, UC: SCCAI) → `active`, pre-flare labels |
| `assemble` | no | `yes` = assemble this shotgun sample for MAG recovery |

Derived automatically (`metadata/`): `active` (score ≥ threshold), `flare_next` (inactive now and active at the
next visit within `design.flare.max_gap_weeks`; 0 if it stays inactive; missing otherwise), flare events with visits
aligned to onset, `is_baseline` (earliest sample per subject and platform), `visit_index`.

Any additional columns are carried through and can be added to `design.covariates_*`.
