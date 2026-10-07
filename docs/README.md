# Documentation index — Part II (GLP-1 outcome analysis)

| Document | Answers |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Where this repo sits in the paper's platform; Mermaid diagrams of the pipeline with file references; runner phases; the one statistical pattern |
| [CODE_REFERENCE.md](CODE_REFERENCE.md) | Every script: inputs, outputs, exact model formulas, constants, shared-module usage, known defects |
| [DATA_SCHEMA.md](DATA_SCHEMA.md) | The source CSV column by column, the analysis-ready CSVs, every intermediate and output artefact, the output tree |
| [PAPER_COVERAGE.md](PAPER_COVERAGE.md) | Methods-to-code and display-item maps; what is upstream and unavailable; modularity assessment; defect list |
| [REUSE_GUIDE.md](REUSE_GUIDE.md) | Running the same design on a different clinical question or on MIMIC; the exact edit points; what you must write yourself |

Companion repository: `glp1_nature2026_extraction_validation` (Part I, extraction validation)
has the same five documents under `docs/`.

Quick facts:

- Python 3.13, `pip install -r requirements.txt` (`lifelines` is needed only by `conf_int/gap_120/step8_survival_plots_and_cox.py`).
- Tests: `python3 tests/test_persistence_and_step0a.py` (no data needed; all pass).
- Nothing else runs without the licensed source CSV.
- The manuscript numbers came from `rerun_conf_int_clean_full.sh`.
