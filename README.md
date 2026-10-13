# melanoma_ICB
Jupyter notebooks and analysis files used to evaluate sex differences in biomarker validity for ICB response.

## Glutamine metabolite impact ranking

`scripts/` + `config/metabolite_panel.tsv` rank glutamine-derived metabolites by
the predicted impact of glutamine supplementation, using Human-GEM constrained
to the Riaz anti-PD-1 melanoma transcriptomes, then test the ranking against
response and survival in the same cohort.

- `docs/METHODS.md` — why the stoichiometric form of the question is not
  answerable, and what replaces it
- `docs/RESULTS.md` — the ranking, and the falsification step it fails
- `scripts/run_all.sh` — end-to-end rerun (`WORK=<dir> scripts/run_all.sh`)
