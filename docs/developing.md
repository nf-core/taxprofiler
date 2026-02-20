# nf-core/taxprofiler developer checklists

This page can act as a reference for new developers who wish to contribute to the nf-core/taxprofiler code base, and act as a reviewing checklist for reviewers.

## Adding a new profiler contribution workflow

> Note: does not have to be in this precise order

> [!WARNING]
> Before starting, always make sure you're on a fork and a new branch!

- [ ] Installed modules `nf-core modules install <tool>/<subtool>`
- [ ] Added profiler to properties in `assets/schema_database.json`
- [ ] Added parameters to `nextflow.config`
  - [ ] Added `--run_<profiler>`
  - [ ] Added other profiler-specific parameters (outside of `meta.db_params`, e.g. saving optional files [read classifications])
- [ ] Added profiler(s) to `subworkflows/local/profiling.nf`
  - [ ] Added relevant modules at the top using `include` statement
  - [ ] Added the tool to the read-database combining section
  - [ ] Added the tool-specific if/else statement
  - [ ] Within if/else added `ch_input_for_profiling.profiler` filtering
  - [ ] Within if/else invoke module(s)
  - [ ] Within if/else raw profile and classification (if necessary) channels mixed
  - [ ] Version and MultiQC (if available) channels mixed
- [ ] Update `modules.conf`
  - [ ] Added `withName:` block
  - [ ] Added `ext.args = { "${meta.db_params}" }` (no other parameters should be added to this for the actual profiling modules!)
  - [ ] Added other args for non-profiling modules (outside of `meta.db_params`, e.g. saving optional files [read classifications])
  - [ ] (Profiling) Added ext.prefix conditional to account for run merging (see other tools)
  - [ ] (Standardisation) Added `ext.prefix = { "ganon_${meta.id}_combined_reports" }` for native multi-sample tables (i.e. not from taxpasta)
  - [ ] Added publish dir for the particular profiler + database, with pattern (see other tools)
- [ ] Added profiler(s) to `subworkflows/local/standardisation_profiles.nf` if necessary
  - [ ] Added relevant modules at the top using `include` statement
  - [ ] Added the tool to the read-database combining section (if necessary)
  - [ ] Invoke module(s)
  - [ ] Version and MultiQC (if available) channels mixed
- [ ] If necessary, added any profiler-specific parameter validation checks within the top of `utils_nfcore_taxprofiler_pipeline/main.nf`
- [ ] Updated Documentation
  - [ ] `nf-core pipelines schema build` has been run and updated
    - [ ] All additional tool specific pipeline parameters have a additional help entry with the `Modifies tool parameter(s)` quote block
  - [ ] Added citation to `citations.md` (citation style: APA 7th edition)
  - [ ] Add citation to the toolCitation/BibliographyText functions in `subworkflows/local/utils_nfcore/taxprofiler_pipeline
    - [ ] Added in-text citation
    - [ ] Added bibliography (citation style: APA 7th edition)
  - [ ] Added relevant documentation to `usage.md`
    - [ ] Entry in example database sheet (full database sheet section)
    - [ ] Entry in 'Expected' paths bullet points (full database sheet section)
    - [ ] Entry in 'Classification and profiling' tips section
    - [ ] Relevant entries in Post Processing section (Visualisation, multi-table generation)
    - [ ] Entry in 'Building custom databases' `docs/usage/tutorial.md`
  - [ ] Described module output in `output.md`
    - [ ] Entry in table of contents
    - [ ] Entry in results section including tool description, expected output files list, and general tips
    - [ ] Entry in taxpasta section
    - [ ] Entry in MultiQC section
  - [ ] Added to pipeline summary list on `README.md`
  - [ ] (Optional) Added to pipeline metro map diagram (can be done just before release)
  - [ ] If it's your first contribution, add or move yourself to the Team list on `README.md`!
- [ ] On nf-core/test-data: updated `database_vX.X.csv` to include the test dataset
- [ ] Updated `test*.config`s
  - [ ] Added `database_vX.X.csv` to all test configs (only once per release with first profiler!)
  - [ ] Added `run_<profiler>` to test configs (where applicable)
- [ ] Test(s) pass

## Outside of pipeline repository

- [ ] Make full test database
  - [ ] See `taxprofiler` branch of test datasets for instructions how to get raw FASTAs
  - [ ] Test database locally
  - [ ] Once built and test, upload to iGenomes s3 (Ask James)
  - [ ] Update `database_full_vX.X.csv` and README to include the new s3 URI and instructions
  - [ ] Open PR against test-datasets, taxprofiler branch
- [ ] Add a MultiQC module
- [ ] Make a Taxpasta module
- [ ] Add the database building module to nf-core/createtaxdb (where possible)
