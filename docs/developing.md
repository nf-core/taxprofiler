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
  - [ ] Add citation to the toolCitation/BibliographyText functions in `subworkflows/local/utils_nfcore_taxprofiler_pipeline/main.nf`
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
- [ ] Add nf-test file and snapshot. See [nf-test specifications](#new-nf-test-procedure-and-specifications) for more information.
- [ ] Test(s) pass

## Outside of pipeline repository

- [ ] Make full test database
  - [ ] See `taxprofiler` branch of test datasets for instructions how to get raw FASTAs
  - [ ] Test database locally
  - [ ] Once built and test, upload to iGenomes s3 (Ask James)
  - [ ] Update `database_full_vX.X.csv` and README to include the new s3 URI and instructions
  - [ ] Open PR against test-datasets, taxprofiler branch
- [ ] Add a [MultiQC](https://github.com/multiqc/multiqc) module
- [ ] Make a [Taxpasta](https://github.com/taxprofiler/taxpasta) module
- [ ] Add the database building module to nf-core/createtaxdb (where possible)

## New nf-test procedure and specifications

### Procedure

When writing a new pipeline-level nf-tests for nf-core/taxprofiler, we recommend the following procedure:

1. Run the test profile locally to have a copy of the expected output files and the results directory structure
2. Check the contents of files are expected (one or two files per directory should be enough)
2. Write the base nf-test file structure, assuming _all_ files are stable (following the specifications below)
3. Run `nf-test --tag <test_name> --profile +docker` once to write the first snapshot
4. Run command above to get the `diff` of unstable files
5. Update the assertions in each directories `.match()` snapshot for unstable files

### Specifications

Write the test files following the specifications below.

The necessary files are are follows:

- [ ] New test files should go under `tests/`
- [ ] Test file should be called `<test_config_name>.nf.test`
- [ ] Snapshot file should be called `<test_config_name>.nf.test.snap`

nf-test file contents:

- Test file header
  - [ ] Specify name as: `Test <config name>` 
  - [ ] Specify two tags: `pipeline` and `<config_name>`
  - [ ] Specify profile as `<config_name>`
- Test block
  - [ ] Specify the name of the test block as `test("-profile <config name>")`
- When block
  - [ ] Specify when block with single param, `outdir`
    - All other parameters should be specified in the config itself
- Then block
- [ ] Specify on the first line, a `stable_name_all` variable to list all file names with the nft-utils `getAllFilesFromDir` function
- [ ] For each top-level output directory under `results` (typically, one per tool), specify a comment with the tool name and a `stable_content_<tool name>` variable
- `assertAll` block
  - [ ] Use the `removeNextflowVersion` function
  - [ ] Check existance of `nf_core_taxprofiler_software_mqc_versions.yml` file
  - [ ] Check existance of  `multiqc_report.html` file
  - [ ] Snapshot `stable_name_all` with a `.match()` name of `all_files`
  - [ ] For each results directory, make a snapshot with a `.match()` name of `directory_name` (typically the tool name all lower case). Inside this:
    - [ ] Specify the `stable_content_<tool_name>` variable
    - [ ] For each unstable file, specify specific path in `.nftignore`
    - [ ] For each unstable file, specify an alternative method of file checks (e.g. sorted file, file size check, contains string, nft-plugin function)
    - [ ] For each unstable file assertion, include a string before the assertion itself with the file name and type of check