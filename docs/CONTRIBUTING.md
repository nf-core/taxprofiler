---
title: Contributing
markdownPlugin: checklist
---

# `nf-core/taxprofiler`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving nf-core/taxprofiler.

This page describes the recommended nf-core way to contribute to both nf-core/taxprofiler and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of nf-core/taxprofiler.

> [!NOTE]
> If you need help using or modifying nf-core/taxprofiler, ask on the nf-core Slack [#taxprofiler](https://nfcore.slack.com/channels/taxprofiler) channel ([join our Slack here](https://nf-co.re/join/slack)).

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/nf-core/taxprofiler/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/taxprofiler repository](https://github.com/nf-core/taxprofiler) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

## Use of AI and LLMs

The nf-core stance on the use of AI and LLMs is that humans are still ultimately responsible for their submitted code, regardless of the tools they use.

If you’re using AI tools, try to stick by these guidelines:

- Keep PRs as small and focussed as possible
- Avoid any unnecessary changes, such as moving or refactoring code (unless that is the explicit intention of the PR)
- Review all generated code yourself before opening a PR, and ensure that you understand it
- Engage with the community review process and expect to make revisions

For more detail, see the the [blog post](https://nf-co.re/blog/2026/statement-on-ai) for a statement from the nf-core/core team.

### Getting help

For further information and help, see the [nf-core/taxprofiler documentation](https://nf-co.re/taxprofiler/usage) or ask on the nf-core [#taxprofiler](https://nfcore.slack.com/channels/taxprofiler) Slack channel ([join our Slack here](https://nf-co.re/join/slack)).

### GitHub Codespaces

You can contribute to nf-core/taxprofiler without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/nf-core/taxprofiler/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test --tag test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Update the snapshots with the following command:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the nf-core/taxprofiler code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](lint) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

### Adding a new profiler contribution workflow

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
  - [ ] (Profiling) Added `ext.prefix` conditional to account for run merging (see other tools)
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

#### Outside of pipeline repository

- [ ] Make full test database
  - [ ] See `taxprofiler` branch of test datasets for instructions how to get raw FASTAs
  - [ ] Test database locally
  - [ ] Once built and test, upload to iGenomes s3 (Ask James)
  - [ ] Update `database_full_vX.X.csv` and README to include the new s3 URI and instructions
  - [ ] Open PR against test-datasets, taxprofiler branch
- [ ] Add a [MultiQC](https://github.com/multiqc/multiqc) module
- [ ] Make a [Taxpasta](https://github.com/taxprofiler/taxpasta) module
- [ ] Add the database building module to nf-core/createtaxdb (where possible)

### New nf-test procedure and specifications

#### Procedure

When writing a new pipeline-level nf-tests for nf-core/taxprofiler, we recommend the following procedure:

1. Run the test profile locally to have a copy of the expected output files and the results directory structure
2. Check the expected results directories and contents of files are expected

- Check that the results directory reflects parameters specified in the test config itself
- One or two files per directory should be enough

2. Write the base nf-test file structure, assuming _all_ files are stable (following the specifications below)
3. Run `nf-test --tag <test_name> --profile +docker` once to write the first snapshot
4. Run command above to get the `diff` of unstable files
5. Update the assertions in each directories `.match()` snapshot for unstable files

#### Specifications

Write the test files following the specifications below.

The necessary files are follows:

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
    - All other parameters should be specified in the config `.conf` file itself
- Then block
  - [ ] Specify on the first line, a `stable_name_all` variable to list all file names with the nft-utils `getAllFilesFromDir` function
  - [ ] For each top-level output directory under `results` (typically, one per tool), specify a `stable_content_<dir name>` variable (exceptions: `multiqc` and `pipeline_info`) in alphabetical order
    - Use syntax: `def stable_content_<dir name> = getAllFilesFromDir(params.outdir, relative: false, includeDir: false, include: ["<dir name>/**"], ignoreFile: 'tests/.nftignore')`
    - [ ] If no stable files, leave comment for that directory
      - All files unstable `// <dir name>: all unstable files, see stable_name_all`
      - Partly unstable (using custom assertions): `// <dir name>: partly unstable files, see custom assertions`
    - [ ] Make sure `relative: false` in function, to capture md5sums
- `assertAll` block
  - [ ] Use the `removeNextflowVersion` function
  - [ ] Check existance of `nf_core_taxprofiler_software_mqc_versions.yml` file
  - [ ] Check existance of `multiqc_report.html` file
  - [ ] Snapshot `stable_name_all` with a `.match()` name of `all_files`
  - [ ] For each results directory:
    - [ ] Add a comment of the directory name
    - [ ] Specify an assert closure with a snapshot and `.match()` name of `directory_name` (typically the tool name all lower case). Exceptions: MultiQC, Pipeline info.
    - [ ] Each `.match()` should be in alphabetical order (same as order in the `--outdir`)
    - [ ] Inside `.match()` assertion, secify the `stable_content_<dir name>` variable (if available), then:
      - [ ] For each unstable file, specify specific path in `.nftignore` including the tool directory prefix
      - [ ] For each unstable file, if the contents are partly stable, specify an alternative method of file checks (e.g. sorted file, file size check, contains string, nft-plugin function)
      - [ ] For each unstable file assertion, include a string before the assertion itself with the file name and type of check (with what being checked for)
      - [ ] If multiple assertions in snapshot, ensure closing `.match()` and closing `{}` are one and two less indents as the assertions

Reviewing:

- `*.conf` local test run
  - [ ] File contents of all files are as expected
- `*.nf.test`
  - [ ] All test names, tags, profiles correct
  - [ ] Structure matches structure described above
  - [ ] All output directories of the `-profile test_<name>` are covered via a `stable_contents_*` variable and an assertions
  - [ ] Non-stable directories not using `stable_contents_*` replaced with a comment
  - [ ] All unstable files covered in custom snapshots
  - [ ] Custom assertion specifies right file name in both sentence and in the file check itself
- `*.nf.test.snap`
  - [ ] All `.match()` sections defined in `*.nf.test` represented
  - [ ] No empty `match()` blocks
  - [ ] All files in `--outdir` listed in each `.match()` section (except if explicitly excluded due to completely unstable files)
    - Comparing with the output of `tree <--outdir>` can be helpful!
  - [ ] No empty `md5sums` (`d41d8cd98f00b204e9800998ecf8427e`)
  - [ ] No custom boolean assertions set as `false`
