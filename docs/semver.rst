===================
Semantic versioning
===================


BALSAMIC is following `Semantic Versioning <https://semver.org/>`_.

Since October 24, 2018 the following changes were added in addition to SemVer also to cover Bioinformatic and data analysis aspect of it:

- **major**:

  - Structural changes to the BALSAMIC workflow. This includes reordering of annotation softwares or sources, variant
    callers, aligners, quality trimmers, and/or anything other than QC reporting and rule `all`.

  - Addition of annotation softwares or sources, variant callers, aligners, quality trimmers, and/or anything other than
    QC reporting.

- **minor**:

  - Under the hood changes to rules that won't affect output results of workflow.

  - Addition of new bioinfo tools for QC reporting.

  - Updating version of a Bioinformatic software or data resource (including annotation sources)

- **patch**:

  - Any bug fix and under the hood changes that won't impact end-users run.

  - Changes to resource allocation of Scheduler job submission

The rational for versioning is heavily inspired from BACTpipe: DOI: 10.5281/zenodo.1254248 and https://github.com/ctmrbio/BACTpipe)


