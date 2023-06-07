---
name: Verification Report - Minor and Patch Releases
about: Report for validating minor and patch releases in Balsamic
title: "Verification Report [vX.X.X]"
labels: ["Verification Report"]
assignees: "@Clinical-Genomics/cancer-bioinfo"
---

# Balsamic Verification Report

## Release Information

- Version: [Version Number]
- Release Type: [Minor/Patch]
- Release Commit: [Release Commit]
- Release Date: [Release Date]

## Person(s) involved in the development of new release
- Developer/s: [Name]
- Verification operator/s: [Name]
- Verification reviewer/s: [Name]

## Summary

<!--
Provide a brief summary of the verification performed for the release, including the objectives and scope.
-->

[Summary]

These are the goals of this verification:
1. [Goal 1]
2. [Goal 2]
3. [Goal 3]

## Changelog

<!--
Provide a summary of the changes and updates included in this minor/patch release.
-->

| Type    | Pull Request | Performed tests              |
|---------|--------------|------------------------------|
| Added   | #PR001       | Test [Test 1], Test [Test 2] |
| Changed | #PR002       | Test [Test 3]                |
| Fixed   | #PR003       | Test [Test 2], Test [Test 4] |

## Pre-Verification Checklist

Before proceeding with the verification process, ensure that the following tasks have been completed:

- [ ] [Install Balsamic in stage and production environments in hasta and build its cache](https://atlas.scilifelab.se/infrastructure/BALSAMIC/balsamic/#instructions-for-installation).
- [ ] Confirm the availability of necessary resources, such as test cases.
- [ ] Review the changelog and ensure all changes and updates are documented:

    | Document               | Sections to Be Updated | Pull Request |
    |------------------------|------------------------|--------------|
    | Balsamic Documentation | Sections A, B, and C   | #PR001       |
    | Atlas Documentation    | Section D and E        | #PR002       |

- [ ] Set up the stage environment with the necessary software and configurations:

    | Software | Current Version | Pull Request with the Required Updates |
    |----------|-----------------|----------------------------------------|
    | CG       | Version X.X.X   | #PR001                                 |
    | Servers  | Version Y.Y.Y   | #PR002                                 |
    | Hermes   | Version Z.Z.Z   | #PR002                                 |


## Verification Results

<!--
List the specific test cases that were executed during the verification process. Include the test case ID,
description, and status (Pass/Fail).
-->

### Workflow Integrity Verification Cases

| Case ID        | Analysis type      | Expected QC                                           | Status    |
|----------------|--------------------|-------------------------------------------------------|-----------|
| `acetuna`      | QC, tumor-only     | Pass                                                  | Pass/Fail |
| `civilsole`    | WGS, tumor-only    | Fail (`PCT_60X=0.004532`)                             | Pass/Fail |
| `fleetjay`     | WGS, tumor-normal  | Fail (`PCT_60X=0.00836`)                              | Pass/Fail |
| `setamoeba`    | TGA, tumor-only    | Pass                                                  | Pass/Fail |
| `unitedbeagle` | TGA, tumor-normal  | Pass                                                  | Pass/Fail |
| `uphippo`      | UMI, tumor-only    | Fail (`GC dropout=1.650392`)                          | Pass/Fail |
| `equalbug`     | UMI, tumor-normal  | Fail (`GC_DROPOUT=1.087173` and `RELATEDNESS=-0.524`) | Pass/Fail |

### Version Specific Verification Cases

| Case ID  | Analysis type     | Expected QC | Status    |
|----------|-------------------|-------------|-----------|
| `TC001`  | [Analysis type]   | Pass/Fail   | Pass/Fail |
| `TC002`  | [Analysis type]   | Pass/Fail   | Pass/Fail |


### Test [Test 1]

<!--
Provide detailed results for the specific test, including observations, passing criteria, and any relevant metrics.
-->

| Case ID  | Description   | Expected outcome | Status    |
|----------|---------------|------------------|-----------|
| `TC001`  | [Description] | [Outcome]        | Pass/Fail |
| `TC002`  | [Description] | [Outcome]        | Pass/Fail |
| `TC003`  | [Description] | [Outcome]        | Pass/Fail |

### Issues Identified

<!--
Document any issues or defects identified during the verification process. Include the issue ID, description, severity, 
and status (Open/Closed).
-->

| Failed test   | Description   | Severity        | Status      |
|---------------|---------------|-----------------|-------------|
| Test [Test 1] | [Description] | Low/Medium/High | Open/Closed |
| Test [Test 2] | [Description] | Low/Medium/High | Open/Closed |

### Recommendations

<!--
Provide any recommendations or suggestions for improvement based on the verification results and observations.
-->

### Conclusion

<!--
Summarize the overall outcome of the verification for the minor or patch release. Include any significant findings, 
achievements, or areas requiring further attention.
-->

---

## Approval
- Verification review: [Name]
- Date: [Date]

## Deployment
- Responsible: [Name]
- Date: [Date]
