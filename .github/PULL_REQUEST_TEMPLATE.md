## Description

<!-- Provide a brief overview of your PR and link any relevant user stories or issues. -->

[PR description]

#### Added

- [Description]

#### Changed

- [Description]

#### Fixed

- [Description]

#### Removed

- [Description]

## Documentation

<!-- Link all added or updated documents. -->

- [ ] N/A
- [ ] Updated Balsamic documentation to reflect the changes as needed for this PR.
    - [Document Name]

## Tests

<!-- Describe in detail how you tested your changes to help reviewers validate the code. -->
<!-- Include screenshots or visual representations of your changes. -->

### Feature Tests

<!-- Include tests relevant to the changes in this PR. -->

- [ ] N/A
- [ ] Test [Description]
    - [Screenshot]

### Pipeline Integrity Tests

<!-- Include tests to verify the integrity of the different Balsamic workflows. -->

- **Report deliver (generation of the `.hk` file)**
    - [x] N/A
    - [ ] Verified
- **TGA T/O Workflow**
    - [x] N/A
    - [ ] Verified
- **TGA T/N Workflow**
    - [x] N/A
    - [ ] Verified
- **UMI T/O Workflow**
    - [x] N/A
    - [ ] Verified
- **UMI T/N Workflow**
    - [x] N/A
    - [ ] Verified
- **WGS T/O Workflow**
    - [x] N/A
    - [ ] Verified
- **WGS T/N Workflow**
    - [x] N/A
    - [ ] Verified
- **QC Workflow**
    - [x] N/A
    - [ ] Verified
- **PON Workflow**
    - [x] N/A
    - [ ] Verified

## Clinical Genomics Stockholm

<details>
<summary></summary>

<!-- Do not reveal clinical data, and if applicable, place it within the internal Google Drive directory. -->

### Documentation

<!-- Link related issues or PRs for necessary changes. -->

- **Atlas documentation**
    - [x] N/A
    - [ ] Updated: [Link]
- **Web portal for Clinical Genomics**
    - [x] N/A
    - [ ] Updated: [Link]

### User Changes

<!-- Please provide justification if you select N/A and the changes affect the output or results. -->

- [x] N/A
- [ ] This PR affects the output files or results.
    - [ ] User feedback is considered unnecessary because [Justification].
    - [ ] Affected users have been included in the development process and given a chance to provide feedback.

### Infrastructure Changes

<!-- Link related issues or PRs for necessary changes. -->
<!-- Make sure that the linked PRs include relevant tests. -->

- **Stored files in Housekeeper**
    - [x] N/A
    - [ ] Updated: [Link]
- **CG (CLI and delivered/uploaded files)**
    - [x] N/A
    - [ ] Updated: [Link]
- **Servers (configuration files on Hasta)**
    - [x] N/A
    - [ ] Updated: [Link]
- **Scout interface**
    - [x] N/A
    - [ ] Updated: [Link]

</details>

## Checklist

> [!IMPORTANT]  
> Ensure that all checkboxes below are ticked before merging.

### For Developers

- **PR Description**
    - [ ] Provided a comprehensive description of the PR.
    - [ ] Linked relevant user stories or issues to the PR.
- **Documentation**
    - [ ] Verified and updated documentation if necessary.
- **Tests**
    - [ ] Described and tested the functionality addressed in the PR.
    - [ ] Ensured integration of the new code with existing workflows.
    - [ ] Confirmed that meaningful unit tests were added for the changes introduced.
    - [ ] Checked that the PR has successfully passed all relevant code smells and coverage checks.
- **Review**
    - [ ] Addressed and resolved all the feedback provided during the code review process.
    - [ ] Obtained final approval from designated reviewers.

### For Reviewers

- **Code**
    - [ ] Code implements the intended features or fixes the reported issue.
    - [ ] Code follows the project's coding standards and style guide.
- **Documentation**
    - [ ] Pipeline changes are well-documented in the CHANGELOG and relevant documentation.
- **Tests**
    - [ ] The author provided a description of their manual testing, including consideration of edge cases and boundary
      conditions where applicable, with satisfactory results.
- **Review**
    - [ ] Confirmed that the developer has addressed all the comments during the code review.

<!-- Add any other relevant information or specific checks necessary for your PR. -->
