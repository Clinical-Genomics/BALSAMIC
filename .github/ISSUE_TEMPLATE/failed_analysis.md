---
name: Report failed analysis
about: Troubleshooting errors in Balsamic
case-id:
lims-id:
Amount read-pairs:
Sequence type (WES, TGA, WGS):

---

**Error message**
Copy either the error message from the e-mail or copy **only the relevant** line(s) in the error log showing the error message(s). Please do not paste the whole error text.

**Link to the error message**
Copy the path for the error file in Hasta which includes the error message(s). In order to get path to the file, run the following command on Hasta: `ejobinfoprint jobid` (replace job id with actual job id).

**Reason for failure**
Provide any theories regarding why the error message arises.

**If workflow, which rules**
If possible, and using the Snakemake workflows, the name of the affected rules and workflows.

**Version (please complete the following information):**
`balsamic --version` output

**Additional context**
Add any other context about the problem here. For e.g. amount sequencing reads or other parameters that are of interest.
