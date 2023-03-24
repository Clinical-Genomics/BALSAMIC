---
name: Release Candidate
about: Tracking the next BALSAMIC release
title: ''
---

## Intended Release TAG 
_e.g. 15.x.x_

## Tests carried out:
_All cases to be used in testing. Should cover all workflows. In Bullet points_
- [ ] Balsamic cache is successfully built
- [ ] Verification cases are passed
- [ ] Specific cases are passed

## Testing

### Prepare for test
_All these should be checked of when creating this issue._
- [ ] Login into `hasta`
- [ ] Reserve testing time: `us; paxa -u <user> -s hasta -r balsamic-stage` 
- [ ] Switch user to HC
- [ ] `us`
- [ ] `conda activate S_balsamic`
- [ ] `pip uninstall balsamic`
- [ ] `pip install --no-cache-dir -U git+https://github.com/Clinical-Genomics/BALSAMIC`
- [ ] Release `paxa`

### Generate balsamic cache
- [ ] Download containers: `balsamic init --outdir /home/proj/stage/cancer/balsamic_cache --cosmic-key ${COSMIC_KEY} --genome-version hg19 --run-mode local --snakemake-opt "--cores 1 --until download_container" -r`
- [ ] Download reference files for hg19: `balsamic init --outdir /home/proj/stage/cancer/balsamic_cache --cosmic-key ${COSMIC_KEY} --genome-version hg19 --run-mode local --snakemake-opt "--cores 30" -r`
- [ ] Download reference files for hg38: `balsamic init --outdir /home/proj/stage/cancer/balsamic_cache --cosmic-key ${COSMIC_KEY} --genome-version hg38 --run-mode local --snakemake-opt "--cores 30" -r`
- [ ] Download reference files for canfam3: `balsamic init --outdir /home/proj/stage/cancer/balsamic_cache --cosmic-key ${COSMIC_KEY} --genome-version canfam3 --run-mode local --snakemake-opt "--cores 30" -r`

### Verification cases
_Common verification cases to be used in testing. Should cover all workflows. In Bullet points_
- [ ] fleetjay (TN, WGS): `cg workflow balsamic start fleetjay -r`
  - Expected QC fail: QC metric PCT_60X: 0.00836 
- [ ] civilsole (TO, WGS): `cg workflow balsamic start civilsole -r`
  - Expected QC fail: QC metric PCT_60X: 0.004532
- [ ] ADD_CASE (TN, WES): 
- [ ] livedeer (TO, WES): `cg workflow balsamic start livedeer -r`
- [ ] unitedbeagle (TN, TGA): `cg workflow balsamic start unitedbeagle -r`
- [ ] setamoeba (TO, TGA): `cg workflow balsamic start setamoeba -r`
- [ ] equalbug (TN, UMI): `cg workflow balsamic-umi start equalbug -r --panel-bed gi_cfdna_3.1`
  - Expected QC fail: QC metric GC_DROPOUT: 1.087643 
  - Expected QC fail: QC metric RELATEDNESS: -0.524
- [ ] uphippo (TO, UMI): `cg workflow balsamic-umi start uphippo -r --panel-bed gi_cfdna_3.1`
  - Expected QC fail: QC metric PCT_TARGET_BASES_1000X: 0.861996 
  - Expected QC fail: QC metric GC_DROPOUT: 1.650392 

### Specific cases
_Specific cases to be used in testing. This section should be updated. In Bullet points_
- [ ] 
    
### Creation ToDo
_All these should be checked of when creating this issue. After that this section
can be removed_
- [ ] Add to release milestone
- [ ] Test cases has been identified and added

