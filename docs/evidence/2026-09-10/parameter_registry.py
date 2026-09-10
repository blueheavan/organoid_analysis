from pathlib import Path
from dataclasses import asdict
import json
from organoid_analysis.config import DEFAULTS
from organoid_analysis.segmentation.cellpose_inference import SegmentationConfig
from organoid_analysis.quantification.multilevel_relationships.config import Multilevel3DConfig

rows=[]
def add(name,value,unit,location,impact,origin='HEURISTIC',basis='HEURISTIC',configurable='Yes',status='INSUFFICIENT EVIDENCE'):
 rows.append(dict(parameter=name,value=value,unit=unit,location=location,origin=origin,basis_strength=basis,rationale_and_consequence=impact,calibration_artifact='None',calibration_data_role='NOT ASSESSED',independence_level='NOT ASSESSED',held_out_validation='NOT ASSESSED',sensitivity_range='NOT ASSESSED',uncertainty='NOT ASSESSED',user_configurable=configurable,evidence_status=status))
impacts={
'channels':'Acquisition channel identity; biological stain identity must be supplied, never inferred from index.',
'segmentation':'Changes foreground, splitting, retention, or segmentation review. No project calibration supports the numerical default.',
'quality':'Changes eligibility or marker background/saturation measurements; acquisition-specific support is missing.',
'viability':'Changes calibration availability or signal-state calls; no independent assay validates this default.',
'report':'Changes resampling precision/reproducibility or display/export scope.',
'stats':'Changes inferential inclusion, outcomes or whether inference runs; no predeclared study design validates this default.'}
for section, values in DEFAULTS.items():
 for name,value in values.items():
  unit='dimensionless'
  if name.endswith('_um3'): unit='um^3'
  elif name.endswith('_um'): unit='um'
  elif 'channel' in section: unit='zero-based channel index'
  elif name=='threshold': unit='algorithm or corrected raw intensity'
  elif 'saturation_value' in name: unit='raw detector intensity'
  elif name in ['min_background_voxels']:unit='voxels'
  elif name in ['min_z_slices']:unit='Z slices'
  elif 'replicates' in name:unit='biological replicates'
  elif name in ['bootstrap_iterations']:unit='resamples'
  engineering=section=='report' and name in ('save_meshes','max_meshes_in_preview','seed')
  add(f'{section}.{name}',value,unit,'src/organoid_analysis/config.py::DEFAULTS',impacts[section],origin='ARBITRARY_OR_UNKNOWN' if engineering else 'HEURISTIC',status='NOT APPLICABLE' if engineering else 'INSUFFICIENT EVIDENCE')
for name,value in asdict(Multilevel3DConfig()).items():
 add('multilevel.'+name,value,'voxels' if name=='minimum_voxels' else 'dimensionless','src/organoid_analysis/quantification/multilevel_relationships/config.py','QC review or core/periphery fraction; no validated biological decision boundary.')
for name,value in asdict(SegmentationConfig()).items():
 if 'source' in name or name.startswith('metadata_'):continue
 unit='full-resolution XY pixels' if 'diameter' in name else 'um/pixel' if name=='xy_spacing_um' else 'working-grid pixels (Gaussian sigma)' if name=='flow3d_smooth' else 'dimensionless'
 impact='Alters model scale, foreground or sampling; project calibration missing.'
 if name in ('nuclei_flow_threshold','cell_flow_threshold'):impact='Ignored by do_3D in installed Cellpose; disabled in UI.'
 if name=='xy_spacing_um':impact='Fallback physical calibration, not measured metadata; scales lengths/areas/volumes. No acquisition evidence supports applying 0.414 universally.'
 if name=='anisotropy':impact='Fallback Z/XY calibration 2.9; changes segmentation and downstream Z spacing.'
 if name=='batch_size':impact='Throughput/backend batch choice; numerical equivalence across batch sizes NOT ASSESSED.'
 add('cellpose.'+name,value,unit,'src/organoid_analysis/segmentation/cellpose_inference.py::SegmentationConfig',impact,status='NOT APPLICABLE' if 'flow_threshold' in name else 'INSUFFICIENT EVIDENCE')
extras=[
('cellpose.eval.normalize','True; whole-stack 1st/99th percentiles','percentile','Network-only normalization; raw intensities retained'),
('cellpose.eval.min_size',15,'working voxels','Small-mask removal'),
('cellpose.eval.max_size_fraction',0.4,'volume fraction','Large-mask removal'),
('cellpose.eval.niter','int(200 / (30 / working_diameter))','iterations','Dynamics convergence/runtime'),
('cellpose.eval.resample',True,'boolean','Flows resampled before mask formation'),
('cellpose.eval.augment',False,'boolean','No tile augmentation'),
('cellpose.eval.tile_overlap',0.1,'fraction','Tile overlap changes seam behavior'),
('cellpose.eval.bsize','384 for DINO; 256 for SAM','pixels','Network tile size'),
('cellpose.eval.stitch_threshold',0.,'IoU','2D stitching inactive in this 3D route'),
]
for n,v,u,i in extras:add(n,v,u,'installed cellpose/models.py, transforms.py, dynamics.py; evidence/2026-09-10/installed-methods.txt',i,'SOFTWARE_OR_MODEL_DEFAULT','CONTEXT_DEPENDENT','No (wrapper uses library default)')
extras=[
('diameter_estimate.max_slices',6,'slices'),('diameter_estimate.downsample',0.5,'XY fraction'),('diameter_estimate.quantile',0.5,'quantile'),('diameter_estimate.accepted_range',[5.,250.],'full-resolution pixels'),('diameter_estimate.flow_threshold',0.,'flow error'),('diameter_estimate.cellprob_threshold',0.,'network score'),('diameter_estimate.min_size',5,'pixels'),('diameter_estimate.batch_size',8,'images')]
for n,v,u in extras:add(n,v,u,'src/organoid_analysis/segmentation/parameter_estimation.py','Coarse model-derived size suggestion can affect final chosen scale. Not an independent size measurement.',configurable='API for slices/downsample; otherwise No')
for n,v,u in [('min_cell_volume_um3',200.,'um^3'),('min_nucleus_volume_um3',25.,'um^3'),('require_nucleus',True,'boolean'),('max_nc_ratio',1.,'volume ratio'),('min_nucleus_containment',0.5,'overlap fraction'),('neighbor_radius_um',25.,'um')]:
 add('cells.'+n,v,u,'src/organoid_analysis/quantification/cellular_measurements.py','Changes selected cells/nuclei, pairing, or spatial density; one-to-one biological assumption unvalidated.')
for n,v,u,loc,impact in [
('spacing_agreement.rtol',0.01,'relative','microscopy_io/metadata.py','Acceptance of grid differences; instrument tolerance not calibrated'),
('spacing_agreement.atol',1e-5,'um','microscopy_io/metadata.py','Acceptance of grid differences'),
('square_xy.rtol_atol',[1e-6,1e-9],'relative; um','microscopy_io/voxel_spacing.py','Reject unsupported rectangular XY for Cellpose'),
('high_anisotropy_ratio',5,'max/min spacing','segmentation/watershed_instances.py','Review flag, not a demonstrated reliability boundary'),
('sphericity_review',1.05,'dimensionless','quantification/features.py','Classical geometry exclusion; not present as a multilevel QC rule'),
('MAD_consistency',[1.4826,0.67448975],'normal scale','phenotyping/viability.py; quantification/features.py; quantification/multilevel_relationships/qc.py','Gaussian scaling identity; does not justify biological cutoff'),
('MAD_zero','finite values not isclose to median; NumPy default rtol=1e-5, atol=1e-8','volume units','quantification/multilevel_relationships/qc.py','Uncalibrated fallback outlier flags'),
('viability.noise_floor',1e-9,'intensity','phenotyping/viability.py','Can make SNR arbitrarily large on constant controls; not instrument noise qualification'),
('bootstrap.minimum_N',3,'replicates','statistics/aggregation.py','Availability of percentile intervals, not precision assurance'),
('bootstrap.quantiles',[0.025,0.975],'quantile','statistics/aggregation.py','Nominal 95% interval; coverage NOT ASSESSED'),
('inference.log10_features',['volume_um3'],'feature','statistics/inference.py','Changes estimand to log10 volume'),
('inference.random_variance_ratio',1e-6,'ratio to residual variance','statistics/inference.py','Switches LMM to clustered OLS'),
('inference.reml_optimizers','True; lbfgs then cg','rule','statistics/inference.py','Conditional fit/fallback policy'),
('inference.contrast_df','G-1','replicate clusters','statistics/inference.py','Heuristic MixedLM small-sample reference; not Satterthwaite/Kenward-Roger'),
('inference.alpha_family','0.05; BH per feature; omnibus unadjusted','p value','statistics/inference.py','Does not control multiplicity across all features/omnibus tests'),
('validation.iou_threshold',0.5,'IoU','validation/segmentation_metrics.py','Determines matches; not calibrated to biological error cost'),
('validation.split_merge_overlap',0.1,'object fraction','validation/segmentation_metrics.py','Hints only, not a validated split/merge accuracy measure'),
('exploration.shapiro_cap_seed',[5000,42],'observations; seed','statistics/exploration.py','Deterministic subsampling; does not cure normality-selection bias'),
('exploration.test_selection','Shapiro + Levene p>0.05: pooled t; otherwise Mann-Whitney; N<3 abstains','rule','statistics/exploration.py','Data-driven method selection remains exploratory'),
('exploration.cohen_bands',[0.2,0.5,0.8],'standard deviations','statistics/exploration.py','Conventional descriptive bands, no assay-specific meaning'),
('exploration.iqr_multiplier',1.5,'IQR','statistics/exploration.py','Descriptive outlier flag'),
('exploration.split','test_size=0.3; stratify=y; random_state=42','objects','statistics/exploration.py','No donor/well grouping; biological leakage remains'),
('exploration.CV','min(5, smallest training-class count); shuffle=True; seed=42','folds','statistics/exploration.py','Train-fold preprocessing; object-level accuracy chooses model'),
('exploration.LogisticRegression','max_iter=1000; class_weight=balanced; seed=42','rule','statistics/exploration.py','Untuned classification'),
('exploration.RandomForest_binary','100 trees; max_depth=10; balanced; seed=42','rule','statistics/exploration.py','Untuned classification'),
('exploration.XGBoost','100 trees; max_depth=5; learning_rate=0.1; scale_pos_weight=training negative/positive; seed=42','rule','statistics/exploration.py','Untuned classification; optional backend may be unavailable'),
('exploration.RandomForest_multi_cluster','200 trees; max_depth=15; seed=42; multiclass balanced','rule','statistics/exploration.py','Predicts object/cluster labels, not biological validity'),
('exploration.KMeans','n_init=10; seed=42; k from silhouette max in 2..8 by helper','rule','statistics/exploration.py','Data-derived clustering; circular ANOVA not independent inference'),
('exploration.shape_quantiles',[25,75],'percentiles','statistics/exploration.py','Dataset-relative Rod/Disk/Sphere labels; not calibrated morphotypes'),
('mask_features.rounding',4,'decimal places','quantification/mask_features.py','Rounds exported quantitative values, so not bit-identical to canonical geometry'),
('surface.level_padding_step',[0.5,1,1],'level; voxel; step','quantification/features.py','Defines discrete surface; existing <5% specification fails'),
('preview.mesh_step',2,'voxels','quantification/features.py','Display only; measurement step remains 1'),
('label_export.max_id',2**32-1,'ID','microscopy_io/tiff_contract.py','OME uint32 export range; overflow rejects')]:
 add(n,v,u,'src/organoid_analysis/'+loc,impact,configurable='No unless API/config exposed')
header='''# Scientific parameter registry — 2026-09-10

Values below are recovered from current code, not justified by prior reports. The machine-readable [registry](evidence/2026-09-10/parameters.json) contains every row's calibration, data role, independence, held-out validation, sensitivity and uncertainty fields. Shared defaults for those fields: **no project calibration artifact; calibration role, independent unit, held-out support, scientifically justified sensitivity range and uncertainty are NOT ASSESSED**. This is a limitation, not evidence of adequacy. Numerical regression tests verify behavior, not suitability of a default.

Origin, basis strength and evidence status are separate. `HEURISTIC` denotes a project rule with no calibration record; historical attribution of an exact value to a paper/manufacturer is not invented. Model-library defaults have `SOFTWARE_OR_MODEL_DEFAULT` origin and likewise lack project validation. User configuration does not establish scientific validity. `NOT APPLICABLE` applies only to inactive 3D flow thresholds and display/persistence/seed settings as scientific cutoffs; stochastic reproducibility remains separately assessed.

All scientific values/defaults are unchanged. Input-domain validation now rejects invalid numeric QC values and undefined statistics abstain. Sensitivity probes, when present in the report, are descriptive perturbations only and do not establish a scientifically plausible range.

| Parameter | Current value/rule | Unit | Origin / basis | Rationale and scientific consequence | Configurable | Evidence status |
|---|---|---|---|---|---|---|
'''
lines=[]
for r in rows:
 val=json.dumps(r['value'],ensure_ascii=False)
 lines.append('| `'+r['parameter']+'` | `'+val+'` | '+r['unit']+' | '+r['origin']+' / '+r['basis_strength']+' | '+r['rationale_and_consequence']+' | '+r['user_configurable']+' | '+r['evidence_status']+' |')
footer='''
## Scientific constants and derived values

Unit conversions (m→um 10^6, cm→um 10^4, mm→um 10^3, nm→um 10^-3, angstrom→um 10^-4, inch→um 25400) are exact definitions. Voxel volume is sz*sy*sx; voxel intrinsic covariance is diag(s²/12); ellipsoid full axes are 2*sqrt(5λ); sphericity is π^(1/3)*(6V)^(2/3)/A. Basis ESTABLISHED by dimensional/analytical derivation; controlled arithmetic PASS. The formulas do not validate the input calibration or segmentation boundary.

Runtime control endpoints are CONTROL_DERIVED, using within-batch eligible control-well/replicate medians. Their artifact is `calibration.csv`, role is calibration, and independence is whatever the supplied biological_replicate IDs actually represent. Held-out assay validation and uncertainty are NOT ASSESSED. Image metadata/explicit spacing is acquisition calibration only if independently qualified; defaults 0.414/2.9 are not such evidence. Model-based diameter suggestions are calibration proposals derived from that image and model, never held-out validation.

## Presets, display and omitted dependency settings

`configs/structural_fluorescence.yaml` selects channels 0/1/2 and controls mode; `brightfield_exploratory.yaml` chooses dark polarity; `brightfield_probability.yaml` uses 0.5 probability threshold/controls; `imported_instances.yaml` uses labels and a watershed agreement check. Their remaining values inherit `config.py`; none is independently calibrated. The demo changes smoothing to 0.6 um, closing to 0, minimum volume to 1000 um^3 and seed distance to 12 um, with generator seed 1729. It is not a validation of production defaults or biological gates.

Display percentiles (typically 1/99), uint8 packing, stride limits and renderer opacity/LUT settings affect inspection but not the raw arrays measured by the core pipeline. Display downsampling must not be confused with the scientifically consequential `cellpose.xy_downsample`. Additional estimator/library constructor defaults are version-bound by `pixi.lock`; their complete estimator parameter dictionaries are captured in the dependency-default evidence. No claim is made that those implicit defaults are scientifically calibrated.

See [algorithm decisions](ALGORITHM_DECISIONS.md) for authoritative method/API sources and their claim boundaries; [validation report](VALIDATION_REPORT.md) records readiness. Numeric range checks in configuration are software domain constraints, not experiment-derived acceptance criteria.
'''
Path('docs/PARAMETERS.md').write_text(header+'\n'.join(lines)+'\n'+footer)
Path('docs/evidence/2026-09-10/parameters.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False))
print(len(rows),'parameter records')
