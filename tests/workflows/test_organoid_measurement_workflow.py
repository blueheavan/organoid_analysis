import json
import numpy as np
import pandas as pd
import pytest
import tifffile
import yaml
from organoid_analysis.workflows.organoid_measurement_workflow import analyze, complete_unit_objects
from organoid_analysis.config import load_config


def write_inputs(tmp_path,empty=False):
    z,y,x=np.indices((20,36,36))
    mask=((z-10)*2)**2+(y-18)**2+(x-18)**2<10**2
    data=np.full((3,20,36,36),100,np.uint16)
    if not empty:
        data[0,mask]=1500
        data[1,mask]=1000
        data[2,mask]=130
    image=tmp_path/'stack.ome.tif'
    tifffile.imwrite(image,data,ome=True,photometric='minisblack',metadata={'axes':'CZYX','PhysicalSizeZ':2.,'PhysicalSizeY':1.,'PhysicalSizeX':1.})
    manifest=tmp_path/'manifest.csv'
    pd.DataFrame([{'sample_id':'F1','condition':'Vehicle','biological_replicate':'R1','image_path':image.name}]).to_csv(manifest,index=False)
    cfg=load_config()
    cfg['channels']={'structure':0,'calcein':1,'pi':2}
    cfg['report']['save_meshes']=False
    config=tmp_path/'config.yaml'
    config.write_text(yaml.safe_dump(cfg))
    return manifest,config


def test_end_to_end_outputs_and_unknown_viability_without_controls(tmp_path):
    manifest,config=write_inputs(tmp_path)
    out=tmp_path/'result'
    result=analyze(manifest,out,config)
    assert result['status']=='complete'
    table=pd.read_csv(out/'organoids.csv')
    assert len(table)==1
    assert table.iloc[0].volume_um3>3000
    assert table.iloc[0].viability_state=='indeterminate'
    assert not (out/'RUN_INCOMPLETE.txt').exists()
    for name in ['report.html','size_comparison.svg','morphology_viability.pdf','replicate_summary.csv','calibration.csv','qc/F1.png','labels/F1.labels.ome.tif']:
        assert (out/name).stat().st_size>100
    report=(out/'report.html').read_text()
    assert 'data:image/png;base64,' in report
    provenance=json.loads((out/'provenance.json').read_text())
    assert len(provenance['input_files'][0]['sha256'])==64
    with pytest.raises(FileExistsError):
        analyze(manifest,out,config)


def test_registered_truth_labels_export_instance_validation_tables(tmp_path):
    manifest,config=write_inputs(tmp_path)
    z,y,x=np.indices((20,36,36))
    truth=((((z-10)*2)**2+(y-18)**2+(x-18)**2)<10**2).astype(np.uint16)
    truth_path=tmp_path/'truth.ome.tif'
    tifffile.imwrite(truth_path,truth,ome=True,photometric='minisblack',
                     metadata={'axes':'ZYX','PhysicalSizeZ':2.,'PhysicalSizeY':1.,'PhysicalSizeX':1.})
    design=pd.read_csv(manifest)
    design['truth_labels_path']=truth_path.name
    design.to_csv(manifest,index=False)
    out=tmp_path/'validated_result'
    result=analyze(manifest,out,config)
    assert result['status']=='complete'
    metrics=pd.read_csv(out/'segmentation_validation_metrics.csv')
    assert metrics.iloc[0].sample_id=='F1'
    assert metrics.iloc[0].true_positives==1
    assert metrics.iloc[0].recall==pytest.approx(1.)
    matches=pd.read_csv(out/'segmentation_validation_matches.csv')
    assert set(matches.match_status)=={'TP'}
    report=(out/'report.html').read_text()
    assert 'Annotated segmentation validation' in report


def test_imported_labels_can_be_compared_with_watershed_reference_without_replacement(tmp_path):
    manifest,config=write_inputs(tmp_path)
    z,y,x=np.indices((20,36,36))
    labels=((((z-10)*2)**2+(y-18)**2+(x-18)**2)<10**2).astype(np.uint16)*42
    labels_path=tmp_path/'imported_labels.ome.tif'
    tifffile.imwrite(labels_path,labels,ome=True,photometric='minisblack',
                     metadata={'axes':'ZYX','PhysicalSizeZ':2.,'PhysicalSizeY':1.,'PhysicalSizeX':1.})
    design=pd.read_csv(manifest)
    design['labels_path']=labels_path.name
    design.to_csv(manifest,index=False)
    cfg=yaml.safe_load(config.read_text())
    cfg['segmentation'].update(method='labels',qc_reference_method='watershed')
    config.write_text(yaml.safe_dump(cfg))
    out=tmp_path/'reference_qc_result'
    analyze(manifest,out,config)
    qc=pd.read_csv(out/'segmentation_qc.csv')
    assert qc.iloc[0].reference_method=='watershed'
    assert qc.iloc[0].primary_instances==1
    assert qc.iloc[0].reference_instances>=0
    measured=tifffile.imread(out/'labels/F1.labels.ome.tif')
    assert set(np.unique(measured))=={0,1}


def test_all_empty_field_has_consistent_tables_and_report(tmp_path):
    manifest,config=write_inputs(tmp_path,empty=True)
    out=tmp_path/'empty_result'
    result=analyze(manifest,out,config)
    assert result['n_detected_organoids']==0
    objects=pd.read_csv(out/'organoids.csv')
    assert objects.empty and 'volume_um3' in objects
    samples=pd.read_csv(out/'sample_summary.csv')
    assert samples.iloc[0].n_detected==0
    assert pd.isna(samples.iloc[0].median_volume_um3)
    assert (out/'report.html').is_file()


def test_complete_unit_objects_excludes_incomplete_units():
    """Regression test for the P2 audit finding: statistical testing must not
    silently include organoids from a partially-failed acquisition unit, the
    same exclusion summary.py already applies to descriptive tables."""
    objects = pd.DataFrame([
        {"batch_id": "B1", "unit_id": "W1", "organoid_id": 1},
        {"batch_id": "B1", "unit_id": "W1", "organoid_id": 2},
        {"batch_id": "B1", "unit_id": "W2", "organoid_id": 3},
    ])
    units = pd.DataFrame([
        {"batch_id": "B1", "unit_id": "W1", "unit_complete": True},
        {"batch_id": "B1", "unit_id": "W2", "unit_complete": False},
    ])
    kept = complete_unit_objects(objects, units)
    assert sorted(kept.organoid_id) == [1, 2]


def test_partial_run_marks_failed_field_and_excludes_its_entire_well(tmp_path):
    manifest,config=write_inputs(tmp_path)
    design=pd.read_csv(manifest)
    design['unit_id']='WELL1'
    bad=tmp_path/'projection.tif'
    tifffile.imwrite(bad,np.zeros((36,36),np.uint16),photometric='minisblack',metadata={'axes':'YX'})
    design=pd.concat([design,pd.DataFrame([{'sample_id':'F2','condition':'Vehicle','biological_replicate':'R1',
                                         'unit_id':'WELL1','image_path':bad.name}])],ignore_index=True)
    design.to_csv(manifest,index=False)
    out=tmp_path/'partial_result'
    result=analyze(manifest,out,config,keep_going=True)
    assert result['status']=='partial' and result['n_failed_samples']==1
    assert (out/'PARTIAL_RUN.txt').is_file()
    assert not (out/'RUN_INCOMPLETE.txt').exists()
    reps=pd.read_csv(out/'replicate_summary.csv')
    assert reps.iloc[0].n_complete_units==0
    assert pd.isna(reps.iloc[0].median_of_unit_medians_volume_um3)
