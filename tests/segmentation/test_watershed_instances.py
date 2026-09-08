import numpy as np
import pytest

from organoid_analysis.config import load_config
from organoid_analysis.experimental_statistics.evaluation import match_instances
from organoid_analysis.segmentation.watershed_instances import compare_instance_qc, segment


def test_touching_spheres_split_into_two_instances_in_anisotropic_grid():
    z,y,x=np.indices((40,64,96))
    mask=(((z-19)*2)**2+(y-31)**2+(x-34)**2<16**2)|(((z-19)*2)**2+(y-31)**2+(x-59)**2<16**2)
    cfg=load_config()['segmentation']
    cfg.update(threshold=500.,gaussian_sigma_um=0.,closing_radius_um=0.,min_volume_um3=1000.,seed_min_distance_um=15.,seed_h_um=2.)
    result=segment(mask.astype(np.uint16)*1000,(2.,1.,1.),cfg)
    assert result.labels.max()==2
    cfg['split_touching']=False
    assert segment(mask.astype(np.uint16)*1000,(2.,1.,1.),cfg).labels.max()==1


def test_empty_field_does_not_create_objects():
    result=segment(np.zeros((10,20,20),np.uint16),(2.,1.,1.),load_config()['segmentation'])
    assert not result.labels.any()
    assert 'no_organoids_detected' in result.flags


def test_probability_and_imported_label_paths():
    probability=np.zeros((15,25,25),float)
    probability[3:12,5:20,5:20]=.9
    cfg=load_config()['segmentation']
    cfg.update(method='probability',closing_radius_um=0.,min_volume_um3=1000.)
    result=segment(np.zeros_like(probability),(2.,1.,1.),cfg,probability=probability)
    assert result.labels.max()==1
    cfg['method']='labels'
    imported=(probability>.5).astype(np.uint32)*999999
    result=segment(np.zeros_like(probability),(2.,1.,1.),cfg,imported_labels=imported)
    assert result.original_ids=={1:999999}
    assert result.labels.max()==1


def test_logit_scores_cannot_be_mistaken_for_probabilities():
    cfg=load_config()['segmentation']
    cfg['method']='probability'
    with pytest.raises(ValueError,match=r'\[0,1\]'):
        segment(np.zeros((7,10,10)),(1.,1.,1.),cfg,probability=np.full((7,10,10),2.))


def test_disconnected_imported_instance_is_rejected():
    labels=np.zeros((12,20,20),np.uint16)
    labels[2:5,3:6,3:6]=3
    labels[7:10,13:16,13:16]=3
    cfg=load_config()['segmentation']
    cfg['method']='labels'
    with pytest.raises(ValueError,match='disconnected'):
        segment(labels,(1.,1.,1.),cfg,imported_labels=labels)


def test_instance_matching_detects_merges_despite_perfect_foreground():
    truth=np.zeros((10,20,30),np.uint16)
    truth[2:8,3:15,2:12]=1
    truth[2:8,3:15,17:27]=2
    merged=(truth>0).astype(np.uint16)
    matches,stats=match_instances(truth,merged)
    assert stats['true_positives']==1
    assert stats['false_negatives']==1
    assert stats['recall']==.5
    assert stats['possible_merges']==1


def test_instance_matching_counts_extra_objects_and_is_id_invariant():
    truth=np.zeros((12,24,36),np.uint16)
    truth[2:8,3:12,3:12]=1
    truth[2:8,3:12,20:29]=2
    predicted=truth.astype(np.uint32)*1001
    predicted[2:8,16:22,15:22]=99
    matches,stats=match_instances(truth,predicted)
    assert stats['precision']==pytest.approx(2/3)
    assert stats['recall']==1.
    assert stats['mean_dice_matched']==1.


def test_method_qc_flags_large_count_and_z_extent_disagreement_without_selecting_mask():
    primary=np.zeros((12,24,24),np.uint16)
    reference=np.zeros_like(primary)
    # The primary imitates a fragmented 3D result: three shallow labels rather
    # than one object extending through the reference stack.
    primary[2:4,3:10,3:10]=1
    primary[5:7,3:10,3:10]=2
    primary[8:10,3:10,3:10]=3
    reference[1:11,3:10,3:10]=17
    qc=compare_instance_qc(primary,reference,(2.,1.,1.),.4,.65)
    assert qc['requires_review']
    assert 'object_count_disagreement_review' in qc['qc_flags']
    assert 'primary_z_fragmentation_review' in qc['qc_flags']
    assert qc['primary_instances']==3
    assert qc['reference_instances']==1


def test_instance_qc_handles_sparse_large_ids_without_max_label_allocation():
    from organoid_analysis.segmentation.watershed_instances import instance_qc_summary
    labels=np.zeros((10,12,14),np.uint32)
    labels[1:4,2:5,3:6]=7
    labels[5:10,7:11,8:13]=1_000_000
    summary=instance_qc_summary(labels,(2.,1.,1.))
    assert summary['instances']==2
    assert summary['median_z_extent_um']==pytest.approx(8.)


def test_instance_qc_rejects_negative_labels_and_invalid_spacing():
    from organoid_analysis.segmentation.watershed_instances import instance_qc_summary
    labels=np.zeros((3,4,5),np.int16)
    labels[0,0,0]=-1
    with pytest.raises(ValueError,match='nonnegative'):
        instance_qc_summary(labels,(1.,1.,1.))
    with pytest.raises(ValueError,match='spacing'):
        instance_qc_summary(np.zeros_like(labels),(0.,1.,1.))
