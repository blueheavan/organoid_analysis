import numpy as np
import pandas as pd

from organoid_analysis.config import load_config
from organoid_analysis.phenotyping.viability import calibrate, classify
from organoid_analysis.statistics.aggregation import make_summaries


def control_rows():
    rows=[]
    for name,c,p in [('live',1050.,50.),('dead',50.,1050.)]:
        for rep in ['R1','R2']:
            rows.append({'sample_id':f'{name}_{rep}','batch_id':'B1','control':name,
                         'condition':f'{name.capitalize()}_control','biological_replicate':rep,
                         'unit_id':f'{name}_{rep}','morphology_eligible':True,'viability_measurement_eligible':True,
                         'viability_measurement_flags':'','calcein_mean_bg_corrected':c,'pi_mean_bg_corrected':p,
                         'calcein_background_noise_mad':5.,'pi_background_noise_mad':5.})
    return pd.DataFrame(rows)


def test_control_scaled_labels_and_both_low_indeterminate():
    cfg=load_config()['viability']
    cfg['mode']='controls'
    data=control_rows()
    reference=data.iloc[0].to_dict()
    tests=[]
    for sid,c,p in [('live_test',950,80),('dead_test',80,950),('mixed_test',550,550),('dim_test',50,50)]:
        tests.append({**reference,'sample_id':sid,'control':'sample','calcein_mean_bg_corrected':c,'pi_mean_bg_corrected':p})
    data=pd.concat([data,pd.DataFrame(tests)],ignore_index=True)
    result=classify(data,calibrate(data,['B1'],cfg),cfg).set_index('sample_id')
    assert result.loc['live_test','viability_state']=='viable_like'
    assert result.loc['dead_test','viability_state']=='compromised_like'
    assert result.loc['mixed_test','viability_state']=='mixed_signal'
    assert result.loc['dim_test','viability_state']=='indeterminate'
    assert result.loc['dim_test','viability_reason']=='both_markers_low'


def test_missing_or_uncalibrated_controls_do_not_generate_confident_states():
    data=control_rows()
    cfg=load_config()['viability']
    result=classify(data,calibrate(data,['B1'],cfg),cfg)
    assert set(result.viability_state)=={'indeterminate'}
    cfg['mode']='controls'
    data=data[data.control=='live']
    assert calibrate(data,['B1'],cfg).iloc[0].status=='unavailable'


def test_multiple_control_images_from_one_biological_replicate_are_not_independent():
    data=control_rows()
    data['biological_replicate']='SAME_DONOR'
    cfg=load_config()['viability']
    cfg['mode']='controls'
    result=calibrate(data,['B1'],cfg).iloc[0]
    assert result.status=='unavailable'
    assert result.live_control_replicates==1


def test_controls_spanning_multiple_conditions_are_rejected():
    """calibrate() must not pool two conditions' same-labeled replicates.

    biological_replicate labels (e.g. "R1") are only unique within a
    condition (see stats.py); a batch whose live (or dead) controls span more
    than one condition must be rejected rather than silently pooled by
    (control, biological_replicate) alone. Regression test for the P1 finding
    from the scientific-software-development-validation audit.
    """
    data = control_rows()
    data.loc[(data.control == 'live') & (data.biological_replicate == 'R1'), 'condition'] = 'LineB_live_control'
    cfg = load_config()['viability']
    cfg['mode'] = 'controls'
    result = calibrate(data, ['B1'], cfg).iloc[0]
    assert result.status == 'unavailable'
    assert result.reason == 'controls_span_multiple_conditions'


def test_badly_separated_control_channels_fail_calibration():
    data=control_rows()
    data['pi_mean_bg_corrected']=55.
    cfg=load_config()['viability']
    cfg['mode']='controls'
    result=calibrate(data,['B1'],cfg).iloc[0]
    assert result.status=='unavailable'
    assert 'pi_controls_not_separated' in result.reason


def test_measurement_qc_blocks_classification():
    data=control_rows()
    test={**data.iloc[0].to_dict(),'control':'sample','sample_id':'saturated','viability_measurement_eligible':False,'viability_measurement_flags':'pi_saturated'}
    data=pd.concat([data,pd.DataFrame([test])],ignore_index=True)
    cfg=load_config()['viability']
    cfg['mode']='controls'
    result=classify(data,calibrate(data,['B1'],cfg),cfg)
    assert result.iloc[-1].viability_state=='indeterminate'
    assert result.iloc[-1].viability_reason=='pi_saturated'


def summary_data():
    metadata=[]
    objects=[]
    # Same donor has many small organoids in one well and one large organoid in another.
    for sid,unit,batch,volume,n in [('F1','W1','B1',100.,90),('F2','W2','B1',10000.,1),('F3','W3','B2',5050.,1)]:
        row={'sample_id':sid,'condition':'A','biological_replicate':'R1','unit_id':unit,'batch_id':batch,'control':'sample','status':'ok'}
        metadata.append(row)
        for i in range(n):
            objects.append({**row,'organoid_id':i+1,'morphology_eligible':True,'volume_um3':volume,
                            'surface_area_um2':100.,'sphericity':.9,'equivalent_diameter_um':10.,'viability_state':'viable_like'})
    metadata.append({'sample_id':'EMPTY','condition':'A','biological_replicate':'R2','unit_id':'W4','batch_id':'B2','control':'sample','status':'ok'})
    return pd.DataFrame(objects),pd.DataFrame(metadata)


def test_wells_are_weighted_equally_and_empty_replicates_are_not_deleted():
    objects,metadata=summary_data()
    samples,units,reps,conditions=make_summaries(objects,metadata,load_config()['report'])
    r1=reps[reps.biological_replicate=='R1'].iloc[0]
    assert r1.median_of_unit_medians_volume_um3==5050.
    assert r1.n_included==92
    assert len(reps)==2  # repeated batches do not create a third independent replicate
    assert conditions.iloc[0].n_biological_replicates_total==2
    assert conditions.iloc[0].n_biological_replicates_with_size_data==1
    assert np.isnan(conditions.iloc[0].ci95_low_volume_um3)
    assert samples[samples.sample_id=='EMPTY'].iloc[0].n_detected==0


def test_incomplete_wells_are_not_silently_included_in_replicate_summaries():
    objects,metadata=summary_data()
    metadata.loc[metadata.sample_id=='F2','status']='failed'
    objects=objects[objects.sample_id!='F2']
    samples,units,reps,conditions=make_summaries(objects,metadata,load_config()['report'])
    r1=reps[reps.biological_replicate=='R1'].iloc[0]
    assert r1.n_complete_units==2
    assert r1.median_of_unit_medians_volume_um3==2575.
