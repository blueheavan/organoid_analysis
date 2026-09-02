from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import tifffile
import analysis.io as analysis_io
from analysis.io import canonical_czyx, ome_spacing, read_manifest, load_sample
from analysis.config import load_config


def test_axis_order_is_explicit_and_exact():
    a=np.arange(2*7*9*11).reshape(2,7,9,11)
    np.testing.assert_array_equal(canonical_czyx(a,'CZYX'),a)
    np.testing.assert_array_equal(canonical_czyx(a.transpose(1,2,0,3),'ZYCX'),a)


@pytest.mark.parametrize('axes,shape',[('QYX',(10,20,20)),('YXS',(20,20,3)),('YX',(20,20)),('ZYX',(1,20,20))])
def test_ambiguous_or_2d_images_are_not_silently_treated_as_volumes(axes,shape):
    with pytest.raises(ValueError):
        canonical_czyx(np.zeros(shape),axes)


def test_multitimepoint_requires_explicit_selection():
    a=np.zeros((2,3,7,9,11))
    a[1]=17
    with pytest.raises(ValueError,match='time_index'):
        canonical_czyx(a,'TCZYX')
    assert (canonical_czyx(a,'TCZYX',1)==17).all()


def test_ome_units_convert_to_micrometres():
    xml='<OME><Image><Pixels PhysicalSizeZ="2000" PhysicalSizeZUnit="nm" PhysicalSizeY="0.001" PhysicalSizeYUnit="mm" PhysicalSizeX="1"/></Image></OME>'
    assert ome_spacing(xml)==(2.,1.,1.)


def test_nonuniform_z_positions_are_rejected():
    xml='<OME><Image><Pixels PhysicalSizeZ="2" PhysicalSizeY="1" PhysicalSizeX="1"><Plane TheZ="0" PositionZ="0"/><Plane TheZ="1" PositionZ="2"/><Plane TheZ="2" PositionZ="7"/></Pixels></Image></OME>'
    with pytest.raises(ValueError,match='uniformly'):
        ome_spacing(xml)


def test_spacing_conflict_is_not_silently_overridden(tmp_path):
    p=tmp_path/'test.ome.tif'
    tifffile.imwrite(p,np.zeros((7,12,12),np.uint16),ome=True,photometric='minisblack',metadata={'axes':'ZYX','PhysicalSizeZ':2.,'PhysicalSizeY':1.,'PhysicalSizeX':1.})
    row={'image_path':str(p),'time_index':'','series_index':'','spacing_z_um':'3','spacing_y_um':'1','spacing_x_um':'1'}
    with pytest.raises(ValueError,match='conflicts'):
        load_sample(row,load_config())


def test_manifest_does_not_inflate_n_by_duplicate_images(tmp_path):
    p=tmp_path/'test.tif'
    tifffile.imwrite(p,np.zeros((7,12,12),np.uint16),photometric='minisblack',metadata={'axes':'ZYX'})
    rows=[{'sample_id':f's{i}','condition':'A','biological_replicate':'R1','image_path':p.name} for i in range(2)]
    manifest=tmp_path/'samples.csv'
    pd.DataFrame(rows).to_csv(manifest,index=False)
    with pytest.raises(ValueError,match='appears twice'):
        read_manifest(manifest)


def test_configuration_typos_fail_instead_of_silently_changing_analysis(tmp_path):
    p=tmp_path/'config.yaml'
    p.write_text('segmentation:\n  min_volum_um3: 500\n')
    with pytest.raises(ValueError,match='Unknown'):
        load_config(p)


def test_separate_paths_cannot_reuse_the_segmentation_channel_for_viability(tmp_path):
    p=tmp_path/'test.ome.tif'
    tifffile.imwrite(p,np.zeros((7,12,12),np.uint16),ome=True,photometric='minisblack',metadata={'axes':'ZYX','PhysicalSizeZ':2.,'PhysicalSizeY':1.,'PhysicalSizeX':1.})
    row={'image_path':str(p),'calcein_path':str(p),'time_index':'','series_index':''}
    with pytest.raises(ValueError,match='same channel'):
        load_sample(row,load_config())


def test_shared_multichannel_source_is_decoded_once_per_sample(tmp_path, monkeypatch):
    p=tmp_path/'channels.ome.tif'
    image=np.zeros((3,7,12,12),np.uint16)
    image[1]=11
    image[2]=22
    tifffile.imwrite(p,image,ome=True,photometric='minisblack',metadata={
        'axes':'CZYX','PhysicalSizeZ':2.,'PhysicalSizeY':1.,'PhysicalSizeX':1.})
    row={'image_path':str(p),'calcein_path':str(p),'pi_path':str(p),
         'calcein_channel':'1','pi_channel':'2','time_index':'','series_index':''}
    original=analysis_io.read_tiff
    calls=[]

    def counted(*args):
        calls.append(args)
        return original(*args)

    monkeypatch.setattr(analysis_io,'read_tiff',counted)
    sample=load_sample(row,load_config())
    assert len(calls)==1
    assert (sample.calcein==11).all()
    assert (sample.pi==22).all()
