import json, inspect, importlib.metadata
from pathlib import Path
from dataclasses import asdict
import numpy as np
import tifffile
from organoid_analysis.microscopy_io.tiff_contract import sha256
from organoid_analysis.segmentation.cellpose_inference import read_stack_multichannel
from organoid_analysis.quantification.features import geometry
out=Path('docs/evidence/2026-09-10')
rows=[]
for p in sorted(Path('data/images').glob('*.tif')):
    row={'file':str(p),'sha256':sha256(p),'size_bytes':p.stat().st_size}
    try:
        s=read_stack_multichannel(p)
        row.update(shape=list(s.shape),axes=s.axes,dtype=str(s.dtype),spacing=asdict(s.spacing),finite=bool(np.isfinite(s.volume).all()),min=float(s.volume.min()),max=float(s.volume.max()))
    except Exception as e: row['error']=repr(e)
    rows.append(row)
(out/'real-data-inventory.json').write_text(json.dumps(rows,indent=2))
from cellpose import models, dynamics
from skimage.measure import marching_cubes
(out/'installed-methods.txt').write_text('\n'.join([f'{p}={importlib.metadata.version(p)}' for p in ['numpy','scipy','scikit-image','tifffile','statsmodels','cellpose','torch','scikit-learn']])+ '\nCellpose.eval '+str(inspect.signature(models.CellposeModel.eval))+'\nnormalize_default='+repr(models.normalize_default)+'\nmarching_cubes '+str(inspect.signature(marching_cubes))+'\n'+inspect.getsource(models.CellposeModel._compute_masks)+'\n'+inspect.getsource(dynamics.compute_masks))
z,y,x=np.indices((25,49,49));r=18.;mask=((z-12)*2)**2+(y-24)**2+(x-24)**2<=r*r
m,_=geometry(mask,(2.,1.,1.))
(out/'analytical-geometry.json').write_text(json.dumps({'oracle':'continuous sphere radius 18 um, rasterized at (2,1,1) um','volume_relative_error':abs(m['volume_um3']/(4*np.pi*r**3/3)-1),'area_relative_error':abs(m['surface_area_um2']/(4*np.pi*r*r)-1),'spec_volume_limit':0.01,'spec_area_limit':0.05,'measurements':m},indent=2))
print('Inventory and analytical probe complete')
