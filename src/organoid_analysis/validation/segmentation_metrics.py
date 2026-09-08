"""One-to-one instance evaluation; merged/extra objects cannot hide in foreground Dice."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def match_instances(truth: np.ndarray, predicted: np.ndarray, iou_threshold: float = .5) -> tuple[pd.DataFrame, dict]:
    """One-to-one instance matching by Hungarian (Kuhn-Munkres) assignment.

    Kuhn, H. W. (1955). The Hungarian method for the assignment problem.
    Naval Research Logistics Quarterly, 2(1-2), 83-97.
    https://doi.org/10.1002/nav.3800020109. See docs/ALGORITHM_DECISIONS.md D9.

    ``iou_threshold`` (default 0.5) is a conventional instance-segmentation
    matching threshold, not independently calibrated for this pipeline; see
    docs/PARAMETERS.md ("Segmentation-validation matching parameter").
    """
    if truth.shape != predicted.shape or truth.ndim != 3:
        raise ValueError("Ground truth and prediction must be matching 3D instance masks")
    for image in [truth, predicted]:
        if not np.issubdtype(image.dtype, np.integer) or (image < 0).any():
            raise ValueError("Evaluation masks must contain nonnegative integer labels")
    if not 0 < iou_threshold <= 1:
        raise ValueError("IoU matching threshold must be in (0,1]")
    tids, ti = np.unique(truth,return_inverse=True)
    pids, pi = np.unique(predicted,return_inverse=True)
    counts = np.bincount(ti.ravel()*len(pids)+pi.ravel(),minlength=len(tids)*len(pids)).reshape(len(tids),len(pids))
    tsel,psel = tids != 0,pids != 0
    intersection = counts[np.ix_(tsel,psel)].astype(float)
    tvol = counts[tsel,:].sum(axis=1).astype(float)
    pvol = counts[:,psel].sum(axis=0).astype(float)
    tids,pids = tids[tsel],pids[psel]
    unions = tvol[:,None]+pvol[None,:]-intersection
    iou = np.divide(intersection,unions,out=np.zeros_like(intersection),where=unions>0)
    dice = np.divide(2*intersection,tvol[:,None]+pvol[None,:],out=np.zeros_like(intersection),where=unions>0)
    matches=[]
    if len(tids) and len(pids):
        # Maximize the number of valid matches first, summed IoU second.
        reward = (iou >= iou_threshold)*(min(len(tids),len(pids))+1)+iou
        tr,pr = linear_sum_assignment(reward,maximize=True)
        matches = [(t,p) for t,p in zip(tr,pr) if iou[t,p] >= iou_threshold]
    matched_t,matched_p = {t for t,p in matches},{p for t,p in matches}
    rows = [{"true_id":int(tids[t]),"predicted_id":int(pids[p]),"match_status":"TP",
             "iou":float(iou[t,p]),"dice":float(dice[t,p])} for t,p in matches]
    rows += [{"true_id":int(tids[t]),"predicted_id":0,"match_status":"FN","iou":0.,"dice":0.}
             for t in range(len(tids)) if t not in matched_t]
    rows += [{"true_id":0,"predicted_id":int(pids[p]),"match_status":"FP","iou":0.,"dice":0.}
             for p in range(len(pids)) if p not in matched_p]
    tp,fn,fp = len(matches),len(tids)-len(matches),len(pids)-len(matches)
    stats = {"iou_threshold":iou_threshold,"n_true":len(tids),"n_predicted":len(pids),"true_positives":tp,
             "false_negatives":fn,"false_positives":fp,
             "precision":tp/(tp+fp) if tp+fp else np.nan,"recall":tp/(tp+fn) if tp+fn else np.nan,
             "detection_f1":2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan,
             "mean_dice_matched":float(np.mean([dice[t,p] for t,p in matches])) if matches else np.nan,
             "mean_iou_matched":float(np.mean([iou[t,p] for t,p in matches])) if matches else np.nan,
             "panoptic_quality":float(sum(iou[t,p] for t,p in matches)/(tp+.5*fp+.5*fn)) if tp+fp+fn else np.nan,
             "possible_splits":int((((intersection / np.maximum(tvol[:,None],1)) >= .1).sum(axis=1)>1).sum()),
             "possible_merges":int((((intersection / np.maximum(pvol[None,:],1)) >= .1).sum(axis=0)>1).sum())}
    return pd.DataFrame(rows,columns=["true_id","predicted_id","match_status","iou","dice"]),stats
