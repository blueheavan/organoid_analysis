from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from organoid_analysis.statistics import exploration, inference


@pytest.mark.parametrize('value', [0., -1., np.inf, -np.inf])
def test_invalid_log_measurement_rejected_before_fitting(value, monkeypatch):
    rows = pd.DataFrame({'condition': ['A'] * 3 + ['B'] * 3,
                         'biological_replicate': ['1', '2', '3'] * 2,
                         'morphology_eligible': True, 'volume_um3': [value, 2, 3, 4, 5, 6]})
    def forbidden(*args):
        pytest.fail('Invalid measurements reached the statistical model')
    monkeypatch.setattr(inference, 'fit_model', forbidden)
    with pytest.raises(ValueError):
        inference.condition_pairwise_tests(rows, features=('volume_um3',))


def test_replicate_identity_does_not_collide_on_separator(monkeypatch):
    rows = pd.DataFrame({'condition': ['A::B'] * 3 + ['A'] * 3,
                         'biological_replicate': ['1', '2', '3', 'B::1', 'B::2', 'B::3'],
                         'morphology_eligible': True, 'volume_um3': [1., 2, 3, 4, 5, 6]})
    class InspectedError(Exception):
        pass
    def inspect(data, feature, group_col):
        assert data[group_col].nunique() == 6
        raise InspectedError
    monkeypatch.setattr(inference, 'fit_model', inspect)
    with pytest.raises(InspectedError):
        inference.condition_pairwise_tests(rows, features=('volume_um3',))


def binary_data():
    rng = np.random.default_rng(123)
    X = rng.normal(size=(60, 2))
    X[::5, 1] = np.nan
    return X, np.repeat([0, 1], 30)


def test_prepare_retains_missingness_until_training_split():
    frame = pd.DataFrame({'condition': ['A', 'A', 'B'], 'v': [1., np.nan, 1000.]})
    X, _, _ = exploration.prepare_binary_data(frame, ['v'], 'condition', 'A', 'B')
    assert np.isnan(X[1, 0])


def test_imputation_and_scaling_are_fit_on_training_rows_and_each_cv_fold(monkeypatch):
    from sklearn import model_selection
    from sklearn.pipeline import Pipeline
    X, y = binary_data()
    X_train, _, _, _ = model_selection.train_test_split(X, y, test_size=.3, random_state=42, stratify=y)
    real_cv = model_selection.cross_val_score
    observed = []
    def audit_cv(estimator, data, labels, **kwargs):
        assert isinstance(estimator, Pipeline)
        assert 'simpleimputer' in estimator.named_steps
        np.testing.assert_array_equal(data, X_train)
        observed.append(True)
        return real_cv(estimator, data, labels, **kwargs)
    monkeypatch.setattr(model_selection, 'cross_val_score', audit_cv)
    summary, _, preprocessing, _, _ = exploration.train_binary_classifiers(X, y, {0: 'A', 1: 'B'})
    np.testing.assert_allclose(preprocessing.named_steps['simpleimputer'].statistics_, np.nanmean(X_train, axis=0))
    assert len(observed) == len(summary)


def test_model_selection_does_not_consult_holdout_performance(monkeypatch):
    from sklearn import metrics, model_selection
    X, y = binary_data()
    counter = iter([.2, .9, .3])
    monkeypatch.setattr(model_selection, 'cross_val_score', lambda *a, **kw: np.array([next(counter)]))
    auc = iter([.99, .1, .2])
    monkeypatch.setattr(metrics, 'roc_auc_score', lambda *a, **kw: next(auc))
    summary, _, _, _, _ = exploration.train_binary_classifiers(X, y, {0: 'A', 1: 'B'})
    assert summary.attrs['best_model'] == 'Random Forest'


def test_entirely_missing_training_feature_is_not_imputed_as_zero():
    X, y = binary_data()
    X[:, 1] = np.nan
    with pytest.raises(ValueError, match='entirely missing'):
        exploration.train_binary_classifiers(X, y, {0: 'A', 1: 'B'})


def test_undefined_effect_sizes_and_assumption_tests_are_not_success():
    from organoid_analysis.statistics.exploration import (
        check_equal_variance,
        check_normality,
        cohens_d,
        effect_size_label,
        field_coefficient_of_variation,
    )
    assert np.isnan(cohens_d(pd.Series([1., 1.]), pd.Series([9., 9.])))
    assert effect_size_label(np.nan) == 'Not estimable'
    small = pd.DataFrame({'group': ['a', 'b'], 'value': [1., 2.]})
    assert check_normality(small, ['value'], 'group').P_value.isna().all()
    assert check_equal_variance(small, ['value'], 'group').P_value.isna().all()
    zero = pd.DataFrame({'group': ['a', 'a'], 'field': [1, 2], 'value': [-1., 1.]})
    assert field_coefficient_of_variation(zero, ['value'], 'group', 'field')['Field_CV_%'].isna().all()


def test_small_group_comparison_and_empty_cluster_abstain():
    from organoid_analysis.statistics.exploration import (
        cluster_characterization,
        compare_two_groups,
    )
    small = pd.DataFrame({'group': ['a', 'a', 'b', 'b'], 'value': [1., 2., 3., 4.]})
    row = compare_two_groups(small, ['value'], 'group').iloc[0]
    assert row.Test_Used == 'NOT ASSESSED'
    assert np.isnan(row.P_value)
    assert pd.isna(row.Significant_Bonferroni)
    clusters = pd.DataFrame({'Cluster': [0, 0, 1], 'value': [1., 2., 3.]})
    row = cluster_characterization(clusters, ['value'], 2).iloc[0]
    assert np.isnan(row.p_value)
    assert np.isnan(row.eta_squared)
