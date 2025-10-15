#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# =============================================================================
#### Modules
# =============================================================================

import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from scipy.stats import f, kruskal, norm, chi2, t
from scipy import stats
from pandas.core.window.rolling import Rolling

proj_dir = os.path.dirname( os.path.abspath('__file__') )
func_dir = os.path.join( proj_dir, 'func' )
sys.path.append( func_dir )

from util import (get_all_cases_to_process_for_regional_maps_or_plumes_or_X11,
                  path_to_fill_to_where_to_save_satellite_files,
                  fill_the_sat_paths, load_csv_files_in_the_package_folder, define_parameters)


# =============================================================================
#### Utility functions
# =============================================================================

def apply_X11_method_and_save_results(values, variable_name, dates, info, X11_dir_out):

    filtered_values, filtered_dates = keep_the_dates_within_full_year(values, pd.to_datetime(dates))

    if info.Temporal_resolution == 'WEEKLY':
        filtered_values = pd.Series(filtered_values).rolling(3, center=True, min_periods=1).median().to_numpy()

    results = temporal_decomp_V2_7_x11(values=filtered_values, dates=filtered_dates,
                                       time_frequency=info.Temporal_resolution,
                                       filter_outlier=False, overall_cutoff=50,
                                       out_limit=3, perc_month_limit=50,
                                       var_stationary=False, lin_interpol=False,
                                       cutoff_fill=30, season_test=True)

    fig, axs = plt.subplots(4, 1, figsize=(24, 14))

    axs[0].plot(results['7_dates'], results['8_values_ini'])
    axs[0].set_title('Initial values')
    axs[1].plot(results['7_dates'], results['10_Interannual_signal'])
    axs[1].set_title(f'Inter-annual signal ({round(results["1_variance_due_to_Interannual_signal"], 1)}%)')
    axs[2].plot(results['7_dates'], results['9_Seasonal_signal'])
    axs[2].set_title(f'Seasonal signal ({round(results["0_variance_due_to_Seasonal_signal"], 1)}%)')
    axs[3].plot(results['7_dates'], results['11_Residual_signal'])
    axs[3].set_title(f'Residual signal ({round(results["2_variance_due_to_Residual_signal"], 1)}%)')

    # Adjust layout
    plt.tight_layout()

    folder_where_to_store_the_plot = os.path.join(X11_dir_out, info.Zone, 'X11_ANALYSIS', variable_name)
    os.makedirs(folder_where_to_store_the_plot, exist_ok=True)
    file_name = "_".join(info.drop(["Zone", "Year", "Satellite_variable"]).astype(str).values)

    plt.savefig(folder_where_to_store_the_plot + "/" + file_name + '.png')

    plt.close(fig)

    pd.DataFrame({'dates': results['7_dates'],
                  'Raw_signal': results['8_values_ini'],
                  'Interannual_signal': results['10_Interannual_signal'],
                  'Seasonal_signal': results['9_Seasonal_signal'],
                  'Residual_signal': results['11_Residual_signal'],
                  'Variation_coefficient': results['5_var_coeff'],
                  'Variance_due_to_Interannual_signal': results['1_variance_due_to_Interannual_signal'],
                  'Variance_due_to_Seasonal_signal': results['0_variance_due_to_Seasonal_signal'],
                  'Variance_due_to_Residual_signal': results['2_variance_due_to_Residual_signal'],
                  'Monotonic_change': results['18_Kendall_Sen_analyses_on_Interannual_signal'][
                      'Is_the_change_of_annual_values_monotonic'],
                  'Rate_of_Change': results['18_Kendall_Sen_analyses_on_Interannual_signal'][
                      'Rate_of_change_of_annual_values_in_percentage_per_time']
                  }).to_csv(folder_where_to_store_the_plot + "/" + file_name + '.csv', index=False)


def keep_the_dates_within_full_year(values, dates):

    min_date = dates.min()
    max_date = dates.max()

    if any(dates <= min_date.replace(month=1, day=1)):
        start_date = min_date.replace(month=1, day=1)
    else:
        start_date = (min_date + pd.DateOffset(years=1)).replace(month=1, day=1)

    if any(dates >= max_date.replace(month=12, day=31)):
        end_date = max_date.replace(month=12, day=31)
    else:
        end_date = (max_date - pd.DateOffset(years=1)).replace(month=12, day=31)

    # Filter values within the determined date range
    index_to_keep = [i for i, d in enumerate(dates) if start_date <= d <= end_date]

    filtered_values = np.array(values)[index_to_keep]
    filtered_dates = np.array(dates)[index_to_keep]

    return filtered_values, filtered_dates


def generate_lagged_matrix(series, p):
    n = len(series)
    lagged_matrix = np.zeros((n - p, p))
    for i in range(p):
        lagged_matrix[:, i] = series[p - i - 1:n - i - 1]
    return lagged_matrix


def autocorrelation(Z, lags, Wt0=None):
    Autocov = np.zeros_like(lags, dtype=float)
    # rk = np.zeros_like(lags, dtype=float)
    if isinstance(Wt0, np.ndarray) == False:
        Wt0 = np.ones(len(Z))

    # Autocovariance calculation for the FIRST ITERATION
    for i, lag in enumerate(lags[:-1]):

        if i == 0:
            Zt = Z
            Zt_K = Z
            Wt_lag = Wt0
            Wt_lag_k = Wt0
        else:
            Zt = Z[:-lag]
            Zt_K = Z[lag:]
            Wt_lag = Wt0[:-lag]
            Wt_lag_k = Wt0[lag:]

        Wt_global = Wt_lag + Wt_lag_k - 1
        # Z_mean_k = np.nansum(Zt * Wt_global) / np.nansum(Wt_lag)
        Z_mean = np.nansum(Zt_K * Wt_global) / np.nansum(Wt_lag_k)

        Autocov[i] = (1.0 / np.nansum(Wt_lag)) * np.nansum(Wt_lag * (Zt - Z_mean) * (Zt_K - Z_mean))
        # rk[i] = 1.96 / np.sqrt(N - lag - 1)
        # Autocov[i] = 1.0 / (N - 1) * np.sum(Wt_lag * (Zt - Z_mean) * (Zt_K - Z_mean_k))

    # Auto Correlation Calculation (considering missing values)
    Autocor = Autocov / Autocov[0]

    return Autocor


def perform_EVF_computation(Z_ini, Z_act, lags, ind_MV, Wt0, ACP_cutoff, level_SD_perc, max_iter, it, test_end,
                            first_iteration):
    N = len(Z_ini)
    dt = 1.0
    ind_Z = np.arange(N) + 1

    if first_iteration == False:
        Wt0 = np.ones(N)

    Autocor = autocorrelation(Z_act, lags, Wt0)

    # Confidence interval of the autocorrelation function
    rk = 1.96 / np.sqrt(N - lags - 1)
    ind_p = np.where(abs(Autocor) < rk)[0]

    # Number of lagged series for PCA
    p = ind_p[0] + 1 if len(ind_p) > 0 else 0

    # Exception to ignore the case of P = 0
    if p == 0:

        if first_iteration:

            Z_ini[ind_MV] = np.nan
            return Z_ini

        else:

            Z_it = np.copy(Z_ini)
            Z_it[ind_MV] = Z_act[ind_MV]

            test_val = np.sqrt(np.sum((Z_it[ind_MV] - Z_act[ind_MV]) ** 2)) / len(ind_MV)

            if test_val < np.std(Z_ini[np.isfinite(Z_ini)]) / level_SD_perc:
                test_end = True

            if it >= max_iter:
                Z_act[ind_MV] = np.nan  # NaN for values where convergence did not occur
                test_end = True

            return [Z_it, it, test_end]

    # P-lagged Matrix construction
    matrix = np.full((p, N - (p - 1) * int(dt)), np.nan)
    matrix_ind_Z = np.copy(matrix)

    # Fill lagged matrices
    index_max = N - (p - 1) * int(dt)
    matrix[0, :] = Z_act[:index_max]
    matrix_ind_Z[0, :] = ind_Z[:index_max]

    M_weight = np.copy(matrix)
    M_weight[0, :] = Wt0[:index_max]

    for i in np.arange(1, p):
        index_max = N - (p - 1 - i) * int(dt)
        matrix[i, :] = Z_act[i: index_max]
        matrix_ind_Z[i, :] = ind_Z[i: index_max]
        M_weight[i, :] = Wt0[i: index_max]

    if (first_iteration == False) and np.isnan(matrix).any():
        sys.exit('There should be no missing data in the second iteration')

    matrix[np.isnan(matrix)] = np.nanmean(Z_act)

    matrix = np.flip(matrix, axis=1)
    matrix_ind_Z = np.flip(matrix_ind_Z, axis=1)
    M_weight = np.flip(M_weight, axis=1)

    # Calculation of the centered matrix
    matrix_std = np.full_like(matrix, np.nan)
    corr_mean = np.full(p, np.nan)

    for i in range(p):
        corr_mean[i] = np.nansum(matrix[i, :] * M_weight[i, :]) / np.nansum(M_weight[i, :])
        matrix_std[i, :] = matrix[i, :] - corr_mean[i]

    # Covariance matrix for PCA
    covar = np.full((p, p), np.nan)

    for i in range(p):
        for j in range(p):
            x1 = matrix_std[i, :]
            x2 = matrix_std[j, :]
            valid_indices = np.where((np.isfinite(x1)) & (np.isfinite(x2)))[0]
            covar[i, j] = np.cov(x1[valid_indices], x2[valid_indices])[0, 1]

    # Perform PCA
    eigen_val, eigen_vector = np.linalg.eigh(covar)
    eigen_val = np.flip(eigen_val)
    eigen_vector = np.flip(np.transpose(eigen_vector), axis=1)
    exp_var = eigen_val / np.sum(eigen_val) * 100.0

    # (np.dot(covar, eigen_vector[:, 0]) - eigen_val[0] * eigen_vector[:, 0]) # Check if eigen_val and eigen_vectors objects are correct. The results hsould be zero.

    # Cumulative variance calculation
    cumul = np.cumsum(exp_var)
    Test_axes = np.where(cumul > ACP_cutoff)[0]
    N_axes = Test_axes[0] + 1 if len(Test_axes) > 0 else 1

    # Components calculation

    # Ctk = np.full((matrix_std.shape[1], len(eigen_vector)), -9999.0)
    # for t in np.arange(matrix_std.shape[1]) :
    #     for k in np.arange(len(eigen_val)) :
    #         a = []
    #         for the_p in np.arange(p) :
    #             a.append( M_weight[the_p, t] * eigen_vector[the_p, k] * matrix_std[the_p, t] )
    #         Ctk[t,k] = np.sum(a)

    Ctk = np.dot(np.transpose(matrix_std * M_weight), eigen_vector)
    Y_p_n_axes = np.full((p, Ctk.shape[1], N_axes), np.nan)

    for i in range(N_axes):
        Y_p_n_axes[:, :, i] = np.transpose(np.dot(Ctk[i, :].reshape(-1, 1), eigen_vector[i, :].reshape(1, -1)))

    Y_p = np.nansum(Y_p_n_axes, axis=2) if N_axes > 1 else Y_p_n_axes[:, :, 0]
    Y_p += np.repeat(corr_mean.reshape(-1, 1), Y_p.shape[1], axis=1)

    if first_iteration:
        # New vector creation after the first iteration
        Z_i = np.full(N, np.nan)
    else:
        Z_i = np.arange(N).astype(float)

    for i in range(N):

        ind = np.where(matrix_ind_Z == (ind_Z)[i])[0]

        if first_iteration == False:

            if len(ind) <= 0:
                sys.exit('Error: no valid indices found.')
                # continue

            if np.isnan(np.mean(Y_p[ind])):
                sys.exit('Error: NaN values present.')

            Z_i[i] = np.nanmean(Y_p[ind])

        else:

            Z_i[i] = np.nanmean(Y_p[ind]) if len(ind) > 0 else np.nan

    Z_it = np.copy(Z_ini)
    Z_it[ind_MV] = Z_i[ind_MV]

    if first_iteration:

        return Z_it

    else:

        it += 1
        test_val = np.sqrt(np.nansum((Z_it[ind_MV] - Z_act[ind_MV]) ** 2)) / len(ind_MV)

        if test_val < np.nanstd(Z_ini) / level_SD_perc:
            test_end = True

        if it >= max_iter:
            Z_it[ind_MV] = np.nan  # NaN for values where convergence did not occur
            test_end = True

        return [Z_it, it, test_end]


def F_EVF_V1_2(X, level_SD_perc, ACP_cutoff, max_iter):
    """
    PURPOSE:
    Fill the gaps present in a time series using the method defined in Ibanez and Conversi (2002).

    INPUTS:
    X: Original time series with gaps
    level_SD_perc: Threshold decision value for deciding the end of the iteration
    ACP_cutoff: Cutoff value of the cumulated variance used for selecting the number of axes of the PCA
    max_iter: Maximum number of iterations

    OUTPUTS:
    The original time series with gaps replaced by estimated values
    """

    Z = np.copy(X)
    ind_MV = np.where(np.isnan(Z))[0]  # Indices of missing values

    if len(ind_MV) == 0:
        return Z

    # Initial weight setup: 1 where not missing, 0 where missing
    Wt0 = np.ones(len(Z))
    Wt0[ind_MV] = 0

    # Definition of the lag for autocorrelation
    lags = np.arange(len(Z) // 3)

    Z_act = perform_EVF_computation(Z_ini=Z, Z_act=Z,
                                    lags=lags, ind_MV=ind_MV, Wt0=Wt0,
                                    ACP_cutoff=ACP_cutoff, level_SD_perc=level_SD_perc, max_iter=max_iter,
                                    it=0, test_end=False, first_iteration=True)

    test_end = False
    it = 0

    while test_end == False:
        [Z_act, it, test_end] = perform_EVF_computation(Z_ini=Z, Z_act=Z_act,
                                                        lags=lags, ind_MV=ind_MV, Wt0=np.ones(len(Z)),
                                                        ACP_cutoff=ACP_cutoff, level_SD_perc=level_SD_perc,
                                                        max_iter=max_iter,
                                                        it=it, test_end=test_end, first_iteration=False)

    return Z_act


def F_census_X_11_pezzulli_V1_2(Xt, period):
    # Test if the period is even or odd
    if period % 2 != 0:
        test_period = 0  # period is odd
    else:
        test_period = 1  # period is even

    Xt = np.array(Xt)  # Reform the input into a NumPy array
    N = len(Xt)

    # Step 1
    if test_period == 1:
        Tt_1 = F_moving_average_v1_0(Xt, 2, period, period)
    else:
        Tt_1 = F_moving_average_v1_0(Xt, 1, period, period)

    Zt_1 = Xt - Tt_1

    SMA_22_1 = np.full(N, np.nan)

    for i in range(period):
        t_SMA_22_1 = F_moving_average_v1_0(Zt_1[i::period], 2, 2, 1)
        SMA_22_1[np.arange(len(t_SMA_22_1)) * period + i] = t_SMA_22_1

    # Normalisation
    if test_period == 1:
        St_1 = SMA_22_1 - F_moving_average_v1_0(SMA_22_1, 2, period, period)
    else:
        St_1 = SMA_22_1 - F_moving_average_v1_0(SMA_22_1, 1, period, period)

    # Step 2
    Yt_2 = Xt - St_1
    Tt_2 = F_henderson_v1_0(X=Yt_2, order=2 * period - 1, period=period)

    Zt_2 = Xt - Tt_2

    SMA_22_2 = np.full(N, np.nan)

    for i in np.arange(period):
        t_SMA_22_2 = F_moving_average_v1_0(Zt_2[i::period], 2, 2, 1)
        # index_to_fill = np.array( [x for x in np.arange(len(t_SMA_22_2)) * period + i if x <= len(SMA_22_2)] )
        index_to_fill = np.array([x for x in np.arange(len(t_SMA_22_2)) * period + i])
        SMA_22_2[index_to_fill] = t_SMA_22_2

    # Normalisation
    if test_period == 1:
        St_2 = SMA_22_2 - F_moving_average_v1_0(SMA_22_2, 2, period, period)
    else:
        St_2 = SMA_22_2 - F_moving_average_v1_0(SMA_22_2, 1, period, period)

    # Step 3
    Yt_3 = Xt - St_2
    Tt_3 = F_henderson_v1_0(Yt_3, 2 * period - 1, period)

    Yt = Xt - Tt_3 - St_2

    Tt = Tt_3
    St = St_2

    decomposed = {'Tt': Tt, 'St': St, 'Yt': Yt}

    return decomposed


def F_moving_average_v1_0(X, m, p, period):
    period_calc = period
    p_calc = p

    if p_calc % 2 == 0:
        p_calc += 1
    half_window = (p_calc - 1) / 2

    X = np.array(X)
    N = len(X)

    if period_calc == 1:
        if p == 2:
            period_calc = 1
        elif p == 3:
            period_calc = 2
        elif p == 5:
            period_calc = 3

    # Calculate indices for slicing
    start1 = int((period_calc - 1) - half_window + 1)
    end1 = int(period_calc)

    start2 = int(N - 1 - period_calc + 1)
    end2 = int(start2 + half_window)

    # Slice the array and concatenate
    X_calc = np.concatenate([
        X[start1:end1],
        X,
        X[start2:end2]
    ])

    N_calc = len(X_calc)
    MA_X = np.arange(N_calc).astype(float)

    if m % 2 == 0:

        for i in np.arange(half_window, N_calc - (half_window)):
            MA_X[int(i)] = (m / 2 * X_calc[int(i - half_window)] +
                            np.nansum(m * X_calc[int((i - half_window) + 1): int(i - 1 + half_window + 1)]) +
                            m / 2 * X_calc[int(i + half_window)]) / ((p_calc - 1) * m)

        MA_X = MA_X[int(half_window): int(N + half_window)]

    else:

        if m == 1:

            MA_X = pd.Series(X_calc).rolling(window=p_calc, center=True, min_periods=1).mean().to_numpy()
            MA_X = MA_X[int(half_window): int(N + half_window)]

        else:

            weight = np.concatenate([np.arange(2) + 1, np.full(p_calc - 2, m), np.arange(2)[::-1] + 1])
            X_calc = np.concatenate([X[int((period_calc - 1) - (p_calc + 2 - 1) / 2 + 1): int((period_calc - 1) + 1)],
                                     X,
                                     X[int(N - 1 - period_calc + 1): int(
                                         (N - 1 - period_calc + 1) + ((p_calc + 2 - 1) / 2 - 1)) + 1]])
            N_calc = len(X_calc)
            MA_X = np.zeros(N_calc)
            for i in range((p_calc + 1) / 2, N_calc - ((p_calc + 1) / 2)):
                MA_X[i] = np.sum(X_calc[int(i - (p_calc + 1) / 2): int(i + (p_calc + 1) / 2 + 1)] * weight) / np.sum(
                    weight)
            MA_X = MA_X[int(half_window + 1): int(N + half_window + 1)]

    return MA_X


def F_henderson_v1_0(X, order, period):
    half_order = (order - 1) / 2

    N = len(X)
    X_calc = np.concatenate([X[int((period - 1) - half_order + 1): int(period)],
                             X,
                             X[int(N - 1 - period + 1): int((N - 1 - period + 1) + (half_order - 1) + 1)]])
    N_calc = len(X_calc)
    MA_X = np.zeros(N_calc)

    if order % 2 == 0:
        raise ValueError(f"ORDER MUST BE ODD {order}")

    i_H = np.concatenate([np.arange(-half_order, 0), np.array([0]), np.arange(half_order) + 1])
    n_H = float((len(i_H) - 1) / 2 + 2)
    h = (315 * ((n_H - 1) ** 2 - i_H ** 2) * (n_H ** 2 - i_H ** 2) * ((n_H + 1) ** 2 - i_H ** 2) *
         (3 * n_H ** 2 - 16 - 11 * i_H ** 2) / (8 * n_H * (n_H ** 2 - 1) * (4 * n_H ** 2 - 1) *
                                                (4 * n_H ** 2 - 9) * (4 * n_H ** 2 - 25)))

    for i in range(int(half_order), int(N_calc - (half_order))):
        MA_X[i] = np.sum(X_calc[int(i - half_order): int(i + half_order + 1)] * h) / np.sum(h)

    MA_X = MA_X[int(half_order): int(N + half_order)]

    return MA_X


def f_test_season_v1_1(X, T, S, period, test_on_X=False, test_on_S=False, n_years=None):
    N_years = len(X) // period

    RES_TEST = np.nan
    S_evolutive = 0.0

    S_I = X - T
    if test_on_S:
        S_I = S

    # Create arrays for the calculations
    M_S_I = np.full((period, N_years), np.nan)
    M_X = np.full((period, N_years), np.nan)

    # Populate the arrays
    for i in range(period):
        M_S_I[i, :len(X[i::period])] = S_I[i::period]
        M_X[i, :len(X[i::period])] = X[i::period]

    ind_bad = np.where(np.isnan(M_S_I))
    ind_valid = np.where(np.isfinite(M_S_I))

    N_col_M_S_I, N_lign_M_S_I = M_S_I.shape

    N_col_M_X, N_lign_M_X = M_X.shape

    if test_on_X:

        # Test on X
        S_A2_X_i = np.full(N_col_M_X, -999.0)
        Mean_X_i = np.full(N_col_M_X, -999.0)

        for i in np.arange(N_col_M_X):
            count_n_i = np.sum(np.isfinite(M_X[i, :]))
            if count_n_i > 0:
                Mean_X_i[i] = np.nanmean(M_X[i, :])
                S_A2_X_i[i] = count_n_i * (Mean_X_i[i] - np.nanmean(M_X[ind_valid])) ** 2

        S_R2_X = np.nansum((M_X - np.tile(Mean_X_i[i], (N_lign_M_X, 1)).T) ** 2)
        S_A2_X = np.nansum(S_A2_X_i)

        Fs = (S_A2_X / (N_col_M_X - 1)) / (S_R2_X / (len(ind_valid[0]) - N_col_M_X))

        dfn = N_col_M_X - 1
        dfd = len(ind_valid[0]) - N_col_M_X

        Limit_FS = 0.1 / 100
        test_Fs = 1 - f.cdf(Fs, dfn, dfd)

    else:

        # Test on M_S_I
        S_A2_X_i = np.full(N_col_M_S_I, -999.0)
        Mean_MSI_i = np.full(N_col_M_S_I, -999.0)

        for i in range(N_col_M_S_I):

            count_n_i = np.where(np.isfinite(M_S_I[i, :]))[0]
            if len(count_n_i) > 0:
                mean_MSI_i = np.nanmean(M_S_I[i, :])
                S_A2_X_i[i] = len(count_n_i) * (mean_MSI_i - np.nanmean(M_S_I[ind_valid])) ** 2

        S_R2_X = np.nansum((M_S_I - np.tile(Mean_MSI_i, (N_lign_M_S_I, 1)).T) ** 2)
        S_A2_X = np.nansum(S_A2_X_i)

        Fs = (S_A2_X / (N_col_M_S_I - 1)) / (S_R2_X / (len(ind_valid[0]) - N_col_M_S_I))
        dfn = N_col_M_S_I - 1
        dfd = len(ind_valid[0]) - N_col_M_S_I

        Limit_FS = 0.1 / 100
        test_Fs = 1 - f.cdf(Fs, dfn, dfd)

    if test_Fs < Limit_FS:

        abs_M_S_I = np.abs(M_S_I - np.nanmean(M_S_I[ind_valid]))

        if len(ind_bad[0]) != 0:
            abs_M_S_I[ind_bad] = np.nan

        # S_2 = np.nanvar(abs_M_S_I) * len(ind_valid[0])
        X_p_p = np.nanmean(abs_M_S_I)

        X_i_p = np.arange(N_col_M_S_I)
        X_p_j = np.arange(N_lign_M_S_I)

        for i in range(N_col_M_S_I):
            valid_indices = np.isfinite(abs_M_S_I[i, :])
            if np.any(valid_indices):
                X_i_p[i] = np.nanmean(abs_M_S_I[i, valid_indices])

        for j in range(N_lign_M_S_I):
            valid_indices = np.isfinite(abs_M_S_I[:, j])
            if np.any(valid_indices):
                X_p_j[j] = np.nanmean(abs_M_S_I[valid_indices, j])

        # S_A_2 = N_lign_M_S_I * np.sum((X_i_p - X_p_p)**2)
        S_B_2 = N_col_M_S_I * np.sum((X_p_j - X_p_p) ** 2)
        S_R_2 = np.sum(
            (abs_M_S_I - np.tile(X_i_p, (N_lign_M_S_I, 1)).T - np.tile(X_p_j, (N_col_M_S_I, 1)) + X_p_p) ** 2)

        FM = (S_B_2 / (N_lign_M_S_I - 1)) / (S_R_2 / ((N_lign_M_S_I - 1) * (N_col_M_S_I - 1)))
        dfn = N_lign_M_S_I - 1
        dfd = (N_col_M_S_I - 1) * (N_lign_M_S_I - 1)

        Limit_Fm = 5.0 / 100
        test_Fm = 1 - f.cdf(FM, dfn, dfd)

        T1 = 7.0 / Fs
        T2 = 3.0 * FM / Fs
        T = np.sqrt(0.5 * (T1 + T2))

        if test_Fm > Limit_Fm:
            if T1 >= 1 or T2 >= 1:
                RES_TEST = 1.0
            else:
                Limit_KW = 0.001
                test_KW_M_X = kw_TEST(M_S_I, df=period - 1, MISSING=-1.0)

                if test_KW_M_X[1] > Limit_KW:
                    RES_TEST = 1.0
                else:
                    RES_TEST = 2.0
        else:
            if T >= 1:
                RES_TEST = 0.0
            else:
                S_evolutive = 1.0

                if T1 >= 1 or T2 >= 1:
                    RES_TEST = 1.0
                else:
                    Limit_KW = 0.1 / 100
                    test_KW_M_X = kw_TEST(M_S_I, df=period - 1, MISSING=-1.0)

                    if test_KW_M_X[1] > Limit_KW:
                        RES_TEST = 1.0
                    else:
                        RES_TEST = 2.0

    else:
        RES_TEST = 0.0

    if RES_TEST != 2.0:
        res_test_season = RES_TEST
    else:
        if S_evolutive == 0:
            res_test_season = 2.0
        else:
            res_test_season = 3.0

    # return res_test_season, S_evolutive
    return res_test_season == 0.0


def kw_TEST(data, df, MISSING):
    """
    Perform the Kruskal-Wallis H-test for independent samples.

    Parameters:
    - data: 2D array-like, where each row represents a group (season).
    - df: Degrees of freedom (number of groups - 1).
    - MISSING: Value used to denote missing data.

    Returns:
    - H: Kruskal-Wallis H statistic.
    - p_value: p-value for the test.
    """

    # Replace missing values with NaN
    data = np.where(data == MISSING, np.nan, data)

    # Split data into groups
    groups = [data[i, ~np.isnan(data[i, :])] for i in range(data.shape[0])]

    # Perform Kruskal-Wallis test
    H, p_value = kruskal(*groups)

    return [H, p_value]


def f_kendall_sen_v1_0(X, alpha=0.05, correction_tie=False):
    """
    Test for monotonic change in a time series using the non-parametric Kendall statistics
    and compute the Sen slope estimator.

    Parameters:
        X (array-like): Input time series vector
        hbad (float, optional): Bad value definition, default is -9999
        alpha (float, optional): Significance level, default is 0.05
        correction_tie (int, optional): If set, correct the calculation for the presence of tied data, default is 0

    Returns:
        tuple: (sen_slope, prob, res_kendall_test, b, RC)
    """

    output = {
        'Is_the_change_of_annual_values_monotonic': np.nan,
        'Sen_slope': np.nan,
        'Sen_intercept': np.nan,
        'Rate_of_change_of_annual_values_in_percentage_per_time': np.nan
    }

    # Prepare data
    dat_in = np.array(X)
    v_indice = np.arange(1, len(dat_in) + 1)
    ind_dat_in_valid = np.where(np.isfinite(dat_in))[0]
    # count_dat_in_valid = len(ind_dat_in_valid)

    n = len(dat_in[ind_dat_in_valid])
    m_dat_in = np.full((n, n), np.nan).astype(float)

    # Compute the matrix of slopes
    for i in range(n):
        index_i = ind_dat_in_valid[i]
        for j in range(n):
            index_j = ind_dat_in_valid[j]
            m_dat_in[i, j] = (dat_in[index_i] - dat_in[ind_dat_in_valid[index_j]]) / \
                             (v_indice[index_i] - v_indice[ind_dat_in_valid[index_j]])

    # Handle invalid data by setting diagonals and upper triangle to hbad
    for i in range(n):
        m_dat_in[i, i:n] = np.nan

    r_m_dat_in = m_dat_in[1:, :-1]

    # Identify valid data
    valid_data = np.where(np.isfinite(r_m_dat_in))
    sign_m_dat_in = np.sign(r_m_dat_in)

    # Handle ties
    if correction_tie:
        test_tie = f_count_tie(dat_in[ind_dat_in_valid])
        if len(test_tie) > 1:
            corr_factor_tie = 0
            for i in range(len(test_tie[0])):
                corr_factor_tie += test_tie[1][i] * (test_tie[0][i] * (test_tie[0][i] - 1) * (2 * test_tie[0][i] + 5))
        else:
            corr_factor_tie = 0
    else:
        corr_factor_tie = 0

    # Sum of signs
    Sum_sign = np.nansum(sign_m_dat_in[valid_data])
    var_S = 1. / 18. * ((n * (n - 1) * (2 * n + 5)) - corr_factor_tie)

    # Compute Sen's slope
    sen_slope = np.nanmedian(r_m_dat_in[valid_data])

    # Compute intercept
    b_mat = dat_in[ind_dat_in_valid] - (sen_slope * v_indice[ind_dat_in_valid])
    b = np.nanmedian(b_mat)

    # Rate of Change (%)
    RC = (sen_slope / abs(b)) * 100

    # Z-values calculation
    # Z_tab = norm.ppf(1 - alpha / 2)

    if Sum_sign > 0:
        Z_calc = (Sum_sign - 1) / np.sqrt(var_S)
    elif Sum_sign < 0:
        Z_calc = (Sum_sign + 1) / np.sqrt(var_S)
    else:
        Z_calc = 0

    # C_alpha = Z_tab * np.sqrt(var_S)

    # M1 = (len(valid_data[0]) - C_alpha) / 2
    # M2 = (len(valid_data[0]) + C_alpha) / 2

    # sorted_slope = np.sort(r_m_dat_in[valid_data])

    # if M1 < 1 :
    #     Lower_limit = sorted_slope[0]
    # else :
    #     Lower_limit = sorted_slope[int(np.round(M1 -1))]

    # if M2 > (len(sorted_slope)-1) :
    #     Upper_limit = sorted_slope[-1]
    # else :
    #     Upper_limit = sorted_slope[int(np.round(M2))]

    prob = abs((norm.cdf(abs(Z_calc)) - 1) * 2)

    if prob < alpha:
        res_kendall_test = 1
    else:
        res_kendall_test = 0

    # return sen_slope, prob, res_kendall_test, b, RC

    output['Is_the_change_of_annual_values_monotonic'] = res_kendall_test == 1
    output['Sen_slope'] = sen_slope
    output['Sen_intercept'] = b
    output['Rate_of_change_of_annual_values_in_percentage_per_time'] = RC

    return output


def f_count_tie(Z):
    """
    Counts the number of ties in the input series.

    Parameters:
        Z (array-like): Input vector

    Returns:
        numpy.ndarray: If tie values exist, returns a matrix with types of ties (double, triple, etc.) and their occurrences.
                       Else, returns 0.
    """

    X = np.array(Z)
    N = len(X)
    tie_count = np.zeros(N + 1, dtype=int)

    # Loop through each element in the array
    for i in range(N):
        # Find indices where elements are equal to X[i]
        i_tie = np.where(X == X[i])[0]
        count_a = len(i_tie)

        # If there are ties (more than one occurrence)
        if count_a > 1:
            if np.isfinite(X[i]):
                tie_count[count_a] += 1

            # Mark the counted ties with -9999
            X[i_tie] = np.nan

    # If there are any ties, prepare the result array
    if np.sum(tie_count) != 0:
        result_tie = np.transpose([
            np.where(tie_count >= 1)[0],  # Types of ties
            tie_count[tie_count >= 1]  # Corresponding occurrences
        ])
    else:
        result_tie = np.array([0])

    return result_tie


def f_seasonal_kendall_v2_0(X, period, alpha=0.05, correction_tie=0):
    """
    Performs the non-parametrical Seasonal Kendall test for trends.
    Computes the Sen slope for the whole time series and the monthly slope.
    Tests the homogeneity of the seasonal trends.
    Provides the significance level of the seasonal and global trends.

    Parameters:
        X (array): Input time series.
        period (int): The periodicity of the time series.
        alpha (float, optional): Significance level to be used. Default is 0.05.
        hbad (float, optional): Bad value. Default is -9999.
        correction_tie (int, optional): Correct the calculation for ties in the data. Default is 0.

    Returns:
        dict: Results containing the seasonal Kendall test information.
    """

    # Initialize the output structure
    # v_bad = -9999.
    # out_seasonal_kendall = {
    #     'sen_slope': np.nan,
    #     'prob_sen_slope': np.nan,
    #     'prob_chi_homog': np.nan,
    #     'prob_chi_trend': np.nan,
    #     'conf_limit': [np.nan, np.nan],
    #     'monthly_sen_slope': np.full(period, np.nan),
    #     'prob_monthly': np.full(period, np.nan),
    #     'intercept': np.nan,
    #     'RC': np.nan
    # }

    out_seasonal_kendall = {
        'Sen_slope': np.nan,
        'prob_Sen_slope': np.nan,
        'Sen_intercept': np.nan,

        'Sen_slope_confidence_limit': [np.nan, np.nan],
        'Is_Sen_slope_significant': np.nan,

        'Sen_monthly_slopes': np.full(period, np.nan),
        'prob_Sen_monthly_slopes': np.full(period, np.nan),
        'prob_test_Homogeneity_of_the_seasonal_signal_across_years': np.nan,
        'prob_test_Trend_in_the_seasonal_signal_across_years': np.nan,
        'Homogeneity_of_the_seasonal_signal_across_years': True,
        'Significant_Trend_in_the_seasonal_signal_across_years': False,

        'Rate_of_change_in_percentage_per_time': np.nan
    }

    # Prepare data
    X = np.array(X)
    ind_dat_in_valid = np.where(np.isfinite(X))[0]
    n_val = len(ind_dat_in_valid)
    dat_in = np.copy(X)

    # Handle missing values
    ind_bad_val_in = np.where(np.isnan(X))[0]
    if len(ind_bad_val_in) > 0:
        dat_in[ind_bad_val_in] = np.nan

    # Reshape data
    n_lign = int(np.ceil(len(dat_in) / period))
    if n_lign * period != len(dat_in):
        complete_y = np.full((n_lign * period) - len(dat_in), np.nan)
        m_dat_in = np.reshape(np.concatenate((dat_in, complete_y)), (period, n_lign), order='F')
    else:
        m_dat_in = np.reshape(dat_in, (period, n_lign), order='F')

    # Initialize arrays
    S_var = np.full(period, np.nan)
    n_month = np.full(period, np.nan)
    Sign_Diff = np.full(period, np.nan)
    sen_slope_m = np.full((period, int(np.sum(np.arange(n_lign)))), np.nan)

    # Calculate Sen's slope for each month
    for i_month in range(period):

        i_valid = np.where(np.isfinite(m_dat_in[i_month, :]))[0]
        count_valid = len(i_valid)

        if count_valid > 0:

            ind = 0
            sub_diff_k = np.full(int(np.sum(np.arange(count_valid))), np.nan)
            sub_slope_s = np.full(int(np.sum(np.arange(count_valid))), np.nan)

            for i in range(count_valid - 1):
                for j in range(i + 1, count_valid):
                    sub_diff_k[ind] = (m_dat_in[i_month, i_valid[j]] - m_dat_in[i_month, i_valid[i]]) / \
                                      abs(m_dat_in[i_month, i_valid[j]] - m_dat_in[i_month, i_valid[i]])
                    sub_slope_s[ind] = (m_dat_in[i_month, i_valid[j]] - m_dat_in[i_month, i_valid[i]]) / (j - i)
                    if sub_slope_s[ind] == 0:
                        sub_diff_k[ind] = 0  # Correct for differences of 0
                    ind += 1

            sen_slope_m[i_month, :len(sub_slope_s)] = sub_slope_s

            # Correction for ties in the data
            if correction_tie == 1:
                test_tie = np.unique(m_dat_in[i_month, i_valid], return_counts=True)
                if len(test_tie[0]) > 1:
                    corr_factor_tie = 0
                    for i in range(len(test_tie[0])):
                        corr_factor_tie += test_tie[1][i] * (
                                    test_tie[0][i] * (test_tie[0][i] - 1) * (2 * test_tie[0][i] + 5))
                    S_var[i_month] = (count_valid * (count_valid - 1) * (2 * count_valid + 5) - corr_factor_tie) / 18
                else:
                    S_var[i_month] = count_valid * (count_valid - 1) * (2 * count_valid + 5) / 18
            else:
                S_var[i_month] = count_valid * (count_valid - 1) * (2 * count_valid + 5) / 18

            Sign_Diff[i_month] = np.nansum(sub_diff_k)
            n_month[i_month] = count_valid

    # Test trend homogeneity over seasons
    Z_calc_i_month = np.full(period, np.nan)

    for i in range(period):
        if np.isfinite(Sign_Diff[i]):
            if Sign_Diff[i] > 0:
                Z_calc_i_month[i] = Sign_Diff[i] / np.sqrt(S_var[i])
            elif Sign_Diff[i] < 0:
                Z_calc_i_month[i] = Sign_Diff[i] / np.sqrt(S_var[i])
            else:
                Z_calc_i_month[i] = 0

    i_m_valid = np.where(np.isfinite(Z_calc_i_month))[0]
    count_n_valid_months = len(i_m_valid)

    Chi_total = np.nansum(Z_calc_i_month ** 2)
    Chi_trend = count_n_valid_months * (np.nanmean(Z_calc_i_month)) ** 2
    Chi_homo = Chi_total - Chi_trend
    Chi_theo_homog = chi2.ppf(0.95, count_n_valid_months - 1)

    out_seasonal_kendall['prob_test_Homogeneity_of_the_seasonal_signal_across_years'] = 1 - chi2.cdf(Chi_homo,
                                                                                                     count_n_valid_months - 1)

    if abs(Chi_homo) < Chi_theo_homog:
        # homog = 1
        out_seasonal_kendall['Homogeneity_of_the_seasonal_signal_across_years'] = False

        if abs(Chi_trend) < chi2.ppf(0.95, 1):
            # overall_chi_test = 1
            out_seasonal_kendall['Significant_Trend_in_the_seasonal_signal_across_years'] = True

        out_seasonal_kendall['prob_test_Trend_in_the_seasonal_signal_across_years'] = 1 - chi2.cdf(Chi_trend, 1)

    # Monthly trend analyses
    out_seasonal_kendall['Sen_monthly_slopes'] = np.full(period, np.nan)

    for i in range(period):
        i_val_i_month = np.where(np.isfinite(sen_slope_m[i, :]))[0]
        i_count_valid_slope = len(i_val_i_month)

        if i_count_valid_slope > 0:
            i_sel_slope = sen_slope_m[i, i_val_i_month]
            if i_count_valid_slope / 2 - round(i_count_valid_slope / 2) == 0:
                out_seasonal_kendall['Sen_monthly_slopes'][i] = np.nanmedian(i_sel_slope)
            else:
                out_seasonal_kendall['Sen_monthly_slopes'][i] = np.nanmedian(i_sel_slope)

    Z_calc_monthly = np.full(period, np.nan)
    out_seasonal_kendall['prob_Sen_monthly_slopes'] = np.full(period, np.nan)

    for i in range(period):
        if np.isfinite(Sign_Diff[i]):
            if Sign_Diff[i] > 0:
                Z_calc_monthly[i] = (Sign_Diff[i] - 1) / np.sqrt(S_var[i])
            elif Sign_Diff[i] < 0:
                Z_calc_monthly[i] = (Sign_Diff[i] + 1) / np.sqrt(S_var[i])
            else:
                Z_calc_monthly[i] = 0

    out_seasonal_kendall['prob_Sen_monthly_slopes'][i_m_valid] = np.abs(
        (norm.cdf(np.abs(Z_calc_monthly[i_m_valid])) - 1) * 2)

    # Sen's slope
    valid_sen_slope = np.where(np.isfinite(sen_slope_m))[0]
    count_valid_slope = len(valid_sen_slope)
    selec_slope = sen_slope_m.flatten()[valid_sen_slope]

    if count_valid_slope / 2 - round(count_valid_slope / 2) == 0:
        sen_slope = np.nanmedian(selec_slope)
    else:
        sen_slope = np.nanmedian(selec_slope)

    out_seasonal_kendall['Sen_slope'] = sen_slope / period

    # Intercept
    b_mat = dat_in[ind_dat_in_valid] - (sen_slope / period * ind_dat_in_valid)
    if n_val / 2 - round(n_val / 2) == 0:
        out_seasonal_kendall['Sen_intercept'] = np.nanmedian(b_mat)
    else:
        out_seasonal_kendall['Sen_intercept'] = np.nanmedian(b_mat)

    # Rate of Change (%)
    out_seasonal_kendall['Rate_of_change_in_percentage_per_time'] = (sen_slope / np.abs(
        out_seasonal_kendall['Sen_intercept'])) * 100

    # Significance of overall trend
    total_sign_diff = np.sum(Sign_Diff[i_m_valid])
    total_s_var = np.sum(S_var[i_m_valid])

    if total_sign_diff > 0:
        Z_calc = (total_sign_diff - 1) / np.sqrt(total_s_var)
    elif total_sign_diff < 0:
        Z_calc = (total_sign_diff + 1) / np.sqrt(total_s_var)
    else:
        Z_calc = 0

    Z_tab = norm.ppf(1 - alpha / 2)
    C_alpha = Z_tab * np.sqrt(total_s_var)

    M1 = (count_valid_slope - C_alpha) / 2
    M2 = (count_valid_slope + C_alpha) / 2

    sorted_slope = np.sort(selec_slope)

    Lower_limit = sorted_slope[int(np.round(M1 - 1))]
    Upper_limit = sorted_slope[int(np.round(M2))]

    out_seasonal_kendall['Sen_slope_confidence_limit'] = [Lower_limit, Upper_limit]
    out_seasonal_kendall['prob_Sen_slope'] = np.abs((norm.cdf(np.abs(Z_calc)) - 1) * 2)

    # Kendall test result
    if np.abs(Z_calc) > np.abs(Z_tab):
        season_kendall_test = 1
    else:
        season_kendall_test = 0

    out_seasonal_kendall['Is_Sen_slope_significant'] = season_kendall_test == 1

    return out_seasonal_kendall


def F_test_pente(x, y, slope):
    """
    Function to test if the slope is significantly different from zero.

    Parameters:
    x : array-like
        Input array for the x-values.
    y : array-like
        Input array for the y-values.
    slope : float
        The slope to test.

    Returns:
    tuple
        decision : int
            1 if slope is significantly different from zero, otherwise 0.
        t_a : float
            The t-value calculated for the slope.
    """

    N = len(x)  # Number of elements in x

    # Calculate the variance of the slope
    s_2_b = ((np.std(y) / np.std(x)) ** 2 - slope ** 2) / (N - 2)

    if s_2_b >= 0:
        # Calculate the t-value
        t_a = slope / np.sqrt(s_2_b)

        # Threshold value from t-distribution for alpha = 0.05 (two-tailed test)
        thresh = t.ppf(1 - 0.025, N - 2)

        if abs(t_a) >= thresh:
            # Slope is significantly different from 0
            decision = 1
        else:
            # Slope is not significantly different from 0
            decision = 0
    else:
        # If s_2_b is negative, set default values
        decision = 0
        t_a = 0

    return decision, t_a


def linfit(x, y):
    """
    Perform a linear fit.

    Parameters:
    x : array-like
        Independent variable.
    y : array-like
        Dependent variable.

    Returns:
    tuple
        Slope and intercept of the linear fit.
    """

    slope, intercept, _, _, _ = stats.linregress(x, y)
    return intercept, slope


# =============================================================================
#### Main functions
# =============================================================================

def temporal_decomp_V2_7_x11(values, dates, time_frequency,
                             filter_outlier=False, overall_cutoff=50.0, 
                             out_limit=3.0, perc_month_limit=50.0, 
                             var_stationary=False, lin_interpol=False, 
                             cutoff_fill=30.0, season_test=False):
    
    """
    ; PURPOSE:
    ; Decompose a time series var into 3 components: T: TREND, S: SEASON, Y: Irregular
    ;
    ; CATEGORY:
    ;
    ; CALLING SEQUENCE:
    ; temporal_decomp_V2_6_x11, Var, month, hbad
    ;
    ; INPUTS:
    ; Var:  vector representing the ANNUAL TIME SERIES to be decomposed (12 months cycle including bad values)
    ; month: vector of month
    ; hbad: bad value
    ; data: ouput structure
    
    ; KEYWORD PARAMETERS:
    ;
    ;filter_outlier : if set filter of the outlier assuming a normal distribution (mean +/- (out_lim*std))
    ;out_limit : limit considered for outlier detection ...default: 3
    ;overall_cutoff: cutoff value for the test on the initial % of valid data in the TS. If %age of valid data < overall cutoff -> flagged. Default 50%
    ;perc_month_limit : cutoff value for selecting the valid months for the "shortened years". in % default 50%
    ;X11_pezzulli: if set use the X11 method defined in Pezzulli et al., J of Climate, 2005
    ;X11_user : if set use the X11 method defined in the X11 user guide manual
    ;var_stationary : if stat compute the contribution of the X11 components to the variance of the stationary part of the initial TS
    ;lin_interpol : if set to 1 fill the gaps in the TS using linear interpolation instead of the evf method (if set to 0)
    ;cutoff_fill : cutoff value for the maximal %age of missing data acceptable for performing the gap filling procedure. Default 30% of missing data max
    ;season_test : if set then the test on the presence of seasonality in the data is performed
    ;
    
    ;========================================================================================
    ;OUTPUT STRUCTURE DEFINITION
    
    ;PtStr = {	var_S:fbad , $ ; Percentage of variance of X due to the Seasonal component
    ;            var_T:fbad , $; Percentage of variance of X due to the Trend component
    ;            var_Y:fbad ,  $; Percentage of variance of X due to the Irregular component
    ;			covar: fbad, $; Percentage of variance of X related to the covariance terms
    ;
    ;			var_X: fbad, $;  variance of the original time series
    ;			var_coeff: fbad, $;  variation coefficient of the original time series
    ;
    ;			code_pres_month: 0, $ ; Integer coding the m_test_year vector of monthly data validity
    ;
    ;			S : make_array(N_months, value =fbad), $ ; vector of N months containing the Seasonal component
    ;			T : make_array(N_months, value =fbad), $ ; vector of N months containing the Trend component
    ;			Y : make_array(N_months, value =fbad), $  ; vector of N months containing the Irregular component
    ;
    ;			year_n_months :fbad ,$ ; Number of months in 1 synthetic year
    ;			N_outlier: 0. ,$ ; Number of outlier
    ;			perc_valid_months: fbad,  $ ; Percentage of valid initial data (Number of valid data in the initial time series/ 87) *100.
    ;			nb_valid_month:  fbad,  $ ;Total number of valid data considered in the analysis (after short_year and gap_filling)
    ;			N_missing_data: 0. ,$  ; Number of data which have been filled (including missing values, outlier...)
    ;			test_season: make_array(2, value = fbad) ,$ ; test the homogeneity of the season
    ;
    ;			Kendall_sen: make_array(5, value = fbad), $
    ;			Seasonal_Kendall_sen: make_array(5, value = fbad), $
    ;
    ;			slope_trend: fbad, $
    ;           intercept_trend : fbad, $               .
    ;           prob_Test_trend:fbad , $
    ;
    ;			flag: fbad  $; TEST the efficiency of the gap filling method
    ;		}
    
    ; SIDE EFFECTS:
    ; none
    ;
    ; COMMON BLOCKS:
    ; none
    ;
    ; MODIFICATION HISTORY:
    ;
    ; Vantrepotte V.
    ;Written  20.10.07
    ;Last Update: 26-09-2008: RC rate of Change + comment
    
    """
    
    # missing_values_are = np.nan # -99999
    # fbad = np.nan
    
    if isinstance(dates, str) : 
        dates = pd.to_datetime(dates)
    
    var = pd.Series(values, index = dates)

    if time_frequency == "ANNUAL" : 
        full_date_range = [pd.Timestamp(f'{date.year}-07-01') for date in 
                           pd.date_range(start=dates.min(), end=dates.max() + pd.offsets.MonthEnd(0), freq = '1Y')]
    
    elif time_frequency == "MONTHLY" : 
        full_date_range = [pd.Timestamp(f'{date.year}-{date.month:02d}-15') for date in 
                           pd.date_range(start=dates.min(), end=dates.max() + pd.offsets.MonthEnd(0), freq = '1M')]
        
    elif time_frequency == "WEEKLY" : 
        full_date_range = [pd.Timestamp(f'{date.year}-{date.month:02d}-{the_day:02d}')
                             for date in pd.to_datetime( np.unique( pd.date_range(start=dates.min(), end=dates.max(), freq = '1D').strftime('%Y-%m') ) )
                             for the_day in [ 4, 12, 20, 28]]
    
    # TODO: The logic gate was left open, nor was there a value for 'DAILY' so I changed that
    elif time_frequency == "DAILY" : 
        full_date_range = [pd.Timestamp(dates) for dates in 
                           pd.date_range(start=dates.min(), end=dates.max(), freq = '1D')]

    # TODO: There does not appear to be any use for 'INITIAL' in the existing code
    # I leave this here for now in case it comes up again later
    # elif time_frequency == "INITIAL" : 
        # full_date_range = dates

    else : 
        raise ValueError("time_frequency must be one of 'MONTHLY', 'ANNUAL', 'WEEKLY', or 'DAILY'")

    var_reindexed = var.reindex(full_date_range, fill_value=np.nan)
    # var_reindexed[np.isnan(var_reindexed)] = missing_values_are
    
    N_values = len(var_reindexed)
    data = {
            "0_variance_due_to_Seasonal_signal": np.nan,  # Percentage of variance of X due to the Seasonal component
            "1_variance_due_to_Interannual_signal": np.nan,  # Percentage of variance of X due to the Trend component
            "2_variance_due_to_Residual_signal": np.nan,  # Percentage of variance of X due to the Irregular component
            "3_covar": np.nan,  # Percentage of variance of X related to the covariance terms
    
            "4_var_X": np.nan,  # Variance of the original time series
            "5_var_coeff": np.nan,  # Variation coefficient of the original time series
    
            "6_code_pres_month": 0,  # Integer coding the m_test_year vector of monthly data validity
    
            "7_dates": np.full(N_values, np.nan).astype(float),
            "8_values_ini": np.full(N_values, np.nan).astype(float),
            "8_values_ini_interpolated": np.full(N_values, np.nan).astype(float),
            
            "9_Seasonal_signal": np.full(N_values, np.nan).astype(float),  # Vector of N months containing the Seasonal component
            "10_Interannual_signal": np.full(N_values, np.nan).astype(float),  # Vector of N months containing the Trend component
            "11_Residual_signal": np.full(N_values, np.nan).astype(float),  # Vector of N months containing the Irregular component
    
            "12_year_n_months": np.nan,  # Number of months in 1 synthetic year
            "13_N_outlier": 0.0,  # Number of outliers
            "14_perc_valid_data": np.nan,  # Percentage of valid initial data (Number of valid data in the initial time series / 87) * 100
            "15_nb_valid_month": np.nan,  # Total number of valid data considered in the analysis (after short_year and gap_filling)
            "16_N_missing_data": 0.0,  # Number of data which have been filled (including missing values, outliers, etc.)
            "17_Homogeneity_of_the_Seasonal_signal": np.full(2, np.nan),  # Test the homogeneity of the season
    
            # "Kendall_sen": np.full(5, np.nan),  # Kendall's tau statistics
            "18_Kendall_Sen_analyses_on_Interannual_signal": np.nan,  # Kendall's tau statistics
            "19_Kendall_Sen_analyses_on_Seasonal_signal": np.nan,  # Seasonal Kendall's tau statistics
    
            "20_slope_trend_per_time_step": np.nan,  # Slope of the linear trend derived from Tt
            
            "22_intercept_trend": np.nan,  # Intercept of the linear trend derived from Tt
            "23_prob_Test_trend": np.nan,  # Probability associated with the linear trend derived from Tt
    
            "24_flag": np.nan,  # Test the efficiency of the gap filling method
            
            "25_Is_the_slope_trend_per_time_step_significant": np.nan,
        }
    
    data['7_dates'] = var_reindexed.index
    data['8_values_ini'] = var_reindexed.values
    
    # Calculate the percentage of valid months
    perc_valid_data = 100 - np.sum(np.isnan(var_reindexed)) / len(var_reindexed) * 100.0
    data['14_perc_valid_data'] = perc_valid_data

    if perc_valid_data < overall_cutoff :        
        data['24_flag'] = -100
        return data

    if time_frequency in ["MONTHLY", "ANNUAL"] : 
        months = np.array([x.month for x in var_reindexed.index])
    elif time_frequency == "WEEKLY" :
        months = np.array([x.strftime('%m-%d') for x in var_reindexed.index])
    elif time_frequency == "DAILY" :
        months = np.array([x.strftime('%m-%d') for x in var_reindexed.index])
    else : 
        raise ValueError("time_frequency must be one of 'MONTHLY', 'ANNUAL', 'WEEKLY', or 'DAILY'")
        
    values_to_use = np.zeros(len(var_reindexed)).astype(bool)
    
    n_time_step_per_year = 12
    if time_frequency == "WEEKLY" : 
        n_time_step_per_year = int(n_time_step_per_year * 4)
    elif time_frequency == "ANNUAL" : 
        n_time_step_per_year = int(n_time_step_per_year / 12)
        # TODO: I added this logic gate, not certain it is correct
    elif time_frequency == "DAILY" :
        n_time_step_per_year = int(len(np.unique(months)))
    else : 
        raise ValueError("time_frequency must be one of 'MONTHLY', 'ANNUAL', 'WEEKLY', or 'DAILY'")
        
    index_month_to_use = np.zeros(n_time_step_per_year).astype(bool)
    
    for i_m in np.unique(months) :
        
        # i_m=np.unique(months)[0]

        ind_month_i = np.where(months == i_m)[0]
        ind_bad_val_month_i = np.where( np.isnan(var_reindexed[ind_month_i]) )[0]
        n_bad_val_month_i = len(ind_bad_val_month_i)
        count_n_month_i = len(ind_month_i)

        if (100 * float(n_bad_val_month_i) / float(count_n_month_i)) < perc_month_limit :
            index_month_to_use[ind_month_i[0]] = True
            values_to_use[ind_month_i] = True
    
    if all(index_month_to_use == False) : 
        data['24_flag'] = -101
        return data
    
    data['12_year_n_months'] = np.sum(index_month_to_use)
    data['15_nb_valid_month'] = np.sum(values_to_use)
    if time_frequency == "MONTHLY" : 
        data['14_perc_valid_month'] = 100 * data['15_nb_valid_month'] / len(np.unique([x.strftime('%Y-%m') for x in var_reindexed.index]))
    if time_frequency == "WEEKLY" :     
        data['14_perc_valid_month'] = 100 * data['15_nb_valid_month'] / len(np.unique([x.strftime('%Y-%m-%d') for x in var_reindexed.index]))

    index_var_to_use = np.where(values_to_use)[0]
    var_to_use = var_reindexed.iloc[ index_var_to_use ]
    
    if time_frequency in ["MONTHLY", "ANNUAL"] : 
        months = np.array([x.month for x in var_reindexed.index])
    if time_frequency == "WEEKLY" :     
        months = np.array([x.strftime('%m-%d') for x in var_reindexed.index])
    months = np.unique(months)

    n_outlier = 0
    
    if filter_outlier : 
    
        for month in months :
            
            var_month = var_to_use.iloc[ np.where([x == month for x in months])[0] ]
            
            # var_month_to_use_for_mean_std = var_month[np.isfinite(var_month)]
            mean_month = np.nanmean(var_month)
            std_month = np.nanstd(var_month)
            
            index_outlier = np.where((var_month > mean_month + out_limit * std_month) | 
                                     (var_month < mean_month - out_limit * std_month))[0]
    
            n_outlier += len(index_outlier)
    
            if len(index_outlier) != 0:                
                dates_outliers = var_month.iloc[index_outlier].index
                var_to_use.loc[dates_outliers] = np.nan
    
    # GAPS FILLING using EVF METHOD
    test_missing = np.where( np.isnan(var_to_use) )[0]
    count_missing = len(test_missing)

    var_interpolated = var_to_use
    if count_missing != 0:
        
        if (100 * float(count_missing) / len(var_to_use)) < cutoff_fill :
            
            if lin_interpol == False :
                
                # Use EVF method
                max_iter = 20
                level_SD_perc = 10
                ACP_cutoff = 80

                var_interpolated = F_EVF_V1_2(X = var_to_use, 
                                              level_SD_perc = level_SD_perc, ACP_cutoff = ACP_cutoff, max_iter = max_iter)
                data['24_flag'] = 1

                if (count_missing > 1) and (var_interpolated[test_missing[0]] == var_interpolated[test_missing[1]]) :
                                                
                    data['24_flag'] = 2
                    level_SD_perc = 10
                    ACP_cutoff = 50
                    var_interpolated[test_missing] = np.nan
                    var_interpolated = F_EVF_V1_2(X = var_interpolated,
                                                  level_SD_perc = level_SD_perc, ACP_cutoff = ACP_cutoff, max_iter = max_iter)

                    if var_interpolated[test_missing[0]] - var_interpolated[test_missing[1]] == 0:
                        data['24_flag'] = 3

            else:
                # Use linear interpolation
                var_interpolated = np.interp(np.arange(len(var_to_use)), 
                                    np.arange(len(var_to_use))[ np.isfinite(var_to_use) ], 
                                    var_to_use[ np.isfinite(var_to_use) ])
                data['24_flag'] = 4
        else:
            data['24_flag'] = -10
    else:
        data['24_flag'] = 0

    data['8_values_ini_interpolated'] = var_interpolated

    test_NAN = np.isnan(var_interpolated).sum()
    if test_NAN != 0:
        data['24_flag'] = -9

    # TIME SERIES analyses
    if test_NAN != 0 or data['24_flag'] == -10 or (time_frequency != "ANNUAL" and data['12_year_n_months'] <= 3) :
        
        return data
                
    Census_X11 = F_census_X_11_pezzulli_V1_2(Xt = var_interpolated, period = data['12_year_n_months'])

    data['9_Seasonal_signal'][index_var_to_use] = Census_X11['St']
    data['11_Residual_signal'][index_var_to_use] = Census_X11['Yt']
    data['10_Interannual_signal'][index_var_to_use] = Census_X11['Tt']

    if season_test:
        data['17_Homogeneity_of_the_Seasonal_signal'] = f_test_season_v1_1(
            X = var_interpolated[:data['12_year_n_months'] * (len(var_interpolated) // data['12_year_n_months'])],
            T = Census_X11['Tt'][:data['12_year_n_months'] * (len(var_interpolated) // data['12_year_n_months'])],
            S = Census_X11['St'][:data['12_year_n_months'] * (len(var_interpolated) // data['12_year_n_months'])],
            period = data['12_year_n_months'])

    data['4_var_X'] = np.var(var_interpolated)
    data['5_var_coeff'] = np.std(var_interpolated) / np.mean(var_interpolated) * 100.0

    if var_stationary:
        
        fit_trend_LINEAR = np.polyfit(np.arange(len(Census_X11['Tt'])), Census_X11['Tt'], 1)
        yfit = np.polyval(fit_trend_LINEAR, np.arange(len(Census_X11['Tt'])))
        data['0_variance_due_to_Seasonal_signal'] = np.var(Census_X11['St']) / np.var(var_interpolated - yfit) * 100.0
        data['1_variance_due_to_Interannual_signal'] = np.var(Census_X11['Tt'] - yfit) / np.var(var_interpolated - yfit) * 100.0
        data['2_variance_due_to_Residual_signal'] = np.var(Census_X11['Yt']) / np.var(var_interpolated - yfit) * 100.0
        data['3_covar'] = 100.0 * (1.0 - np.sum([np.var(Census_X11['St']) / np.var(var_interpolated - yfit), 
                                                np.var(Census_X11['Tt'] - yfit) / np.var(var_interpolated - yfit), 
                                                np.var(Census_X11['Yt']) / np.var(var_interpolated - yfit)]))
        
    else:
        
        var_St = np.var(Census_X11['St'])
        var_Tt = np.var(Census_X11['Tt'])
        var_Yt =np.var(Census_X11['Yt'])
        sum_vars = var_St + var_Tt + var_Yt
        
        data['0_variance_due_to_Seasonal_signal'] = var_St / sum_vars * 100.0
        data['1_variance_due_to_Interannual_signal'] = var_Tt / sum_vars * 100.0
        data['2_variance_due_to_Residual_signal'] = var_Yt / sum_vars * 100.0
        data['3_covar'] = 100.0 * (1.0 - np.sum([np.var(Census_X11['St']) / np.var(var_interpolated), 
                                                np.var(Census_X11['Tt']) / np.var(var_interpolated), 
                                                np.var(Census_X11['Yt']) / np.var(var_interpolated)]))
       
    if time_frequency != 'ANNUAL' :
        v_mean_year = [np.mean(var_interpolated[ int(i * data['12_year_n_months']) : int( (i + 1) * data['12_year_n_months'] -1) ]) 
                       for i in np.arange(len(var_interpolated) / data['12_year_n_months'])
                       if len(var_interpolated[ int(i * data['12_year_n_months']) : int( (i + 1) * data['12_year_n_months'] ) ]) == data['12_year_n_months']]
    else : 
        v_mean_year = var_interpolated
    
    if len( np.where(np.isfinite(v_mean_year))[0] ) > 2 :
        data['18_Kendall_Sen_analyses_on_Interannual_signal'] = f_kendall_sen_v1_0(X = v_mean_year, alpha=0.05, correction_tie = True)
        
    data['19_Kendall_Sen_analyses_on_Seasonal_signal'] = f_seasonal_kendall_v2_0(var_interpolated, data['12_year_n_months'], alpha=0.05)
    
    # Generate x-values corresponding to the length of census_X11.Tt
    x_values = np.arange(len(Census_X11['Tt']))
    
    # Perform a linear fit
    fit_trend_test = linfit(x_values, Census_X11['Tt'])
    
    # Test the significance of the slope
    test_slope = F_test_pente(x_values, Census_X11['Tt'], fit_trend_test[1])
    
    # Store the slope and intercept in a data object (assuming 'data' is a pre-defined object with the necessary attributes)
    data['20_slope_trend_per_time_step'] = fit_trend_test[1]
    data['22_intercept_trend'] = fit_trend_test[0]
    
    # Calculate the significance level
    data['15_nb_valid_month'] = len(Census_X11['Tt'])  # Assuming this represents the number of valid months
    data['21_Significance_level_of_the_slope'] = 1.0 - abs(1.0 - 2.0 * (1.0 - t.cdf(abs(test_slope[1]), data['15_nb_valid_month'] - 2)))
        
    data['25_Is_the_slope_trend_per_time_step_significant'] = data['21_Significance_level_of_the_slope'] < 0.05
    
    return data


def Apply_X11_method_on_time_series(core_arguments, Zones, nb_cores,
                                    plume_time_step,
                                    plume_dir_in, X11_dir_out,
                                    include_river_flow=False):
    
    core_arguments.update({'Zones': Zones,
                           'Years': "*",
                           'Satellite_variables': ['SPM'],
                           'Temporal_resolution': ([plume_time_step]
                                                   if isinstance(plume_time_step, str)
                                                   else plume_time_step)})

    cases_to_process = get_all_cases_to_process_for_regional_maps_or_plumes_or_X11(core_arguments)

    var_to_use = 'area_of_the_plume_mask_in_km2'

    for i, info in cases_to_process.iterrows():

        # info = cases_to_process.iloc[i]

        file_names_pattern = (fill_the_sat_paths(info,
                                                 path_to_fill_to_where_to_save_satellite_files(
                                                     plume_dir_in + "/" + info.Zone),
                                                 local_path=True)
                              .replace(info.atmospheric_correction, f'{info.atmospheric_correction}/PLUME_DETECTION/')
                              .replace('/*/*/*', '/Time_series_of_plume_area_and_SPM_threshold.csv'))

        if os.path.exists(file_names_pattern) == False:
            print(f'File does not exists : {file_names_pattern}')
            continue

        ts_data = pd.read_csv(file_names_pattern)

        apply_X11_method_and_save_results(values=ts_data[var_to_use].tolist(), variable_name=var_to_use,
                                          dates=ts_data.date, info=info,
                                          X11_dir_out=X11_dir_out)

        if include_river_flow:
            river_flow_data = load_csv_files_in_the_package_folder(RIVER_FLOW=True, Zone_of_river_flow=info.Zone,
                                                                   RIVER_FLOW_time_resolution=info.Temporal_resolution)

            apply_X11_method_and_save_results(values=river_flow_data.Values.tolist(),
                                              variable_name='river_flow',
                                              dates=river_flow_data.Date,
                                              info=info.replace({info.Satellite_variable: "River_flow",
                                                                 info.atmospheric_correction: "",
                                                                 info.sensor_name: "",
                                                                 info.Data_source: "River_flow"}),
                                              X11_dir_out=X11_dir_out)

            # Source the R script
            X11_R_path = os.path.join(func_dir, 'X11.R')
            robjects.r['source'](X11_R_path)
            # robjects.r['source']("myRIOMAR_dev/_4_X11_analysis/utils.R")

            r_function = robjects.r['plot_time_series_of_plume_area_and_river_flow']

            # Call the R function
            r_function(
                where_are_saved_X11_results=robjects.StrVector([X11_dir_out]),
                Zone=robjects.StrVector([info.Zone]),
                Data_source=robjects.StrVector([info.Data_source]),
                sensor_name=robjects.StrVector([info.sensor_name]),
                atmospheric_correction=robjects.StrVector([info.atmospheric_correction]),
                Temporal_resolution=robjects.StrVector([info.Temporal_resolution])
            )

