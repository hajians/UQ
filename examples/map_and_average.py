#! /usr/bin/env python2.7

from examples.utils_map_and_average import get_pipe_true, plot_lambda, plot_pressure, get_data_frame, \
    get_pipe_from_file_df


def get_setting(expan_coef):
    return dict(
        c_sound=1.0,
        t_final=5.0,
        x_l=0.0,
        x_r=1.0,
        dx=0.005,
        boundary_eps=0.05,
        expan_coef=expan_coef
    )


def main(filepath, save_fig_lambda, save_fig_pressure):
    df_coef, df_proba = get_data_frame(filepath)

    expan_coef = df_coef.shape[1]
    setting = get_setting(expan_coef=expan_coef)

    pipe = get_pipe_from_file_df(setting=setting)

    pipe_map = get_pipe_from_file_df(setting=setting)

    pipe_true = get_pipe_true(setting=setting)

    plot_lambda(pipe=pipe, pipe_true=pipe_true, pipe_map=pipe_map, df_coef=df_coef, df_proba=df_proba,
                save_fig_path=save_fig_lambda)

    plot_pressure(pipe=pipe, pipe_true=pipe_true, pipe_map=pipe_map, df_coef=df_coef, df_proba=df_proba,
                  save_fig_path=save_fig_pressure)


if __name__ == '__main__':
    filepath = 'results/uq_pcn_n_4_sample_10000_burn_1000.dat'
    save_fig_lambda = 'results/uq_pcn_n_4_lambda.pdf'
    save_fig_pressure = 'results/uq_pcn_n_4_pressure.pdf'

    main(filepath=filepath, save_fig_lambda=save_fig_lambda, save_fig_pressure=save_fig_pressure)
