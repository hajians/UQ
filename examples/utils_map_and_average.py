import pandas as pd
from matplotlib import pyplot as plt

from UQuant.SemilinearSystem import SemiLinSystem


def get_true_friction():
    return [0.185938, -0.0519335, 0., 0., -0.0696583, 0.0336323, 0., 0.,
            0.0348292, -0.0121076, 0., 0., -0.00773981, 0.00105987, 0., 0., 0.,
            -0.000641154, 0., 0., -0.00278633, 0.00250158, 0., 0., 0.00386991,
            -0.00179107, 0., 0., -0.0014216, 0.000230816, 0., 0., 0.,
            -0.000179701, 0., 0., -0.000859979, 0.000838478, 0., 0., 0.00139317]


def get_pipe_true(setting):
    pipe = SemiLinSystem(c_sound=setting.get('c_sound'),
                         t_final=setting.get('t_final'),
                         x_l=setting.get('x_l'),
                         x_r=setting.get('x_r'),
                         dx=setting.get('dx'),
                         lambda_len=len(get_true_friction()),
                         eps=setting.get('boundary_eps')
                         )

    return pipe


def plot_lambda(pipe, pipe_true, pipe_map, df_coef, df_proba, save_fig_path):
    pipe.get_lambda_average(df_coef.mean().tolist())

    pipe_true.get_lambda_average(get_true_friction())

    coef_max = df_coef[df_proba.max() == df_proba].values.tolist()[0]

    pipe_map.get_lambda_average(coef_max)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)

    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18, rotation=0)
    plt.xlabel("$x$", fontsize=24)

    ax.plot(pipe.mesh, pipe.lambda_avg,
            label="$E(\lambda)(x)$")

    ax.plot(pipe_true.mesh, pipe_true.lambda_avg, "--",
            label="$\lambda_{true}(x)$")

    ax.plot(pipe_map.mesh, pipe_map.lambda_avg,
            label="$\lambda_{MAP}(x)$")

    ax.legend(loc='upper left', fontsize=12)

    plt.tight_layout()
    plt.savefig(save_fig_path)
    plt.show()


def plot_pressure(pipe, pipe_true, pipe_map, df_coef, df_proba, save_fig_path, time_ins=20):
    pipe.run(df_coef.mean().tolist())
    pipe.get_presure_drop(time_instance=time_ins)

    pipe_true.run(get_true_friction())
    pipe_true.get_presure_drop(time_instance=time_ins)

    coef_max = df_coef[df_proba.max() == df_proba].values.tolist()[0]
    pipe_map.run(coef_max)
    pipe_map.get_presure_drop(time_instance=time_ins)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)

    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18, rotation=0)
    plt.xlabel("$x$", fontsize=24)

    ax.plot(pipe.timeslices, pipe.pressure_drop, 'x', label="Mean", markersize=10)

    ax.plot(pipe_true.timeslices, pipe_true.pressure_drop, "--", label="$\delta p_{true}^n$")

    ax.plot(pipe_map.timeslices, pipe_map.pressure_drop, '^', label="MAP", markersize=10)

    ax.legend(loc='lower right', fontsize=12)

    plt.ylabel("$\delta p_h^n$", fontsize=24)

    plt.tight_layout()
    plt.savefig(save_fig_path)
    plt.show()


def get_data_frame(filepath):
    raw_df = pd.read_csv(filepath, header=None)
    return raw_df.iloc[:, 0:-1], raw_df.iloc[:, -1]


def get_pipe_from_file_df(setting):
    pipe = SemiLinSystem(c_sound=setting.get('c_sound'),
                         t_final=setting.get('t_final'),
                         x_l=setting.get('x_l'),
                         x_r=setting.get('x_r'),
                         dx=setting.get('dx'),
                         lambda_len=setting.get('expan_coef'),
                         eps=setting.get('boundary_eps')
                         )

    return pipe
