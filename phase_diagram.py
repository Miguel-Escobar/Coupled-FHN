from instantaneous_neurons import *
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import multiprocessing

def calculate_phase_diagram_point(params):
    i, j, a1, a2, eps, c = params
    initial_state = [-a1*3, 0.0, a2*1.0, 0.0]
    t_span = (0, 1000)
    t_values, x1_values, y1_values, x2_values, y2_values, norm_difference, peak_times, peak_values = system_observables(a1, a2, eps, c, initial_state, t_span, max_step=10)
    print(f"Hey! I'm done with point ({i}, {j})")
    return i, j, peak_values[-1]/np.max(peak_values)

def main():
    a1_range = np.linspace(-1.0, 1.0, 500)
    a2_range = np.linspace(-1.0, 1.0, 500)
    c = 1.0
    eps = 0.01

    # Create a list of parameter values for parallel processing
    params_list = [(i, j, a1, a2, eps, c) for i, a1 in enumerate(a1_range) for j, a2 in enumerate(a2_range)]

    # Use multiprocessing to calculate phase diagram
    with multiprocessing.Pool() as pool:
        results = pool.map(calculate_phase_diagram_point, params_list)

    # Convert the results into the phase_diagram array
    phase_diagram = np.zeros((len(a1_range), len(a2_range)))
    for i, j, value in results:
        phase_diagram[i, j] = value
    return a1_range, a2_range, phase_diagram

if __name__ == "__main__":
    # a1_range = np.linspace(-1.0, 1.0, 500)
    # a2_range = np.linspace(-1.0, 1.0, 500)
    # c = 1.0
    # eps = 0.01
    # phase_diagram = np.zeros((len(a1_range), len(a2_range)))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    # for i, a1 in tqdm(enumerate(a1_range)):
    #     for j, a2 in enumerate(a2_range):
    #         initial_state = [-a1*3, 0.0, a2*1.0, 0.0]
    #         t_span = (0, 1000)
    #         t_values, x1_values, y1_values, x2_values, y2_values, norm_difference, peak_times, peak_values = system_observables(a1, a2, eps, c, initial_state, t_span, max_step=10)
    #         phase_diagram[i, j] = peak_values[-1]/np.max(peak_values)



    # a1_range = np.linspace(-1.0, 1.0, 500)
    # a2_range = np.linspace(-1.0, 1.0, 500)
    # c = 1.0
    # eps = 0.01

    # # Create a list of parameter values for parallel processing
    # params_list = [(i, j, a1, a2, eps, c) for i, a1 in enumerate(a1_range) for j, a2 in enumerate(a2_range)]

    # # Use multiprocessing to calculate phase diagram
    # with multiprocessing.Pool() as pool:
    #     results = pool.map(calculate_phase_diagram_point, params_list)

    # # Convert the results into the phase_diagram array
    # phase_diagram = np.zeros((len(a1_range), len(a2_range)))
    # for i, j, value in results:
    #     phase_diagram[i, j] = value

    a1_range, a2_range, phase_diagram = main()
    
    fig, ax = plt.subplots()
    im = ax.pcolormesh(a1_range, a2_range, phase_diagram, shading='auto')
    ax.set_title('Phase diagram')
    ax.set_xlabel('a1')
    ax.set_ylabel('a2')

    fig.colorbar(im, ax=ax, label="Final amplitude of the difference")

    plt.show()
