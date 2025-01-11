import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

# Function to calculate RMSE
def calculate_rmse(benchmark, comparison):
    rmse_x = np.sqrt(mean_squared_error(benchmark['x'], comparison['x']))
    rmse_y = np.sqrt(mean_squared_error(benchmark['y'], comparison['y']))
    rmse_vx = np.sqrt(mean_squared_error(benchmark['vx'], comparison['vx']))
    rmse_vy = np.sqrt(mean_squared_error(benchmark['vy'], comparison['vy']))
    return rmse_x, rmse_y, rmse_vx, rmse_vy


benchmark_orbit = pd.read_csv('benchmark_orbit.csv')
barnes_hut_orbit = pd.read_csv('barnes_hut_orbit.csv')
brute_force_orbit = pd.read_csv('brute_force_orbit.csv')
fmm_orbit = pd.read_csv('fmm_orbit.csv')
parallel_barnes_hut_orbit = pd.read_csv('parallel_barnes_hut_orbit.csv')
parallel_brute_force_orbit = pd.read_csv('parallel_brute_force_orbit.csv')
parallel_fmm_orbit = pd.read_csv('parallel_fmm_orbit.csv')


data_files = {
    "Benchmark": benchmark_orbit,
    "Barnes Hut": barnes_hut_orbit,
    "Brute Force": brute_force_orbit,
    "FMM": fmm_orbit,
    "Parallel Barnes Hut": parallel_barnes_hut_orbit,
    "Parallel Brute Force": parallel_brute_force_orbit,
    "Parallel FMM": parallel_fmm_orbit
}

# Calculate RMSE for each method
results = {}
for key, df in data_files.items():
    if key != "Benchmark":
        rmse = calculate_rmse(benchmark_orbit, df)
        results[key] = {"RMSE": rmse}


rmse_data = {method: values['RMSE'] for method, values in results.items()}
labels = ['x', 'y', 'vx', 'vy']
variable_titles = ['x Position', 'y Position', 'vx Velocity', 'vy Velocity']

# Creating individual plots for each variable
for i, var_title in enumerate(variable_titles):
    fig, ax = plt.subplots(figsize=(8, 4))
    for method, rmses in rmse_data.items():
        ax.bar(method, rmses[i], width=0.4, label=f'{method} {labels[i]}')
    ax.set_title(f'RMSE for {var_title}')
    ax.set_ylabel('RMSE')
    ax.set_xlabel('Method')
    
    plt.xticks(fontsize=7)
    
    plt.show()
