import numpy as np
import matplotlib.pyplot as plt
import re

from scipy import ndimage
from scipy.ndimage import gaussian_filter, median_filter, uniform_filter, generic_filter

#==================================
# PARAMETERS
#==================================
factor = 2
N = factor * 7000
L = 336
L_x = factor * L
L_y =  L
delta = 200
x_wave = 20 + 30

BINS = 20
DELTA_FRAME = 30
DIR_PATH = "data/"
			
DIR_PATH = "data/delta 100/N=14000, vel=0.1, tau=1-10, delta=100/"
FILES_PATH = [DIR_PATH + f"animacion_{i}.txt" for i in range(0, 19)]

#==================================
# FUNCTIONS
#==================================
def generate_histogram(num): 
	simulations_data = []
	for file in FILES_PATH:
		try:
			simulation  = np.loadtxt(file, usecols=[0, 1, 3], skiprows=num * N, max_rows=N)
			mask = np.where(simulation[:, 2]==1) 
			simulations_data.append(simulation[mask])
		except Exception as e:
			print(f"Error: {e}")
			pass
	
	# Define the bin edges
	bins = BINS
	xedges = np.linspace(0, L_x, factor *  bins + 1)
	yedges = np.linspace(0, L_y,  bins + 1)
	heatmap = np.zeros((factor * bins,  bins))
	
	# Loop through each simulation's positions and accumulate counts in the histogram
	for sim_data in simulations_data:
		hist, _, _ = np.histogram2d(sim_data[:, 0], sim_data[:, 1], bins=[xedges, yedges])
		heatmap += hist
	
	# Save result
	np.save(DIR_PATH + f"heatmap_{num}", heatmap)
	
	# Create extent to set limits for the heatmap
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	
	x_r = np.arange(0, L_x, 5)
	y_r = lambda x : (L / delta) * (x - x_wave)

	plt.plot(x_r, y_r(x_r), color='white', linestyle="dashed", linewidth=2)
	plt.plot(x_r + 250, y_r(x_r), color='white', linestyle="dashed", linewidth=2)

	filtered_heatmap = gaussian_filter(heatmap.T, sigma=2.5)  
	#filtered_heatmap = median_filter(filtered_heatmap, size=5)  
	#filtered_heatmap = uniform_filter(heatmap.T, size=10)  
	#filtered_heatmap = generic_filter(heatmap.T, np.mean, size=10)  
	#filtered_heatmap = gaussian_filter(filtered_heatmap, sigma=1.5)  
	
	# Plot the heatmap
	c_map = 'viridis'
	c_map = 'hot'
	#plt.imshow(filtered_heatmap, extent=extent, origin='lower', cmap=c_map)  
	plt.imshow(heatmap.T, extent=extent, origin='lower', cmap=c_map, vmin=30) 
		
	plt.colorbar(label='Counts')
	plt.xlabel('X-axis')
	plt.ylabel('Y-axis')
	plt.title('Particle Distribution Heatmap')
	plt.pause(0.05)
	plt.savefig(DIR_PATH + f"heatmap_{num}.png")
	plt.clf()



if __name__ == '__main__':
	print("Angle: ", round(np.arctan(L / delta) * 180 / np.pi, 3))
	
	# Read the simulation.txt file
	file_path = DIR_PATH + 'simulation_data.txt'  # Adjust the path accordingly
	with open(file_path, 'r') as file:
		content = file.read()

	# Find the seed number using regex
	anim_step_time_match = re.search(r'anim step time \s*=\s*(\d+)', content)
	if anim_step_time_match:
		time_step = int(anim_step_time_match.group(1))

	file_path = DIR_PATH + 'evolution.txt'  # Adjust the path accordingly
	durations = np.loadtxt(file_path, usecols=3)
	final_time = int(durations.min() / time_step)
	

	FINAL_TIME = final_time
	STEP = final_time // DELTA_FRAME	
	print("Final Time: ", FINAL_TIME)
	print("STEP: ", STEP) 

	for i in range(0, FINAL_TIME, STEP):
		generate_histogram(i)