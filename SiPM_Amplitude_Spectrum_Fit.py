import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.stats as stats
from scipy.signal import find_peaks
import math
from sklearn.mixture import GaussianMixture

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 15
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = 5.0
plt.rcParams['xtick.minor.size'] = 3.0
plt.rcParams['ytick.major.size'] = 5.0
plt.rcParams['ytick.minor.size'] = 3.0
plt.rcParams['figure.constrained_layout.use'] = True

## Input parameters

SiPM = "KETEK 3350-WB"

bin_0 = 20  # First Bin to use for the Fit (must be before 0 pe Peak)
NoP = 10 	# Number of Peaks
y_min = 0.5 # Maximum y-Value in output Plots for Amplitude Spectrum
y_max = 2e6 # Minimum y-Value in output Plots for Amplitude Spectrum

files_to_fit = []
out_dict = {}

## Find all histograms in directory 
for root, dirs, files in os.walk("./data"):
	 
		for file in files:
			
			if file.endswith(".csv"):
				files_to_fit.append([file, os.path.join(root, file)])


## Print files.
print("Files to fit:")
for i, (file, _) in enumerate(files_to_fit):
	print('File {0}: {1}'.format(i, file))

## Loop thorough files
for i, (file, path) in enumerate(files_to_fit):
	items = os.path.splitext(file)[0].split('_')

	I = items[2].replace('muA', '_$\\mu$A').replace('nA', '_nA').replace('dark', '0_nA')
	I = I.split('_')
	I[0] = float(I[0])



	print("\n\n")
	print("===============================================================")
	print("FIT {:d} - {:s}".format(i, file))
	print("===============================================================")
	print("\n\n")

	data = pd.read_csv(path).to_numpy()
	energy, counts = data.transpose()

	fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10.8), dpi=72, constrained_layout=True)
	transform = transforms.blended_transform_factory(ax1.transAxes, fig.transFigure)

	## For effective gain
	y_f = np.abs(np.fft.fft(counts)) ** 2
	x_f = np.fft.fftfreq(energy.size)

	g = scipy.interpolate.CubicSpline(x_f[:x_f.size//2], y_f[:x_f.size//2], bc_type='natural')
	x = np.linspace(np.min(x_f[:x_f.size//2]), np.max(x_f[:x_f.size//2]), 1000)
	
	dg = g.derivative(1)
	r = dg.roots()
	Qeff = r[np.argmax(g(r))]
	M = 1/Qeff
	print(f'effective Gain: {M}')

	ax1.set_title('\\textbf{Fourier Transform of the Amplitude Spectrum}')
	ax1.set_ylabel('Power Spectral density')
	ax1.set_xlabel('ADC counts in Fourier Space')
	ax1.grid()
	ax1.plot(x, g(x), label="Absolute Square of FFT")
	ax1.axvline(x=Qeff, color='black', ls='--', linewidth=0.8, label=f'$G_{{FFT}}^* = {M:.3f}$ ADC counts')
	ax1.legend(edgecolor='black', facecolor='white', fancybox=False, framealpha=1).get_frame().set_linewidth(0.8)
	ax1.set_yscale('log')
	ax1.set_xscale('log')

	## Primary peak position
	peaks, _ = find_peaks(counts, height=6e5)
	mu_1 = peaks[0]

	## For local minima
	invert = np.array([counts[mu_1]-counts[k] for k in energy])
	peaks, _ = find_peaks(invert[bin_0:], distance=3*M/4)
	e_min = np.insert(peaks+bin_0, 0, 0)
	c_min = counts[e_min]
	max_adc = e_min[NoP]

	f = scipy.interpolate.CubicSpline(e_min, c_min, bc_type='natural')
	x = np.linspace(0, max_adc, 1000)

	ax2.set_title('\\textbf{Estimation of 1 pe Peak and Background}')
	ax2.set_ylabel('Counts (in 30 seconds)')
	ax2.set_xlabel('ADC counts')
	ax2.set_ylim([y_min, y_max])
	ax2.fill_between(energy, counts, 0, step='mid', alpha=0.4, label="Spectrum")
	ax2.set_yscale('log')
	ax2.set_xlim([0, max_adc])
	sub_area = [f(k) if f(k)<counts[k] else counts[k] for k in energy[:max_adc+1]]
	ax2.fill_between(energy[:max_adc+1], sub_area, 0, step='mid', color='red', facecolor='None', hatch='////', label="Background")
	ax2.axvline(x=mu_1, color='black', ls='--', linewidth=0.8, label="1 pe Peak")
	ax2.legend(edgecolor='black', facecolor='white', fancybox=False, framealpha=1).get_frame().set_linewidth(0.8)

	print(f'e_min: {e_min}')
	print(f'max_adc: {max_adc}')

	## subtraction of backgroundnoise
	sub_counts = np.array([np.trunc(counts[n] - f(n)) for n in energy[:max_adc]], dtype=int)
	sub_counts[sub_counts<0] = 0

	sub_data = np.array([energy[:max_adc], sub_counts]).transpose()
	X = np.array([], dtype=int)

	for [e, c] in sub_data:
		if c>0:
			X = np.append(X, np.full(c, e, dtype=int))

	means_init = [[mu_1 + (n-1)*M] for n in range(NoP)]

	## Gaussian Mixture Model
	j = 0
	while j<2:
		print("GMM iteration {0} of 3".format(j+1))
		model = GaussianMixture(NoP, init_params='k-means++', covariance_type='spherical',
								means_init=means_init).fit(X.reshape(-1,1))

		means = model.means_
		M = means.flatten()[2] - means.flatten()[1]
		means_init = [[mu_1 + (n-1)*M] for n in range(NoP)]

		j += 1

	
	
	means_init = model.means_
	precisions_init = 1/model.covariances_
	print("GMM iteration 3 of 3")
	model = GaussianMixture(NoP, init_params='k-means++', covariance_type='spherical',
							means_init=means_init, precisions_init=precisions_init, tol=1e-6, max_iter=1000).fit(X.reshape(-1,1))


	N = len(X)
	x = np.linspace(bin_0, max_adc, 1000)
	logprob = model.score_samples(x.reshape(-1, 1))
	pdf = np.exp(logprob)

	fit = pdf*N + f(x)
	ax3.set_title('\\textbf{Fit result (GMM + Background)}')
	ax3.set_ylabel('Counts (in 30 seconds)')
	ax3.set_xlabel('ADC counts')
	ax3.set_ylim([y_min, y_max])
	ax3.fill_between(energy, counts, 0, step='mid', alpha=0.4, label="Spectrum")
	ax3.plot(x, fit, color='orange', ls='--', label="GMM + Background")
	ax3.set_xlim(0, max_adc)
	ax3.set_yscale('log')
	ax3.legend(edgecolor='black', facecolor='white', fancybox=False, framealpha=1).get_frame().set_linewidth(0.8)

	fig.suptitle("\\textbf{{{0} Amplitude Spectrum Fit ($I_{{LED}} =$ {1} {2})}}\n".format(SiPM, I[0], I[1]), transform=transform, x=0.5)
	plt.savefig("./results/{0}_fit.pdf".format(os.path.splitext(file)[0]))

	## For total number of detected photons
	peaks = model.weights_
	total_counts = sum([peaks[n]*n for n in range(len(peaks))])*N

	## For output CSV
	I[1] = I[1].replace('nA', '1').replace('$\\mu$A', '1000')
	I[1] = float(I[1])

	fit_out = {}
	fit_out["SiPM"] = SiPM
	fit_out["I(nA)"] = I[0]*I[1]
	fit_out["NoP"] = NoP
	fit_out["NoE"] = N

	mean = model.means_
	var = model.covariances_
	weight = model.weights_

	for j in range(NoP):
		fit_out["mu_{0}".format(j)] = mean.flatten()[j]
		fit_out["sigma2_{0}".format(j)] = var[j]
		fit_out["w_{0}".format(j)] = weight[j]

	fit_out["photon_counts"] = total_counts

	
	if out_dict == {}:
		for key in fit_out.keys():
			out_dict[key] = []

	for key in fit_out.keys():
		out_dict[key].append(fit_out[key])


df = pd.DataFrame.from_dict(out_dict)
df.to_csv("./fit_results_{:s}.csv".format(SiPM))



