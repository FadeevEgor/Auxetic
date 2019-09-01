#TODO: errors
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

prediction_filename = "pred_0.3.txt"
experiment_filename = "exp1.txt"

l0 = 28.90000
S = 0.0000252
deviation_error = 0.030*5
defformation_error = deviation_error/l0

def read_experiment(filename):
	f = open(filename, 'r')
	load = []
	delta_l = []
	f.readline()
	for line in f:
		data = line.split(sep = ",")
		load.append(float(data[1]))
		delta_l.append(float(data[2]))
	return delta_l, load
	

def read_predictons(filename):
	f = open(filename, 'r')
	load = []
	delta_l = []
	for line in f:
		data = line.split(sep = "\t")
		delta_l.append(float(data[0]))
		load.append(float(data[1]))
	return delta_l, load
	

exp_delta_l, exp_load = read_experiment(experiment_filename)
exp_delta_l = np.array(exp_delta_l)
exp_load = np.array(exp_load)
exp_load -= exp_load[0]

pred_delta_l, pred_load = read_predictons(prediction_filename)
pred_delta_l = np.array(pred_delta_l)
pred_load = np.array(pred_load)
#pred_load /= 1.45


exp_deformation = exp_delta_l / l0
exp_stress = exp_load / S
exp_load_spl = UnivariateSpline(exp_delta_l, exp_load)
exp_load_derivative_spl = exp_load_spl.derivative()
exp_young_modulii = exp_load_derivative_spl(exp_delta_l) / (S * l0 / 1000)
#exp_stress_spl = UnivariateSpline(exp_deformation, exp_stress)
#exp_young_modulii_spl = exp_stress_spl.derivative()
#exp_young_modulii = exp_young_modulii_spl(exp_deformation)



pred_deformation = pred_delta_l / l0
pred_stress = pred_load / S
pred_load_spl = UnivariateSpline(pred_delta_l, pred_load)
pred_load_derivative_spl = pred_load_spl.derivative()
pred_young_modulii = pred_load_derivative_spl(pred_delta_l) / (S * l0 / 1000)
#pred_stress_spl = UnivariateSpline(pred_deformation, pred_stress)
#pred_young_modulii_spl = pred_stress_spl.derivative()
#pred_young_modulii = pred_young_modulii_spl(pred_deformation)

plt.subplot(1, 2, 1)
plt.plot(exp_deformation*100, exp_stress, pred_deformation*100, pred_stress)
#plt.plot(exp_delta_l, exp_load_spl(exp_delta_l))
#plt.plot(pred_delta_l, pred_load_spl(pred_delta_l))
plt.xlim(0, 20)
plt.ylim(0, max(exp_stress))
plt.plot((exp_deformation + defformation_error)*100, exp_stress, color="green")
plt.plot((exp_deformation - defformation_error)*100, exp_stress, color="green")
plt.xlabel("deformation (%)")
plt.ylabel("stress (Pa)")
plt.legend(["experiment", "calculation (constant)", "error"])



	
plt.subplot(1, 2, 2)
plt.plot(exp_deformation*100, exp_young_modulii)
plt.plot(pred_deformation*100, pred_young_modulii)
plt.xlim(0, 20)
plt.ylim(0., 1.3 * max(pred_young_modulii))
plt.xlabel("deformation(%)")
plt.ylabel("Young's modulii (Pa)")
plt.legend(["experiment", "calculation (constant)"])
plt.show()