import dnest4.classic as dn4
import matplotlib.pyplot as plt
import numpy as np

# Run postprocess from DNest4
dn4.postprocess()

# Load posterior samples
posterior_sample = dn4.my_loadtxt("posterior_sample.txt")
indices = dn4.load_column_names("posterior_sample.txt")["indices"]


# Extract amplitudes
start = indices["amplitude[0]"]
end   = indices["period[0]"]
all_amplitudes = posterior_sample[:, start:end].flatten()
all_amplitudes = all_amplitudes[all_amplitudes != 0.0]

# Periods
start = indices["period[0]"]
end   = indices["quality[0]"]
all_periods = posterior_sample[:, start:end].flatten()
all_periods = all_periods[all_periods != 0.0]

# Extract quality factors
start = indices["quality[0]"]
end   = indices["sigma_boost_factor"]
all_qualities = posterior_sample[:, start:end].flatten()
all_qualities = all_qualities[all_qualities != 0.0]

# Trim periods
temp = all_periods.copy()
temp.sort()
use = (all_periods >= temp[int(0.01*len(temp))]) & \
      (all_periods <= temp[int(0.99*len(temp))])

# Histogram of inferred periods
plt.hist(all_periods[use], 500, color=[0.2, 0.2, 0.4])
plt.xlabel("Period")
plt.ylabel("Relative probability")
plt.show()

# Histogram of inferred periods, weighted by amplitude
plt.hist(all_periods[use], bins=500,
         weights=all_amplitudes[use], alpha=0.5)
plt.xlabel("Period")
plt.ylabel("Relative expected amplitude")
plt.show()

# Plot period vs. quality factor
plt.plot(all_periods[use],
         all_qualities[use],
         ".", markersize=1, alpha=0.1)
plt.xlabel("Period")
plt.ylabel("Quality factor")
plt.show()

# Histogram of number of modes
width = 0.7
bins = np.arange(0, posterior_sample[0, indices["max_num_components"]]+1)\
        - 0.5*width
plt.hist(posterior_sample[:, indices["num_components"]],
         bins,
         width=width,
         color=[0.2, 0.2, 0.2],
         normed=True)
plt.xlabel("num_components")
plt.show()

# Plot the marginal posterior for the error bar boost parameter
plt.hist(posterior_sample[:, indices["sigma_boost_factor"]], 100, color=[0.2, 0.2, 0.4])
plt.xlabel("sigma_boost_factor")
plt.show()

