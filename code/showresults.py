import dnest4.classic as dn4
import matplotlib.pyplot as plt

# Run postprocess from DNest4
dn4.postprocess()

# Load posterior samples
posterior_sample = dn4.my_loadtxt("posterior_sample.txt")
indices = dn4.load_column_names("posterior_sample.txt")["indices"]

# Where time delays start and end in the file
start = indices["time_delays[1]"]
end = indices["mu_amplitude"]

# Histogram of time delays
k = 1
for i in range(start, end):
    plt.hist(posterior_sample[:,i], 100,
             alpha=0.2, label="$\\tau_{k}$".format(k=k),
             normed=True)
    k += 1
plt.legend()
plt.show()

