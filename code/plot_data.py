import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("j1131.txt")

def plot_subset(image):
    subset = data[:,3] == image
    plt.errorbar(data[subset,0], data[subset,1], yerr=data[subset,2],
                 fmt="o", markersize=3, alpha=0.3,
                 label="Image {image}".format(image=image))

for image in range(0, 4):
    plot_subset(image)
plt.legend()
plt.show()

