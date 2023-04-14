import numpy as np
from emcpy.plots.plots import LinePlot, HorizontalLine, VerticalLine
from emcpy.plots.create_plots import CreatePlot, CreateFigure
import matplotlib

#filter_cutoffs = [600, 900, 960, 1200, 1800]  # km
#colors = ['black', 'black', 'red', 'black', 'black']
filter_cutoffs = [960]  # km
colors = ['black']

wavelengths = np.arange(10,3510,10)  # km
pi = 3.1415

dxs = [3, 3, 3, 3, 3]  # km

nwavelengths = len(wavelengths)
ncutoffs = len(filter_cutoffs)
h = np.zeros(shape=(ncutoffs,nwavelengths))*np.nan

for filter_cutoff in filter_cutoffs:
    i = filter_cutoffs.index(filter_cutoff)
    dx = dxs[i]
    e = (np.tan(pi*dx/filter_cutoffs[i]))**-6
    for wavelength in wavelengths:
        j = wavelengths.tolist().index(wavelength)
        h[i,j] = ( 1.0 + e*(np.tan(pi*dx/wavelength))**6 )**-1

plot = CreatePlot()
plt_list = []
x = wavelengths

for filter_cutoff in filter_cutoffs:
    i = filter_cutoffs.index(filter_cutoff)
    y = h[i,:]
    lp = LinePlot(x, y)
    lp.color = colors[i] 
    lp.linestyle = "-"
    lp.linewidth = 3
    lp.marker = None
    lp.markersize = 3
    lp.alpha = None
    lp.label = f"{filter_cutoffs[i]}"
    plt_list.append(lp)

    # Vertical Line at filter_cutoff
    lp = VerticalLine(filter_cutoffs[i])
    lp.color = colors[i]
    lp.linestyle = "--"
    lp.linewidth = 1
    lp.label = None
    plt_list.append(lp)

# Horizontal and Vertical Lines
lp = HorizontalLine(0.5)
lp.color = "black"
lp.linestyle = "-"
lp.linewidth = 1
lp.label = None
plt_list.append(lp)



yticks=np.arange(0,1.1,.1)
xticks=np.arange(0,3600,100)
x_str = [str(item) for item in xticks]
plot.plot_layers = plt_list
plot.add_ylabel
plot.add_grid()
plot.set_yticks(yticks)
plot.set_xticks(xticks)
plot.set_xticklabels(x_str, rotation=60)
plot.set_ylim(0, 1)
plot.set_xlim(0, 3500)
plot.add_xlabel("Wavelength (km)")
plot.add_ylabel("Amplitude Response")
plot.add_legend(loc="upper left", bbox_to_anchor=(1, 1), fancybox=True, framealpha=0.80)


fig = CreateFigure()
fig.plot_list = [plot]
fig.create_figure()
fig.tight_layout()
fig.save_figure("./blending_tangent_filter.png")



