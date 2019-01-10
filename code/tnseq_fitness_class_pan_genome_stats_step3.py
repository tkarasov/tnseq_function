#!/usr/bin/env python
import pickle, scipy, pandas, numpy
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem  
  
# This function takes an array of numbers and smoothes them out.  
# Smoothing is useful for making plots a little easier to read.  
def sliding_mean(data_array, window=5):  
    data_array = numpy.array(data_array)  
    new_list = []  
    for i in range(len(data_array)):  
        indices = range(max(i - window + 1, 0),  
                        min(i + window + 1, len(data_array)))  
        avg = 0  
        for j in indices:  
            avg += data_array[j]  
        avg /= float(len(indices))  
        new_list.append(avg)  
          
    return numpy.array(new_list)  


fit = pickle.load(open("tnseq_panx_fitness_data.cpk"))
fit_reduced = fit[fit.columns[3:142]]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
fit_sort = fit_reduced.ix[fit_reduced.mean(axis=1).sort_values().index]
#graph the relationships together
x = numpy.array(fit['ttree'])
y = numpy.array(fit['set5IT028 Minimal media pH9'])
#plt.scatter(x, y)
#ys = lowess(x, y)
m9 = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"})) #, x_bins = 10, x_ci = 'sd')
m9.set_xlim(0,)
m9_fig = m9.get_figure()
m9_fig.savefig("m9_fitness.pdf")

#histogram of ages
plt.clf()
age_dist = sns.distplot(fit['ttree'])
age_dist.set_xlim(0,)
fig_dist = age_dist.get_figure()
fig_dist.savefig("tree_age_dist.pdf")


#Do a heatmap of all loci and fitness over all conditions
my_heatmap = sns.heatmap(fit_sort, xticklabels=False, yticklabels=False)
fig_heat = my_heatmap.get_figure()
#fig_heat.set(xlabel = "Fitness Experiment")
fig_heat.savefig("fitness_heatmap.pdf")

#Now graph the relationship between average fitness effect per locus and depth of tree 
plt.clf()
sns.set_style("white")
y = fit_reduced.mean(1)
x = numpy.array(fit['ttree'])
mean_fit = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"})
mean_fit.set_xlim(0,)
mean_fit.set(xlabel = "Time in Tree", ylabel = "Average fitness across environments")
mean_fit.set_xlim(0,)
fig_mean = mean_fit.get_figure()
fig_mean.savefig("ave_fitness.pdf")

#Now graph the relationship between variance in fitness effect per locus and depth of tree 
plt.clf()
sns.set_style("white")
y = fit_reduced.var(1)
x = numpy.array(fit['ttree'])
mean_fit = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"})
mean_fit.set(xlabel = "Time in Tree", ylabel = "Variance in fitness")
mean_fit.set_xlim(0,)
fig_mean = mean_fit.get_figure()
fig_mean.savefig("var_fitness.pdf")


plt.clf()
sns.set_style("white")
y = fit_reduced.mean(1)
x = numpy.array([int(rec) for rec in fit['gene_ev']])
mean_fit = sns.pointplot(x, y, scatter_kws={"color": "black"}, ci=95)
plt.plot(np.linspace(-20,120,1000), [0]*1000, 'r')
mean_fit.set(xlabel = "Gene Loss/Gain Events", ylabel = "Average fitness across environments")
mean_fit.set_xlim(-0.1,)
fig_mean = mean_fit.get_figure()
fig_mean.savefig("ave_fitness_events.pdf")

plt.clf()
sns.set_style("white")
y = fit_reduced.median(1)
x = numpy.array([int(rec) for rec in fit['gene_ev']])
mean_fit = sns.pointplot(x, y, scatter_kws={"color": "black"}, ci=95)
plt.plot(np.linspace(-20,120,1000), [0]*1000, 'r')
mean_fit.set(xlabel = "Gene Loss/Gain Events", ylabel = "Min fitness across environments")
mean_fit.set_xlim(-0.1,)
fig_mean = mean_fit.get_figure()
fig_mean.savefig("min_fitness_events.pdf")


