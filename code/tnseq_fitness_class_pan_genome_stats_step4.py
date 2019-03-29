#!/usr/bin/env python
import pickle
import scipy
import pandas
import numpy
import matplotlib
import pandas as pd
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


full_tag_pd_with_fitness = pd.read_csv("/ebio/abt6_projects9/tnseq/tnseq_function/data/full_tag_pd_with_fitness.csv", index_col=1)
# fit = pickle.load(open("tnseq_panx_fitness_data.cpk"))
fit_reduced = full_tag_pd_with_fitness[full_tag_pd_with_fitness.columns[1:106010]]
fit_reduced.index = full_tag_pd_with_fitness.index + "_" + full_tag_pd_with_fitness['Unnamed: 0']
fit_reduced = full_tag_pd_with_fitness.iloc[0:1451]
fit_mean = fit_reduced.mean()
fit_max = fit_reduced.max()
fit_min = fit_reduced.min()
fit_var = fit_reduced.var()
time_in_tree = full_tag_pd_with_fitness.iloc[1451][1:]
Genetic_Diversity = full_tag_pd_with_fitness.iloc[1452][1:]
num_gene_events = full_tag_pd_with_fitness.iloc[1453][1:]
fit_sort = fit_reduced.ix[fit_reduced.mean(axis=1).sort_values().index]


# graph the relationship between time in tree and average fitness across all environment in every genome
x = numpy.array(time_in_tree)
y = numpy.array(fit_mean)
# plt.scatter(x, y)
# ys = lowess(x, y)
time = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"})  # , x_bins = 10, x_ci = 'sd')
time.set_xlim(0,)
time_fig = time.get_figure()
time_fig.savefig("/ebio/abt6_projects9/tnseq/tnseq_function/figs/ave_fitness_ttree.pdf")
time_fig.clf()

# graph the relationship between number of gene gains and losses and average fitness across all environment in every genome
x = numpy.array(num_gene_events)
y = numpy.array(fit_mean)
# plt.scatter(x, y)
# ys = lowess(x, y)
gain = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"}, x_bins=10, x_ci='sd')
gain.set_xlim(0,)
gain_fig = gain.get_figure()
gain_fig.savefig("/ebio/abt6_projects9/tnseq/tnseq_function/figs/ave_fitness_num_events.pdf")


# graph the relationship between number of gene gains and losses and average fitness across all environment in every genome
x = numpy.array(fit_reduced['Genetic_Diversity'])
y = numpy.array(fit_mean)
# plt.scatter(x, y)
# ys = lowess(x, y)
div = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"}))  # , x_bins = 10, x_ci = 'sd')
div.set_xlim(0,)
div_fig=div.get_figure()
div_fig.savefig("ave_fitness_gene_diversity.pdf")








# histogram of ages
plt.clf()
age_dist=sns.distplot(fit['ttree'])
age_dist.set_xlim(0,)
fig_dist=age_dist.get_figure()
fig_dist.savefig("/ebio/abt6_projects9/tnseq/tnseq_function/figs/tree_age_dist.pdf")


# Do a heatmap of all loci and fitness over all conditions
my_heatmap=sns.heatmap(fit_sort, xticklabels = False, yticklabels = False)
fig_heat=my_heatmap.get_figure()
# fig_heat.set(xlabel = "Fitness Experiment")
fig_heat.savefig("fitness_heatmap.pdf")

# Now graph the relationship between average fitness effect per locus and depth of tree
plt.clf()
sns.set_style("white")
y=fit_reduced.mean(1)
x=numpy.array(fit['ttree'])
mean_fit=sns.regplot(x, y, fit_reg = False, scatter_kws = {"color": "black"})
mean_fit.set_xlim(0,)
mean_fit.set(xlabel = "Time in Tree", ylabel = "Average fitness across environments")
mean_fit.set_xlim(0,)
fig_mean=mean_fit.get_figure()
fig_mean.savefig("ave_fitness.pdf")

# Now graph the relationship between variance in fitness effect per locus and depth of tree
plt.clf()
sns.set_style("white")
y=fit_reduced.var(1)
x=numpy.array(fit['ttree'])
mean_fit=sns.regplot(x, y, fit_reg = False, scatter_kws = {"color": "black"})
mean_fit.set(xlabel = "Time in Tree", ylabel = "Variance in fitness")
mean_fit.set_xlim(0,)
fig_mean=mean_fit.get_figure()
fig_mean.savefig("var_fitness.pdf")


plt.clf()
sns.set_style("white")
y=fit_reduced.mean(1)
x=numpy.array([int(rec) for rec in fit['gene_ev']])
mean_fit=sns.pointplot(x, y, scatter_kws = {"color": "black"}, ci = 95)
plt.plot(np.linspace(-20, 120, 1000), [0] * 1000, 'r')
mean_fit.set(xlabel = "Gene Loss/Gain Events", ylabel = "Average fitness across environments")
mean_fit.set_xlim(-0.1,)
fig_mean=mean_fit.get_figure()
fig_mean.savefig("ave_fitness_events.pdf")

plt.clf()
sns.set_style("white")
y=fit_reduced.median(1)
x=numpy.array([int(rec) for rec in fit['gene_ev']])
mean_fit=sns.pointplot(x, y, scatter_kws = {"color": "black"}, ci = 95)
plt.plot(np.linspace(-20, 120, 1000), [0] * 1000, 'r')
mean_fit.set(xlabel = "Gene Loss/Gain Events", ylabel = "Min fitness across environments")
mean_fit.set_xlim(-0.1,)
fig_mean=mean_fit.get_figure()
fig_mean.savefig("min_fitness_events.pdf")
