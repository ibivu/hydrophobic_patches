import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

def create_scatter(x, y, xlab, ylab, title):
    plt.scatter(x,y)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    max_point = min(max(list(x)),max(list(y)))
    plt.plot([0,max_point],[0,max_point],'--') # identity line
    plt.show()

def create_heatmap(x, y, xlab, ylab, title):
    plt.hist2d(x,y, bins=100)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.ylim(0,10000)
    plt.xlim(0,10000)
    plt.title(title)
    max_point = min(max(list(x)),max(list(y)))
    plt.plot([0,10000],[0,10000],'--') # identity line
    plt.show()

def perc_corr(y_test, pred, x):
    error = abs(y_test-pred)
    y = [sum(error/y_test < i)/len(error) for i in x]
    return(y)

def perc_perc_corr(y_test, pred, x):
    error = abs(y_test-pred)
    y = [sum(error < i)/len(error) for i in x]
    return(y)

def plot_curve(y_test, pred_dict, step, xlab, ylab, title):
    x = [x/step for x in range(step+1)]
    corr_dict = {pred:perc_corr(y_test, pred_dict[pred], x) for pred in pred_dict}

    plt.figure(figsize=(12,8))

    for corr in corr_dict:
        plt.plot([i*100 for i in x],[i*100 for i in corr_dict[corr]], label=corr)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_perc_curve(y_test, pred_dict, step, xlab, ylab, title):
    x = [x/step for x in range(step+1)]
    corr_dict = {pred:perc_perc_corr(y_test, pred_dict[pred], x) for pred in pred_dict}

    plt.figure(figsize=(12,8))

    for corr in corr_dict:
        plt.plot([i*100 for i in x],[i*100 for i in corr_dict[corr]], label=corr)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()
