#!/usr/bin/env python

import numpy as np
import pandas as pd
from matplotlib import pyplot as pl
import sys
sys.path.append("/home/mlparker/programs/python_modules/plotting_functions/")

C1='#D81B60'
C2='#1E88E5'
C3='#FFC107'
C4='#135247'
colors=[C1,C2,C3,C4]


from step_plots import *


def spin_plot(data):
    fig=pl.figure(figsize=(7,7),facecolor='w')

    ax=fig.add_axes([0.1,0.1,0.55,0.55])
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    ax.set_xlabel('Simulated $a$')
    ax.set_ylabel('Fit $a$')

    pl.plot([0,1],[0,1],c=C1)
    pl.plot(data['a'],data['fit_a'],marker='.',ls='none',alpha=1,color='k')

    ax2=fig.add_axes([0.1,0.7,0.55,0.2])
    ax2.set_xlim(0,1)
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('N')
    pl.hist(data['a'],bins=np.linspace(0,1,11),alpha=0.3,color=C2)

    ax3=fig.add_axes([0.7,0.1,0.2,0.55])
    ax3.set_ylim(0,1)
    ax3.yaxis.set_ticklabels([])
    ax3.set_xlabel('N')
    pl.hist(data['fit_a'],bins=np.linspace(0,1,11),orientation='horizontal',alpha=0.3,color=C2)
    # fig.close()


def inc_plot(data,alpha=1):
    fig=pl.figure(figsize=(7,7),facecolor='w')

    low=15
    high=80

    ax=fig.add_axes([0.1,0.1,0.55,0.55])
    ax.set_ylim(low,high)
    ax.set_xlim(low,high)
    ax.set_xlabel('Simulated $i$')
    ax.set_ylabel('Fit $i$')

    pl.plot([low,high],[low,high],c=C1)
    pl.plot(data['i'],data['fit_i'],marker='.',ls='none',alpha=alpha,color='k')

    ax2=fig.add_axes([0.1,0.7,0.55,0.2])
    ax2.set_xlim(low,high)
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('N')
    pl.hist(data['i'],bins=np.linspace(low,high,11),alpha=0.3,color=C2)

    ax3=fig.add_axes([0.7,0.1,0.2,0.55])
    ax3.set_ylim(low,high)
    ax3.yaxis.set_ticklabels([])
    ax3.set_xlabel('N')
    pl.hist(data['fit_i'],bins=np.linspace(low,high,11),orientation='horizontal',alpha=0.3,color=C2)
    # fig.close()


def Afe_plot(data):
    fig=pl.figure(figsize=(7,7),facecolor='w')

    low=1
    high=5

    ax=fig.add_axes([0.1,0.1,0.55,0.55])
    ax.set_ylim(low,high)
    ax.set_xlim(low,high)
    ax.set_xlabel('Simulated $i$')
    ax.set_ylabel('Fit $i$')

    pl.plot([low,high],[low,high],c='r')
    pl.plot(data['Afe'],data['fit_Afe'],marker='.',ls='none',alpha=1,color='k')

    ax2=fig.add_axes([0.1,0.7,0.55,0.2])
    ax2.set_xlim(low,high)
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('N')
    pl.hist(data['Afe'],bins=np.linspace(low,high,11),alpha=0.3,color=C2)

    ax3=fig.add_axes([0.7,0.1,0.2,0.55])
    ax3.set_ylim(low,high)
    ax3.yaxis.set_ticklabels([])
    ax3.set_xlabel('N')
    pl.hist(data['fit_Afe'],bins=np.linspace(low,high,11),orientation='horizontal',alpha=0.3,color=C2)



def mdot_plot(data):
    fig=pl.figure(figsize=(7,7),facecolor='w')

    ax=fig.add_axes([0.1,0.1,0.55,0.55])
    ax.set_ylim(0.1,0.5)
    ax.set_xlim(0.1,0.5)
    ax.set_xlabel('Simulated $\dot{M}$')
    ax.set_ylabel('Fit $\dot{M}$')

    pl.plot([0,1],[0,1],c=C1)
    pl.plot(data['mdot'],data['fit_mdot'],marker='.',ls='none',alpha=1,color='k')

    ax2=fig.add_axes([0.1,0.7,0.55,0.2])
    ax2.set_xlim(0.1,0.5)
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('N')
    pl.hist(data['mdot'],bins=np.linspace(0.1,0.5,11),alpha=0.3,color=C2)

    ax3=fig.add_axes([0.7,0.1,0.2,0.55])
    ax3.set_ylim(0.1,0.5)
    ax3.yaxis.set_ticklabels([])
    ax3.set_xlabel('N')
    pl.hist(data['fit_mdot'],bins=np.linspace(0.1,0.5,11),orientation='horizontal',alpha=0.3,color=C2)
    # fig.close()

def fv_plot(data):
    length=data.shape[0]
    # noise1=np.random.randn(length)*0.25-0.125
    # noise2=np.random.randn(length)*0.25-0.125

    fig=pl.figure(figsize=(7,7),facecolor='w')

    ax=fig.add_axes([0.1,0.1,0.55,0.55])
    ax.set_ylim(0.25,1.5)
    ax.set_xlim(0.25,1.5)
    ax.set_xlabel('Simulated $f_v$')
    ax.set_ylabel('Fit $f_v$')

    pl.plot([0,2],[0,2],c=C1)
    pl.plot(data['fv'],data['fit_fv'],marker='.',ls='none',color='k',alpha=0.1)

    ax2=fig.add_axes([0.1,0.7,0.55,0.2])
    ax2.set_xlim(0.25,1.5)
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('N')
    pl.hist(data['fv'],bins=np.linspace(0.25,1.5,11),alpha=0.3,color=C2)

    ax3=fig.add_axes([0.7,0.1,0.2,0.55])
    ax3.set_ylim(0.25,1.5)
    ax3.yaxis.set_ticklabels([])
    ax3.set_xlabel('N')
    pl.hist(data['fit_fv'],bins=np.linspace(0.25,1.5,11),orientation='horizontal',alpha=0.3,color=C2)
    # fig.close()


def lx_plot(data):
    fig=pl.figure(figsize=(7,7),facecolor='w')

    ax=fig.add_axes([0.1,0.1,0.55,0.55])
    ax.set_ylim(0.5,2)
    ax.set_xlim(0.5,2)
    ax.set_xlabel('Simulated $\dot{M}$')
    ax.set_ylabel('Fit $\dot{M}$')

    pl.plot([0,2],[0,2],c=C1)
    pl.plot(data['lx'],data['fit_lx'],marker='.',ls='none',alpha=1,color='k')

    ax2=fig.add_axes([0.1,0.7,0.55,0.2])
    ax2.set_xlim(0.5,2)
    ax2.xaxis.set_ticklabels([])
    ax2.set_ylabel('N')
    pl.hist(data['lx'],bins=np.linspace(0.5,2,11),alpha=0.3,color=C2)

    ax3=fig.add_axes([0.7,0.1,0.2,0.55])
    ax3.set_ylim(0.5,2)
    ax3.yaxis.set_ticklabels([])
    ax3.set_xlabel('N')
    pl.hist(data['fit_lx'],bins=np.linspace(0.5,2,11),orientation='horizontal',alpha=0.3,color=C2)
    # fig.close()
