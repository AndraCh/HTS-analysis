# -*- coding: utf-8 -*-
import pyodbc
import pandas as pd
import seaborn as sns
from ptitprince import PtitPrince as pt
import matplotlib.pyplot as plt

colour_palette = 'plasma'

Connstr_icd = ... # constructor for SQL connection
PATH_SAVE_IMG = ... # path to save the final image
TARGET = ... # target
TABLE = ... # SQL table
# Source: https://github.com/pog87/PtitPrince/blob/master/RainCloud_Plot.ipynb
   
def get_nmlogEC50_data(table, target):
    db_query = "SELECT mlogEC50, hits FROM "+ table
    df = pd.read_sql(db_query, pyodbc.connect(Connstr_icd)) 
    df['Target'] = [target for i in range(0, len(df))]
    return df

def process_data_nmlogEC50_targets(targets):
    frames = []
    for target in targets:
        table = TABLE
        frame = get_nmlogEC50_data(table, target)
        frames.append (frame)
    df = pd.concat(frames)
    df = df[df['mlogEC50']>=0]
    return df

def create_rain_plot (targets, img_name):
    if type(targets) != list:
        targets = [targets]
    df = process_data_nmlogEC50_targets(targets)
    f, ax = plt.subplots(figsize=(int(1.5*len(targets)), 10))

    dy = "mlogEC50"; dx = "Target"; ort = "v"
    # Draw a violinplot with a narrower bandwidth than the default
    ax=pt.half_violinplot(data = df, palette=colour_palette, bw=.2,  linewidth=1,cut=0.,\
                       scale="area", width=.8, inner=None,orient=ort,x=dx,y=dy)
    ax=sns.stripplot(data=df[df['hits']=='positive'], palette=colour_palette, edgecolor="white",size=4,orient=ort,\
                     x=dx,y=dy,jitter=1,zorder=0)
    ax=sns.stripplot(data=df[df['hits']=='negative'],palette=colour_palette,  edgecolor="white",size=2, alpha=0.5, orient=ort,\
                     x=dx,y=dy,jitter=1,zorder=0)
    ax=sns.boxplot(data=df, color="black",orient=ort,width=.15,  fliersize=0,x=dx,y=dy,zorder=10,\
                  showcaps=True,boxprops={'facecolor':'none', "zorder":10},\
                   showfliers=True,whiskerprops={'linewidth':2, "zorder":10},saturation=1)
    ax.set_ylabel ('-logEC50').set_fontsize('12')
    ax.set_xlabel ('', weight='bold').set_fontsize('12')
    #ax = sns.pointplot(x=dx, y=dy, data=ddf,color='red')
    # Finalize the figure
    #ax.set(ylim=(3.5, -.7))
    sns.despine(left=True)
    plt.savefig(PATH_SAVE_IMG + img_name + '.pdf', dpi=300, bbox_inches='tight') 
    plt.show()
create_rain_plot (TARGET,'Reactivity_boxplot_target' )
