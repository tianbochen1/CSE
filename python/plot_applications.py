"""Plot the reported lower-upper quantile spectra from numerical summaries."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('figures/applications_spectra.pdf'))
    args=parser.parse_args()
    directory=Path(__file__).resolve().parents[1]/'results'/'reported'
    plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.titlesize':11,
                         'axes.labelsize':10,'legend.fontsize':8,'pdf.fonttype':42})
    fig,axs=plt.subplots(2,2,figsize=(9,5.6),sharex=True)
    for row,(app,title) in enumerate([('EURUSD','EURUSD'),('SP500','S&P 500')]):
        for col,part in enumerate(['real','imaginary']):
            d=pd.read_csv(directory/f'{app}_{part}_curve.csv'); ax=axs[row,col]
            x=d['lambda']/np.pi
            if app=='SP500':
                ax.fill_between(x,1000*d.lower,1000*d.upper,color='#a9c6da',alpha=.55,
                                label='95% pointwise reference envelope')
            ax.plot(x,1000*d.observed,color='#174e73',lw=1.65,label='Empirical spectrum')
            ax.plot(x,1000*d.model,color='#c35732',ls='--',lw=1.5,label='Fitted SV spectrum')
            ax.axhline(0,color='.65',lw=.5,zorder=0)
            ax.set_title(f'{title}: {part} component')
            ax.set_ylabel(r'$10^3 K(\lambda;0.1,0.9)$')
            ax.set_xlim(0,1); ax.grid(alpha=.15)
            if row==1:ax.set_xlabel(r'Integration endpoint $\lambda/\pi$')
            ax.spines[['top','right']].set_visible(False)
    handles,labels=axs[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.005))
    fig.tight_layout(rect=(0,.055,1,1))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':main()
