import numpy as np
import pandas as pd
import streamlit as st
#import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.decomposition import PCA
from bokeh.layouts import column, row
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, CustomJS, LabelSet, Select
#import re
#import json
#from pathlib import Path

class MathOps:
    """
    The class for mathematical operations of dataframes.
    """
    def __init__(self):
        pass
    
    @staticmethod
    def logtrans_noerr(df, base):
        """
        Log-transforms the data frame
        Only numeric values are transformed
        base should be a number
        Values == 0 are converted to NaN
        Non-numeric values are returned as is
        """
        def local_transformer(value_in):
            try:
                return(
                    np.log2(value_in) / np.log2(base)
                )
            except:
                return value_in
                
        return df.replace(0, np.nan).transform(local_transformer)
                
    
    @staticmethod 
    def pca_on_columns_new(
        in_df, plotname, nComp=5, compToPlot=(0,1), figWidth = 8
        ):
        """
        Takes an expression table with samples as columns and proteins in rows.
        Returns the pyplot fig object, pandas dataframe of principal components,
        pandas dataframe of PCA loadings, variance by component and the plot title.
        """
        if len(in_df.columns) <= nComp:
            nComp = len(in_df.columns) - 1
        if len(in_df.index) <= nComp:
            nComp = len(in_df.index) - 1
            
        samplenames = list(in_df.columns)
        logdata = in_df.to_numpy()

        X = logdata.transpose()

        pca = PCA(n_components=nComp, svd_solver='arpack')
        principalComponents = pca.fit_transform(X)
        
        print('Variance explained by components:')
        print(pca.explained_variance_ratio_)
        rounded_var = [
            np.round(x,3) for x in pca.explained_variance_ratio_
        ]
        plotname = plotname + f'\nVar by comps: {str(rounded_var)}'
        
        pcaLoadings = pd.DataFrame(
            pca.components_.T, columns=[ f'PC{n+1}' for n in range(nComp) ], index=in_df.index
        )
        
        pcColnames = []
        for n in range(nComp):
            pcColnames.append(
                f'PC {n+1} ({pca.explained_variance_ratio_[n]*100:.2f}%)'
            )

        principalDf = pd.DataFrame(
            principalComponents, columns = pcColnames
        )
        principalDf['sample'] = samplenames
        principalDf.set_index('sample', inplace = True)
        
        #Find min and max values on x and y axes in order to detrmine the figsize
        arrX = principalDf[ pcColnames[ compToPlot[0] ] ]
        arrY = principalDf[ pcColnames[ compToPlot[1] ] ]
        spanX = max(arrX) - min(arrX)
        spanY = max(arrY) - min(arrY)
        plotH = round(figWidth * spanY / spanX, 1)
        print(f'Figure width set to {figWidth} and height to {plotH}')
        
        fig, ax1 = plt.subplots(figsize=(figWidth, plotH))
        fig.patch.set_facecolor('white')
        principalDf.plot.scatter(
            x=pcColnames[ compToPlot[0] ], y=pcColnames[ compToPlot[1] ],
            ax = ax1, color='#213A8F',
            s=50, alpha=0.7
        )
        ax1.axhline(color='grey', alpha=0.4, linestyle='--')
        ax1.axvline(color='grey', alpha=0.4, linestyle='--')
        
        for i, row in principalDf.iloc[:, list(compToPlot) ].iterrows():
            ax1.annotate(i, row,
                        xytext=(10,-5), textcoords='offset points',
                        fontsize=12, color='#213A8F')
        ax1.set_title(plotname, fontsize=16)
        ax1.grid(False)
        
        return fig, principalDf, pcaLoadings, rounded_var, plotname

        
    @staticmethod 
    def ratio_var(df, ratio_columns, accession_col):
        """
        Calculates peptide Ratio variability.
        Takes the ratio columns and the protein accession/identifier column
        """
                
        dfN = df[ ([accession_col,] + ratio_columns) ]
            
        dfMeans = dfN.groupby([accession_col]).mean()
        dfDev = dfN.groupby([accession_col]).std()
        dfVar = np.divide(dfDev, dfMeans)
        dfVar = dfVar*100
        dfVar.replace(0, np.nan, inplace=True)
        dfVar = dfVar.round(2)            
        return dfVar            
    
    @staticmethod 
    def scale_row_mean(df):
        """
        Takes numerical dataframe.
        Scales each row on it's mean value, NaNs excluded.
        """
        
        dfOut = df.T.copy()
        dfOut = dfOut / dfOut.mean()
        dfOut = dfOut.T
        return dfOut
        
class SequenceOps:
    """
    The class for operations with peptide sequences and modifications.
    """
    def __init__(self):
        pass

    @classmethod    
    def find_cys(
        cls, df, seqc, modc,
        dict_mods: dict = {
            'carbamidomethyl': 'CAM', 'methylthio': 'MeS'
        }
    ):
        """
        Takes either the peptide sequence as it is,
        Or reads the annotated sequence with the flanking residues from PD.
        It is possible to extend the list of cysteine modifications
        by supplying a custom dictionary to the optional argument list_mods
        """
        print('Received dict_mods:')
        print(dict_mods)
        cys_cname = 'Contains Cys'
        cys_mod_cname = 'Modified Cys'
        #Check if is an annotated sequence, expressed as [K]. and .[A] in PD
        def seq_read(str_in):
            if type(str_in)== str:
                if '].' in str_in:
                    return str_in.split('.')[1]
                else:
                    return str_in
            else:
                return ''
        df[cys_cname] = [
            'Cys' if 'C' in seq_read(x) else 'None' for x in df[modc]
        ]
        
        def match_mod(str_mod, dict_mods):
            found_mods = []
            if type(str_mod) == str:
                for i in dict_mods.keys():
                    if i in str_mod:
                        found_mods.append( dict_mods[i] )
                    
            if len(found_mods) == 0:
                return 'None'
            elif len(found_mods) == 1:
                return found_mods[0]
            elif len(found_mods) > 1:
                str_out = found_mods[0]
                for i in found_mods[1:]:
                    str_out = str_out + ' and ' + i
                return str_out
        df[cys_mod_cname] = [
            match_mod(cls.to_lowercase(x), dict_mods) for x in df[modc]
        ]
                    
        return (df, cys_cname, cys_mod_cname)        

    @classmethod    
    def find_deam(cls, df, modc):
        d_cname = 'Deamidation'
        def match_deam(str_in):
            str_in = cls.to_lowercase(str_in)
            if 'deamida' in str_in:
                return 'Deam'
            else:
                return 'None'
        df[d_cname] = [ match_deam(x) for x in df[modc] ]
        return (df, d_cname)
    
    @classmethod    
    def find_ox(cls, df, modc):
        ox_cname = 'Oxidation'
        def match_ox(str_in):
            str_in = cls.to_lowercase(str_in)
            if 'oxidation [m' in str_in:
                return 'M-Ox'
            elif 'oxidatiion' in str_in:
                return 'Ox'
            else:
                return 'None'
        df[ox_cname] = [ match_ox(x) for x in df[modc] ]
        return (df, ox_cname)
    
    @classmethod    
    def find_phospho(cls, df, modc):
        phtype_cname = 'Phosphorylation'
        is_ph_cname = 'Is Phospho'
        def match_ph(str_in):
            str_in = cls.to_lowercase(str_in)
            if 'phospho (sty)' in str_in:
                return ('Ph(STY)', 1)
            elif ('phospho [s/t' in str_in) or ('phospho [t/s' in str_in) or ('phospho (st)' in str_in):
                return ('Ph(ST)', 1)
            elif ('phospho [s' in str_in) and ('phospho [t' in str_in):
                return ('Ph(ST)', 1)
            elif 'phospho [s' in str_in:
                return ('Ph(S)', 1)
            elif 'phospho [t' in str_in:
                return ('Ph(T)', 1)
            elif 'phospho [y' in str_in:
                return ('Ph(Y)', 1)
            elif 'phospho' in str_in:
                return ('Ph', 1)
            elif 'phospho' not in str_in:
                return ('None', 0)
        df[phtype_cname] = [ match_ph(x)[0] for x in df[modc] ]
        df[is_ph_cname] = [ match_ph(x)[1] for x in df[modc] ]            
        return (df, phtype_cname, is_ph_cname)
        
    @classmethod    
    def find_tmt(cls, df, modc):
        tmt_cname = 'Labeling'
        numtmt_cname = 'Is TMT/TMTpro'
        def match_tmt(str_in):
            str_in = cls.to_lowercase(str_in)
            if 'tmtpr' in str_in:
                return ('TMTpro', 1)
#            elif 'tmt pro' in str_in:
#                return ('TMTpro', 1)
            elif 'tmt pr' in str_in:
                return ('TMTpro', 1)
            elif 'tmt' in str_in:
                return ('TMT', 1)
            else:
                return ('Unlabeled', 0)
        df[tmt_cname] = [ match_tmt(x)[0] for x in df[modc] ]
        df[numtmt_cname] = [ match_tmt(x)[1] for x in df[modc] ]
        return (df, tmt_cname, numtmt_cname)
    
    @staticmethod    
    def to_lowercase(str_in):
        if type(str_in) == str:
            return str_in.lower()
        else:
            return ''

class PlottingOps(SequenceOps, MathOps):
    """
    The class for plotting operations on a dataframe.
    """
    def __init__(self):
        pass

    @staticmethod
    def bokeh_scatter(df, x_cname, y_cname, title, fig_w = 800, fig_h=600):
        """
        Interactive scatter plot in Bokeh.
        Returns the Bokeh Figure object.
        """
        cds = ColumnDataSource(data=df)
        row_name = df.index.name
        p = figure(
            title=title,
            x_axis_label = x_cname, y_axis_label = y_cname,
            tooltips = [('Row', '@' + row_name)],
            active_scroll = 'wheel_zoom',
            width=fig_w, height=fig_h
        )
        p.title.text_font_size = '14pt'
        p.xaxis.axis_label_text_font_size = '12pt'
        p.yaxis.axis_label_text_font_size = '12pt'
        p.circle(
            x_cname, y_cname, source = cds,
            size = 10, color = '#1f77b4', alpha = 0.5
        )
        labels = LabelSet(
            x = x_cname, y = y_cname, text = row_name, source = cds,
            x_offset = 3, y_offset = 3, text_font_size = '10pt' #render_mode='canvas'
        )
        p.add_layout(labels)
        return p
    
    @staticmethod
    def bokeh_scatter_select(
        df, annot_cols, abund_cols, title, fig_w = 800, fig_h=600
        ):
        """
        Interactive scatter plot in Bokeh.
        The items on axes x and y can be selected in an interactive fashion
        from the list of abund_cols.
        annot_cols: columns that are shown on the hover tooltips.
        abund_cols: abundance columns with quantitative values
        Returns the Bokeh Figure object.
        """
        annot_cols = [ x for x in annot_cols if type(x) == str ]
        cds = ColumnDataSource(
            data = df[ annot_cols + abund_cols ]
        )
        def create_annotation(x):
            if ' ' in x:
                return (x, "@{" + x + "}")
            else:
                return (x, "@" + x)
        hover_display = [ create_annotation(x) for x in annot_cols ]
        p = figure(
            title=title,
            #x_axis_label = abund_cols[0], y_axis_label = abund_cols[1],
            tooltips = hover_display,
            active_scroll = 'wheel_zoom',
            width=fig_w, height=fig_h
        )
        minval = df[abund_cols].min().mean()
        maxval = df[abund_cols].max().mean()
        
        p.line([minval, maxval], [minval, maxval], line_width = 2, alpha = 0.5)
        p.title.text_font_size = '14pt'
        p.xaxis.axis_label_text_font_size = '12pt'
        p.yaxis.axis_label_text_font_size = '12pt'
        sc = p.circle(
            abund_cols[0], abund_cols[1], source = cds,
            size = 10, color = '#1f77b4', alpha = 0.2
        )
        selectX = Select(
            title = 'Choose column to plot on the x-axis:',
            value = abund_cols[0], options = abund_cols,
        )
        selectY = Select(
            title = 'Choose column to plot on the y-axis:',
            value = abund_cols[1], options = abund_cols,
        )
        callbackX = CustomJS(
            args=dict(sc=sc, select=selectX),
            code = """sc.glyph.x.field = select.value;
            sc.glyph.x.axis_label = select.value;
                        sc.glyph.change.emit()"""
        )
        callbackY = CustomJS(
            args=dict(sc=sc, select=selectY),
            code = """sc.glyph.y.field = select.value;
                        sc.glyph.change.emit();"""
        )
        selectX.js_on_change('value', callbackX)
        selectY.js_on_change('value', callbackY)
#        p.add_layout(labels)
        return column(row(selectX, selectY), p)

    @staticmethod
    def box_plot_intensities(
        df, ab_cols, abn_cols, title_main,
        title_1, xlab_1, title_2 = '', xlab_2 = ''
    ):
        """
        Box plots for normalized and (if present) non-normalized intensities.
        ab_cols should contain alid columns, abn_cols can be empty.
        Returns the pyplot fig object.
        """
        if len(abn_cols) > 0:
            fig, (ax1, ax2) = plt.subplots(
                nrows = 2, ncols = 1,
                figsize = ( 10, 0.5 + 1*len(ab_cols) ),
                constrained_layout = True, sharex = True
            )            
        elif len(abn_cols) == 0:
            fig, ax1 = plt.subplots(
                nrows = 1, ncols = 1,
                figsize = ( 10, 0.5 + 0.5*len(ab_cols) ),
                constrained_layout = True, sharex = True
            ) 
        
        fig.suptitle(title_main, fontsize = 18)
        np.log10(
            df[ ab_cols ].replace(0, np.nan).dropna(
                axis='rows', thresh=int( len(ab_cols)/2 )
            )
        ).boxplot(
            ax=ax1, notch=True, showmeans=True, vert=False,
            boxprops= dict(linewidth=1, color='#1214AD'),
            whiskerprops= dict(color='#1214AD'),
            medianprops= dict(linewidth=2), fontsize=12
        )
        ax1.set_title(title_1, fontsize = 16)
        ax1.set_xlabel(xlab_1)
        ax1.set_ylabel('')
        if len(abn_cols) > 0:
            np.log10(
                df[ abn_cols ].replace(0, np.nan).dropna(
                    axis='rows', thresh=int( len(abn_cols)/2 )
                )
            ).boxplot(
                ax=ax2, notch=True, showmeans=True, vert=False,
                boxprops= dict(linewidth=1, color='#1214AD'),
                whiskerprops= dict(color='#1214AD'),
                medianprops= dict(linewidth=2), fontsize=12
            )
            ax2.set_title(title_2, fontsize = 16)
            ax2.set_xlabel(xlab_2)
            ax2.set_ylabel('')
        try:
            st.pyplot(fig)
        except:
            print('Could not display the streamlit plot Boxplots.')
        return fig
    
    @classmethod
    def box_plot_variability(
        cls, df, abund_cols, acc, title, xlabel='Log10 Ratio Variability(%)'
    ):
        """
        Calculates variability on the supplied abundances or abundance ratios.
        Plots the log10-transformed variability in %
        """
        #Calculate variability for each protein,
        #defined by the same identifier in the acc column. 
        dfVar = cls.ratio_var(df, abund_cols, acc)
              
        fig, ax1 = plt.subplots(
            nrows = 1, ncols = 1, figsize = ( 10, 0.5 + 0.5*len(abund_cols) ),
            constrained_layout = True
        ) 
        fig.suptitle( title, fontsize = 18)
        np.log10(dfVar).boxplot(
            ax=ax1, notch=True, showmeans=True, vert=False,
            boxprops= dict(linewidth=1, color='#6e2287'),#'#1214AD'),
            whiskerprops= dict(color='#6e2287'),#'#1214AD'),
            medianprops= dict(linewidth=2), fontsize=12
        )
        ax1.set_xlabel(xlabel, fontsize=16)
        ax1.set_ylabel('')
        try:
            st.pyplot(fig)
        except:
            print('Could not display the streamlit plot Boxplots.')
        return fig

    @staticmethod
    def convert_cvs(val):
        """
        Converts FAIMS CV values to string:
        the string representation of the number if CV is a number,
        "None" if CV is np.nan
        """
        if np.isnan(val):
            return 'None'
        else:
            return str(val)
    
    @staticmethod
    def find_sets(df, fidc):
        """
        Find separate TMT sets based on the FIle ID from Proteome DIscoverer.
        Fractions look like F1.1, F1.2 and so on, where F1 is one set.
        """
        new_cname = 'Set or File'
        if type(fidc) == str:
            #Remove the fraction numbers, if present
            finder_func = (
                lambda x: x.split('.')[0] if ('.' in x) else x
            )
            df[new_cname] = [ finder_func(x) for x in df[fidc] ]
            all_files = df[new_cname].unique()
            all_files = sorted( list( set(all_files) ) )
        else:
            df[new_cname] = ['Unknown'] * len(df.index)
            all_files = ['Unknown']
            
        return (df, all_files, new_cname)
        

    @staticmethod
    def hist_by_category(
        df, plot_c: str, iter_c: str, bins: int = 60,
        x_lab: str = '', y_lab: str = '',
        iter_list = None, title = None,
        color: str = '#1f77b4', n_cols: int = 3, plot_width = 10
    ):
        """
        Plots the content of one column as the multiple histograms,
        grouped by the variable in another column. 
        plot_c: column name to plot
        iter_c: column name to group by
        bins: number of bins
        iter_list: iterable or None. If None: iterates over unique values in the iter_c column.
        Alternatively, iterates over items in the provided list.
        title: string or None.
        color, n_cols, plot_width.
        Returns the pyplot fig object.
        """
        if iter_list is None:
            iter_list = sorted( list( df[iter_c].unique() ) )
                
        #Set number of rows for the multi-file display
        n_rows = int( np.ceil( len(iter_list)/n_cols ) )
        starting_row = 0
        if len(iter_list) > 1:
            if n_rows == 1:
                n_cols = 2
                n_rows = int( np.ceil(n_cols/2) )
                if n_rows == 1:
                    n_rows = 2
                    starting_row = 1
        elif len(iter_list) == 1:
            n_rows = 1
            n_cols = 1
        
        fig, axes = plt.subplots(
            nrows=n_rows, ncols=n_cols, figsize = (plot_width, 2*n_rows),
            constrained_layout = True, sharex = True
        )
        if title is not None:
            fig.suptitle(title, fontsize = 14)
        i = starting_row
        j = 0
        for f in iter_list:
            if len(iter_list) > 1:
                ax = axes[i][j]
            elif len(iter_list) == 1:
                ax = axes
            df[ df[iter_c] == f ][plot_c].hist(
                bins = bins, ax = ax, color = color
            )
            ax.set_title(f)
            if j == 0:
                ax.set_ylabel(y_lab)
            if i == (n_rows - 1):
                ax.set_xlabel(x_lab)
                
            j += 1
            if j == n_cols:
                j = 0
                i += 1
        try:
            st.pyplot(fig)
        except:
            print('Could not display the streamlit plot Hist by Category.')
        return fig
    
    @staticmethod
    def hist_iso_injt(df, isoc, itc, title_iso, lab_iso, title_inj, lab_inj):
        """
        Returns the pyplot fig object.
        """
        if (type(isoc) == str) or (type(itc) == str):
            fig, (ax1, ax2) = plt.subplots(
                nrows=1, ncols=2, figsize = (10, 3),
                constrained_layout = True
            )
        if type(isoc) == str:
            df[isoc].hist(
                bins = 50, ax = ax1
            )
            ax1.set_title(title_iso, fontsize = 14)
            ax1.set_xlabel(lab_iso)
        if type(itc) == str:
            df[itc].hist(
                bins = 50, ax = ax2
            )
            ax2.set_title(title_inj, fontsize = 14)
            #ax1.legend('')
            ax2.set_xlabel(lab_inj)
        try:
            st.pyplot(fig)
        except:
            print('Could not display the streamlit plot Hist Twi Items.')
        return fig

    @staticmethod
    def rename_abund(df, abc, abnc):
        """
        Used for the file QuanSpectra for TMT data
        abc: base name for abundances, abnc: base name for normalized abundances
        """
        renaming_dict = {}
        abc_updated = []
        abnc_updated = []
        for i in df.columns:
            for j in abc:
                if j in i:
                    new_name = i.replace('(', '').replace(')', '')
                    new_name = new_name.split(' ')[-1]
                    new_name = new_name.lstrip(':')
                    renaming_dict[i] = new_name
                    abc_updated.append(new_name)
            for j in abnc:
                if j in i:
                    new_name = i.replace('(', '').replace(')', '')
                    new_name = new_name.split(' ')[-1]
                    new_name = new_name.lstrip(':')
                    new_name = 'Norm_' + new_name
                    renaming_dict[i] = new_name
                    abnc_updated.append(new_name)
        return (renaming_dict, abc_updated, abnc_updated)
    
    @staticmethod
    def pie_chart(df, ax1, plot_col, aux_col, title):
        df[
            [plot_col, aux_col]
        ].groupby( [plot_col] ).count().plot(
            kind = 'pie', ax = ax1,
            y = aux_col, autopct = '%1.1f%%', fontsize = 12
        )
        ax1.set_title(title)
        ax1.legend('')
        ax1.set_xlabel('')
        ax1.set_ylabel('')
    
    @staticmethod
    def rename_ratios(
        df, abund_only = False,
        abc = [], abnc = [], abgrc = [], abscc = [], abratc = [],
        exclstr = []
    ):
        """
        Used for the file PeptideGroups for TMT data
        abund_only: bool, if True - skip looking for ratio columns,
        if False, look for ratio columns in the first place.
        abc: base name for abundances,
        abnc: base name for normalized abundances,
        abgrc: base name for grouped abundances,
        abscc: base name for scaled abundances,
        abratc: base name for abundance ratios,
        exclstr: substrings to exclude from the quantitative column search
        (if the substring is present, do not count the column as quantitative).
        """
        renaming_dict = {}
        ab_type = None
        ab_cols_new = []
        #Look for the abundance columns in the order of priority
        found_abundances = False
        
        def repetitive_elements(
            df, basenames, renaming_dict, ab_cols_new, ab_type, found_abundances, exclstr
        ):
            all_columns = []
            #Check for the exclusion strings
            for i in df.columns:
                if np.sum([ x in i for x in exclstr ]) == 0:
                    all_columns.append(i)
            if type(basenames) == str:
                basenames = [basenames]
            for i in basenames:
                if not found_abundances:
                    for j in all_columns:
                        if i in j:
                            new_name = j.split(i)[-1]
                            new_name = new_name.replace('(', '').replace(')', '').replace('  ', ' ')
                            new_name = new_name.replace(':', '')
                            new_name = new_name.lstrip(' ')
                            renaming_dict[j] = new_name
                            ab_cols_new.append(new_name)
                            found_abundances = True
            return (renaming_dict, ab_cols_new, ab_type, found_abundances)
        
        #Ratio has priority 1
        renaming_dict, ab_cols_new, ab_type, found_abundances = repetitive_elements(
            df, abratc, renaming_dict, ab_cols_new, ab_type, found_abundances, exclstr
        )
        if found_abundances:
            ab_type = 'Ratios'
        #Grouped abundance has priority 2
        if not found_abundances:
            renaming_dict, ab_cols_new, ab_type, found_abundances = repetitive_elements(
                df, abgrc, renaming_dict, ab_cols_new, ab_type, found_abundances, exclstr
            )
            if found_abundances:
                ab_type = 'Ab Grouped'
        #Scaled abundance has priority 3
        if not found_abundances:
            renaming_dict, ab_cols_new, ab_type, found_abundances = repetitive_elements(
                df, abscc, renaming_dict, ab_cols_new, ab_type, found_abundances, exclstr
            )
            if found_abundances:
                ab_type = 'Ab Scaled'
        #Normalized abundance has priority 4
        if not found_abundances:
            renaming_dict, ab_cols_new, ab_type, found_abundances = repetitive_elements(
                df, abnc, renaming_dict, ab_cols_new, ab_type, found_abundances, exclstr
            )
            if found_abundances:
                ab_type = 'Ab Norm'
        #Abundance has priority 5
        if not found_abundances:
            renaming_dict, ab_cols_new, ab_type, found_abundances = repetitive_elements(
                df, abc, renaming_dict, ab_cols_new, ab_type, found_abundances, exclstr
            )
            if found_abundances:
                ab_type = 'Ab'   
        return (renaming_dict, ab_cols_new, ab_type)
    
    @staticmethod
    def scatter_hexbin_by_category(
        df, plot_type: str,
        varx_c: str, vary_c: str, iter_c: str,
        x_lab: str = '', y_lab: str = '',
        iter_list = None, title = None,
        n_cols: int = 2, plot_width = 10
    ):
        """
        Plots the content of one column multiple scatter plots or hexbin plots;
        plot_type = 'scatter' or 'hexbin';
        grouped by the variable in another column; 
        varx_c, vary_c: column names to plot against one another;
        iter_c: column name to group by;
        bins: number of bins;
        iter_list: iterable or None. If None: iterates over unique values in the iter_c column;
        Alternatively, iterates over items in the provided list;
        title: string or None; plot_width.
        Returns the pyplot fig object.
        """
        if (plot_type != 'scatter') and (plot_type != 'hexbin'):
            print('Plot type should be scatter or hexbin.')
            
        if iter_list is None:
            iter_list = sorted( list( df[iter_c].unique() ) )
                
        #Set number of rows for the multi-file display

        n_rows = int( np.ceil( len(iter_list)/n_cols ) )
        starting_row = 0
        if len(iter_list) > 1:
            if n_rows == 1:
                n_cols = 2
                n_rows = int( np.ceil(n_cols/2) )
                if n_rows == 1:
                    n_rows = 2
                    starting_row = 1
        elif len(iter_list) == 1:
            n_rows = 1
            n_cols = 1
        
        fig, axes = plt.subplots(
            nrows=n_rows, ncols=n_cols, figsize = (plot_width, 1 + 2.5*n_rows),
            constrained_layout = True, sharex = True, sharey = True
        )
        if title is not None:
            fig.suptitle(title, fontsize = 14)
        i = starting_row
        j = 0
        for f in iter_list:
            if len(iter_list) > 1:
                ax = axes[i][j]
            elif len(iter_list) == 1:
                ax = axes
            
            if plot_type == 'scatter':
                df[ df[iter_c] == f ].plot.scatter(
                    x = varx_c, y = vary_c, ax = ax,
                    s = 3, color='navy', alpha = 0.3
                )
            elif plot_type == 'hexbin':
                dfLocal = df[ df[iter_c] == f ]
                hb = ax.hexbin(
                    x = dfLocal[varx_c], y = dfLocal[vary_c],
                    gridsize = 60,
                    cmap = 'coolwarm'
                )
                cb = fig.colorbar(hb, ax=ax)
                cb.set_label('Count')
                
            ax.set_title(f)
            if j == 0:
                ax.set_ylabel(y_lab)
            if i == (n_rows - 1):
                ax.set_xlabel(x_lab)
                
            j += 1
            if j == n_cols:
                j = 0
                i += 1

        try:
            st.pyplot(fig)
        except:
            print('Could not display the streamlit plot Hist/Hexbin by File.')
        return fig
    
    @classmethod
    def plot(cls, df, ftype: str, params_dict: dict):
        """
        Takes the dataframe, format and the column name format dictionary.
        Chooses which table-specific plotting function to use.
        Select which column names are present in the actual table,
        and passes them through to the plotting function 
        """
        if "MaxQuant-" in ftype:
            pass
        else:
            cnames = params_dict['formats']['cnames']['pd']
            cols = set(df.columns)
            #Find which version of each column name is actually used in this table
            #If the column name is found, the value of the dictionary turns into a string
            #If no column name has been matched, the corresponding values in cnames remains a list
            for i in cnames:
                for j in cnames[i]:
                    if j in cols:
                        cnames[i] = j
            
            if ftype == 'MSMS':
                cls.plot_msms(
                    df, cnames['psm_num'], cnames['raw_fname'],
                    cnames['rtime'], cnames['prec_z'],
                    cnames['prec_mz'], cnames['prec_mh'],
                    cnames['iso'], cnames['injtime'], cnames['faims_cv']
                )
            elif ftype == 'QuanSpectra':
                cls.plot_quan(
                    df, cnames['psm_num'], cnames['raw_fname'],
                    cnames['avg_sn'], cnames['ab'], cnames['ab_norm'],
                    cnames['fileid'], cnames['faims_cv']
                )
            elif ftype == 'PSMs':
                cls.plot_psms(
                    df, cnames['rtime'], cnames['mass_err'],
                    cnames['raw_fname'], cnames['se_score'], cnames['sps_match'],
                    cnames['iso'], cnames['injtime'],
                    cnames['prec_mz'], cnames['prec_mh'],
                    cnames['mods'], cnames['missed_cleav'],
                    cnames['faims_cv'], cnames['sequence'], cnames['ab'],
                    params_dict['formats']['mod_abbreviations']
                )
            elif ftype == 'Peptides':
                cls.plot_peptides(
                    df, cnames['accession'], cnames['sequence'],
                    cnames['mods'], cnames['psm_num'], cnames['missed_cleav'],
                    cnames['ab'], cnames['ab_norm'], cnames['ab_gr'],
                    cnames['ab_sc'], cnames['ab_rat'], cnames['excl_strings'],
                    params_dict['formats']['mod_abbreviations']
                )
            elif ftype == 'Proteins':
                cls.plot_proteins(
                    df, cnames['accession'], cnames['descr'], cnames['psm_num'],
                    cnames['pepts_all'], cnames['pepts_un'],
                    cnames['master_prot'], cnames['prot_id_qval'],
                    cnames['mw'], cnames['coverage'],
                    cnames['ab'], cnames['ab_norm'], cnames['ab_gr'],
                    cnames['ab_sc'], cnames['ab_rat'], cnames['excl_strings']
                )
            else:
                st.warning(f'File format {ftype} is not recognized for plotting.')

    @classmethod
    def plot_msms(
        cls, df, psmc, rawfc, rtimec,
        zc, mzc, mhc, isoc, itc, faimsc
    ):
        """
        Takes the dataframe and the column names:
        df, psmc: PSMs, rawc: RAW file name,
        rtimec: retention time in min, zc: charge,
        mzc: precursor m/z , mhc: precursor MH+
        isoc: isolation interference, itc: injection time in milliseconds
        faimsc: FAIMS compensation voltage in V
        """
        #Assert the most important column names, as we want to stop the run if they have not been supplied
        assert type(psmc) == str, f'Could not find the column for "Number of PSMs", got {psmc}'
        assert type(rawfc) == str, f'Could not find the column with RAW file names, got {rawfc}'
        assert type(mzc) == str, f'Could not find the column with precursor m/z, got {mzc}'
        assert type(rtimec) == str, f'Could not find the column with retention times, got {rtimec}'
        
        st.write(
            f'ID rate across all files is {df[psmc].mean():.03f} PSM per MSMS'
        )
        rawfiles = sorted( list(df[rawfc].unique()) )
        n_rawfiles = len( rawfiles )
        st.write(f'Found raw files {rawfiles}')
        
        #-------------------------------------------------------------------
        #st.write('ID rate by File')
        ax = df[ [rawfc, psmc] ].groupby(
            [rawfc]
        ).mean().plot.barh( figsize=( 8, 0.5 + 0.3*n_rawfiles ) )
        ax.set_title('ID Rate by File', fontsize = 14)
        ax.legend('')
        ax.set_xlabel('Average PSM per MSMS')
        ax.set_ylabel('')
        st.pyplot(ax.get_figure())
        
        #-------------------------------------------------------------------
        #Find max RT in mins and use it as the number of bins
        max_RT = int( max( df[rtimec] ) )
        #st.write('MSMS Spectra per Min by File')
        cls.hist_by_category(
            df, rtimec, rawfc, max_RT,
            'RT, min', 'MSMS per minute',
            None, 'MSMS Spectra per Min'
        )

        #-------------------------------------------------------------------
        #st.write('PSM per Min by File')
        cls.hist_by_category(
            df[ df[psmc] > 0 ],
            rtimec, rawfc, max_RT,
            'RT, min', 'PSM per minute',
            None, 'Peptide-to-Spectrum Matches per Min',
            color = '#1ca641'
        )
        
        #-------------------------------------------------------------------
        #Distributions of precursor MH+ with and without ID
        if type(mhc) == str:
            fig, (ax1, ax2) = plt.subplots(
                nrows=2, ncols=1, figsize = (10, 4),
                sharex = True, sharey = False,
                constrained_layout = True
            )
            fig.suptitle('Precursor MH+ Distribution', fontsize = 14)
            df[ df[psmc] == 0 ][mhc].hist(bins = 100, ax = ax1)
            df[ df[psmc] > 0 ][mhc].hist(bins = 100, ax = ax2, color = '#1ca641')
            ax1.set_title('MSMS Without Identifications')
            ax2.set_title('Peptide-to-Spectrum Matches')
            ax2.set_xlabel('Precursor MH+ (Da)')
            st.pyplot(fig)
        
        #-------------------------------------------------------------------
        #Historgams of ion injection times and isolation interference
        cls.hist_iso_injt(
            df, isoc, itc,
            'Precursor Isolation Interference', 'Isolation interference, %',
            'PSM Injection Time', 'Injection time, ms'
        )
        
        #-------------------------------------------------------------------
        #Historgams of ion injection time by file
        if (type(itc) == str) and (len(df[rawfc].unique()) > 1):
            cls.hist_by_category(
                df, itc, rawfc, 50,
                'Injection time, ms', 'Count',
                None, 'MSMS Injection Time by File',
                color = '#1f77b4'
            )
        
        #-------------------------------------------------------------------
        #Pie charts of precursor charge with and without ID
        if type(zc) == str:
            fig, (ax1, ax2) = plt.subplots(
                nrows=1, ncols=2, figsize = (10, 3),
                constrained_layout = True
            )
            fig.suptitle('Precursor Charge States')
            cls.pie_chart(
                df[ df[psmc] == 0 ], ax1, zc, psmc, 'MSMS Without ID'
            )
            cls.pie_chart(
                df[ df[psmc] > 0 ], ax2, zc, psmc, 'PSMs'
            )        
            st.pyplot(fig)
        
        #-------------------------------------------------------------------
        #Hexbin plots of precursor m/z vs retention time by raw file
        cls.scatter_hexbin_by_category(
            df, 'hexbin', rtimec, mzc, rawfc,
            'RT, min', 'Precursor m/z',
            None, 'Precursor m/z vs Retention Time for All MSMS'
        )
        
        #-------------------------------------------------------------------
        #Plots by FAIMS CV, if there is more than one value
        if type(faimsc) == str:
            df[faimsc] = [ cls.convert_cvs(x) for x in df[faimsc] ]
            faims_voltages = sorted( list( df[faimsc].unique() ) )
            print(f'Detected FAIMS CVs {faims_voltages}')
            if len(faims_voltages) == 1:
                if faims_voltages[0] != 0:
                    st.write(f'Detected FAIMS CV {faims_voltages[0]} V')
            elif len(faims_voltages) > 1:
                st.write(f'Detected FAIMS CVs {faims_voltages} V')
                #------------------------------------------------------------
                #ID rate by CV
                ax = df[ [faimsc, psmc] ].groupby(
                    [faimsc]
                ).mean().plot.barh( figsize=( 10, 0.2 + 0.4*len(faims_voltages) ) )
                ax.set_title('ID Rate by FAIMS CV', fontsize = 16)
                ax.legend('')
                ax.set_xlabel('Average PSM per MSMS')
                ax.set_ylabel('CV, V')
                st.pyplot(ax.get_figure())                
                #------------------------------------------------------------
                #ID rate by File and CV
                if n_rawfiles > 1:
                    ax = df[ [rawfc, faimsc, psmc] ].groupby(
                        [rawfc, faimsc]
                    ).mean().plot.barh( figsize=( 10, 0.5 + 0.3*n_rawfiles*len(faims_voltages) ) )
                    ax.set_title('ID Rate by File and FAIMS CV', fontsize = 16)
                    ax.legend('')
                    ax.set_xlabel('Average PSM per MSMS')
                    ax.set_ylabel('')
                    st.pyplot(ax.get_figure())
                
                #------------------------------------------------------------
                #MSMS or PSMs vs RT by CV
                cls.hist_by_category(
                    df,
                    rtimec, faimsc, max_RT,
                    'RT, min', 'PSM per minute',
                    None, 'MSMS per Min by FAIMS CV',
                    color = '#1f77b4', n_cols = 2
                )
                
                cls.hist_by_category(
                    df[ df[psmc] > 0 ],
                    rtimec, faimsc, max_RT,
                    'RT, min', 'PSM per minute',
                    None, 'Peptide-to-Spectrum Matches per Min by FAIMS CV',
                    color = '#1ca641', n_cols = 2
                )
                #------------------------------------------------------------
                #MSMS or PSMs histrogram vs Precursor MH+ by CV
                if type(mhc) == str:
                    cls.hist_by_category(
                        df,
                        mhc, faimsc, 50,
                        'MH+, Da', 'Count',
                        None, 'MSMS Precursor MH+ by FAIMS CV',
                        color = '#1f77b4', n_cols = 2
                    )
                    
                    cls.hist_by_category(
                        df[ df[psmc] > 0 ],
                        mhc, faimsc, 50,
                        'MH+, Da', 'Count',
                        None, 'PSM Precursor MH+ by FAIMS CV',
                        color = '#1ca641', n_cols = 2
                    )
                #------------------------------------------------------------
                #MSMS or PSMs histogram vs m/z by CV
                cls.hist_by_category(
                    df,
                    mzc, faimsc, 50,
                    'm/z', 'Count',
                    None, 'MSMS Precursor m/z by FAIMS CV',
                    color = '#1f77b4', n_cols = 2
                )
                
                cls.hist_by_category(
                    df[ df[psmc] > 0 ],
                    mzc, faimsc, 50,
                    'm/z', 'Count',
                    None, 'PSM Precursor m/z by FAIMS CV',
                    color = '#1ca641', n_cols = 2
                )
                #--------------------------------------------------------------
                #Hexbin plots of precursor m/z vs retention time by FAIMS CV
                cls.scatter_hexbin_by_category(
                    df, 'hexbin', rtimec, mzc, faimsc,
                    'RT, min', 'Precursor m/z',
                    None, 'Precursor m/z vs Retention Time for All MSMS by FAIMS CV',
                    n_cols = 2
                )
                #Hexbin plots of PSM m/z vs retention time by FAIMS CV
                cls.scatter_hexbin_by_category(
                    df[ df[psmc] > 0 ], 'hexbin', rtimec, mzc, faimsc,
                    'RT, min', 'Precursor m/z',
                    None, 'PSM Precursor m/z vs Retention Time by FAIMS CV',
                    n_cols = 2
                )
                

    @classmethod
    def plot_quan(
        cls, df, psmc, rawfc, snc,
        abc, abnc, fidc, faimsc
    ):
        """
        Takes the dataframe and the column names:
        df, psmc: PSMs, rawfc: RAW file name,
        snc: average reporter S/N
        abc: base name for abundance columns
        abnc: base name for normalized abundance columns
        fidc: File ID column from PD
        faimsc: FAIMS compensation voltage in V
        Abundances and normalized abundances are found automatically
        """
        #Assert the most important column names, as we want to stop the run if they have not been supplied
        assert type(psmc) == str, f'Could not find the column for "Number of PSMs", got {psmc}'
        assert type(rawfc) == str, f'Could not find the column with RAW file names, got {rawfc}'
        
        rawfiles = sorted( list(df[rawfc].unique()) )
        n_rawfiles = len( rawfiles )
            
        rename_dict, ab_cols, abn_cols = cls.rename_abund(df, abc, abnc)
        print(f'Detected abundance cols:\n{ab_cols}\nNorm abundance cols:\n{abn_cols}')
        df.rename(columns = rename_dict, inplace = True)
        print(f'Renamed columns:\n{df.columns}')
        
        #-------------------------------------------------------------------
        #Find if there was more than one labeling set        
        df, sets, sets_cname = cls.find_sets(df, fidc)
        print('Found separate sets or files:')
        print(sets)
        st.write(
            f'Found the following files or fractionation sets: {str(sets)}'
        )
        #Add column with set/file information, but without fraction numbers
        
        #-------------------------------------------------------------------
        #Print normalization coefficients and abundances per file/set
        for i in sets:
            #Use quan spectra with PSMs for the normalization coeffs
            #If the norm abundance is present, print the norm coefficients
            #Also check the lists of column names are of the same length
            dfLocal = df[ df[sets_cname] == i ]
            if len(ab_cols) == len(abn_cols):
                dfLocal2 = dfLocal[ df[psmc] > 0 ]
                norm_coeff = dfLocal2[abn_cols].mean().to_numpy() / dfLocal2[ab_cols].mean().to_numpy()
                for idx, j in enumerate(abn_cols):
                    st.write(f'Normalization coefficient for set {i} sample {j}/{ab_cols[idx]}' +
                            f' is {norm_coeff[idx]:.02f}')
        
            #----------------------------------------------------------------
            #Box plots of abundance and (if present), normalized abundance
            cls.box_plot_intensities(
                dfLocal, ab_cols, abn_cols,
                f'All Reporter Abundances in Set {i}',
                'Reporter Abundance', 'Log10 Reporter SN/Int',
                'Normalized Reporter Abundance', 'Log10 Norm Reporter SN/Int'
            )
            
            #----------------------------------------------------------------
            #Box plots of abundance and (if present), normalized abundance
            #This time, only for the quan spectra with PSMs
            cls.box_plot_intensities(
                dfLocal[ df[psmc] > 0 ],
                ab_cols, abn_cols,
                f'Reporter Abundances for Spectra with PSMs in Set {i}',
                'Reporter Abundance', 'Log10 Reporter SN/Int',
                'Normalized Reporter Abundance', 'Log10 Norm Reporter SN/Int'
            )
        
        #-------------------------------------------------------------------
        #Check if the SN column has been found, it will be a string in this case
        #If the snc column is not provided, calculate the mean abundance instead
        #Use the non-normalzed abundances ab_cols for the calculattion
        if type(snc) == str:
            print(f'Found Average SN in the column {snc}')
        else:
            print('Did not get Average SN column name from the config file. Calculating average sn instead.')
            snc = 'AVG Intensity Recalc'
            print(f'Reassigned the Average Reporter SN Column to {snc}')
            df[snc] =  np.nanmean(
                df[ab_cols].replace(np.nan, 0).to_numpy(),
                axis = 1
                )
            st.write(f'Did not find the Average SN column. Calculating {snc} instead.')
            
        #Reporter SN by file
        ax = df[ [rawfc, snc] ].groupby(
            [rawfc]
        ).mean().plot.barh( figsize=( 8, 0.5 + 0.3*n_rawfiles ) )
        ax.set_title('Reporter Signal-to-Noise by File', fontsize = 14)
        ax.legend('')
        ax.set_xlabel('Average Reporter SN')
        ax.set_ylabel('')
        st.pyplot(ax.get_figure())
        #----------------------------------------------------------------
        #Square root of SN with ID and without ID
        st.write('Plot the square root instead of log10 to preserve eventual 0 abundances:')
        sqrtsnc = 'Sqrt ' + snc
        df[sqrtsnc] = np.sqrt(df[snc])
        fig, (ax1, ax2) = plt.subplots(
            nrows=2, ncols=1, figsize = (10, 5),
            constrained_layout = True, sharex = True
        )
        fig.suptitle('Square Root of Average Reporter SN', fontsize = 16)
        df[ df[psmc] == 0 ][sqrtsnc].hist(bins = 100, ax = ax1)
        df[ df[psmc] > 0 ][sqrtsnc].hist(bins = 100, ax = ax2, color = '#1ca641')
        ax1.set_title('Without Identifications', fontsize = 14)
        ax2.set_title('With Peptide-to-Spectrum Matches', fontsize = 14)
        ax2.set_xlabel('Sqrt Average Reporter SN')
        st.pyplot(fig)
        #----------------------------------------------------------------
        #Plots by FAIMS CV, if there is more than one value
        if type(faimsc) == str:
            df[faimsc] = [ cls.convert_cvs(x) for x in df[faimsc] ]
            faims_voltages = sorted( list( df[faimsc].unique() ) )
            print(f'Detected FAIMS CVs {faims_voltages}')
            if len(faims_voltages) == 1:
                if faims_voltages[0] != 0:
                    st.write(f'Detected FAIMS CV {faims_voltages[0]} V')
            elif len(faims_voltages) > 1:
                st.write(f'Detected FAIMS CVs {faims_voltages} V')
                #------------------------------------------------------------
                #SN by CV
                ax = df[ [faimsc, snc] ].groupby(
                    [faimsc]
                ).mean().plot.barh( figsize=( 10, 0.2 + 0.4*len(faims_voltages) ) )
                ax.set_title('Average Reporter SN (All Spectra) by FAIMS CV', fontsize = 16)
                ax.legend('')
                ax.set_xlabel('Average Reporter SN')
                ax.set_ylabel('')
                st.pyplot(ax.get_figure())
                #------------------------------------------------------------
                #Histograms of Sqrt SN by CV
                #Quan spectra without ID
                cls.hist_by_category(
                    df[ df[psmc] == 0 ],
                    sqrtsnc, faimsc, 50,
                    'Sqrt Reporter SN', 'Count',
                    None, 'Reporter SN by FAIMS CV for Spectra Without ID',
                    color = '#1f77b4', n_cols = 2
                )                
                #Quan spectra with PSMs
                cls.hist_by_category(
                    df[ df[psmc] > 0 ],
                    sqrtsnc, faimsc, 50,
                    'Sqrt Reporter SN', 'Count',
                    None, 'Reporter SN by FAIMS CV for Spectra With PSM',
                    color = '#1ca641', n_cols = 2
                )
                
    @classmethod
    def plot_psms(
        cls, df, rtimec, merrc,
        rawfc, sescorec, spsc, isoc, itc, mzc, mhc,
        modc, mcleavc, faimsc, seqc, abc,
        dict_mods
    ):
        """
        Takes the dataframe and the column names:
        df, spsc: SPS match percentage
        rtimec: retention time in min, merrc: precursor mass error in ppm
        rawfc: RAW file name, sescorec: seach engine score,
        abc: base name for abundance columns
        isoc: isolation interference, itc: injection time in milliseconds
        mzc: precursor m/z, mhc: precursor MH+, modc: modifications
        mcleavc: number of missed cleavages
        faimsc: FAIMS compensation voltage in V
        seqc: peptide sequence
        Abundances are found automatically
        dict_mods: dictionary {type: {string to search: name of ptm}}
        """
        #Assert the most important column names, as we want to stop the run if they have not been supplied
        assert type(merrc) == str, f'Could not find the column with precursor mass errors, got {merrc}'
        assert type(mzc) == str, f'Could not find the column with precursor m/z, got {mzc}'
        assert type(rtimec) == str, f'Could not find the column with retention times, got {rtimec}'
        
        #Check if there are different values of FAIMS CV
        #If there is more than one value, use them for plotting
        to_plot_faims = False
        if type(faimsc) == str:
            df[faimsc] = [ cls.convert_cvs(x) for x in df[faimsc] ]
            faims_voltages = sorted( list( df[faimsc].unique() ) )
            print(f'Detected FAIMS CVs {faims_voltages}')
            if len(faims_voltages) > 1:
                to_plot_faims = True
        #-------------------------------------------------------------------
        #PSMs per file  
        n_rawfiles = 1
        if (type(rawfc) == str):
            st.write(f'Found raw files {sorted(list(df[rawfc].unique()))}')
            n_rawfiles = len( list(df[rawfc].unique()) )
            ax = df[ [rawfc, merrc] ].groupby(
                [rawfc]
            ).count().plot.barh( figsize=( 8, 0.5 + 0.4*n_rawfiles ) )
            ax.set_title('Number of PSMs per File', fontsize = 16)
            ax.legend('')
            ax.set_xlabel('PSMs')
            ax.set_ylabel('')
            st.pyplot(ax.get_figure())
        else:
            st.write('Could not identify the column with the RAW file names!')
        #-------------------------------------------------------------------
        if to_plot_faims:
            #PSMs by FAIMS CV
            ax = df[ [faimsc, merrc] ].groupby(
                [faimsc]
            ).count().plot.barh( figsize=( 8, 0.5 + 0.3*len(faims_voltages) ) )
            ax.set_title('Number of PSMs by FAIMS CV', fontsize = 16)
            ax.legend('')
            ax.set_xlabel('PSMs')
            ax.set_ylabel('CV, V')
            st.pyplot(ax.get_figure())
            #PSMs per file and FAIMS CV 
            if n_rawfiles > 1:
                ax = df[ [rawfc, faimsc, merrc] ].groupby(
                    [rawfc, faimsc]
                ).count().plot.barh( figsize=( 10, 0.5 + 0.3*n_rawfiles*len(faims_voltages) ) )
                ax.set_title('Number of PSMs by File and FAIMS CV', fontsize = 16)
                ax.legend('')
                ax.set_xlabel('PSMs')
                ax.set_ylabel('')
                st.pyplot(ax.get_figure())
        #-------------------------------------------------------------------
        #If there were missed cleavages, show the average number of MCs per file
        if (type(mcleavc) == str) and (type(rawfc) == str):
            m_cleav_values = list(df[mcleavc].unique())
            if len(m_cleav_values) == 1:
                st.write(f'Only PSMs with {m_cleav_values[0]} missed cleavages were detected.')
            elif len(m_cleav_values) > 1:
                st.write(f'PSMs with {str(m_cleav_values)} missed cleavages were detected.')
                ax = df[ [rawfc, mcleavc] ].groupby(
                    [rawfc]
                ).mean().plot.barh(
                    color = '#87158f', figsize=( 8, 0.5 + 0.4*n_rawfiles )
                )
                ax.set_title('Average Missed Cleavages per PSM by File', fontsize = 14)
                ax.legend('')
                #ax.set_xlim(0,1)
                ax.set_xlabel('Missed Cleavages per PSM')
                ax.set_ylabel('')
                st.pyplot(ax.get_figure())
                #------------------------------------------------------------
                #Average missed cleavages by FAIMS CV
                if to_plot_faims:
                    ax = df[ [faimsc, mcleavc] ].groupby(
                        [faimsc]
                    ).mean().plot.barh(
                        color = '#87158f', figsize=( 8, 0.5 + 0.4*len(faims_voltages) )
                    )
                    ax.set_title('Average Missed Cleavages per PSM by FAIMS CV', fontsize = 14)
                    ax.legend('')
                    #ax.set_xlim(0,1)
                    ax.set_xlabel('Missed Cleavages per PSM')
                    ax.set_ylabel('CV, V')
                    st.pyplot(ax.get_figure())
        #-------------------------------------------------------------------
        #Find if TMT-labeled, plot percentage of labeled PSMs by file          
        if type(modc) == str:
            df, tmt_cname, tmt_n_cname = cls.find_tmt(df, modc)
            all_lab_values = list(df[tmt_cname].unique())
            if len(all_lab_values) == 1:
                st.write(f'Only PSMs with {all_lab_values[0]} labeling were detected.')
            elif len(all_lab_values) > 1:
                st.write(f'PSMs with {all_lab_values} labeling were detected.')
                st.write(f'Average fraction of labeled PSMs was {df[tmt_n_cname].mean():.03f} across the whole data set.')
                ax = df[ [rawfc, tmt_n_cname] ].groupby(
                    [rawfc]
                ).mean().plot.barh(
                    color = '#ab5411', figsize=( 8, 0.5 + 0.4*n_rawfiles )
                )
                ax.set_title('Average Fraction of Labeled PSMs by File', fontsize = 14)
                ax.legend('')
                ax.set_xlim(0,1)
                ax.set_xlabel('Labeled PSMs / All PSMs')
                ax.set_ylabel('')
                st.pyplot(ax.get_figure())
            #-------------------------------------------------------------------
            #Cysteine-containing PSMs per file
            print('Checking Cys')
            if (type(seqc) == str) and (type(rawfc) == str):
                df, cys_cname, cys_mod_cname = cls.find_cys(
                    df,seqc,modc,dict_mods['cys']
                )
                all_cys_values = list(df[cys_mod_cname].unique())
                if len(all_cys_values) == 1:
                    st.write(f'Only PSMs with {all_cys_values[0]} modification on Cys were detected.')
                elif len(all_cys_values) > 1:
                    #Calculate the fraction of the modified Cys PSMs
                    dfGrCys = df[
                        df[cys_mod_cname] != 'None'
                    ][
                        [rawfc, cys_mod_cname, modc]
                    ].groupby( [rawfc, cys_mod_cname] ).count()
                    dfGrCys.columns = ['Cys Modified']
                    dfGrAll = df[[rawfc, modc]].groupby( [rawfc] ).count()
                    dfGrAll.columns = ['All']
                    dfGrCys = dfGrCys.join(dfGrAll, how='left')
                    dfGrCys['Fraction Mod'] = dfGrCys['Cys Modified'] / dfGrCys['All']
                    ax =  dfGrCys[['Fraction Mod']].plot.barh(
                        #color = '#ab5411',
                        figsize=( 8, 0.5 + 0.4*n_rawfiles )
                    )
                    ax.set_title('PSM Fraction by Cys Modification and File', fontsize = 14)
                    ax.legend('')
                    #ax.set_xlim(0,1)
                    ax.set_xlabel('Cys Modified PSMs / All PSMs')
                    ax.set_ylabel('')
                    st.pyplot(ax.get_figure())
                #st.dataframe(df[
                #    [cys_mod_cname,rawfc,modc]
                #].groupby(
                #    [rawfc,cys_cname,cys_mod_cname]
                #).count())
                #st.write('Checked Cys')
                #st.dataframe(df[[cys_cname,cys_mod_cname,rawfc,modc]])
                #st.dataframe(df[[cys_mod_cname,rawfc]])
            #---------------------------------------------------------------
            #Find deamidations. If there are deamidated peptides, plot scatters Mass error vs m/z
            #For deamidated peptides and peptides without deamidation
            df, deam_cname = cls.find_deam(df, modc)
            all_deam_values = list(df[deam_cname].unique())
            if 'Deam' in all_deam_values:
                if len(all_deam_values) == 1:
                    fig, ax1 = plt.subplots(
                        nrows=1, ncols=1, figsize = (10, 5),
                        constrained_layout = True
                    )
                    fig.suptitle('Mass Accuracy vs Precursor m/z for Deamidated PSMs', fontsize = 14)
                    df[ df[deam_cname] == 'Deam' ].plot.scatter(
                        x = mzc, y = merrc, ax = ax1,
                        s = 4, color='navy', alpha = 0.2
                    )
                    ax1.set_xlabel('Precursor m/z')
                    ax1.set_ylabel('Mass error, ppm')
                    st.pyplot(fig)
                else:
                    fig, (ax1, ax2) = plt.subplots(
                        nrows=1, ncols=2, figsize = (10, 5),
                        constrained_layout = True, sharex = True, sharey = True
                    )
                    fig.suptitle('Mass Accuracy vs Precursor m/z', fontsize = 16)
                    df[ df[deam_cname] == 'Deam' ].plot.scatter(
                        x = mzc, y = merrc, ax = ax1,
                        s = 4, color='navy', alpha = 0.2
                    )
                    ax1.set_title('Deamidated PSMs', fontsize = 16)
                    ax1.set_xlabel('Precursor m/z')
                    ax1.set_ylabel('Mass error, ppm')
                    df[ df[deam_cname] != 'Deam' ].plot.scatter(
                        x = mzc, y = merrc, ax = ax2,
                        s = 4, color='blue', alpha = 0.2
                    )
                    ax2.set_title('Other PSMs', fontsize = 16)
                    ax2.set_xlabel('Precursor m/z')
                    ax2.set_ylabel('Mass error, ppm')
                    st.pyplot(fig)
        #-------------------------------------------------------------------
        # Is the spsc is defined, plot the distribution of SPS matches:
        if type(spsc) == str:
            ax = df[
                [spsc, rtimec]
            ].groupby([spsc]).count().plot.bar(
                figsize = (10, 3), rot = 0
            )
            ax.set_title('SPS Match % Counts', fontsize = 16)
            ax.set_xlabel('SPS Match (%)')
            ax.set_ylabel('')
            ax.legend('')
            st.pyplot(ax.get_figure())
        #-------------------------------------------------------------------
        #Historgams of ion injection times and isolation interference
        cls.hist_iso_injt(
            df, isoc, itc,
            'Precursor Isolation Interference', 'Isolation interference, %',
            'PSM Injection Time', 'Injection time, ms'
        )
        #-------------------------------------------------------------------
        #Historgams of the search engine score and precursor MH+
        cls.hist_iso_injt(
            df, sescorec, mhc,
            'Search Engine Score', sescorec,
            'PSM Precursor MH+ Masses', 'Precursor MH+ (Da)'
        )
        #-------------------------------------------------------------------
        #Scatter search engine score vs precursor MH+
        if (type(sescorec) == str) and (type(mhc) == str):
            ax = df.plot.scatter(
                x = mhc, y = sescorec, s = 2, color = 'navy', alpha = 0.3,
                figsize = (10, 3)#, logy = True
            )
            ax.set_title('Search Engine Score by Precursor MH+', fontsize = 16)
            ax.set_xlabel('Precursor MH+, Da')
            ax.set_ylabel(sescorec)
            ax.legend('')
            st.pyplot(ax.get_figure())
            
        #-------------------------------------------------------------------
        #Box plot of PSM abundances
        rename_dict, ab_cols, _ = cls.rename_abund(df, abc, [])
        print(f'Detected abundance cols:\n{ab_cols}')
        df.rename(columns = rename_dict, inplace = True)
        print(f'Renamed columns:\n{df.columns}')
        if len(ab_cols) > 0:
            cls.box_plot_intensities(
                df,
                ab_cols, [],
                'Reporter Abundances', '', 'Log10 Reporter SN/Int',
                None, None
            )
        
        #-------------------------------------------------------------------
        #Scatter plots of mass error vs retention time by file        
        cls.scatter_hexbin_by_category(
            df, 'scatter', rtimec, merrc, rawfc,
            'RT, min', 'Mass error, ppm',
            None, 'Mass Error vs Retention Time'
        )
        
        #-------------------------------------------------------------------
        #Hexbin plots of mass error vs retention time by file        
        cls.scatter_hexbin_by_category(
            df, 'hexbin', rtimec, merrc, rawfc,
            'RT, min', 'Mass error, ppm',
            None, 'Mass Error vs Retention Time'
        )
        
        #-------------------------------------------------------------------
        #Hexbin plots of precursor m/z vs retention time by raw file
        if type(mzc) == str:
            cls.scatter_hexbin_by_category(
                df, 'hexbin', rtimec, mzc, rawfc,
                'RT, min', 'Precursor m/z',
                None, 'PSM Precursor m/z vs Retention Time'
            )

    @classmethod
    def plot_peptides(
        cls, df, acc, seqc, modc, psmc,
        mcleavc,
        abc, abnc, abgrc, abscc, abratc, exclstr,
        dict_mods
    ):
        """
        df, acc: Uniprot accession or other protein identifier,
        seqc: peptide sequence, modc: modifications, psmc: PSMs,
        mcleavc: number of missed cleavages,
        dict_mods: dictionary {type: {string to search: name of ptm}}
            type can be "cys" or "other"
        """
        #Assert the most important column names, as we want to stop the run if they have not been supplied
        assert type(acc) == str, f'Could not find the column with unique protein identifiers, got {acc}'
        assert type(seqc) == str, f'Could not find the column with peptide sequences, got {seqc}'
        assert type(modc) == str, f'Could not find the column with peptide modifications, got {modc}'

        #Annotate important modifications with separate columns
        df, tmt_cname, tmt_n_cname = cls.find_tmt(df, modc)
        df, cys_cname, cys_mod_cname = cls.find_cys(df,seqc,modc,dict_mods['cys'])
        df, ox_cname = cls.find_ox(df, modc)
        df, phtype_cname, is_ph_cname = cls.find_phospho(df, modc)
        df, deam_cname = cls.find_deam(df, modc)

        #Find abundance columns, with the Ratio as the first priority
        renaming_dict, abund_cols, abund_type = cls.rename_ratios(
            df, False, abc, abnc, abgrc, abscc, abratc, exclstr
        )
        print('Abundance columns:')
        print(abund_cols)
        print('Abundance type:')
        print(abund_type)
        df.rename(columns=renaming_dict, inplace=True)
        if abund_type is None:
            st.write('Did not find abundance columns in the table.')
        else:
            n_pepts = len(df.index)
            print(df.shape)
            dfNum = df[abund_cols].replace(0, np.nan).dropna(
                axis = 'columns', thresh = int(n_pepts/3)
            ).dropna(
                axis = 'rows'
            ).copy()
            filtered_abund_cols = [ x for x in abund_cols if x in dfNum.columns ]
            print(f'Numeric column after filtering for PCA has dimensions {str(dfNum.shape)}')
            #Principal components to be plotted on the plt figure
            x_PC_index = 0
            y_PC_index = 1
            fig, df_PCs, _, var_byPC, pca_title = cls.pca_on_columns_new(
                np.log2(dfNum[filtered_abund_cols]),
                f'PCA on Log2 {abund_type}',
                nComp=5, compToPlot=(x_PC_index,y_PC_index), figWidth = 8
            )
            

            fig_w = 800
            fig_h = 500
#            print(f'PCA plot height is {fig_h}')
            st.write(f'Interactive PCA plot PC1 vs PC2 on {abund_type}. Use mouse scroll to zoom in and out:')
            fig_bokeh = cls.bokeh_scatter(
                df_PCs, df_PCs.columns[x_PC_index], df_PCs.columns[y_PC_index],
                pca_title,
                fig_w = fig_w, fig_h=fig_h
                )
            st.bokeh_chart(fig_bokeh, use_container_width = True)
            #Plot components 3 and 4 if the component 3 has more than 10% of variation
            if len(var_byPC) >= 4:
                if var_byPC[2] >= 0.1:
                    x_PC_index = 2
                    y_PC_index = 3
                    st.write(f'Interactive PCA plot PC3 vs PC4 on {abund_type}. Use mouse scroll to zoom in and out:')
                    fig_bokeh = cls.bokeh_scatter(
                        df_PCs, df_PCs.columns[x_PC_index], df_PCs.columns[y_PC_index],
                        pca_title,
                        fig_w = fig_w, fig_h=fig_h
                        )
                    st.bokeh_chart(fig_bokeh, use_container_width = True)
#            st.write('The same PCA plot as above, but in the static format:')
#            st.pyplot(fig)
            
            if abund_type == 'Ratios':
                #-------------------------------------------------------------------
                #If the quan columns are ratios, calculate the peptide ratio variabilites
                #And put them on the box plot
                cls.box_plot_variability(
                    df, abund_cols, acc,
                    'Peptide Ratio Variability (%) within Proteins\n1.0=10%, 1.25=17.8%, 1.5=31.6%, 2.0 = 100%',
                    'Log10 (%) Ratio Variability'
                )
                #----------------------------------------------------------------------
                #Pearson correlations on ratios
                st.write(f'Pearson correlation on Peptide Log2 {abund_type}')
                fig, ax1 = plt.subplots( figsize=( 10, 9 ) )
                sns.heatmap(
                    np.log2(dfNum).corr(method='pearson').round(2),
                    ax = ax1, square=True, cmap='coolwarm', vmin=-1, vmax=1,
                    annot=False, cbar=True, linewidth=2, linecolor='white'
                )
                fig.suptitle(f'Pearson correlation on Peptide Log2 {abund_type}', fontsize=16)
                st.pyplot(fig)
                
                #-------------------------------------------------------------------
                #Look for abundance columns again, skipping Ratio columns this time
                renaming_dict, new_abund_cols, abund_type = cls.rename_ratios(
                    df, True, abc, abnc, abgrc, abscc, abratc, exclstr
                )
                df.rename(columns=renaming_dict, inplace=True)
                if abund_type is None:
                    st.write('Did not find abundance columns in the table.')
                    #If there were no abundance values, revert to the ratios
                    abund_type == 'Ratios'
                else:
                    n_pepts = len(df.index)
                    abund_cols = new_abund_cols
                    dfNum = df[abund_cols].replace(0, np.nan).dropna(
                        axis = 'columns', thresh = int(n_pepts/3)
                    ).dropna(
                        axis = 'rows'
                    ).copy()
                    filtered_abund_cols = [ x for x in abund_cols if x in dfNum.columns ]
                    st.write(f'Numeric value array has size {str(dfNum.shape)}')
                    x_PC_index = 0
                    y_PC_index = 1
                    fig, df_PCs, _, var_byPC, pca_title = cls.pca_on_columns_new(
                        np.log2(dfNum[filtered_abund_cols]),
                        f'PCA on Log2 {abund_type}',
                        nComp=5, compToPlot=(x_PC_index,y_PC_index), figWidth = 8
                    )
                    st.write(f'Interactive PCA plot PC1 vs PC2 on {abund_type}. Use mouse scroll to zoom in and out:')
                    fig_bokeh = cls.bokeh_scatter(
                        df_PCs, df_PCs.columns[x_PC_index], df_PCs.columns[y_PC_index],
                        pca_title,
                        fig_w = fig_w, fig_h=fig_h
                    )
                    st.bokeh_chart(fig_bokeh, use_container_width = True)
    #                st.pyplot(fig)
                    
            if abund_type != 'Ratios':
                #----------------------------------------------------------------
                #Re-scale abundances on Mean, if they are not Ratios 
                dfScaled = cls.scale_row_mean(df[abund_cols])
                print('Abundances after scaling:')
                print(dfScaled.dropna(axis='rows').head(3))
                dfScaled[acc] = df[acc]
                #Boxplot variabilities on the re-scaled abundances
                cls.box_plot_variability(
                    dfScaled, abund_cols, acc,
                    f'Peptide Variability (%) within Proteins on Re-Scaled {abund_type}\n1.0=10%, 1.25=17.8%, 1.5=31.6%, 2.0 = 100%',
                    f'Log10 (%) Variability on Re-Scaled {abund_type}'
                )
            
                #----------------------------------------------------------------------
                #Pairplot re-scaled abundances
                n_pepts = len(df.index)
                dfNum = dfScaled[abund_cols].replace(0, np.nan).dropna(
                    axis = 'columns', thresh = int(n_pepts/3)
                ).dropna( axis = 'rows' ).copy()
                #Pearson correlations on re-scaled abundances
                st.write(f'Pearson correlation on Peptide Log2 Re-Scaled {abund_type}')
                fig, ax1 = plt.subplots( figsize=( 10, 9 ) )
                sns.heatmap(
                    np.log2(dfNum).corr(method='pearson').round(2),
                    ax = ax1, square=True, cmap='coolwarm', vmin=-1, vmax=1,
                    annot=False, cbar=True, linewidth=2, linecolor='white'
                )
                fig.suptitle(f'Pearson correlation on Peptide Log2 Re-Scaled {abund_type}', fontsize=16)
                st.pyplot(fig)
        
        #----------------------------------------------------------------------
        #Pie charts of numbers of modified peptides
        fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(
            nrows=2, ncols=3, figsize = (10, 8),
            constrained_layout = True
        )
        fig.suptitle('Peptide Numbers by Modification/Cleavage', fontsize = 16)
        cls.pie_chart(df, ax1, ox_cname, seqc, 'Oxidation')
        cls.pie_chart(df, ax2, cys_mod_cname, seqc, 'Cys Modifications')
        cls.pie_chart(df, ax3, phtype_cname, seqc, 'Phosphopeptides')
        cls.pie_chart(df, ax4, deam_cname, seqc, 'Deamidated')
        cls.pie_chart(df, ax5, tmt_cname, seqc, 'TMT/TMTpro Peptides')
        if type(mcleavc) == str:
            cls.pie_chart(df, ax6, mcleavc, seqc, 'Missed Cleavages')
        
        st.pyplot(fig)
        if abund_type is not None:
            dfNum = df[
                [ox_cname, cys_mod_cname, phtype_cname] + abund_cols
            ].replace(0, np.nan).copy()
            print(dfNum.shape)
            dfNum[mcleavc] = df[mcleavc]
            print(dfNum.shape)
            dfNum = dfNum.dropna(
                axis = 'columns', thresh = int(n_pepts/3)
            ).dropna(
                axis = 'rows'
            ).copy()
            print(dfNum.shape)
            filtered_abund_cols = [ x for x in abund_cols if x in dfNum.columns ]
            st.write(f'Numeric value array has size {str(dfNum.shape)}')
            
            #----------------------------------------------------------------
            #Peptide abundance plots by modification state
            for i in (ox_cname, cys_mod_cname, phtype_cname, mcleavc):
                #Plot only if there is more than one type of modification per category
                if len(dfNum[i].unique()) > 1:
                    #Boxplots grouped by sample and modification category
                    num_rows = len(filtered_abund_cols) * len(dfNum[i].unique())
                    fig, ax1 = plt.subplots(
                        nrows = 1, ncols = 1,
                        figsize = ( 10, 0.5 + 0.4*num_rows ),
                        constrained_layout = True, sharex = True
                    ) 
        
                    fig.suptitle(
                        f'Peptide {abund_type} by Sample and {i} State', fontsize = 18
                    )
                    dfLocal = df.copy()
                    #Transform numeric missed cleavages into strings
                    if i == mcleavc:
                        dfLocal[i] = [ str(int(x)) for x in dfLocal[i] ]
                        
                    dfLog = pd.melt(
                        dfLocal[ [i,] + filtered_abund_cols ],
                        id_vars=[i], value_vars=filtered_abund_cols,
                        var_name='Sample', value_name=' '
                    )
                    dfLog[' '] = np.log2(dfLog[' '].replace(0, np.nan))
                    print('Log-trans')
                    print(dfLog)
                    print(dfLog.groupby(by=['Sample',i]).mean())
                        
                    dfLog.groupby(['Sample',i]).boxplot(
                        ax=ax1, subplots=False,
                        notch=True, showmeans=True, vert=False,
                        boxprops= dict(linewidth=1, color='#174e57'),
                        whiskerprops= dict(color='#174e57'),
                        medianprops= dict(linewidth=2), fontsize=12
                    )
#                    ax1.set_title(, fontsize = 16)
                    ax1.set_xlabel(f'Log2 Peptide {abund_type}',fontsize=16)
                    ax1.set_ylabel('')
                    #Re-scale, as abundance usually has huge outliers
                    center_scale = dfLog[' '].mean()
                    old_span = dfLog[' '].max() - dfLog[' '].min()
                    new_half_span = old_span / 6
                    ax1.set_xlim(
                        center_scale - new_half_span, center_scale + new_half_span
                    )
                    st.pyplot(fig)
                    del fig
                        
            #----------------------------------------------------------------
            #Relative abundances of peptides with missed cleavages
            mc_states = sorted( list( dfNum[mcleavc].unique() ) )
            if (len(mc_states) > 1) and (mc_states[0] == 0):
                mc_abunds = {'Description': [], 'Value': []} 
                for i in filtered_abund_cols:
                    mean_0mc = dfNum[
                        dfNum[mcleavc] == 0
                    ][i].mean()
                    for j in mc_states[1:]:
                        mean_local = dfNum[
                            dfNum[mcleavc] == j
                        ][i].mean()
                        mc_abunds['Description'].append(
                            (f'Ratio of average {abund_type} {j}mc/0mc for {i} ')
                        )
                        mc_abunds['Value'].append(np.round(mean_local/mean_0mc, 2))
#                        st.write(
#                            f'Ratio of average {abund_type} {j}mc/0mc for {i} is {mean_local/mean_0mc:.2f}'
#                        )
                df_mc_abund = pd.DataFrame(
                    {'Value': mc_abunds['Value'] },
                    index = mc_abunds['Description']
                )
                ax = df_mc_abund.plot.barh( figsize=( 8, 0.5*len(mc_abunds['Description']) ) )
                ax.set_title('Average Abundances by Missed Cleavage', fontsize = 14)
                ax.legend('')
                ax.set_xlabel('Average Abundance Ratio')
                ax.set_ylabel('')
                ax.grid(visible = True)
                st.pyplot(ax.get_figure())
                    
                    
                

    @classmethod
    def plot_proteins(
        cls, df, acc: str, descc, psmc, peptalls, peptunc,
        masterc, idqvalc, mwc, covc,
        abc, abnc, abgrc, abscc, abratc, exclstr
    ):
        """
        acc: protein accession or unique identifier, descc: protein description,
        psmc: number of psms, peptalls: number of peptides,
        peptunc: number of unique peptides,
        masterc: Master protein column,
        idqvalc: identification q-value column,
        mwc: molecular weight, covc: sequence coverage,
        exclstr: strings to exclude as quantitative columns.
        """
        #Assert the most important column names, as we want to stop the run if they have not been supplied
        assert type(acc) == str, f'Could not find the column with unique protein identifiers, got {acc}'
        #assert type(covc) == str, f'Could not find the column with peptide sequences, got {covc}'
        #assert type(modc) == str, f'Could not find the column with peptide modifications, got {modc}'
        
        if type(masterc) == str:
            print('Values in the Master column:', df[masterc].unique())
            if df[masterc].isin(['IsMasterProtein']).sum() > 0:
                df = df[ df[masterc].isin(['IsMasterProtein']) ]
                print(f'Master protein column {masterc}')
                print(f"Protein table after filtering for Master proteins has size {str(df.shape)}")
#            print(df[masterc].unique())
        if type(idqvalc) == str:
            df = df[df[idqvalc] <= 0.05]
            print(f"Protein table after filtering for q-value has size {str(df.shape)}")
        st.write(f"Protein table after filtering for q-value and Master proteins has size {str(df.shape)}")
        #Find abundance columns, with the Ratio as the first priority
        renaming_dict, abund_cols, abund_type = cls.rename_ratios(
            df, False, abc, abnc, abgrc, abscc, abratc, exclstr
        )
        print('Abundance columns:')
        print(abund_cols)
        print('Abundance type:')
        print(abund_type)
        df.rename(columns=renaming_dict, inplace=True)
        #----------------------------------------------------------------
        #Interactive table of proteins with most PSMs
        st.write('Proteins with highest PSM count:')
        annot_cols = [
            x for x in [acc,descc,psmc,peptalls,peptunc,mwc,covc] if type(x) == str
        ]
        dfLocal = df[annot_cols + abund_cols].sort_values(by=[psmc], ascending=False).copy()
        if len(dfLocal.index) > 10:
            dfLocal = dfLocal.iloc[:10, :]
        st.dataframe(dfLocal)
        if abund_type is None:
            st.write('Did not find abundance columns in the table.')
        else:

            #----------------------------------------------------------------
            #Principal Component Analysis
            dfNum = df[abund_cols].replace(0, np.nan).dropna(
                axis = 'columns', thresh = int(len(df.index)/3)
            ).dropna(
                axis = 'rows'
            ).copy()
            filtered_abund_cols = [ x for x in abund_cols if x in dfNum.columns ]
            print(f'Numeric dataframe after filtering for PCA has dimensions {str(dfNum.shape)}')
            #Principal components to be plotted on the plt figure
            x_PC_index = 0
            y_PC_index = 1
            fig, df_PCs, _, var_byPC, pca_title = cls.pca_on_columns_new(
                np.log2(dfNum[filtered_abund_cols]),
                f'PCA on Log2 {abund_type}',
                nComp=5, compToPlot=(x_PC_index,y_PC_index), figWidth = 8
            )
            fig_w = 800
            fig_h = 500
            st.write(f'Interactive PCA plot PC1 vs PC2 on {abund_type}. Use mouse scroll to zoom in and out:')
            fig_bokeh = cls.bokeh_scatter(
                df_PCs, df_PCs.columns[x_PC_index], df_PCs.columns[y_PC_index],
                pca_title,
                fig_w = fig_w, fig_h=fig_h
                )
            st.bokeh_chart(fig_bokeh, use_container_width = True)
            #Plot components 3 and 4 if the component 3 has more than 10% of variation
            if len(var_byPC) >= 4:
                if var_byPC[2] >= 0.1:
                    x_PC_index = 2
                    y_PC_index = 3
                    st.write(f'Interactive PCA plot PC3 vs PC4 on {abund_type}. Use mouse scroll to zoom in and out:')
                    fig_bokeh = cls.bokeh_scatter(
                        df_PCs, df_PCs.columns[x_PC_index], df_PCs.columns[y_PC_index],
                        pca_title,
                        fig_w = fig_w, fig_h=fig_h
                        )
                    st.bokeh_chart(fig_bokeh, use_container_width = True)
            #----------------------------------------------------------------------
            #Pearson correlations
            st.write(f'Pearson correlation on Protein Log2 {abund_type} without clustering')
            fig, ax1 = plt.subplots( figsize=( 10, 9 ) )
            dfCorr = np.log2(dfNum).corr(method='pearson')
            sns.heatmap(
                dfCorr.round(2),
                ax = ax1, square=True, cmap='coolwarm', vmin=-1, vmax=1,
                annot=False, cbar=True, linewidth=2, linecolor='white'
            )
            fig.suptitle(f'Pearson correlation on Protein Log2 {abund_type}', fontsize=16)
            st.pyplot(fig)
            #----------------------------------------------------------------------
            #Clustered Pearson correlations
            st.write(f'Clustered Pearson correlations on Protein Log2 {abund_type}')
            fig, ax1 = plt.subplots( figsize=( 10, 9 ) )
            clust_indices = leaves_list( linkage(dfCorr, 'ward') )
            clust_cols = [ dfCorr.columns[x] for x in clust_indices ]
            
            sns.heatmap(
                dfCorr.loc[clust_cols, clust_cols].round(2),
                ax = ax1, square=True, cmap='coolwarm', vmin=-1, vmax=1,
                annot=False, cbar=True, linewidth=2, linecolor='white'
            )
            fig.suptitle(f'Clustered Pearson correlations on Protein Log2 {abund_type}', fontsize=16)
            st.pyplot(fig)
            #----------------------------------------------------------------
            #Pairplot on the abundances.
#            dfNum = df[abund_cols].replace(0, np.nan).dropna(
#                axis = 'columns', thresh = int(len(df.index)/3)
#            ).dropna( axis = 'rows' ).copy()
            if len(dfNum.columns) > 20:
                st.write(f'There were more than 20 columns, skipping the pairplot.')
            elif len(dfNum.columns) > 0:
                st.write(f'Pairplot on Protein Log2 {abund_type}')
                pair_plot = sns.pairplot(
                    np.log2( dfNum ),
                    diag_kind='kde', kind='scatter', markers=',', height = 1.5,
                    plot_kws={'s': 8,  'alpha':0.7, 'color': '#1f77b4'}
                )
                
                fig = pair_plot.figure
                fig.suptitle(f'Pairplot on Protein Log2 {abund_type}', fontsize=20)
                pair_plot.tight_layout()
                st.pyplot(fig)
            #----------------------------------------------------------------
            #Interactive scatter plot in Bokeh
            #User can choose which columns to plot
            st.write('Interactive scatter plot. Choose what to display on x-axis and y-axis using the two selector widgets below:')
            dfLocal = np.log2( df[abund_cols].replace(0, np.nan) ).copy()
            for i in annot_cols:
                dfLocal[i] = df[i]
            fig_bokeh = cls.bokeh_scatter_select(
                dfLocal, annot_cols, abund_cols,
                f'Log2 {abund_type} Scatter Plot on Proteins'
            )
            st.bokeh_chart(fig_bokeh, use_container_width = True)
            

        
class FileFormatOps:
    """
    The class for operations with specified formats.
    Static methods only.
    """
    def __init__(self):
        pass
    
    @staticmethod
    def determine_type(fname: str):
       if '_Proteins.txt' in fname:
           return 'Proteins'
       elif '_PeptideGroups.txt' in fname:
           return 'Peptides'
       elif '_PSMs.txt' in fname:
           return 'PSMs'
       elif '_MSMSSpectrumInfo.txt' in fname:
           return 'MSMS'
       elif '_QuanSpectra.txt' in fname:
           return 'QuanSpectra'
       else:
           return 'Unknown'
    
    @staticmethod
    def return_types():
        return (
            'Automatic', 'Unknown', 'Proteins',
            'Peptides', 'PSMs', 'MSMS', 'QuanSpectra'
        )
