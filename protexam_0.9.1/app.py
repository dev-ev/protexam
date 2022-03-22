import numpy as np
import pandas as pd
import streamlit as st
import json
from helpers import FileFormatOps, PlottingOps
from pathlib import Path

if __name__ == '__main__':
    st.set_page_config(page_title = 'protexam: Quality Check for Proteomics')
    st.title('Quality Check App')

    #Load the config files, which contains the column names among other stuff
    cwd = Path(__file__).parent.resolve()
    
    try:
        config_path = cwd / 'config.json'#cwd.parents[0] / 'config.json'
        print(f'Config file in {config_path}')
        with open(config_path, 'r') as j:
            PARAMS = json.loads( j.read() )
            print('Loaded configuration file:')
            print(PARAMS)
    except:
        st.warning('Could not read the config file, stopping the processing.')
        st.stop()
        
    # Upload the file as in ByteIO object
    st.markdown('''Upload a tab-separated text file from Proteome Discoverer (PD) or MaxQuant (MQ).  \nSupported PD tables are QuanSpectra, MSMSSpectrumInfo, PSMs, PeptideGroups and Proteins.  \nSupported MQ tables are msmsScans, modificationSpecificPeptides and proteinGroups:''')
    file_in = st.file_uploader(
        '', type=(['tsv','txt','tab'])
    )
    #Find the table type from the filename
    if file_in is not None:
        fname = file_in.name
        print(fname)
        ftype = FileFormatOps.determine_type(fname)
    else:
        ftype = 'Unknown'
        st.stop()
    st.write(ftype)
    #Give an option to override the type of the incoming table
    ftype_override = st.selectbox(
        'If the file type above is incorrect, choose manually here:',
        FileFormatOps.return_types(),
        0
    )
    #If the user selects "Automatic" or "Unknown", do not correct the automatically determined type
    if ftype_override != 'Automatic' and ftype_override != 'Unknown':
        ftype = ftype_override
        st.write(f'Changed the file type to {ftype}')

    to_run = False
    to_run = st.button('Generate plots')

    if to_run == True:
    #    st.write(ftype)
        with st.spinner('Reading the table.'):
            df = pd.read_csv(file_in, sep='\t')
        st.write(f'The {ftype} table has dimensions {df.shape}')
        print(f'The table {fname} has columns:')
        print(df.columns)
        with st.spinner('Drawing the plots, this may take a few seconds.'):
            PlottingOps.plot(df, ftype, PARAMS)
     
