from pathlib import Path

#Streamlit runs the single python script, app.py in our case.
#Thus, on running the package as a script, we will not start the streamlit app.
#Instead, we will display the instructions on how to run the app,
#With the path to the app.py and additional info

if __name__ == '__main__':
    """
    Displays the information on ghow to run the Streamlit app.
    Displays information on the helper classes.
    """
    
    cwd = Path(__file__).parent.resolve()
    
    print(
        f'''This package contains the Stremlit web-app
and classes for sequence operations and plotting.\n
To run the streamlit app:
    1) Go to the Windows command line or Linux console
    2) Change directory to the app folder with:
        cd {cwd}
    3) Run the app with additional arguments for the max file size
       and port, if desired (default is 8501):
        streamlit run app.py --server.maxUploadSize 1000 --server.port 8501\n
The package contains four classes for import:
    1) FileFormatOps: determining the type of the PD table
    2) MathOps: a few math operations for subsequent plotting
    3) SequenceOps: searching for particular amino acids or modifications 
    4) PlottingOps: the main class for displaying streamlit pltos,
       inherits from MathOps and SequenceOps
    The operations are available as static or class methods,
    depending on the particular method.
    '''
    )
