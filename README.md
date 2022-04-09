# protexam
## Inspect your quantitative proteomics results using this streamlit-powered dashboard. The app is specifically tailored for comprehensive examination of result files from isobaric labeling-based quantitative experiments.

*protexam* has been tested on Python 3.8 and 3.10 on Windows 10 and Python 3.8 on Ubuntu 20.04. *protexam* is dependent on multiple other packages, such as *streamlit*, *bokeh*, *matplotlib*. Because of that, we strongly advise creating a virtual environment that will contain the package and all of it's dependencies. One option is [*venv*](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments), which is a part of Python's standard library.\
First, we need to choose a path and create a virtual environment. On Linux, create the environment using the command line, with your custom path and project name:

```
python3 -m venv /your_project_path/project_name
```

Now we can activate the newly created environment:

```
source /your_project_path/project_name/bin/activate
```

If the environment has been activated correctly, it's name should be displayed in the console.

```
(project_name) user:~$
```

It is a good idea to update pip within the virtual environment before installing other packages:

```
pip3 install --upgrade pip
```

We can now install *protexam*. The installation will take a while due to the long list of dependencies.

```
pip3 install protexam
```

The dashboard is based on *streamlit*, so we will need to run the script *app.py* in order to deploy the dashboard. The module will help us to find the path to the *app.py* file. Run the module as follows:

```
python3 -m protexam
```

The output in the console will contain a few notes, including the path to the source file. We should change the working directory using the suggested path:

```
cd /path_to_app_source
```

Finally, it's time to run the *streamlit* app:

```
streamlit run app.py --server.maxUploadSize 1000 --server.port 8501
```

where:\
    server.maxUploadSize - max size of a single text file in MB\
    server.port - port for the app
    
We would need to modify the above-mentioned commands if we are using a Windows machine:

```
>python -m venv /your_project_path/project_name
>your_project_path/project_name/Scripts/activate
>pip install --upgrade pip
>pip install protexam
>python -m protexam
>cd path_to_app_source
>streamlit run app.py --server.maxUploadSize 1000 --server.port 8501
``` 

Right after launching the app, *streamlit* will give us a hint on how to access the dashboard in a web browser. To view locally, go to:
```
http://localhost:8501
```

Use the IP address to view the dashboard in the local network:
```
http://ip-of-the-server:8501
```
If the app is running correctly, we will see the file upload widget:

<img src="https://github.com/dev-ev/protexam/blob/main/img/app_screenshot_01.png" width="700">

The dashboard has been specifically adapted for the [Proteome Discoverer (PD)](https://www.thermofisher.com/se/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html) tab-delimited output files. However, a few of the MaxQuant tables are also supported. There are two options for creating tab-delimited text output in PD 2.4:
* add the Result Exporter node to the consensus workflow. We usually choose "True" for the "R-Friendly header" option, but the app works with "False" as well. The resulting text files can be found in the study folder.
* if the workflows have been already completed, open the PD result file and go to the menu File -> Export -> To Text (tab delimited).

The files are uploaded and visualized one-by-one. Choose one of the supported tables and upload. The app will try to determine the type of the uploaded file and show that below the upload widget. If the type has been inferred correctly, press the "Generate Plots" button. If we know the type of the table, but it has not been determined automatically, we can override it manually by selecting from the list, an then pressing the "Generate Plots" button.

<img src="https://github.com/dev-ev/protexam/blob/main/img/app_screenshot_02.png" width="700">

Displaying all of the plots can take time, but we can start scrolling down while the images are being rendered:

<img src="https://github.com/dev-ev/protexam/blob/main/img/app_screenshot_03.png" width="700">
