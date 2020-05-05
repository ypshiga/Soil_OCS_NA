# Soil_OCS_NA
Python 2 code to produce soil OCS fluxes over North America using NARR inputs

The main code 'COS_yps_NA_combo_mac.py' uses inputs from NARR https://rda.ucar.edu/#!lfd?nb=y&b=proj&v=NCEP%20North%20American%20Regional%20Reanalysis that you will have to download ahead of time.
The two variables needed are soil temperature and soil volumetric moisture.
You will also need the constant NARR fields file from here: https://rda.ucar.edu/datasets/ds608.0/index.html#!sfol-wl-/data/ds608.0?g=1
These specify the vegetation type mask.


The 'ploting_beta.py' code produces images that can be used to create an animation. I used ffmpeg from my terminal with the following options

ffmpeg -framerate 8 -i Jul/frames_%d.png -c:v h264 -r 30 -s 1920x1080 ./soil_test_jul.mp4

