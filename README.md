# SiPM-Amplitude-Spectrum-Fit
---
## Disclaimer
I do not guarantee, that the Code will work for all SiPM Amplitude Spectra! As of now it is very simpel and could very probably be improved a lot. Bear in mind, that it was used for a very specific task: measuring photon currents of a LED. Thus, using it for your purposes probably will require some changes to the Code.

## How to use
The Program requires a csv as input in the "data" folder with the histogram data of the Amplitude Spectrum you want to fit.
apart from that it only excepts 5 inputs in the Preamble code: The SiPM name, the first Bin to use for the Fit (must be before 0 pe Peak), The number of Peaks with good resolution (inluding 0 pe) and values for y-min and y-max for the output plots it will generate.

## Input files
The input files must be csv files with the data already in hsitogram form (binned). Scince I used them very specificly for measuring photons from a LED, it requires you to input a LED current in the csv name. For examples on how to name files, take a look in the "data" folder.  

## Output files
As output it will generate a csv file with the following columns:

|      "SiPM"      |      "I(nA)"     |      "NoP"      |                 "NoE"                 |        "mu_i"        |    "sigma2_{0}"    |      "w_{0}"     |                "photon_counts"               |
|:----------------:|:----------------:|:---------------:|:-------------------------------------:|:--------------------:|:------------------:|:----------------:|:--------------------------------------------:|
| Name of the SiPM | LED current used | Number of Peaks | Number of Events (- Background noise) | Mean value of Peak i | Variance of Peak i | Weight of Peak i | Total number of photons detected by the SiPM |

It will also generate a plot for each csv file in "data" and save them in the "results" folder.
