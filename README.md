# FINC Code Readme

This is meant as a reference for running all code related to the Freezing Ice Nuclei Counter (FINC) instrument used at the institut fur atmosphare und klima at ETH-Zurich. The original MATLAB software was written and developed by Killian Brennan in 2019, and is currently maintained by Andrew Loeppky (2020)

Legacy Repo - all code developed prior to September 2020 can be found at:

```https://github.com/AndrewLoeppky/FINC/tree/legacy_code_2019```

Main branch of the repository is the working version of the code. Create dev branches as needed.

## Running a Droplet Freezing Experiment

Read the SOP and FINC paper for instructions on using the hardware. Once plugged in/turned on, in the matlab command prompt do:

```matlab
>> FINCrunner('measName')
```
Replace ```measName``` with something sensible that describes your experiment ie "```SAcontrol```", "```lignin5```", "```whispeak19```", etc. Note there are also optional arguments to change the number of trays/name them.

Possible example inputformats:
```matlab
FINCrunner('test3')  % runs autowell with 3 trays
FINCrunner('test3','nTrays',2)  % runs autowell with two trays
FINCrunner('test3','nTrays',2,'names',{'sample1','sample2'})% adds current date to foldername automatically
```

Experiment runs are by default ~30mins long and are completely automatic. Once finished, you should see a folder containing ```ramp.txt``` (time and temp log), a huge pile of .JPG (images of the PCR trays). To then perform the analysis, do:

```matlab
>> autowell(‘foldername’, ‘manual’, 1)
```

(The ```manual``` argument is optional, but is recommended since the automatic well detection is not completely reliable). Follow the prompts, selecting the center of the top left and bottom right wells. 