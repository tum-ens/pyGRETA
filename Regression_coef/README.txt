		###########################
		# Regression Coefficients #
		###########################

!Required Files!

The following folder has been created:
00 Assumptions/Regression_Inputs/


Populate the Inputs folder with the required TS .csv files corresponding to the technologies specified in the config.py file.
If a file is missing, the script will prompt you to provide it. If multiple hub heights are detected, the optimization is done both over the quantiles and hubheights.

!Process!

The intersection between IRENA and EMHIRES region list is first determined. It is, based on that list, that the regions are extracted from the TS files.
If a region is present in both IRENA and EMHIRES data but not in the generated TS file, it is considered as missing data.

The existance of a solution to the constrained least square analysis is first evaluated. In case no solution exists, a place holder is set instead of the solution.
The no-solution place holder can be specified in the config.py file.

!Important Note!

The .csv files must follow the system's decimal separator setting. 
This means that the decimal separator in the .csv file might have to be changed manually from a point ('.') to a comma (',') or visversa.

!OUTPUTS!

The summary .csv files are present in the output folder:

02 Intermediate files/Files *region*/Renewable energy/*time_stamp*/Regression_Ouputs