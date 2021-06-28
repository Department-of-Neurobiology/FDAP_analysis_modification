@ECHO OFF 
cd %cd%
FOR /R %%G IN (CW_Analysis09raw_DecayFit03_mod01.m*) DO "C:\Program Files\MATLAB\R2020b\bin\matlab.exe" -nodisplay -nosplash -nodesktop -r "run('%%G'); exit" & pause 
