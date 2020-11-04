function [] = do_finc(foldername)

addpath 'C:\Users\Owner\UBC_F2020\FINC\outputs';
addpath 'C:\Users\Owner\UBC_F2020\FINC\code';
autowell(strcat('C:\Users\Owner\UBC_F2020\FINC\outputs\',foldername),'manual',1);