install.packages("RGtk2") ## the command to install RGtk2 package.
library(RGtk2) ## the command to check whether GTK-libraries has been installed in your system. 
               ## If not, a popup-menu will show up, please choose "Install GTK+" option to install them. Please also ignore those warnings after installation.
install.packages("gWidgetsRGtk2",dependencies=TRUE) ## the command to install other packages needed: gWidgetsRGtk2, gWidgets and cairoDevice.
install.packages("Bridge",repos="https://www.msu.edu/~qlu/doc/") ## the command to install Bridge.
