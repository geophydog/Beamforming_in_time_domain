:hotel:[Return to Home page](https://github.com/geophydog/geophydog.github.io)
# Beam-forming
- To constrain the back-azimuth and slowness of some seismic phase with seismic array

***

### Example
- This is an example and necessary files in the diretory "Example" 
![EXAMPLE](https://github.com/geophydog/Beamforming_in_time_domain/blob/master/images/630.000-700.000.png)

***

![results](https://github.com/geophydog/Beamforming_in_time_domain/blob/master/images/Results.jpg)
# Dependencies
- Linux or Mac OS platform;  
-  SAC(Seismic Analysis Code, here version 101.6a)  
      - SAC is used to do bandpass filtering;
      - [Request for SAC](http://ds.iris.edu/ds/nodes/dmc/forms/sac/)
-  GMT(the Generic Mapping Tools, here version 5.3.1)  
      - GMT is used to plot the results after executing beam-forming.
      - [Download GMT](http://gmt.soest.hawaii.edu/projects/gmt/wiki/Download)
# Installation
- Just run "make" to compile and get the executable command "beamforming".
# Usage
- beamforming <sacfile.lst> <t1> <t2> <fre_low> <fre_high> <slow_low> <slow_high> <slow_step> <baz_step>  <output_file> 

| parameter |  mean |
| --------- | ----- |
| sacfile.lst| file containing names of SAC format files |
|     t1     |     begin time of doing beam-forming      |
|     t2     |     end time of doing beam-forming        |
|  fre_low   |low limitation of corner frequency of SAC files |
|  fre_high  |high limitation of corner frequency of SAC files |
|  slow_low  |low limitation of scanning horizontal slowness or ray parameter|
|  slow_high |high limitation of scanning horizontal slowness or ray parameter|
|  slow_step | step length of scanning horizontal slowness or ray parameter|
|  baz_step | step length of scanning back-azimuth |
| output_file | file name of outputting results |

- NOTICE!!!
- [x] Unit of slowness: second/degree;
- [x] Unit of back-azimuth: degree.
# Contribution
-  Author: Xuping Feng
- Email: geophydogvon@gmail.com
