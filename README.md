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
-  GMT(the Generic Mapping Tools, here version 5.3.1)  
      - GMT is used to plot the results after executing beam-forming.
      - [Download GMT](http://gmt.soest.hawaii.edu/projects/gmt/wiki/Download)
# Installation
- Just run "make" to compile and get the executable command "beamforming".
# Usage
- beamforming <sacfile.lst> <t1> <t2> <fre_low> <fre_high> <slow_low> <slow_high> <slow_step> <baz_step>  <output_file>

- NOTICE!!!
- [x] Unit of slowness: second/degree;
- [x] Unit of back-azimuth: degree.
# Contribution
-  Author: Xuping Feng
- Email: geophydogvon@gmail.com
