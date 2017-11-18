# Beam-forming
To constrain the back-azimuth and slowness of some seismic phase with seismic array
# Dependencies
- [x] SAC(Seismic Analysis Code, here version 101.6a)  
      - SAC is used to do bandpass filtering;
- [x] GMT(the Generic Mapping Tools, here version 5.3.1)  
      - GMT is used to plot the results after executing beam-forming;
# Installation
Just run "make" to compile and get the executable command "beamforming".
# Usage
beamforming <sacfile.lst> <t1> <t2> <fre_low> <fre_high> <slow_low> <slow_high> <slow_step> <baz_step>  <output_file>
# Contribution
- [x] Author: Xuping Feng
- [x] Email: geophydogvon@gmail.com
