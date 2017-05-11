/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+  data:        array saving reading into SAC files                    +
+  t1:          beginning time of doing beam-forming                   +
+  t2:          ending time of doing beam-forming                      +
+  fre_low:     low limitation frequency of inputint SAC files         +
+  frehigh:     high limitation frequency of inputing SAC files        +
+  slow_low:    low limitation scaning hirizontal slowness             +
+  slow_high:   high limitation scaning hirizontal slowness            +
+  slow_step:   step length of scaning slowness                        +
+  baz_tesp:    step length of scaning backazimuth                     +
+  2017-5-11    Initial coding by Xuping Feng @ NJU                    +
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <time.h>
#include "sacio.h"

/*---------------------------------------------------------------------------------------------------*/
void no_spa(char *ps) {
    char *pt = ps;
    while ( *ps != '\0' ) {
        if ( *ps != ' ' && *ps != '\n' ) {
            *pt++ = *ps;
        }
        ++ps;
    }
    *pt = '\0';
}

/*---------------------------------------------------------------------------------------------------*/
int near_pow2( int n ) {
    int m;
    float f;
    f = log((float)n) / log(2.) + 1;
    m = (int) pow(2.,(int)f);
    return m;
}

/*---------------------------------------------------------------------------------------------------*/
#define pi 3.1415926  // The approximation value of PI
#define g 111.195     // The value of transformation from degree to km in the Earth


/*---------------------------------------------------------------------------------------------------*/
int main( int argc, char *argv[] ) {
    float *data, slow_low, slow_high, slow_step, slow_scan, slow_x, slow_y, baz_scan, baz_step,\
        shifttime, center_lon = 0., center_lat = 0., center_ele = 0., fre_low, fre_high, t1, t2,\
        delta, dx, dy, *sum, time_start, time_end, time_used, **coordi, **amp, peak, tmp = 0.;
    int i, j, size = 256, count = 0, sac_npts, sta_index = 0, shift_index, begin_index, end_index, beam_npts;
    char *ss;
    FILE *fin, *fout, *fp, *fbp;
    SACHEAD hd;

    if ( argc != 11 ) {
        fprintf(stderr,"Usage: beamforming <sacfile.slt> <t1> <t2> <fre_low> <fre_high> <slow_low> <slow_high> <slow_step> <baz_step> <output_file>\n");
        fprintf(stderr,"       <t1>          Beginning time of inputing SAC files;\n");
        fprintf(stderr,"       <t2>          Ending time of inputing SAC files;\n");
        fprintf(stderr,"       <sacfile.lst> File contains these names of SAC format files;\n");
        fprintf(stderr,"       <fre_low>     The low limitation of frequency of bandpass of SAC file;\n");
        fprintf(stderr,"       <fre_high>    The high limitation of frequency of bandpass of SAC file;\n");
        fprintf(stderr,"       <slow_low>    The low limitation of scaning horizontal slowness;\n");
        fprintf(stderr,"       <slow_high>   The high limitation of scaning horizontal slowness;\n");
        fprintf(stderr,"       <slow_step>   THe step length of scaning horizontal slowness;\n");
        fprintf(stderr,"       <baz_step>    The step length of scaning backazimuth;\n");
        fprintf(stderr,"       <output_file> The file name of outputing result;\n");
        exit(1);
    }

    fin = fopen(argv[1], "r");
    t1 = atof(argv[2]);
    t2 = atof(argv[3]);
    fre_low = atof(argv[4]);
    fre_high = atof(argv[5]);
    slow_low = atof(argv[6]);
    slow_high = atof(argv[7]);
    slow_step = atof(argv[8]);
    baz_step = atof(argv[9]);
    fout = fopen(argv[10],"w");
    fp = fopen("plot.sh", "w");

    time_start = clock();
    ss = (char*) malloc(size);
    while ( fgets(ss, size, fin) ) {
        no_spa(ss);
        data = read_sac(ss, &hd);
        delta = hd.delta;
        sac_npts = hd.npts;
        center_lon += hd.stlo;
        center_lat += hd.stla;
        center_ele += hd.stel;
        count +=1;
    }
    free(data); fclose(fin);

    begin_index = (int)(t1/delta);
    end_index = (int)(t2/delta);
    beam_npts = (int)((t2-t1)/delta);

    amp = (float**) malloc(sizeof(float*) * count);
    coordi = (float**) malloc(sizeof(float*) * count);
    for ( i = 0; i < count; i ++ ) {
        amp[i] = (float*) malloc(sizeof(float) * sac_npts);
        coordi[i] = (float*) malloc(sizeof(float) * 3);
    }
    sum = (float*) malloc(sizeof(float) * beam_npts);


/*---------------------bandpass filter inputing seismic data and saving them in array "amp"----------------------*/
    fin = fopen(argv[1], "r");
    while ( fgets(ss, size, fin) ) {
        no_spa(ss);
        fbp = fopen("bandpassfilter","w");
        fprintf(fbp,"sac<<END\n");
        fprintf(fbp,"r %s\n", ss);
        fprintf(fbp,"bp c %f %f n 4 p 2\n", fre_low, fre_high);
        fprintf(fbp,"w over %s\n", ss);
        fprintf(fbp,"q\n");
        fprintf(fbp,"END\n");
        fclose(fbp);
        system("sh bandpassfilter");
        data = read_sac(ss,&hd);
        coordi[sta_index][0] = center_lon - hd.stlo; coordi[sta_index][1] = center_lat - hd.stla; coordi[sta_index][2] = center_ele - hd.stel;
        for ( i = 0; i < sac_npts; i ++ ) amp[sta_index][i] = data[i];
        sta_index += 1;
    }
    free(data);fclose(fin);

/*----------------------------------------scaning slowness and backazimuth-----------------------------------------*/
    slow_scan = slow_low;

    while ( slow_scan <= slow_high ) {
        baz_scan = -185;
        while( baz_scan <= 180 ) {
            for ( i = 0; i < beam_npts; i++ ) sum[i] = 0.; peak = 0.; tmp = 0.;
            for ( i = 0; i < count; i ++ ) {
                slow_x = slow_scan*cos(baz_scan/180.*pi);
                slow_y = slow_scan*sin(baz_scan/180.*pi);
                dx = (center_lon - coordi[i][0])*g; dy = (center_lat - coordi[i][1])*g;
                shifttime = slow_x*dx + slow_y*dy;
                shift_index = (int)(shifttime/delta);
                if ( (shift_index+begin_index) >=0 && (shift_index+begin_index + beam_npts) < sac_npts ) {
                    for ( j = 0; j < beam_npts; j ++ ) {
                        //printf("AMP_INDEX: %d  AMP: %f\n", (j+shift_index+begin_index, amp[i][j+begin_index+shift_index]));
                        sum[j] += amp[i][j+begin_index+shift_index]/count;
                    }
                }
                else continue;
            }
            for ( j = 0; j < beam_npts; j ++ ) {
                if ( peak < fabs(sum[i]) ) peak = fabs(sum[i]);
                //printf("INDEX: %d  SUM: %f\n", j, sum[j]);
            }
            fprintf(fout,"%f %f %f %f %f %d\n", slow_scan, baz_scan, peak, dx, dy, shift_index);
            baz_scan += baz_step;
        }
        slow_scan += slow_step;
    }
    time_end = clock();
    time_used = (time_end - time_start)/CLOCKS_PER_SEC;
    printf("Time used: %f seconds!!!\n", time_used);

    fclose(fp); fclose(fout);
    for ( i = 0; i < count; i++ ) {
        free(coordi[i]); free(amp[i]);
    }
    free(coordi); free(amp); free(sum);
    return 0;
}
