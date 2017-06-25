/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+  data:        array saving reading into SAC files                    +
+  t1:          beginning time of executing beam-forming               +
+  t2:          ending time of executing beam-forming                  +
+  fre_low:     low limitation corner frequency of bandpassfilter      +
+  frehigh:     high limitation corner frequency of bandpass filter    +
+  slow_low:    low limitation of scaning hirizontal slowness          +
+  slow_high:   high limitation of scaning hirizontal slowness         +
+  slow_step:   step length of scaning slowness                        +
+  baz_tesp:    step length of scaning backazimuth                     +
+  2017-5-11    Initially coded by Xuping Feng @ NJU                   +
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <time.h>
#include "sacio.h"

/*-------------------------------------------------to exclude space or new line chracter of a string-------------------------------------------------------*/
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

/*---------------------------------------------------------to get 2^k( k is an integer ), and m <= 2^k ----------------------------------------------------*/
int near_pow2( int n ) {
    int m;
    float f;
    f = log((float)n) / log(2.) + 1;
    m = (int) pow(2.,(int)f);
    return m;
}

/*------------------------------------to calculate cross-coeffient of two  inputing arrays------------------------------------------------------------------*/
float coef( float *data1, float *data2, int npts ) {
    int i;
    float cof, x_y = 0., x_sqa = 0., y_sqa = 0., xy_sqa;
    for ( i = 0; i < npts; i ++ ) {
        x_sqa += data1[i]*data1[i];
        y_sqa += data2[i]*data2[i], x_y += data1[i]*data2[i];
    }
	xy_sqa = sqrt(x_sqa*y_sqa);
	if ( fabs(xy_sqa) == 0. ) {
		cof = 0.;
	}
	else {
		cof = x_y / xy_sqa;
	}
	return fabs(cof);
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------*/
#define pi 3.1415926  // The approximation value of PI
#define g 111.195     // The value of transformation from degree to km in the Earth


/*-----------------------------------------------------MAIN PROGRAM-----------------------------------------------------------------------------------------*/
int main( int argc, char *argv[] ) {
    float *data, slow_low, slow_high, slow_step, slow_scan, slow_x, slow_y, baz_scan, baz_step, grid, grid_step,\
        shifttime, center_lon = 0., center_lat = 0., center_ele = 0., fre_low, fre_high, t1, t2,\
        delta, dx, dy, *sum, time_start, time_end, time_used, **coordi, **amp, cof = 0., **tmp, cof_peak = 0., cof_low = 1.;
    int i, j, size = 256, count = 0, sac_npts, sta_index = 0, shift_index, begin_index, end_index, beam_npts, k = 0;
    char *ss, ch[16];
    FILE *fin, *fout, *fp, *fbp;
    SACHEAD hd;

    if ( argc != 11 ) {
	    fprintf(stderr,"***********************************************************************************************************************************\n");
        fprintf(stderr,"** Usage: beamforming <sacfile.lst> <t1> <t2> <fre_low> <fre_high> <slow_low> <slow_high> <slow_step> <baz_step> <output_file>   **\n");
        fprintf(stderr,"**       <t1>          Beginning time of inputing SAC files;                                                                     **\n");
        fprintf(stderr,"**       <t2>          Ending time of inputing SAC files;                                                                        **\n");
        fprintf(stderr,"**       <sacfile.lst> File contains these names of SAC format files;                                                            **\n");
        fprintf(stderr,"**       <fre_low>     The low limitation of frequency of bandpass of SAC file;                                                  **\n");
        fprintf(stderr,"**       <fre_high>    The high limitation of frequency of bandpass of SAC file;                                                 **\n");
        fprintf(stderr,"**       <slow_low>    The low limitation of scaning horizontal slowness;                                                        **\n");
        fprintf(stderr,"**       <slow_high>   The high limitation of scaning horizontal slowness;                                                       **\n");
        fprintf(stderr,"**       <slow_step>   THe step length of scaning horizontal slowness;                                                           **\n");
        fprintf(stderr,"**       <baz_step>    The step length of scaning backazimuth;                                                                   **\n");
        fprintf(stderr,"**       <output_file> The file name of outputing result; Containing 3 columns Col1: baz  Col2: slowness  Col3: cross-coeffient  **\n");
        fprintf(stderr,"**       ATTENTION!!! plot shell script will be saved in file \"plot.sh\", and just run command \"sh plot.sh!\";                     **\n");
        fprintf(stderr,"**       RUN \"sh plot.sh\" REQUIRES GMT(the Generic Mapping Tools, major version 5).                                              **\n");
	    fprintf(stderr,"***********************************************************************************************************************************\n");
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

	center_lon /= count; center_lat /= count; center_ele /= count;

    begin_index = (int)(t1/delta);
    end_index = (int)(t2/delta);
    beam_npts = (int)((t2-t1)/delta);

    amp = (float**) malloc(sizeof(float*) * count);
    coordi = (float**) malloc(sizeof(float*) * count);
    tmp = (float**) malloc(sizeof(float*) * count);
    for ( i = 0; i < count; i ++ ) {
        amp[i] = (float*) malloc(sizeof(float) * sac_npts);
        coordi[i] = (float*) malloc(sizeof(float) * 3);
        tmp[i] = (float*) malloc(sizeof(float) * beam_npts);
    }
    sum = (float*) malloc(sizeof(float) * beam_npts);


/*---------------------call SAC(Seismic Analysis Code) to bandpass filter inputing seismic data and saving them in array "amp"----------------------*/
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
        system("rm bandpassfilter");
        data = read_sac(ss,&hd);
        coordi[sta_index][0] = hd.stlo; coordi[sta_index][1] = hd.stla; coordi[sta_index][2] = hd.stel;
        for ( i = 0; i < sac_npts; i ++ ) amp[sta_index][i] = data[i];
        sta_index += 1;
    }
    free(data);fclose(fin);

/*-------------------------------------------------------------------scaning slowness and backazimuth------------------------------------------------*/
    slow_scan = slow_low;

    while ( slow_scan <= slow_high ) {
        baz_scan = -50;
        while( baz_scan <= 365 ) {
            for ( i = 0; i < beam_npts; i++ ) {
                sum[i] = 0., cof = 0.;
            }
            for ( i = 0; i < count; i ++ )
                for ( j = 0; j < beam_npts; j ++ )
                    tmp[i][j] = 0.;
            for ( i = 0; i < count; i ++ ) {
                slow_x = slow_scan*cos((90 - baz_scan)/180.*pi);
                slow_y = slow_scan*sin((90 - baz_scan)/180.*pi);
	            dx = (center_lon - coordi[i][0])*g; dy = (center_lat - coordi[i][1])*g;
                //dx = (coordi[i][0] - center_lon)*g; dy = (coordi[i][1] - center_lat)*g;
                shifttime = slow_x*dx + slow_y*dy;
                shift_index = (int)(shifttime/delta);
				//printf("SHIFT_INDEX : %d  dx:%f  dy: %f slow_x: %f  slow_y: %f,  center_lon: %f  center_lat: %f %f %f\n", shift_index, dx, dy, slow_x, slow_y, center_lon, center_lat, coordi[i][0], coordi[i][1]);
                if ( (shift_index+begin_index) >= 0 && (shift_index+begin_index + beam_npts) < sac_npts ) {
                    for ( j = 0; j < beam_npts; j ++ ) {
                        sum[j] += amp[i][j+begin_index+shift_index]/count;
                        tmp[i][j] = amp[i][j+begin_index+shift_index];
                    }
                }
                else if ( (shift_index+begin_index) < 0 ) {
	                for (j = 0; j < beam_npts; j ++) {
                        if ( (j+begin_index+shift_index) < 0 ) {
                        	sum[j] += 0.; tmp[i][j] = 0.;
                        }
                    	else {
                        sum[j] += amp[i][j+begin_index+shift_index];
                        tmp[i][j] = amp[i][j+begin_index+shift_index];
                    	}
              	    }
	            }
                else if ( (shift_index+begin_index+beam_npts) >= sac_npts ) {
                    for ( j = 0; j < beam_npts; j ++ ) {
	                    if( (shift_index+begin_index+j) >= sac_npts ) {
                        	sum[j] += 0.; tmp[i][j] = 0.;
                    	}
                    	else  {
                        	sum[j] += amp[i][shift_index+begin_index+j];
                        	tmp[i][j] = amp[i][j+begin_index+shift_index];
                    	}
                	}
                }
                else continue;
            }
            for ( i = 0; i < count; i ++ ) {
                cof += coef(sum, tmp[i], beam_npts)/count;
            }
			if ( cof_peak < cof ) cof_peak = cof;
			if ( cof_low > cof ) cof_low = cof;
            fprintf(fout,"%f %f %f\n", baz_scan, slow_scan, cof);
			//sprintf(ch,"%d.sac",k);
			//k += 1; hd.b = t1; hd.e = t2; hd.npts = beam_npts;
			//write_sac(ch,hd,sum);
            baz_scan += baz_step;
        }
        slow_scan += slow_step;
    }

/*-------------------------------------------------data virtialization shell script saved in file "plot.sh"-----------------------------------------------*/
	grid = slow_high;
	grid_step = grid/10.;
	fprintf(fp,"gmt gmtset FONT_LABEL 15p,Times-Bold,black\n");
	fprintf(fp,"gmt gmtset FONT_TITLE 20p,27,black\n");
	fprintf(fp,"gmt gmtset MAP_TITLE_OFFSET 40p\n");
	fprintf(fp,"gmt gmtset MAP_GRID_PEN_PRIMARY 1p,white,\"..\"\n");
	fprintf(fp,"gmt gmtset MAP_GRID_PEN_SECONDARY 1p,white,\"..\"\n");
	if ( slow_low >= 0.2*slow_high ) {
		fprintf(fp,"R1=-5/365/%f/%f\n", slow_low, slow_high);
		fprintf(fp,"R2=0/360/%f/%f\n", slow_low, slow_high);
	}
	else {
		fprintf(fp,"R1=-5/365/%f/%f\n", slow_high*0.2, slow_high);
		fprintf(fp,"R2=0/360/%f/%f\n", slow_high*0.2, slow_high);
	}
	fprintf(fp,"J=Pa6i\n");
	fprintf(fp,"PS=%f-%f.ps\n", t1, t2); fprintf(fp,"PDF=%f-%f.pdf\n", t1, t2);
	fprintf(fp,"awk '{print $1,$2,$3/%f}' %s > tmp.file\n", cof_peak, argv[10]);
	fprintf(fp,"gmt surface tmp.file -R$R1 -I%f/%f -Gtmp.grd\n", baz_step/10., slow_step/5.);
	fprintf(fp,"gmt makecpt -Cjet -T%f/1/0.1 -Z >tmp.cpt\n", cof_low/cof_peak);
	fprintf(fp,"gmt psxy -R$R2 -J$J -K -T>$PS\n");
	fprintf(fp,"gmt grdimage tmp.grd -R -J -K -O -Ctmp.cpt -Bx30g15+l\"backazimuth(deg)\" -By%fg%f+l\"slowness(s/km)\" -BwsEN+t\"bandpass: %.3f ~ %.3f Hz\" >>$PS\n", slow_high/5., slow_high/10., fre_low, fre_high);
    fprintf(fp,"gmt psxy -R -J -K -O -W1p,white,\"..\" >>$PS<<EOF\n");
    for ( i = 0; i < 6; i ++ ) {
        fprintf(fp,">\n");
        fprintf(fp,"%d %f\n", i*30, slow_high);
        fprintf(fp,"%d %f\n", i*30+180, slow_high);
    }
    fprintf(fp,"EOF\n");
	while ( grid > (2.8*grid_step) ) {
		grid -= grid_step;
		fprintf(fp,"echo 0 %.3f %.3f | gmt pstext -R -J -K -O -F+f12p>>$PS\n", grid, grid);
	}
	fprintf(fp,"gmt psscale -Ctmp.cpt -D7i/3i/12/0.8 -Ba0.1g0:\"Normalized cross-coeffient\": -K -O >>$PS\n");
	fprintf(fp,"echo 0 %f N | gmt pstext -R -J -K -O -F+f15p,27,red -Y1.9i>>$PS\n", slow_high/2.);
	fprintf(fp,"echo 90 %f E | gmt pstext -R -J -K -O -F+f15p,27,red -X1.9i -Y-1.9i>>$PS\n", slow_high/2.);
	fprintf(fp,"echo 180 %f S | gmt pstext -R -J -K -O -F+f15p,27,red -Y-1.9i -X-1.9i>>$PS\n", slow_high/2.);
	fprintf(fp,"echo 270 %f W | gmt pstext -R -J -K -O -F+f15p,27,red -X-1.9i -Y1.9i>>$PS\n", slow_high/2.);
	fprintf(fp,"gmt psxy -R -J -O -T>>$PS\n");
	fprintf(fp,"ps2pdf $PS $PDF\n");
    fprintf(fp,"gmt psconvert %f-%f.ps -A -P -Tg\n", t1, t2);
	fprintf(fp,"rm gmt.* tmp.*\n");

/*-----------------------------------------------------------------array coordinates------------------------------------------------------------------------*/
	for ( i = 0; i < count; i ++ ) {
		printf("station %4d   lon -> %10.6f  lat -> %10.6f\n", i+1, coordi[i][0], coordi[i][1]);
	}
	    printf("array center:  lon -> %10.6f  lat -> %10.6f\n", center_lon, center_lat);
/*-------------------------------------------release dynamic memory of array "coordi", "amp", "tmp" and "sum"-----------------------------------------------*/
    fclose(fp); fclose(fout);
    for ( i = 0; i < count; i++ ) {
        free(coordi[i]); free(amp[i]); free(tmp[i]);
    }
    free(coordi); free(amp); free(tmp); free(sum);

    time_end = clock();
    time_used = (time_end - time_start)/CLOCKS_PER_SEC;
    printf("Time used: %f seconds!!!\n", time_used);
    return 0;
}
