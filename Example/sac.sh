gmt psxy -R500/1900/65.5/69.5 -JX8i/6i -K -T > SAC.ps
#saclst gcarc f *.SAC | awk '{print $1, "0", $2}' | \
    gmt pssac -R -J -K -O -C500/1900 -M1 \
    -Bx100f10+l"Travel time (s)" \
    -By0.5f0.1+l"Distance (deg)" -BWSen -W1p -Edt *.SAC >> SAC.ps
gmt psxy -R -J -O -T >> SAC.ps
gmt psconvert SAC.ps -A -P -Tg

ps2pdf SAC.ps SAC.pdf
evince SAC.pdf
rm SAC.ps gmt.*
