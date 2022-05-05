#GIFS
#convert -delay 17.5 -loop 0 `ls -v` WT.gif

reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,10' fontscale 1 size 1024,768
set encoding koi8u
set encoding utf8
 
set border lw 4
set tics scale 2
set tit font ',20'
set tics font ',30'
set key font ',25'
set xlabel 'time (h)' font ',35'
set xtics offset 0,graph -0.02
set xlabel offset -1.5,graph -0.035
set ytics offset 0,graph 0
set ylabel offset -5,graph 0
set tmargin 1
set bmargin 6
set rmargin 1
set lmargin 12

set style fill transparent solid 0.3 noborder

Error=sqrt(1)
ExpError=sqrt(1)

unset key
#~ set key at 41,15.8
set xrange [20:80]
set yrange [0:16]
set xtics 24, 12, 78
set ytics 2, 2, 14
set output 'Mean.png'
set ylabel 'mean length of intervals' font ',35'
set zeroaxis linetype 0 linewidth 2.5
plot 'SimData.dat' u 1:'NWTMean' w l lw 6 lc rgb '#228B22' notit,\
	 'SimData.dat' u 1:(column('NWTMean')-column('NWTMeanDev')/Error):(column('NWTMean')+column('NWTMeanDev')/Error) with filledcurves lc rgb '#228B22' notit,\
	 'ExpData/ExpData.dat' u 1:'WTMean':(column('WTMeanDev')/ExpError) with yerrorbars lw 4 ps 2 pt 4 lc rgb '#228B22' tit 'wild type',\
	 'SimData.dat' u 1:'NpatSdelMean' w l lw 6 lc '#0000CD' notit,\
	 'SimData.dat' u 1:(column('NpatSdelMean')-column('NpatSdelMeanDev')/Error):(column('NpatSdelMean')+column('NpatSdelMeanDev')/Error) with filledcurves lc rgb '#0000CD' notit,\
	 'ExpData/ExpData.dat' u 1:'patSdelMean':(column('patSdelMeanDev')/ExpError) with yerrorbars lw 4 ps 2 pt 8 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}',\
	 'SimData.dat' u 1:'NhetNdelMean' w l lw 6 lc '#DC143C' notit,\
	 'SimData.dat' u 1:(column('NhetNdelMean')-column('NhetNdelMeanDev')/Error):(column('NhetNdelMean')+column('NhetNdelMeanDev')/Error) with filledcurves lc rgb '#DC143C' notit,\
	 'ExpData/ExpData.dat' u 1:'hetNdelMean':(column('hetNdelMeanDev')/ExpError) with yerrorbars lw 4 ps 2 pt 6 lc '#DC143C' tit 'Δ{/Euclid-Italic-Bold hetN}';
 
set key at 80,24
set xrange [20:80]
set yrange [0:25]
set xtics 24, 12, 78
set ytics 0, 4, 25
set output 'Percentage.png'
unset tit
set ylabel offset -5,graph -0.02
set xlabel offset -1.5,graph -0.04
set ylabel 'percentage of heterocysts' font ',35'
set zeroaxis linetype 0 linewidth 2.5
plot 'SimData.dat' u 1:'NWTHetPer' w l lw 6 lc rgb '#228B22' notit,\
	'SimData.dat' u 1:(column('NWTHetPer')-column('NWTHetPerDev')/Error):(column('NWTHetPer')+column('NWTHetPerDev')/Error) with filledcurves lc rgb '#228B22' notit,\
	'ExpData/ExpData.dat' u 1:'WTHetPer':(column('WTHetPerDev')/ExpError) with yerrorbars lw 4 ps 2 pt 4 lc rgb '#228B22' tit 'wild type',\
	'SimData.dat' u 1:'NpatSdelHetPer' w l lw 6 lc '#0000CD' notit,\
	'SimData.dat' u 1:(column('NpatSdelHetPer')-column('NpatSdelHetPerDev')/Error):(column('NpatSdelHetPer')+column('NpatSdelHetPerDev')/Error) with filledcurves lc rgb '#0000CD' notit,\
	'ExpData/ExpData.dat' u 1:'patSdelHetPer':(column('patSdelHetPerDev')/ExpError) with yerrorbars lw 4 ps 2 pt 8 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}',\
	'SimData.dat' u 1:'NhetNdelHetPer' w l lw 6 lc '#DC143C' notit,\
	'SimData.dat' u 1:(column('NhetNdelHetPer')-column('NhetNdelHetPerDev')/Error):(column('NhetNdelHetPer')+column('NhetNdelHetPerDev')/Error) with filledcurves lc rgb '#DC143C' notit,\
	'ExpData/ExpData.dat' u 1:'hetNdelHetPer':(column('hetNdelHetPerDev')/ExpError) with yerrorbars lw 4 ps 2 pt 6 lc '#DC143C' tit 'Δ{/Euclid-Italic-Bold hetN}';

set key default
set xrange [20:96]
set yrange [0:5]
set xtics 24, 24, 98
set ytics 0, 1, 5
set output 'PercentageNewMut.png'
unset tit
set ylabel 'percentage of internal heterocysts' font ',32'
set zeroaxis linetype 0 linewidth 2.5
plot 'SimData.dat' u 1:'NpatAdelHetPer':'NpatAdelHetPerDev' w l lw 6 lc rgb '#ffeb63' tit 'Δ{/Euclid-Italic-Bold patA}',\
     'SimData.dat' u 1:(column('NpatAdelHetPer')-column('NpatAdelHetPerDev')/Error):(column('NpatAdelHetPer')+column('NpatAdelHetPerDev')/Error) with filledcurves lc rgb '#ffeb63' notit,\
	 'SimData.dat' u 1:'NpatAhetNdelHetPer':'NpatAhetNdelHetPerDev' w l lw 6 lc '#d18e4a' tit 'ΔpatAΔ{/Euclid-Italic-Bold hetN}',\
	 'SimData.dat' u 1:(column('NpatAhetNdelHetPer')-column('NpatAhetNdelHetPerDev')/Error):(column('NpatAhetNdelHetPer')+column('NpatAhetNdelHetPerDev')/Error) with filledcurves lc rgb '#d18e4a' notit;

	
reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,20' fontscale 1 size 3072,768
unset xlabel
set boxwidth 1. absolute
set style fill solid 1.00 border lt -1
set key fixed right top vertical Right noreverse noenhanced nobox
set style increment default
set style histogram errorbars gap 1 tit textcolor lt -1 lw 1.5
set datafile missing ' '
set style data histograms
set xtics out nomirror autojustify
set xtics norangelimit
set xtics ()
set xrange [ -0.9 : 26.8 ] noreverse writeback
set tics scale 0.8
set tit font ',30'
set tics font ',40'
set key font ',35'
set border lw 3
set ylabel 
set tmargin 0
set bmargin 1
set xtics offset 0,graph -10

LMARGIN = "set lmargin at screen 0.08; set rmargin at screen 0.38"
MMARGIN = "set lmargin at screen 0.385; set rmargin at screen 0.685"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"

set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
set format y
set yrange [ 0 : 0.18 ] noreverse writeback
set ytics 0, 0.03, 0.16

set output 'HistoWT.png'
set multiplot layout 1,3 title "wild type" font ",50" offset 5, graph -0.25
#~ set output 'HistoWT24h.png'
set label "24h" font ",45" at 1,0.169
set key
@LMARGIN 
plot 'Histograms.dat' using 'WTExp24hMean':(column('WTExp24hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#114511' tit 'Experiment','' u 'NWTSim24hMean':(column('NWTSim24hMeanDev')/Error) lc '#90c590' tit 'Simulation';
												
#~ set output 'HistoWT48h.png' 
unset label
set label "48h" font ",45" at 1,0.169
@MMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'WTExp48hMean':(column('WTExp48hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#114511' tit 'Experiment','' u 'NWTSim48hMean':(column('NWTSim48hMeanDev')/Error) lc '#90c590' tit 'Simulation';
											
#~ set output 'HistoWT72h.png' 
unset label
set label "72h" font ",45" at 1,0.169
@RMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'WTExp72hMean':(column('WTExp72hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#114511' tit 'Experiment','' u 'NWTSim72hMean':(column('NWTSim72hMeanDev')/Error) lc '#90c590' tit 'Simulation';

unset multiplot

set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
set format y
set yrange [ 0 : 0.37 ] noreverse writeback
set ytics 0, 0.05, 0.35
set xlabel 'interval length' font ',40'

set output 'HistopatS.png'
#~ set multiplot layout 1,3 title "Δ{/Euclid-Italic-Bold patS} without border diffusion" font ",50" offset 5, graph -0.025
set multiplot layout 1,3 title "Δ{/Euclid-Italic-Bold patS}" font ",50" offset 5, graph -0.15
set xtics offset 0,graph 0.01
#~ set output 'HistopatS24h.png'
unset label
set label "24h" font ",45" at 1,0.348
set key
@LMARGIN 
plot 'Histograms.dat' using 'patSdelExp24hMean':(column('patSdelExp24hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Experiment','' u 'NpatSdelSim24hMean':(column('NpatSdelSim24hMeanDev')/Error) lc '#7f7fe6' tit 'Simulation';
											
#~ set output 'HistopatS48h.png' 
unset label
set label "48h" font ",45" at 1,0.348
@MMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'patSdelExp48hMean':(column('patSdelExp48hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Experiment','' u 'NpatSdelSim48hMean':(column('NpatSdelSim48hMeanDev')/Error) lc '#7f7fe6' tit 'Simulation';
											
#~ set output 'HistopatS72h.png' 
unset label
set label "72h" font ",45" at 1,0.348
@RMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'patSdelExp72hMean':(column('patSdelExp72hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Experiment','' u 'NpatSdelSim72hMean':(column('NpatSdelSim72hMeanDev')/Error) lc '#7f7fe6' tit 'Simulation';

unset multiplot

set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
set format y
set yrange [ 0 : 0.43 ] noreverse writeback
set ytics 0, 0.05, 0.4
set xtics offset 0,graph -10

set output 'HistohetN.png'
set multiplot layout 1,3 title "Δ{/Euclid-Italic-Bold hetN}" font ",50" offset 5, graph -0.15
#~ set output 'HistohetN24h.png'
unset label
set label "24h" font ",45" at 1,0.405
@LMARGIN 
set key
plot 'Histograms.dat' using 'hetNdelExp24hMean':(column('hetNdelExp24hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#6e0a1e' tit 'Experiment','' u 'NhetNdelSim24hMean':(column('NhetNdelSim24hMeanDev')/Error) lc '#ed899d' tit 'Simulation';
											
#~ set output 'HistohetN48h.png' 
unset label
set label "48h" font ",45" at 1,0.405
@MMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'hetNdelExp48hMean':(column('hetNdelExp48hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#6e0a1e' tit 'Experiment','' u 'NhetNdelSim48hMean':(column('NhetNdelSim48hMeanDev')/Error) lc '#ed899d' tit 'Simulation';
												
#~ set output 'HistohetN72h.png' 
unset label
set label "72h" font ",45" at 1,0.405
@RMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'hetNdelExp72hMean':(column('hetNdelExp72hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#6e0a1e' tit 'Experiment','' u 'NhetNdelSim72hMean':(column('NhetNdelSim72hMeanDev')/Error) lc '#ed899d' tit 'Simulation';

unset multiplot

#~ set tmargin 0.5
#~ set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
#~ set format y
#~ set yrange [ 0 : 0.35 ] noreverse writeback
#~ set xtics offset 0,graph 0.01
#~ set tics scale 0.8
#~ set ytics 0, 0.05, 0.34

#~ set output 'HistoMpatS.png'
#~ set multiplot layout 1,3
#~ unset label
#~ set label "24h" font ",45" at 1,0.328
#~ @LMARGIN
#~ set key
#~ plot 'Histograms.dat' using 'patSdelExp24hMean':(column('patSdelExp24hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Exp ΔpatS','' u 'NpatSdelSim24hMean':(column('NpatSdelSim24hMeanDev')/Error) lc '#7f7fe6' tit 'Sim ΔpatS', '' u 'NpatXpatSdelSim24hMean':(column('NpatXpatSdelSim24hMeanDev')/Error) fs pattern 2 lc '#3f3fda' tit 'Sim ΔpatXΔpatS';
								
#~ unset label
#~ set label "48h" font ",45" at 1,0.328
#~ @MMARGIN 
#~ unset key 
#~ set format y ''
#~ unset ylabel 
#~ plot 'Histograms.dat' using 'patSdelExp48hMean':(column('patSdelExp48hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Exp ΔpatS','' u 'NpatSdelSim48hMean':(column('NpatSdelSim48hMeanDev')/Error) lc '#7f7fe6' tit 'Sim ΔpatS', '' u 'NpatXpatSdelSim48hMean':(column('NpatXpatSdelSim48hMeanDev')/Error) fs pattern 2 lc '#3f3fda' tit 'Sim ΔpatXΔpatS';
								
#~ unset label
#~ set label "72h" font ",45" at 1,0.328
#~ @RMARGIN 
#~ unset key 
#~ set format y ''
#~ unset ylabel 
#~ plot 'Histograms.dat' using 'patSdelExp72hMean':(column('patSdelExp72hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Exp ΔpatS','' u 'NpatSdelSim72hMean':(column('NpatSdelSim72hMeanDev')/Error) lc '#7f7fe6' tit 'Sim ΔpatS', '' u 'NpatXpatSdelSim72hMean':(column('NpatXpatSdelSim72hMeanDev')/Error) fs pattern 2 lc '#3f3fda' tit 'Sim ΔpatXΔpatS';

#~ unset multiplot

#~ set key right top vertical Right noreverse noenhanced nobox
#~ set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
#~ set format y
#~ set yrange [ 0 : 0.23 ] noreverse writeback
#~ set ytics 0, 0.03, 0.21

#~ set output 'HistopatApatS.png'
#~ set multiplot layout 1,3
#~ unset label
#~ set label "24h" font ",45" at 22,0.212
#~ @LMARGIN 
#~ set key at 26,0.19
#~ plot 'Histograms.dat' using 'NpatSdelSim24hMean':(column('NpatSdelSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'ΔpatS', '' u 'NpatApatSdelSim24hMean':(column('NpatApatSdelSim24hMeanDev')/Error) lc '#999933' tit 'ΔpatAΔpatS','' u 'NWTSim24hMean':(column('NWTSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type';
									
#~ unset label
#~ set label "48h" font ",45" at 22,0.212
#~ @MMARGIN 
#~ unset key 
#~ set format y ''
#~ unset ylabel 
#~ plot 'Histograms.dat' using 'NpatSdelSim48hMean':(column('NpatSdelSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'ΔpatS', '' u 'NpatApatSdelSim48hMean':(column('NpatApatSdelSim48hMeanDev')/Error) lc '#999933' tit 'ΔpatAΔpatS','' u 'NWTSim48hMean':(column('NWTSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type';
									
#~ unset label
#~ set label "72h" font ",45" at 22,0.212
#~ @RMARGIN
#~ unset key 
#~ set format y ''
#~ unset ylabel 
#~ plot 'Histograms.dat' using 'NpatSdelSim72hMean':(column('NpatSdelSim72hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'ΔpatS', '' u 'NpatApatSdelSim72hMean':(column('NpatApatSdelSim72hMeanDev')/Error) lc '#999933' tit 'ΔpatAΔpatS','' u 'NWTSim72hMean':(column('NWTSim72hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type';

#~ unset multiplot

#~ set key fixed right top vertical Right noreverse noenhanced nobox
#~ set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
#~ set format y
#~ set xrange [ -0.8 : 26.8 ] noreverse writeback
#~ set yrange [ 0 : 0.18 ] noreverse writeback
#~ set ytics 0, 0.03, 0.16

#~ set output 'ExpHistopatX.png'
#~ set multiplot layout 1,3
#~ unset label
#~ set label "24h" font ",45" at 1,0.165
#~ @LMARGIN 
#~ set key
#~ plot 'Histograms.dat' using 'WTExp24hMean':(column('WTExp24hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#114511' tit 'Exp wild type','' u 'NWTSim24hMean':(column('NWTSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'Sim wild type','' u 'NpatXdelSim24hMean':(column('NpatXdelSim24hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'Sim ΔpatX';
								
#~ unset label
#~ set label "48h" font ",45" at 1,0.165
#~ @MMARGIN
#~ unset key 
#~ set format y ''
#~ unset ylabel
#~ plot 'Histograms.dat' using 'WTExp48hMean':(column('WTExp48hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#114511' tit 'Exp wild type','' u 'NWTSim48hMean':(column('NWTSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'Sim wild type','' u 'NpatXdelSim48hMean':(column('NpatXdelSim48hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'Sim ΔpatX';
								
#~ unset label
#~ set label "72h" font ",45" at 1,0.165
#~ @RMARGIN
#~ unset key 
#~ set format y ''
#~ unset ylabel
#~ plot 'Histograms.dat' using 'WTExp72hMean':(column('WTExp72hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#114511' tit 'Exp wild type','' u 'NWTSim72hMean':(column('NWTSim72hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'Sim wild type','' u 'NpatXdelSim72hMean':(column('NpatXdelSim72hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'Sim ΔpatX';

#~ unset multiplot

#~ set key fixed right top vertical Right noreverse noenhanced nobox
#~ set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
#~ set format y
#~ set xrange [ -0.8 : 26.8 ] noreverse writeback
#~ set yrange [ 0 : 0.18 ] noreverse writeback
#~ set ytics 0, 0.03, 0.16

#~ set output 'HistopatX.png'
#~ set multiplot layout 1,3
#~ unset label
#~ set label "24h" font ",45" at 1,0.165
#~ @LMARGIN 
#~ set key
#~ plot 'Histograms.dat' using 'NWTSim24hMean':(column('NWTSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type','' u 'NpatXdelSim24hMean':(column('NpatXdelSim24hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'ΔpatX';
								
#~ unset label
#~ set label "48h" font ",45" at 1,0.165
#~ @MMARGIN
#~ unset key 
#~ set format y ''
#~ unset ylabel
#~ plot 'Histograms.dat' using 'NWTSim48hMean':(column('NWTSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type','' u 'NpatXdelSim48hMean':(column('NpatXdelSim48hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'ΔpatX';
							
#~ unset label
#~ set label "72h" font ",45" at 1,0.165
#~ @RMARGIN
#~ unset key 
#~ set format y ''
#~ unset ylabel
#~ plot 'Histograms.dat' using 'NWTSim72hMean':(column('NWTSim72hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type','' u 'NpatXdelSim72hMean':(column('NpatXdelSim72hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'ΔpatX';

#~ unset multiplot


reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,20' fontscale 1 size 3072,768
unset xlabel
set boxwidth 1. absolute
set style fill solid 1.00 border lt -1
set key fixed right top vertical Right noreverse noenhanced nobox
set style increment default
set style histogram errorbars gap 1 tit textcolor lt -1 lw 1.5 
set datafile missing ' '
set style data histograms
set xtics out nomirror autojustify
set xtics norangelimit
set xtics ()
set xrange [ -0.5 : 4.5 ] noreverse writeback
unset ylabel
set yrange [ * : * ] noreverse 
set tit font ',30'
set tics font ',40'
set key font ',30'
set border lw 3
set ylabel 
set tmargin 0.5
set bmargin 1
set xtics offset 0,graph 0.01
set xlabel 'terminal heterocysts' font ',44'

LMARGIN = "set lmargin at screen 0.0624; set rmargin at screen 0.2924"
MLMARGIN = "set lmargin at screen 0.2974; set rmargin at screen 0.5274"
MRMARGIN = "set lmargin at screen 0.5324; set rmargin at screen 0.7624"
RMARGIN = "set lmargin at screen 0.7674; set rmargin at screen 0.99"

set yrange [ 0 : 1 ] noreverse writeback
set ytics 0, 0.2, 1

set ylabel 'fraction of filaments' font ',48' offset -2.5 ,graph -0.01
set output 'HistopatA.png'
set multiplot layout 1,4

#~ set output 'HistopatAdel24h.png'
set label "24h" font ",45" at 3,0.9
@LMARGIN 
set key at 4.5,0.8
plot 'BorderHistograms.dat' using 'patAdelExp24hMean':(column('patAdelExp24hMeanDev')/ExpError):xtic(1) lc '#b29b00' tit 'Experiment','' u 'NpatAdelSim24hMean':(column('NpatAdelSim24hMeanDev')/Error) lc '#ffeb63' tit 'Simulation';
																																			
#~ set output 'HistopatAdel48h.png' 
unset label
set label "48h" font ",45" at 3,0.9
@MLMARGIN 
unset key 
set format y ''
unset ylabel
plot 'BorderHistograms.dat' using 'patAdelExp48hMean':(column('patAdelExp48hMeanDev')/ExpError):xtic(1) lc '#b29b00' tit 'Experiment','' u 'NpatAdelSim48hMean':(column('NpatAdelSim48hMeanDev')/Error) lc '#ffeb63' tit 'Simulation';
																																			
#~ set output 'HistopatAdel72h.png' 
unset label
set label "72h" font ",45" at 3,0.9
@MRMARGIN 
unset key 
set format y ''
unset ylabel
plot 'BorderHistograms.dat' using 'patAdelExp72hMean':(column('patAdelExp72hMeanDev')/ExpError):xtic(1) lc '#b29b00' tit 'Experiment','' u 'NpatAdelSim72hMean':(column('NpatAdelSim72hMeanDev')/Error) lc '#ffeb63' tit 'Simulation';
																																					
#~ set output 'HistopatAdel96h.png' 
unset label
set label "96h" font ",45" at 3,0.9
@RMARGIN 
unset key 
set format y ''
unset ylabel
plot 'BorderHistograms.dat' using 'patAdelExp96hMean':(column('patAdelExp96hMeanDev')/ExpError):xtic(1) lc '#b29b00' tit 'Experiment','' u 'NpatAdelSim96hMean':(column('NpatAdelSim96hMeanDev')/Error) lc '#ffeb63' tit 'Simulation';
 
unset multiplot 

set format y
set yrange [ 0 : 1 ] noreverse writeback 
set ytics 0, 0.2, 1

set ylabel 'fraction of filaments' font ',48' offset -2.5 ,graph -0.01
set output 'HistopatAhetN.png'
set multiplot layout 1,4

#~ set output 'HistopatAhetNdel24h.png'
unset label
set label "24h" font ",45" at 3,0.9
@LMARGIN 
set key at 4.5,0.8
plot 'BorderHistograms.dat' using 'patAhetNdelExp24hMean':(column('patAhetNdelExp24hMeanDev')/ExpError):xtic(1) lc '#68421a' tit 'Experiment','' u 'NpatAhetNdelSim24hMean':(column('NpatAhetNdelSim24hMeanDev')/Error) lc '#d18e4a' tit 'Simulation';
																																			
#~ set output 'HistopatAhetNdel48h.png' 
unset label
set label "48h" font ",45" at 3,0.9
@MLMARGIN 
unset key 
set format y ''
unset ylabel
plot 'BorderHistograms.dat' using 'patAhetNdelExp48hMean':(column('patAhetNdelExp48hMeanDev')/ExpError):xtic(1) lc '#68421a' tit 'Experiment','' u 'NpatAhetNdelSim48hMean':(column('NpatAhetNdelSim48hMeanDev')/Error) lc '#d18e4a' tit 'Simulation';
																																			
#~ set output 'HistopatAhetNdel72h.png' 
unset label
set label "72h" font ",45" at 3,0.9
@MRMARGIN 
unset key 
set format y ''
unset ylabel
plot 'BorderHistograms.dat' using 'patAhetNdelExp72hMean':(column('patAhetNdelExp72hMeanDev')/ExpError):xtic(1) lc '#68421a' tit 'Experiment','' u 'NpatAhetNdelSim72hMean':(column('NpatAhetNdelSim72hMeanDev')/Error) lc '#d18e4a' tit 'Simulation';
																																					
#~ set output 'HistopatAhetNdel96h.png' 
unset label
set label "96h" font ",45" at 3,0.9
@RMARGIN 
unset key 
set format y ''
unset ylabel
plot 'BorderHistograms.dat' using 'patAhetNdelExp96hMean':(column('patAhetNdelExp96hMeanDev')/ExpError):xtic(1) lc '#68421a' tit 'Experiment','' u 'NpatAhetNdelSim96hMean':(column('NpatAhetNdelSim96hMeanDev')/Error) lc '#d18e4a' tit 'Simulation';
 
unset multiplot


reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,10' fontscale 1 size 1024,768
set encoding koi8u
set encoding utf8
 
set border lw 4
set tics scale 2
set tit font ',20'
set tics font ',30'
set key font ',25'
set xlabel 'time (h)' font ',35'
set xtics offset 0,graph -0.02
set xlabel offset -1.5,graph -0.035
set ytics offset 0,graph 0
set ylabel offset -5,graph 0
set tmargin 1
set bmargin 6
set rmargin 1
set lmargin 12

set style fill transparent solid 0.3 noborder

unset key
#~ set key at 70,16
set xrange [20:80]
set yrange [2:16]
set xtics 24, 12, 78
set ytics 4, 2, 14
set output 'Mean(patS).png'
unset tit
set ylabel 'mean length of intervals' font ',35'
set zeroaxis linetype 0 linewidth 2.5
plot 'SimData.dat' u 1:'NWTMean' w l lw 6 lc rgb '#228B22' notit,\
	 'SimData.dat' u 1:(column('NWTMean')-column('NWTMeanDev')/Error):(column('NWTMean')+column('NWTMeanDev')/Error) with filledcurves lc rgb '#228B22' notit,\
	 'SimData.dat' u 1:'NpatXdelMean' w l dt "-" lt 1 lw 6 lc '#228B22' tit 'Δ{/Euclid-Italic-Bold patX}',\
	 'SimData.dat' u 1:(column('NpatXdelMean')-column('NpatXdelMeanDev')/Error):(column('NpatXdelMean')+column('NpatXdelMeanDev')/Error) with filledcurves lc rgb '#228B22' notit,\
	 'SimData.dat' u 1:'NpatSdelMean' w l lw 6 lc '#0000CD' notit,\
	 'SimData.dat' u 1:(column('NpatSdelMean')-column('NpatSdelMeanDev')/Error):(column('NpatSdelMean')+column('NpatSdelMeanDev')/Error) with filledcurves lc rgb '#0000CD' notit,\
	 'SimData.dat' u 1:'NpatXpatSdelMean' w l dt "-" lt 1 lw 6 lc '#3f3fda' tit 'Δ{/Euclid-Italic-Bold patS}Δ{/Euclid-Italic-Bold patX}',\
	 'SimData.dat' u 1:(column('NpatXpatSdelMean')-column('NpatXpatSdelMeanDev')/Error):(column('NpatXpatSdelMean')+column('NpatXpatSdelMeanDev')/Error) with filledcurves fs transparent solid 0.1 lc rgb '#3f3fda' notit,\
	 'SimData.dat' u 1:'NpatApatSdelMean' w l lw 6 lc '#999933' tit 'ΔpatAΔ{/Euclid-Italic-Bold patS}',\
	 'SimData.dat' u 1:(column('NpatApatSdelMean')-column('NpatApatSdelMeanDev')/Error):(column('NpatApatSdelMean')+column('NpatApatSdelMeanDev')/Error) with filledcurves lc rgb '#999933' notit,\
	 'ExpData/ExpData.dat' u 1:'WTMean':(column('WTMeanDev')/ExpError) with yerrorbars lw 4 ps 2 pt 4 lc rgb '#228B22' tit 'wild type',\
	 'ExpData/ExpData.dat' u 1:'patSdelMean':(column('patSdelMeanDev')/ExpError) with yerrorbars lw 4 ps 2 pt 8 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}';

set key top vertical maxrows 3 
set key at 80,27.3
set xrange [20:80]
set yrange [6:28]
set xtics 24, 12, 78
set ytics 10, 5, 25
set output 'Percentage(patS).png'
unset tit
set ylabel 'percentage of heterocysts' font ',35'
set zeroaxis linetype 0 linewidth 2.5
plot 'ExpData/ExpData.dat' u 1:'WTHetPer':(column('WTHetPerDev')/ExpError) with yerrorbars lw 4 ps 2 pt 4 lc rgb '#228B22' tit 'wild type',\
     'ExpData/ExpData.dat' u 1:'patSdelHetPer':(column('patSdelHetPerDev')/ExpError) with yerrorbars lw 4 ps 2 pt 8 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}',\
     "+" u 1:(NaN) title " " w dots linecolor rgb "white",\
     'SimData.dat' u 1:'NWTHetPer' w l lw 6 lc rgb '#228B22' notit,\
     'SimData.dat' u 1:(column('NWTHetPer')-column('NWTHetPerDev')/Error):(column('NWTHetPer')+column('NWTHetPerDev')/Error) with filledcurves lc rgb '#228B22' notit,\
     'SimData.dat' u 1:'NpatXdelHetPer' w l dt "-" lt 1 lw 6 lc '#228B22' tit 'Δ{/Euclid-Italic-Bold patX}',\
     'SimData.dat' u 1:(column('NpatXdelHetPer')-column('NpatXdelHetPerDev')/Error):(column('NpatXdelHetPer')+column('NpatXdelHetPerDev')/Error) with filledcurves lc rgb '#228B22' notit,\
     'SimData.dat' u 1:'NpatSdelHetPer' w l lw 6 lc '#0000CD' notit,\
     'SimData.dat' u 1:(column('NpatSdelHetPer')-column('NpatSdelHetPerDev')/Error):(column('NpatSdelHetPer')+column('NpatSdelHetPerDev')/Error) with filledcurves lc rgb '#0000CD' notit,\
     'SimData.dat' u 1:'NpatXpatSdelHetPer' w l dt "-" lt 1 lw 6 lc '#3f3fda' tit 'Δ{/Euclid-Italic-Bold patX}Δ{/Euclid-Italic-Bold patS}',\
     'SimData.dat' u 1:(column('NpatXpatSdelHetPer')-column('NpatXpatSdelHetPerDev')/Error):(column('NpatXpatSdelHetPer')+column('NpatXpatSdelHetPerDev')/Error) with filledcurves fs transparent solid 0.1 lc rgb '#3f3fda' notit,\
     'SimData.dat' u 1:'NpatApatSdelHetPer' w l lw 6 lc '#999933' tit 'Δ{/Euclid-Italic-Bold patA}Δ{/Euclid-Italic-Bold patS}',\
     'SimData.dat' u 1:(column('NpatApatSdelHetPer')-column('NpatApatSdelHetPerDev')/Error):(column('NpatApatSdelHetPer')+column('NpatApatSdelHetPerDev')/Error) with filledcurves lc rgb '#999933' notit;


reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,10' fontscale 1 size 1024,768
set encoding koi8u
set encoding utf8

set border lw 4
set tics scale 2
set tit font ',20'
set tics font ',30'
set key font ',25'
set xlabel 'time (h)' font ',35'
set xtics offset -1,graph -0.01
set xlabel offset -1.5,graph -0.03
set ytics offset 0,graph 0
set ylabel offset -5,graph 0
set tmargin 1
set bmargin 5.5
set rmargin 3
set lmargin 13

set style fill transparent solid 0.2 noborder

#~ set key below vertical maxrows 3 spacing 1 #~ width -1
set key vertical maxrows 4 spacing 1 
set key at 95,3.2
set xrange [0:96]
set xtics (0,6,24,48,72,96)
set yrange [0.25:3.25]
set ytics 0.5, 0.5, 3
set output 'Concentrations.png'
unset tit
set ylabel "[HetR] relative to WT" font ',30'
set zeroaxis linetype 0 linewidth 2
plot 'SimData.dat' u 1:(column('NWTHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 1 ps 2 lc '#228B22' tit 'wild type',\
	 'SimData.dat' u 1:(column('NpatSdelHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 10 ps 1.5 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}',\
	 'SimData.dat' u 1:(column('NhetNdelHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 6 ps 1.5 lc '#DC143C' tit 'Δ{/Euclid-Italic-Bold hetN}',\
	 'SimData.dat' u 1:(column('NpatAdelHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 4 ps 1.5 lc '#ffeb63' tit 'Δ{/Euclid-Italic-Bold patA}',\
	 'SimData.dat' u 1:(column('NpatAhetNdelHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 8 ps 1.5 lc '#d18e4a' tit 'Δ{/Euclid-Italic-Bold patA}Δ{/Euclid-Italic-Bold hetN}',\
	 'SimData.dat' u 1:(column('NpatApatSdelHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 12 ps 1.5 lc '#999933' tit 'Δ{/Euclid-Italic-Bold patA}Δ{/Euclid-Italic-Bold patS}',\
	 'SimData.dat' u 1:(column('NhetFdelHtRCon')/column('NWTHtRCon')) w linesp lt 1 lw 2 pt 2 ps 1.5 lc '#CC9900' tit 'Δ{/Euclid-Italic-Bold hetF}';

	
reset
set terminal pngcairo transparent truecolor enhanced font 'Euclid-Bold,10' fontscale 1 size 1024,768
set encoding koi8u
set encoding utf8

set border lw 4
set tics nomirror scale 2
set tit font ',20'
set tics font ',30'
set key font ',20'
set xlabel 'time (h)' font ',35'
set xtics offset -1,graph -0.01
set xlabel offset -1.5,graph -0.03
set ytics offset 0,graph 0
set ylabel offset -5.5,graph 0
set tmargin 1
set bmargin 5.5
set rmargin 3
set lmargin 13.5

set style fill transparent 
stats "SimData.dat" u 1:'NWTHtRCon'
set xrange [0:96]
set xtics (0,6,24,48,72,96)
set yrange [0:3.15]
set ytics 0, 0.5, 3
set output 'ConcentrationN.png'
unset tit
set key vertical maxrows 3 spacing 0.8 width -1
set key at 96,3.1
set ylabel 'normalized [HetR]' font ',33'
set zeroaxis linetype 0 linewidth 2
set multiplot
plot "+" u 1:(NaN) title " " w dots linecolor rgb "white",\
	 'SimData.dat' u 1:(column('NhetFdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 2 ps 1.5 lc '#CC9900' tit 'Δ{/Euclid-Italic-Bold hetF}',\
	 'SimData.dat' u 1:(column('NWTHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 1 ps 2 lc rgb '#228B22' tit 'wild type',\
	 'SimData.dat' u 1:(column('NpatAdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 4 ps 1.5 lc '#ffeb63' tit 'Δ{/Euclid-Italic-Bold patA}',\
	 'SimData.dat' u 1:(column('NhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 6 ps 1.5 lc '#DC143C' tit 'Δ{/Euclid-Italic-Bold hetN}',\
	 'SimData.dat' u 1:(column('NpatAhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 8 ps 1.5 lc '#d18e4a' tit 'Δ{/Euclid-Italic-Bold patAΔ{/Euclid-Italic-Bold patN}',\
	 'SimData.dat' u 1:(column('NpatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 10 ps 1.5 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}',\
	 'SimData.dat' u 1:(column('NpatApatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 12 ps 1.5 lc '#999933' tit 'Δ{/Euclid-Italic-Bold patA}Δ{/Euclid-Italic-Bold patS}',\
	 "+" u 1:(NaN) title " " w dots linecolor rgb "white";


set origin 0.38,0.38
set size 0.53,0.45
set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0
clear
set style fill transparent 
unset tit
unset key
unset xlabel
unset ylabel
set xrange [0:48]
set yrange [0:1.05]
set xtics offset -1,graph -0.02
set ytics offset 0,graph 0
set tics nomirror scale 2
set xtics (0,6,24,48,72,96)
set ytics 0, 0.25, 1
set tics font ',30'
plot 'SimData.dat' u 1:(column('NpatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 10 ps 1.5 lc '#0000CD' notit,\
	 'SimData.dat' u 1:(column('NhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 6 ps 1.5 lc '#DC143C' notit,\
	 'SimData.dat' u 1:(column('NWTHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 1 ps 2 lc rgb '#228B22' notit ,\
	 'SimData.dat' u 1:(column('NpatAdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 4 ps 1.5 lc '#ffeb63' notit,\
	 'SimData.dat' u 1:(column('NpatAhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 8 ps 1.5 lc '#d18e4a' notit,\
	 'SimData.dat' u 1:(column('NpatApatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 12 ps 1.5 lc '#999933' notit,\
	 'SimData.dat' u 1:(column('NhetFdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 2 ps 1.5 lc '#CC9900' notit;
	
unset multiplot

reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,10' fontscale 1 size 1024,768
set encoding koi8u
set encoding utf8

set key vertical maxrows 4 spacing 1 
set key at 95,3
set border lw 4
set tics nomirror scale 2
set tit font ',20'
set tics font ',30'
set key font ',25'
set xlabel 'time (h)' font ',35'
set xrange [0:96]
set xtics (0,6,24,48,72,96)
set xtics offset -1,graph -0.02
set xlabel offset -1.5,graph -0.15
set ytics offset 0,graph 0
set ylabel 'normalized [HetR]' font ',33'
set ylabel offset -8 ,graph 0.55
set output 'BConcentrationN.png'
bm = 0.11
lm = 0.15
rm = 0.97
gap = 0.03
size = 0.85
y1 = 0
y2 = 1.05
y3 = 2.4
y4 = 3.15

set multiplot
set border 1+2+8
set xtics nomirror
set ytics nomirror
set lmargin at screen lm
set rmargin at screen rm
set bmargin at screen bm
set tmargin at screen bm + size * (abs(y2-y1) / (abs(y2-y1) + abs(y4-y3) ) )

set yrange [y1:y2]
set ytics 0, 0.25, 1

plot 'SimData.dat' u 1:(column('NpatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 11 ps 1.5 lc '#0000CD' notit,\
	 'SimData.dat' u 1:(column('NhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 7 ps 1.5 lc '#DC143C' notit,\
	 'SimData.dat' u 1:(column('NWTHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 1 ps 2 lc '#228B22' notit ,\
	 'SimData.dat' u 1:(column('NpatAdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 5 ps 1.5 lc '#ffeb63' notit,\
	 'SimData.dat' u 1:(column('NpatAhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 9 ps 1.5 lc '#d18e4a' notit,\
	 'SimData.dat' u 1:(column('NpatApatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 13 ps 1.5 lc '#999933' notit,\
	 'SimData.dat' u 1:(column('NhetFdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 2 ps 1.5 lc '#CC9900' notit;

unset xtics
unset xlabel
set border 2+4+8
set bmargin at screen bm + size * (abs(y2-y1) / (abs(y2-y1) + abs(y4-y3) ) ) + gap
set tmargin at screen bm + size + gap
set yrange [y3:y4]
set ytics 2, 0.25, 3
unset ylabel

set arrow from screen lm - gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1)+abs(y4-y3) ) ) - gap / 4.0 to screen \
lm + gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1) + abs(y4-y3) ) ) + gap / 4.0 nohead lw 4

set arrow from screen lm - gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1)+abs(y4-y3) ) ) - gap / 4.0  + gap to screen \
lm + gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1) + abs(y4-y3) ) ) + gap / 4.0 + gap nohead lw 4

set arrow from screen rm - gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1)+abs(y4-y3) ) ) - gap / 4.0 to screen \
rm + gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1) + abs(y4-y3) ) ) + gap / 4.0 nohead lw 4

set arrow from screen rm - gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1)+abs(y4-y3) ) ) - gap / 4.0  + gap to screen \
rm + gap / 4.0, bm + size * (abs(y2-y1) / (abs(y2-y1) + abs(y4-y3) ) ) + gap / 4.0 + gap nohead lw 4

plot 'SimData.dat' u 1:(column('NWTHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 1 ps 2 lc '#228B22' tit 'wild type',\
	 'SimData.dat' u 1:(column('NhetFdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 2 ps 1.5 lc '#CC9900' tit 'Δ{/Euclid-Italic-Bold hetF}',\
	 'SimData.dat' u 1:(column('NpatAdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 5 ps 1.5 lc '#ffeb63' tit 'Δ{/Euclid-Italic-Bold patA}',\
	 'SimData.dat' u 1:(column('NhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 7 ps 1.5 lc '#DC143C' tit 'Δ{/Euclid-Italic-Bold hetN}',\
	 'SimData.dat' u 1:(column('NpatAhetNdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 9 ps 1.5 lc '#d18e4a' tit 'Δ{/Euclid-Italic-Bold patAΔ{/Euclid-Italic-Bold patN}',\
	 'SimData.dat' u 1:(column('NpatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 11 ps 1.5 lc '#0000CD' tit 'Δ{/Euclid-Italic-Bold patS}',\
	 'SimData.dat' u 1:(column('NpatApatSdelHtRCon')/STATS_max_y) w linesp lt 1 lw 2 pt 13 ps 1.5 lc '#999933' tit 'Δ{/Euclid-Italic-Bold patA}Δ{/Euclid-Italic-Bold patS}';

unset multiplot

	
reset
set terminal pngcairo transparent enhanced font 'Euclid-Bold,20' fontscale 1 size 3072,768
unset xlabel
set boxwidth 1. absolute
set style fill solid 1.00 border lt -1
set key fixed right top vertical Right noreverse noenhanced nobox
set style increment default
set style histogram errorbars gap 1 tit textcolor lt -1 lw 1.5
set datafile missing ' '
set style data histograms
set xtics out nomirror autojustify
set xtics norangelimit
set xtics ()
set xrange [ -0.9 : 26.8 ] noreverse writeback
set tics scale 0.8
set tit font ',30'
set tics font ',40'
set key font ',35'
set border lw 3
set ylabel 
set tmargin 0
set bmargin 1
set xtics offset 0,graph -10
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.55"
RMARGIN = "set lmargin at screen 0.56; set rmargin at screen 0.96"

set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
set format y
set yrange [ 0 : 0.35 ] noreverse writeback
set tmargin 0.5
set bmargin 4
set xtics offset 0,graph 0.01
set tics scale 0.8
set ytics 0, 0.05, 0.34

set output 'HistoMpatS.png'
set multiplot layout 1,3
unset label
set label "24h" font ",45" at 1,0.328
@LMARGIN
set key
plot 'Histograms.dat' using 'patSdelExp24hMean':(column('patSdelExp24hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Exp ΔpatS','' u 'NpatSdelSim24hMean':(column('NpatSdelSim24hMeanDev')/Error) lc '#7f7fe6' tit 'Sim ΔpatS', '' u 'NpatXpatSdelSim24hMean':(column('NpatXpatSdelSim24hMeanDev')/Error) fs pattern 2 lc '#3f3fda' tit 'Sim ΔpatXΔpatS';
									
unset label
set label "48h" font ",45" at 1,0.328
@RMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'patSdelExp48hMean':(column('patSdelExp48hMeanDev')/ExpError):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'Exp ΔpatS','' u 'NpatSdelSim48hMean':(column('NpatSdelSim48hMeanDev')/Error) lc '#7f7fe6' tit 'Sim ΔpatS', '' u 'NpatXpatSdelSim48hMean':(column('NpatXpatSdelSim48hMeanDev')/Error) fs pattern 2 lc '#3f3fda' tit 'Sim ΔpatXΔpatS';

unset multiplot

set key right top vertical Right noreverse noenhanced nobox
set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
set format y
set yrange [ 0 : 0.23 ] noreverse writeback
set ytics 0, 0.03, 0.21

set output 'HistopatApatS.png'
set multiplot layout 1,3
unset label
set label "24h" font ",45" at 22,0.212
@LMARGIN 
set key at 26,0.19
plot 'Histograms.dat' using 'NpatSdelSim24hMean':(column('NpatSdelSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'ΔpatS', '' u 'NpatApatSdelSim24hMean':(column('NpatApatSdelSim24hMeanDev')/Error) lc '#999933' tit 'ΔpatAΔpatS','' u 'NWTSim24hMean':(column('NWTSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type';
									
unset label
set label "48h" font ",45" at 22,0.212
@RMARGIN 
unset key 
set format y ''
unset ylabel 
plot 'Histograms.dat' using 'NpatSdelSim48hMean':(column('NpatSdelSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#000066' tit 'ΔpatS', '' u 'NpatApatSdelSim48hMean':(column('NpatApatSdelSim48hMeanDev')/Error) lc '#999933' tit 'ΔpatAΔpatS','' u 'NWTSim48hMean':(column('NWTSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type';
		
unset multiplot

set key fixed right top vertical Right noreverse noenhanced nobox
set ylabel 'fraction of intervals' font ',48' offset -3.5 ,graph 0.01
set format y
set xrange [ -0.8 : 26.8 ] noreverse writeback
set yrange [ 0 : 0.18 ] noreverse writeback
set ytics 0, 0.03, 0.16

set output 'HistopatX.png'
set multiplot layout 1,3
unset label
set label "24h" font ",45" at 1,0.165
@LMARGIN 
set key
plot 'Histograms.dat' using 'NWTSim24hMean':(column('NWTSim24hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type','' u 'NpatXdelSim24hMean':(column('NpatXdelSim24hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'ΔpatX';
								
unset label
set label "48h" font ",45" at 1,0.165
@RMARGIN
unset key 
set format y ''
unset ylabel
plot 'Histograms.dat' using 'NWTSim48hMean':(column('NWTSim48hMeanDev')/Error):xtic(int($0)%3 == 0 ? stringcolumn(1) : '') lc '#228B22' tit 'wild type','' u 'NpatXdelSim48hMean':(column('NpatXdelSim48hMeanDev')/Error) fs pattern 2 lc '#228B22' tit 'ΔpatX';

unset multiplot
