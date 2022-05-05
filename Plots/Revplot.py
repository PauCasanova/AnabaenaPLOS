import funcRevplot
import matplotlib.pyplot as plt

plt.style.use('seaborn-dark-palette')
plt.rcParams['font.size'] = 20

OMutantes=['WT','patSdel','hetNdel']
NMutantes=['patAdel','patAhetNdel']

DifpatAdel=['patAdel','patAdel full leakage','patAdel no leakage']
DifpatAhetNdel=['patAhetNdel','patAhetNdel full leakage','patAhetNdel no leakage']

AltpatS=['WT','patXdel','patSdel','patXpatSdel','patXdel2','patSdel2']

Tetra=['WT','WT(tetramers&0.53KR)','patSdel','patSdel(tetramers&0.53KR)','hetNdel','hetNdel(tetramers&0.53KR)']

Leak=['patAdel no leakage','patAdel','patAhetNdel no leakage','patAhetNdel']

Cellnoise=['WT(Cell noise x0.1)','WT','WT(Cell noise x10)']
Genenoise=['WT(Gene noise x0.1)','WT','WT(Gene noise x4)']
Allnoise=['WT(Low noise)','WT','WT(High noise)']
noise=['WT','WT(Cell noise x0.1)','WT(Cell noise x10)','WT(Gene noise x0.1)','WT(Gene noise x4)','WT(Low noise)','WT(High noise)']

colors={
	'WT':'#228B22',
	'patXdel':'#1a684d',
	'patSdel':'#0000CD',
	'hetNdel':'#DC143C',
	'patAdel':'#ffeb63',
	'patAdelExp':'#b29b00',
	'patAdelhet':'#6F4EAE',
	'patAhetNdel':'#d18e4a',
	'patAhetNdelExp':'#68421a',
	'patAhetNdelhet':'#0B4758',
	'hetFdel':'#CC9900',
	'patApatSdel':'#999933',	
	'patXpatSdel':'#0A0A5C',
	'patAdel full leakage':'#9544A8',
	'patAdel no leakage':'#c90076',
	'patAhetNdel full leakage':'#3D4D8E',
	'patAhetNdel no leakage':'#379B55',
	'patXdel2':'#AE682A',
	'patSdel2':'#892160',
	'WT(Cell noise x0.1)':'#000000',
	'WT(Gene noise x0.1)':'#000000',
	'WT(Low noise)':'#000000',
	'WT(Cell noise x10)':'#ff8c00',
	'WT(Gene noise x4)':'#ff8c00',
	'WT(High noise)':'#ff8c00',
	'WT(tetramers&0.53KR)':'#004D00',
	'WT(tetramers&0.5KR)':'#68C468',
	'patSdel(tetramers&0.53KR)':'#00007E',
	'patSdel(tetramers&0.5KR)':'#2A2AFC',
	'hetNdel(tetramers&0.53KR)':'#88001B',
	'hetNdel(tetramers&0.5KR)':'#F16A84'	
}
# ~ funcRevplot.IntHistograms('L3',OMutantes,'OldP',colors)
# ~ funcRevplot.MeanPer('L3',OMutantes,'OldP',colors,0,16,2,0,25,4,'B')
# ~ funcRevplot.BorHistograms('L3',NMutantes,'NewP',colors,2)


# ~ funcRevplot.BorHistograms('L3',DifpatAdel,'DifpatA',colors,2)
# ~ funcRevplot.BorHistograms('L3',DifpatAhetNdel,'DifpatAhetN',colors,2)
# ~ funcRevplot.MeanPer('L3',AltpatS,'AltpatS',colors,2,18.5,2,3,28,4,'S')

# ~ Leak=['patAdel no leakage','patAhetNdel no leakage','patAdel','patAhetNdel']
# ~ funcRevplot.BorHistograms('L3',Leak,'Leak',colors,2)
# ~ Leak=['patAdel no leakage','patAdel','patAdel full leakage','patAhetNdel no leakage','patAhetNdel','patAhetNdel full leakage']
# ~ funcRevplot.BorHistograms('L3',Leak,'AllLeak',colors,2)

# ~ funcRevplot.MeanPer('L3',Tetra,'Tetra',colors,0,16,2,0,25,4,'B')

# ~ funcRevplot.MeanPer('L3',Cellnoise,'CellNoise',colors,7,17,2,2.9,13,4,'L')
# ~ funcRevplot.MeanPer('L3',Genenoise,'GeneNoise',colors,3,19.9,3,0,21,4,'S')
# ~ funcRevplot.MeanPer('L3',Allnoise,'AllNoise',colors,0,24,4,2.5,18,5,'L')
#noise subplot
# ~ funcRevplot.MeanPerSubplots('L3',noise,'Noise',colors,0,24,4,0,21,5,3,1)

funcRevplot.ClustHistograms('L3',OMutantes,'OldClusters',colors,3,25,3,5,24,3)


# ~ Bord=['patAdel','patAhetNdel']
# ~ funcRevplot.BorHistograms('L3',Bord,'Bord',colors,2)
