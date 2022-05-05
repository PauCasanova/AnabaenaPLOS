import os
import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import re
from itertools import combinations
import matplotlib.transforms

expfile=os.path.dirname(os.getcwd()) + "/Plots/ExpData/exphisto.dat"
with open(expfile, newline = '') as _file_2:
	expdata = pd.read_csv(_file_2, sep='\t', engine='python', index_col=0)

	# ~ #print(expdata)	
for col in expdata.columns:
	Total=sum(expdata[col])
	if Total!=1:
		expdata[col]=[a/Total for a in expdata[col]]
	# ~ #print(expdata)
expdata.to_csv(expfile, sep='\t', float_format='%.4f', index=True)

sim='L3'
OMutantes=['WT','patSdel','hetNdel']
DataSets=[]
simdata=pd.DataFrame(index=expdata.index)
for mut in OMutantes:
	if mut=='WT':
		Time=[24,30,48,50,72,96]
	else:
		Time=[24,48,72,96]
	for t in Time:
		simfile=os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/histogram_med_t={}.dat".format(sim,mut,t)
		with open(simfile, newline = '') as _file_2:
			hdata = pd.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
		MValues = hdata.iloc[:,0]
		name='{}{}hSimulation'.format(mut,t)
		simdata[name]=MValues
		# ~ #print(simdata)	
data=pd.concat([simdata,expdata],axis=1)	

for col in data.columns:
    data.rename(columns={col:col.replace("patSdel","ΔpatS")},inplace=True)
    data.rename(columns={col:col.replace("hetNdel","ΔhetN")},inplace=True)
    data.rename(columns={col:col.replace("WT","wild type")},inplace=True)
	
OMutantes=['wild type','ΔpatS','ΔhetN']
ATime=[24,30,48,50,72,96]
KSdistance=pd.DataFrame()
#Datasets
for t in ATime:
	for mut in OMutantes:
		LHistos=[]
		for name in data.columns:
			if re.search(mut,name) and re.search(str(t),name):
				LHistos.append(name)
		if len(LHistos)>1:
			perm = combinations(LHistos, 2)
			for i in list(perm):
				cum1=np.cumsum(data[i[0]])
				cum2=np.cumsum(data[i[1]])
				KolSmr=abs(cum1-cum2)	
				Tlabel='{}h'.format(t)
				common=mut+' '+Tlabel
				diff0=i[0].replace(Tlabel,'')
				diff0=diff0.replace(mut,'')
				diff1=i[1].replace(Tlabel,'')
				diff1=diff1.replace(mut,'')
				if diff0 not in DataSets:
					DataSets.append(diff0)
				if diff1 not in DataSets:
					DataSets.append(diff1)	
				row=pd.DataFrame([{'Dataset':diff0+'-'+diff1,'Mutant':mut,'Time':Tlabel,'Conserved':common,'Class':'Datasets','Different':diff0+'<->'+diff1,'KSdistance':max(KolSmr)}])	
				KSdistance=pd.concat([KSdistance, row], axis=0,ignore_index=True)
#Mutants				
	for dat in DataSets:	
		LHistos=[]
		for name in data.columns:
			if re.search(dat,name) and re.search(str(t),name):
				LHistos.append(name)
		if len(LHistos)>1:
			perm = combinations(LHistos, 2)
			for i in list(perm):
				cum1=np.cumsum(data[i[0]])
				cum2=np.cumsum(data[i[1]])
				KolSmr=abs(cum1-cum2)		
				common=dat+' {}h'.format(t)
				diff0=i[0].replace(Tlabel,'')
				diff0=diff0.replace(dat,'')
				diff1=i[1].replace(Tlabel,'')
				diff1=diff1.replace(dat,'')
				row=pd.DataFrame([{'Dataset':dat,'Mutant':diff0+'-'+diff1,'Time':Tlabel,'Conserved':common,'Class':'Mutants','Different':diff0+'<->'+diff1,'KSdistance':max(KolSmr)}])
				KSdistance=pd.concat([KSdistance, row], axis=0,ignore_index=True)
#Times
for dat in DataSets:
	for mut in OMutantes:	
		LHistos=[]
		for name in data.columns:
			if re.search(dat,name) and re.search(mut,name):
				LHistos.append(name)
		if len(LHistos)>1:
			perm = combinations(LHistos, 2)
			for i in list(perm):
				cum1=np.cumsum(data[i[0]])
				cum2=np.cumsum(data[i[1]])
				KolSmr=abs(cum1-cum2)		
				common=dat+' '+mut
				diff0=i[0].replace(dat,'')
				diff0=diff0.replace(mut,'')
				diff1=i[1].replace(dat,'')
				diff1=diff1.replace(mut,'')
				row=pd.DataFrame([{'Dataset':dat,'Mutant':mut,'Time':diff0+'-'+diff1,'Conserved':common,'Class':'Times','Different':diff0+'<->'+diff1,'KSdistance':max(KolSmr)}])
				KSdistance=pd.concat([KSdistance, row], axis=0,ignore_index=True)
				
KSdistance=KSdistance.sort_values(['Class','KSdistance'])
KSfile=os.path.dirname(os.getcwd()) + "/Plots/KSdistance.dat"
KSdistance.to_csv(KSfile, sep='\t', float_format='%.4f', index=False)

plotfile=os.path.dirname(os.getcwd()) + "/Plots/PlotKSdistance.dat"
KSdistance=KSdistance[~((KSdistance['Different'].str.contains('30|50'))&(KSdistance['Class']!='Datasets'))]
KSdistance.to_csv("PlotKSdistance.dat", sep='\t', float_format='%.4f', index=False)

plotfile=os.path.dirname(os.getcwd()) + "/Plots/PlotKSdistance.dat"
with open(plotfile, newline = '') as _file_2:
	KSdistance = pd.read_csv(_file_2, sep='\t', engine='python', index_col=0)

plt.style.use('seaborn-dark-palette')
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.weight'] ='bold'
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 1

plt.rcParams["figure.figsize"] = (16,6)
KSdistance=KSdistance.sort_values(['Class','Conserved','Different'])
filename="SwarmDifferent.png"


maxD=max(KSdistance['KSdistance'])
for c in np.unique(KSdistance['Class']):

	filename1="Swarm_Conserved_{}.png".format(c)
	filename2="Swarm_Conserved_{}.png".format(c)
	PlotData=KSdistance[KSdistance['Class']==c]
	
	PlotData=PlotData.sort_values(['Conserved','Different'])
	color=sns.color_palette("husl", len(np.unique(PlotData['Conserved'])))
	ax1=sns.swarmplot(x="Different",y="KSdistance",hue="Conserved",data=PlotData, palette=color,size=15,linewidth=1)
	# ~ ax1.grid(color='k', linestyle=':', linewidth=0.5)
	ax1.set_xticklabels(ax1.get_xticklabels(),rotation=10, ha="center")
	ax1.set_ylim((0,maxD+0.02))
	
	handles, labels = ax1.get_legend_handles_labels()
	r = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',visible=False)
	if c=='Datasets':
		ax1.set_title('dataset comparison', y=1.0, pad=+10, fontweight='bold')
		handles.insert(6,r)
		labels.insert(6,'')
		ax1.legend(handles,labels,frameon=True,fontsize=plt.rcParams['font.size']*0.95,loc='best',ncol=2)
	if c=='Mutants':    
		ax1.set_title('mutant comparison', y=1.0, pad=+10, fontweight='bold')                                                        
		ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.95,loc='best',ncol=1)	
	if c=='Times':  
		ax1.set_title('time comparison', y=1.0, pad=+10, fontweight='bold')                                                            
		ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.95,loc='best',ncol=2)
		
	ax1.set_xlabel('')
	ax1.set_ylabel('K-S distance',fontweight='bold',fontsize=plt.rcParams['font.size']*1.3)
	plt.tight_layout()
	plt.savefig(filename1, transparent=True)
	plt.close()
	
	
	PlotData=PlotData.sort_values(['Conserved','Different'])	
	color=sns.color_palette("husl", len(np.unique(PlotData['Different'])))
	fig, ax2 = plt.subplots(figsize=(10,7))
	sns.swarmplot(x="Conserved",y="KSdistance",hue="Different",data=PlotData, palette=color,size=15,linewidth=1,ax=ax2)	
	ax2.grid(axis='x',color='k', linestyle=':', linewidth=0.8)
	# ~ ax2.set_xticklabels(ax2.get_xticklabels(),rotation=25, ha="right")
	labels = [item.get_text() for item in ax2.get_xticklabels()]
	if c=='Datasets':		 
		for m in ['wild type ','ΔhetN ','ΔpatS ']:
			labels=[item.replace(m,'') for item in labels]
		ax2.set_xticklabels(labels,rotation=20, ha="center")
		ax2.text(1.2,-0.12, 'wild type',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(5.3,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(6.05,-0.12, 'ΔhetN',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(8.3,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(9.7,-0.12, 'ΔpatS',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
	if c=='Mutants':
		for m in ['Corrales ','Simulation ','Yoon ']:    
			labels=[item.replace(m,'') for item in labels]
		ax2.set_xticklabels(labels,rotation=20, ha="center")
		ax2.text(-0.3,-0.12,'Corrales',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(2.35,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(2.9,-0.12,'Simulation',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(6.35,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(7.9,-0.12,'Yoon',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
	if c=='Times':  
		for m in ['Corrales ','Simulation ','Yoon ']:    
			labels=[item.replace(m,'') for item in labels]
		ax2.set_xticklabels(labels,rotation=20, ha="center")
		ax2.text(0.3,-0.12,'Corrales',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(2.4,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(3.1,-0.12,'Simulation',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(5.4,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(6.1,-0.12,'Yoon',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		
	dx = 5/72.; dy = 0/72. 
	offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)

	# apply offset transform to all x ticklabels.
	for label in ax2.xaxis.get_majorticklabels():
		label.set_transform(label.get_transform() + offset)

	
	ax2.set_ylim((0,maxD+0.02))
	if len(np.unique(PlotData['Different']))>4:
		plt.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.9,loc='best',title=c,ncol=2)
	else:
		plt.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.9,loc='best',title=c,ncol=2)
		
	ax2.set_xlabel('')
	ax2.set_ylabel('K-S distance',fontweight='bold',fontsize=plt.rcParams['font.size']*1.3)
	plt.tight_layout()
	plt.savefig(filename2, transparent=True)
	plt.close()

#Subplots	
i=0
fig1, axs1 = plt.subplots(1, 3,figsize=(20, 8), sharey='row', gridspec_kw={'width_ratios': [1.,0.9,1]})
fig2, axs2 = plt.subplots(1, 3,figsize=(20, 8), sharey='row', gridspec_kw={'width_ratios': [1,0.9,1.1]})
filename1="SubSwarm_Conserved.png"
filename2="SubSwarm_Different.png"
	
maxD=max(KSdistance['KSdistance'])
for c in np.unique(KSdistance['Class']):
	ax1=axs1[i]
	ax2=axs2[i]
	
	PlotData=KSdistance[KSdistance['Class']==c]
	
	PlotData=PlotData.sort_values(['Conserved','Different'])
	color=sns.color_palette("husl", len(np.unique(PlotData['Conserved'])))
	sns.swarmplot(x="Different",y="KSdistance",hue="Conserved",data=PlotData, palette=color,size=15,linewidth=1,ax=ax1)
	# ~ ax1.grid(color='k', linestyle=':', linewidth=0.5)
	ax1.set_xticklabels(ax1.get_xticklabels(),rotation=10, ha="center")
	ax1.set_ylim((0,maxD+0.02))
	
	handles, labels = ax1.get_legend_handles_labels()
	r = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',visible=False)
	if c=='Datasets':
		ax1.set_title('dataset comparison', y=1.0, pad=+10, fontweight='bold')
		handles.insert(6,r)
		labels.insert(6,'')
		ax1.legend(handles,labels,frameon=True,fontsize=plt.rcParams['font.size']*0.95,loc='best',ncol=2)
	if c=='Mutants':    
		ax1.set_title('mutant comparison', y=1.0, pad=+10, fontweight='bold')                                                        
		ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.95,loc='best',ncol=1)	
	if c=='Times':  
		ax1.set_title('time comparison', y=1.0, pad=+10, fontweight='bold')                                                            
		ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.95,loc='best',ncol=2)
		
	PlotData=PlotData.sort_values(['Conserved','Different'])	
	color=sns.color_palette("husl", len(np.unique(PlotData['Different'])))
	sns.swarmplot(x="Conserved",y="KSdistance",hue="Different",data=PlotData, palette=color,size=15,linewidth=1,ax=ax2)	
	ax2.grid(axis='x',color='k', linestyle=':', linewidth=0.8)
	# ~ ax2.set_xticklabels(ax2.get_xticklabels(),rotation=25, ha="right")
	labels = [item.get_text() for item in ax2.get_xticklabels()]
	if c=='Datasets':		 
		for m in ['wild type ','ΔhetN ','ΔpatS ']:
			labels=[item.replace(m,'') for item in labels]
		ax2.set_xticklabels(labels,rotation=20, ha="center")
		ax2.text(1.2,-0.12, 'wild type',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(5.3,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(6.05,-0.12, 'ΔhetN',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(8.3,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(9.7,-0.12, 'ΔpatS',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
	if c=='Mutants':
		for m in ['Corrales ','Simulation ','Yoon ']:    
			labels=[item.replace(m,'') for item in labels]
		ax2.set_xticklabels(labels,rotation=20, ha="center")
		ax2.text(-0.3,-0.12,'Corrales',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(2.35,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(2.9,-0.12,'Simulation',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(6.35,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(7.9,-0.12,'Yoon',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
	if c=='Times':  
		for m in ['Corrales ','Simulation ','Yoon ']:    
			labels=[item.replace(m,'') for item in labels]
		ax2.set_xticklabels(labels,rotation=20, ha="center")
		ax2.text(0.3,-0.12,'Corrales',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(2.4,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(3.1,-0.12,'Simulation',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		ax2.text(5.4,-0.12, '|',fontsize=plt.rcParams['font.size']*1.5,fontweight='bold')
		ax2.text(6.1,-0.12,'Yoon',fontsize=plt.rcParams['font.size']*1.2,fontweight='bold')
		
	dx = 5/72.; dy = 0/72. 
	offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig2.dpi_scale_trans)

	# apply offset transform to all x ticklabels.
	for label in ax2.xaxis.get_majorticklabels():
		label.set_transform(label.get_transform() + offset)
	
	ax2.set_ylim((0,maxD+0.02))	
		
	handles, labels = ax2.get_legend_handles_labels()
	if c=='Datasets':
		ax2.legend(frameon=True,fontsize=plt.rcParams['font.size']*1.1,loc='upper left',title='dataset comparison',ncol=1)
	if c=='Mutants':
		ax2.legend(frameon=True,fontsize=plt.rcParams['font.size']*1.1,loc='best',title='mutant comparison',ncol=1)	
	if c=='Times':
		order=[0,1,3,2,4,5]
		ax2.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=True,fontsize=plt.rcParams['font.size']*1.1,loc='upper left',title='time comparison',ncol=2)
	
	ax1.set_xlabel('')
	ax2.set_xlabel('')
	if i>0:		
		ax1.set_ylabel('')
		ax2.set_ylabel('')
	else:
		ax1.set_ylabel('K-S distance',fontweight='bold',fontsize=plt.rcParams['font.size']*1.3)
		ax2.set_ylabel('K-S distance',fontweight='bold',fontsize=plt.rcParams['font.size']*1.3)
	i+=1
	
fig1.subplots_adjust(wspace=0.03,left=0.1, bottom=0.11, right=0.99, top=0.94)
fig1.savefig(filename1, transparent=True)
plt.close()
fig2.subplots_adjust(wspace=0.03,left=0.1, bottom=0.17, right=0.99, top=0.99)
fig2.savefig(filename2, transparent=True)
plt.close()

