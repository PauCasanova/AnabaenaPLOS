import pandas
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import math
import re
from matplotlib.ticker import MaxNLocator

OMutantes=['WT','patSdel','hetNdel']
NMutantes=['patAdel','patAhetNdel']

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

lalist=['no leakage','(tetramers&0.53KR)','Gene noise ','Cell noise ', 'WT','patXpatSdel','patAhetNdel','patApatSdel','patSdel','patXdel','hetNdel','patAdel']
labels={'no leakage':'no leak','(tetramers&0.53KR)':' r($R_j,A_j,I_j,G_i$)','Gene noise ':'$\Omega_{\Phi}$','Cell noise ':'$\Omega_{\Lambda}$','WT':'wild type','patSdel':'$\Delta\it{patS}$','patXdel':'$\Delta\it{patX}$','patXpatSdel':'$\Delta\it{patX}\Delta\it{patS}$','hetNdel':'$\Delta\it{hetN}$','patAdel':'$\Delta\it{patA}$','patAhetNdel':'$\Delta\it{patA}\Delta\it{hetN}$','patApatSdel':'$\Delta\it{patA}\Delta\it{patS}$'} 

#Histograms intervals
def IntHistograms(sim,HMutantes,Htag,colors):
	plt.rcParams["figure.figsize"] = (13,5)
	for t in [24,48,72]:	
		filename="{}Histogram{}.png".format(Htag,t)
		bar_width = 1/len(HMutantes)-len(HMutantes)*0.01
		opacity = 0.8
		index = [x for x in np.arange(0, 27, 1)]
		indexT = [x+bar_width*(len(HMutantes)-1)/2 for x in index]
		Maxfrec=0
		Itag=[]
		for i in np.arange(0, 27, 1):
			if not i%2:
				Itag.append(i)
			else:
				Itag.append('')
		
		for mut in HMutantes:
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/histogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
				data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			
			histoMean = data.iloc[:,0].tolist()
			histoDev = data.iloc[:,2]
				
			if max(histoMean)+max(histoDev)>Maxfrec:
				Maxfrec=max(histoMean)+max(histoDev)+0.01
			
			plabel=mut
			for l in lalist:
				plabel=plabel.replace(l,labels[l])
				
			rectsl = plt.bar(index, histoMean, bar_width,alpha=opacity,color=colors[mut],edgecolor='white',label=r'{}'.format(plabel), yerr=histoDev,align='center', ecolor='black', capsize=5)
			index = [x+bar_width for x in index]

		plt.xlabel('vegetative interval length', fontweight='bold')
		plt.ylabel('fraction of intervals', fontweight='bold')
		
		plt.ylim((0, Maxfrec))
		if len(HMutantes)%2:
			plt.xlim((-bar_width,max(index)+bar_width/2))
		else:
			plt.xlim((-1.5*bar_width,max(index)+bar_width/2))
			
		plt.xticks(indexT,Itag)
		
		if len(HMutantes)>3:
			plt.legend(loc='best',ncol=int(len(HMutantes)/2), frameon=False)
		else:
			plt.legend(loc='best', frameon=False)
			
		plt.title('Histogram for {}h'.format(t), y=1.0, pad=-20, fontweight='bold')

		plt.locator_params(axis="y", nbins=5)
		plt.tight_layout()
		plt.savefig(filename, transparent=True)
		plt.close()

#Mean and percentage
def MeanPer(sim,PMutantes,Ptag,colors,mMin,mMax,dm,pMin,pMax,dp,ptype):
	if(ptype=='L'):
		plt.rcParams["figure.figsize"] = (11,5)
	elif(ptype=='B'):
		plt.rcParams["figure.figsize"] = (7,5)
	else:
		plt.rcParams["figure.figsize"] = (7,5)
		
	Meanname="{}Mean.png".format(Ptag)
	Pername="{}Percentage.png".format(Ptag)
	fig1 = plt.figure() #figure object
	fig2 = plt.figure()
	ax1 = fig1.gca() #axis object
	ax2 = fig2.gca()

	with open(os.path.dirname(os.getcwd()) + "/Plots/ExpData/ExpData.dat", newline = '') as _file_2:
		ExpData = pandas.read_csv(_file_2, sep='\t', engine='python', index_col=0)

	for mut in PMutantes:
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/histogram_mean.dat".format(sim,mut), newline = '') as _file_2:
			Mean = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
			HetPer = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,mut), newline = '') as _file_2:
			HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
		
		IValues = Mean.index.tolist()
		MValues = Mean.iloc[:,0].tolist()
		HPValues = HetPer.iloc[:,0].tolist()
		HtCValues = HtRCon.iloc[:,0].tolist()

		MDeviation = Mean.iloc[:,2]
		HPDeviation = HetPer.iloc[:,2]
			
		plabel=mut
		for l in lalist:
			plabel=plabel.replace(l,labels[l])
			
		ax1.plot(IValues,MValues,color=colors[mut],linewidth=2.5,label=r'{}'.format(plabel))	
		ax1.fill_between(IValues, MValues-MDeviation, MValues+MDeviation,color=colors[mut],alpha=0.5)
		
		ax2.plot(IValues,HPValues,color=colors[mut],linewidth=2.5,label=r'{}'.format(plabel))	
		ax2.fill_between(IValues, HPValues-HPDeviation, HPValues+HPDeviation,color=colors[mut],alpha=0.5)
		
		if mut in OMutantes:
			IExp=ExpData["{}Mean".format(mut)].index.tolist()
			MVExp=ExpData["{}Mean".format(mut)].tolist()
			MDExp=ExpData["{}MeanDev".format(mut)].tolist()
			ax1.errorbar(IExp,MVExp, yerr=MDExp,color=colors[mut], fmt="o", capsize=5)
			HVExp=ExpData["{}HetPer".format(mut)].tolist()
			HDExp=ExpData["{}HetPerDev".format(mut)].tolist()
			ax2.errorbar(IExp,HVExp, yerr=HDExp,color=colors[mut], fmt="o", capsize=5)

	ax1.set_ylabel('mean length of intervals', fontweight='bold')
	ax1.set_yticks(np.arange(mMin-1, mMax+1, dm))
	ax1.set_ylim((mMin, mMax))

	ax2.set_ylabel('percentage of heterocysts', fontweight='bold')
	ax2.set_yticks(np.arange(pMin-1, pMax+1, dp))
	ax2.set_ylim((pMin, pMax))

	ax1.set_xlabel('time (h)', fontweight='bold')
	ax1.set_xticks(np.arange(0, 100, 24))
	ax1.set_xlim((18, 80))
	
	if(ptype=='S'):
		ax1.legend().remove()
	elif(ptype=='B'):
		if len(PMutantes)>3:
			ax1.legend(loc='best',ncol=int(len(PMutantes)/3), frameon=False,fontsize=plt.rcParams['font.size']*0.9)
		else:
			ax1.legend(loc='best', frameon=False)		
	else:
		ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*1,loc='center left', bbox_to_anchor=(1, 0.5))
		
	fig1.tight_layout()
	fig1.savefig(Meanname, transparent=True)
	plt.close()

	ax2.set_xlabel('time (h)', fontweight='bold')
	ax2.set_xticks(np.arange(0, 100, 24))
	ax2.set_xlim((18, 80))
		
	if(ptype=='S'):
		ax2.legend().remove()
	elif(ptype=='B'):
		if len(PMutantes)>3:
			ax2.legend(loc='best',ncol=int(len(PMutantes)/3), frameon=False,fontsize=plt.rcParams['font.size']*0.9)
		else:
			ax2.legend(loc='best', frameon=False)
	else:
		ax2.legend(frameon=False,fontsize=plt.rcParams['font.size']*1,loc='center left', bbox_to_anchor=(1, 0.5))
	
	fig2.tight_layout()
	fig2.savefig(Pername, transparent=True)
	plt.close()

#Histograms Border
def BorHistograms(sim,BMutantes,Btag,colors,exp):
	plt.rcParams["figure.figsize"] = (10,5)
	
	with open(os.path.dirname(os.getcwd()) + "/Plots/ExpData/BorderHistograms.dat", newline = '') as _file_2:
		ExpData = pandas.read_csv(_file_2, sep='\t', engine='python', skipfooter=1, index_col=0)
		
	for t in [24,48,72,96]:
		
		filename="{}Histogram{}.png".format(Btag,t)
		bar_width = 1/(len(BMutantes)+exp)-(len(BMutantes)+exp)*0.001
		opacity = 1
		index = [x for x in np.arange(0,len(ExpData.index), 1)]
		indexT = [x+bar_width*(len(BMutantes)-1+exp)/2 for x in index]
		Maxfrec=0
		Itag=[]
		for i in np.arange(0, len(ExpData.index), 1):
			if not i%1:
				Itag.append(i)
			else:
				Itag.append('')
				
		for mut in BMutantes:
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/Bhistogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
				data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=6,index_col=0)
			BhistoMean = data.iloc[:,0].tolist()
			BhistoDev = data.iloc[:,2]
				
			if max(BhistoMean)+max(BhistoDev)>Maxfrec:
				Maxfrec=max(BhistoMean)+max(BhistoDev)+0.01
			
			plabel=mut
			for l in lalist:
				plabel=plabel.replace(l,labels[l])
			if exp>0:
				if mut in NMutantes:
					MVExp=ExpData["{}Exp{}hMean".format(mut,t)].tolist()
					MDExp=ExpData["{}Exp{}hMeanDev".format(mut,t)].tolist()
					plt.bar(index, MVExp, bar_width,alpha=opacity,color=colors['{}Exp'.format(mut)],edgecolor='white',label=r'{} experimental'.format(plabel), yerr=MDExp,align='center', ecolor='black', capsize=10)
					index = [x+bar_width for x in index]
			
			plt.bar(index, BhistoMean, bar_width,alpha=opacity,color=colors[mut],edgecolor='white',label=r'{}'.format(plabel), yerr=BhistoDev,align='center', ecolor='black', capsize=10)
			index = [x+bar_width for x in index]	
			
		Dy=(indexT[2]-indexT[1])/2
		for y in indexT:
			pl=y+Dy
			plt.plot([pl,pl],[0,1], color='black', lw = 0.5, ls='--')
			
		plt.xlabel('terminal heterocysts', fontweight='bold')
		plt.ylabel('fraction of filaments', fontweight='bold')
		plt.yticks(np.arange(0, 1.1, 0.2))
		plt.ylim((0, 1))
		plt.xlim((-bar_width,max(index)+bar_width/2))
		plt.xticks(indexT,Itag)
		
		plt.legend(bbox_to_anchor=(1, 1),loc=1, frameon=False,fontsize=plt.rcParams['font.size']*0.5, title='{}h'.format(t))
		# ~ plt.legend(frameon=False,fontsize=plt.rcParams['font.size']*1,loc='center left', bbox_to_anchor=(1, 0.5), title='{}h'.format(t))

		plt.tight_layout()
		plt.savefig(filename, transparent=True)
		plt.close()
	
	Meanname="{}BHet.png".format(Btag)
	Compname="{}CompHet.png".format(Btag)
	fig1 = plt.figure(figsize=(12, 5)) #figure object
	fig2 = plt.figure(figsize=(12, 5))
	ax1 = fig1.gca() #axis object
	ax2 = fig2.gca()
		
	with open(os.path.dirname(os.getcwd()) + "/Plots/ExpData/ExpDataBord.dat", newline = '') as _file_2:
		ExpData = pandas.read_csv(_file_2, sep='\t', engine='python', index_col=0)	
	with open(os.path.dirname(os.getcwd()) + "/Plots/ExpData/ExpDataBordII.dat", newline = '') as _file_2:
		Exp2Data = pandas.read_csv(_file_2, sep='\t', engine='python', index_col=0)	
		
	for mut in BMutantes:
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/n_heterocists_boundaries.dat".format(sim,mut), newline = '') as _file_2:
			Mean = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/n_heterocists.dat".format(sim,mut), newline = '') as _file_2:
			hetMean = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
		
		IValues = Mean.index.tolist()
		MValues = Mean.iloc[:,0]
		MDeviation = Mean.iloc[:,2]
		AValues = hetMean.iloc[:,0]
		ADeviation = hetMean.iloc[:,2]	
		
		DV=np.array([2*a for a in MValues])
		DDev=np.array([2*a for a in MDeviation])	
		
		PV=np.array([a-2*c for a,c in zip(AValues,MValues)])
		PDev=np.array([math.sqrt(b**2+(2*d)**2) for b,d in zip(ADeviation,MDeviation)])			
		
		if re.search('hetN',mut):
			Palpha=0.4
		else:
			Palpha=0.8			
		
		plabel=mut
		for l in lalist:
			plabel=plabel.replace(l,labels[l])
		if mut in NMutantes:
			IExp=ExpData["{}Mean".format(mut)].index.tolist()
			MVExp=ExpData["{}Mean".format(mut)].tolist()
			MDExp=ExpData["{}MeanDev".format(mut)].tolist()
			ax1.errorbar(IExp,MVExp, yerr=MDExp,color=colors['{}Exp'.format(mut)],label=r'{} experimental'.format(plabel), fmt="o", capsize=5)
			MDExp=Exp2Data["{}Mean".format(mut)].tolist()
			MDDExp=Exp2Data["{}MeanDev".format(mut)].tolist()
			ax2.errorbar(IExp,MDExp, yerr=MDDExp,color=colors['{}Exp'.format(mut)],label=r'Border in {} experimental'.format(plabel), fmt="o", capsize=5)
			
			
		ax1.plot(IValues,MValues,color=colors[mut],linewidth=2.5,label=r'{}'.format(plabel))	
		ax1.fill_between(IValues, MValues-MDeviation, MValues+MDeviation,color=colors[mut],alpha=Palpha)
		
		if mut in NMutantes:
			ax2.plot(IValues,PV,color=colors['{}het'.format(mut)],linewidth=2.5,label=r'Internal in {}'.format(plabel))		
			ax2.fill_between(IValues, PV-PDev, PV+PDev,color=colors['{}het'.format(mut)],alpha=0.5)		
			ax2.plot(IValues,DV,color=colors[mut],linewidth=2.5,label=r'Border in {}'.format(plabel))
			ax2.fill_between(IValues, DV-DDev, DV+DDev,color=colors[mut],alpha=Palpha)	
		else:
			ax2.plot(IValues,DV,color=colors[mut],linewidth=2.5,label=r'Border in {}'.format(plabel))
			ax2.fill_between(IValues, DV-DDev, DV+DDev,color=colors[mut],alpha=Palpha)
				
		
		
	ax1.set_ylabel('number of border heterocysts', fontweight='bold')
	ax1.yaxis.set_major_locator(MaxNLocator(5))
	ax1.set_ylim((0, 2.5))

	ax1.set_xlabel('time (h)', fontweight='bold')
	ax1.set_xticks(np.arange(0, 100, 24))
	ax1.set_xlim((18, 98))
	
	# ~ handles, tags = ax1.get_legend_handles_labels()
	# ~ if len(handles)>3:
		# ~ ax1.legend(loc='upper right',ncol=int(len(handles)/2), frameon=False,fontsize=plt.rcParams['font.size']*0.69, bbox_to_anchor=(1.01, 0.9999))
	# ~ ax1.legend().remove()
	ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*1,loc='center left', bbox_to_anchor=(1, 0.5))
	
	fig1.tight_layout()
	fig1.savefig(Meanname, transparent=True)
	plt.close()
		
	ax2.set_ylabel('number of heterocysts', fontweight='bold')
	ax2.yaxis.set_major_locator(MaxNLocator(5))
	ax2.set_ylim((0, 5))

	ax2.set_xlabel('time (h)', fontweight='bold')
	ax2.set_xticks(np.arange(0, 100, 24))
	ax2.set_xlim((18, 98))
	
	# ~ handles, tags = ax2.get_legend_handles_labels()
	# ~ if len(handles)>3:
		# ~ ax2.legend(loc='upper right',ncol=int(len(handles)/2), frameon=False,fontsize=plt.rcParams['font.size']*0.69, bbox_to_anchor=(1.01, 0.9999))
	# ~ ax2.legend().remove()
	ax2.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.85,loc='center left', bbox_to_anchor=(1, 0.5))
		
	fig2.tight_layout()
	fig2.savefig(Compname, transparent=True)
	plt.close()


#Mean and percentage bigplot
def MeanPerSubplots(sim,PMutantes,Ptag,colors,mMin,mMax,dm,pMin,pMax,dp,nrow,ncol):
	if ncol>1:
		plt.rcParams["figure.figsize"] = (15,17)
		plt.rcParams['font.size'] = 16
	else: 
		plt.rcParams["figure.figsize"] = (20,7)
	
	Name="{}.png".format(Ptag)
	
	if ncol>1:
		fig, axs = plt.subplots(nrow, ncol, sharex='col')
		fig.subplots_adjust(wspace=0.1, hspace=0.06)
	else: 
		fig, axs = plt.subplots(ncol, nrow, sharey='row')
		fig.subplots_adjust(wspace=0.1, hspace=0.06)
	
	with open(os.path.dirname(os.getcwd()) + "/Plots/ExpData/ExpData.dat", newline = '') as _file_2:
		ExpData = pandas.read_csv(_file_2, sep='\t', engine='python', index_col=0)
	
	for i in range(nrow):
		PartMutantes=['WT']
		a=2*i+1
		PartMutantes.append(PMutantes[a])
		PartMutantes.append(PMutantes[a+1])
		# ~ PartMutantes.append(PMutantes[a+2])
		if ncol>1:
			ax=axs[i, j]
		else:
			ax=axs[i]
			if i==0:
				ax.set_title('Cell noise', fontweight='bold')
			elif i==1:
				ax.set_title('Gene noise', fontweight='bold')
			elif i==2:
				ax.set_title('Overall noise', fontweight='bold')

		for j in range(ncol):
			for mut in PartMutantes:
				with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/histogram_mean.dat".format(sim,mut), newline = '') as _file_2:
					Mean = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
				with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
					HetPer = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
				
				IValues = Mean.index.tolist()
				MValues = Mean.iloc[:,0].tolist()
				HPValues = HetPer.iloc[:,0].tolist()

				MDeviation = Mean.iloc[:,2]
				HPDeviation = HetPer.iloc[:,2]
				
				plabel=mut
				for l in lalist:
					plabel=plabel.replace(l,labels[l])
				
				if j==0:
					ax.plot(IValues,MValues,color=colors[mut],linewidth=2.5,label=r'{}'.format(plabel))	
					ax.fill_between(IValues, MValues-MDeviation, MValues+MDeviation,color=colors[mut],alpha=0.5)
				elif j==1:
					axs[i, j].plot(IValues,HPValues,color=colors[mut],linewidth=2.5)	
					axs[i, j].fill_between(IValues, HPValues-HPDeviation, HPValues+HPDeviation,color=colors[mut],alpha=0.5)
				
				if mut in OMutantes:
					IExp=ExpData["{}Mean".format(mut)].index.tolist()
					if j==0:
						MVExp=ExpData["{}Mean".format(mut)].tolist()
						MDExp=ExpData["{}MeanDev".format(mut)].tolist()
						ax.errorbar(IExp,MVExp, yerr=MDExp,color=colors[mut], fmt="o", capsize=5)
					elif j==1:
						HVExp=ExpData["{}HetPer".format(mut)].tolist()
						HDExp=ExpData["{}HetPerDev".format(mut)].tolist()
						axs[i, j].errorbar(IExp,HVExp, yerr=HDExp,color=colors[mut], fmt="o", capsize=5)			
			if i==0:
				if j==0:	
					if ncol>1:
						ax.set_title('mean length of intervals', fontweight='bold')
					else:
						ax.set_ylabel('mean length of intervals', fontweight='bold')
						
					ax.set_yticks(np.arange(mMin-1, mMax+1, dm))
					ax.set_ylim((mMin, mMax))
				elif j==1:
					axs[i, j].set_title('percentage of heterocysts', fontweight='bold')
					axs[i, j].set_yticks(np.arange(mMin-1, pMax+1, dp))	
					axs[i, j].set_ylim((pMin, pMax))

				ax.set_xticks(np.arange(0, 100, 24))
				ax.set_xlim((18, 80))
			else:				
				if ncol>1:
					ax.sharey(axs[0, j])
				else:
					ax.sharex(axs[0])							
				
			if i==int(nrow/2):
				ax.set_xlabel('time (h)', fontweight='bold')
				
			ax.legend(loc='best', frameon=False,fontsize=plt.rcParams['font.size']*0.9)
			
	fig.savefig(Name, transparent=True)
	plt.close()

#Histograms Clusters
def ClustHistograms(sim,CMutantes,Ctag,colors,mMin,mMax,dm,pMin,pMax,dp):
	plt.rcParams["figure.figsize"] = (6,6)
	plt.rcParams['font.size'] = 21
	for t in [24,48,72]:		
		filename1="{}Histogram{}.png".format(Ctag,t)
		filename2="{}HistogramN{}.png".format(Ctag,t)	
		fig1 = plt.figure() #figure object
		fig2 = plt.figure()
		ax1 = fig1.gca() #axis object
		ax2 = fig2.gca()
		
		bar_width = 1/len(CMutantes)-len(CMutantes)*0.0001
		opacity = 1
		index = [x for x in np.arange(0, 5, 1)]
		indexT = [x+bar_width*(len(CMutantes)-1)/2 for x in index]
		Maxfrec=0
		MaxY=0
		Itag=[]
		
		for i in np.arange(1, 6, 1):
			if not i%1:
				Itag.append(i)
			else:
				Itag.append('')
				
		for mut in CMutantes:
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/CHhistogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
				data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=6,index_col=0)
			ChistoMean = data.iloc[:,0].tolist()
			ChistoDev = data.iloc[:,2]
			
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/histogram_nclusters.dat".format(sim,mut), newline = '') as _file_2:
				nclusters = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			
			ChistoN=np.array([a*nclusters.iloc[t,0]  for a in ChistoMean])
			ChistoNDev=np.array([a*(nclusters.iloc[t,0]+nclusters.iloc[t,2])  for a in ChistoDev])
			
			if max(ChistoMean)+max(ChistoDev)>Maxfrec:
				Maxfrec=max(ChistoMean)+max(ChistoDev)+0.01
				
			if nclusters.iloc[t,0]>MaxY:
				MaxY=int(nclusters.iloc[t,0]+0.1)
			
			plabel=mut
			for l in lalist:
				plabel=plabel.replace(l,labels[l])
				
			ax1.bar(index, ChistoMean, bar_width,alpha=opacity,color=colors[mut],edgecolor='white',label=r'{}'.format(plabel), yerr=ChistoDev,align='center', ecolor='black', capsize=10)
			ax2.bar(index, ChistoN, bar_width,alpha=opacity,color=colors[mut],edgecolor='white',label=r'{}'.format(plabel), yerr=ChistoNDev,align='center', ecolor='black', capsize=10)
			index = [x+bar_width for x in index]
						
		ax1.set_xlabel('number of heterocysts', fontweight='bold')
		ax1.set_ylabel('fraction of clusters', fontweight='bold')
		ax1.set_yticks(np.arange(0, 1.1, 0.2))
		ax1.set_ylim((0, 1))
		ax1.set_xlim((-bar_width,max(index)+bar_width/2))
		ax1.set_xticks(indexT,Itag)
		if t==24:
			if len(CMutantes)>3:
				ax1.legend(loc='best',ncol=int(len(CMutantes)/2), frameon=False)
			else:
				ax1.legend(bbox_to_anchor=(1, 0.85),loc=1, frameon=False)
		else:
			ax1.legend().remove()
				
		ax1.set_title('{}h'.format(t), y=1.0, pad=-30, fontweight='bold')
		fig1.tight_layout()
		fig1.savefig(filename1, transparent=True)
		plt.close()
		
		ax2.set_xlabel('number of heterocysts', fontweight='bold')
		ax2.set_ylabel('number of clusters', fontweight='bold')
		ax2.yaxis.set_major_locator(MaxNLocator(5,integer=True))
		ax2.set_ylim((0, MaxY+0.01))
		ax2.set_xlim((-bar_width,max(index)+bar_width/2))
		ax2.set_xticks(indexT,Itag)
		
		if t==48:
			if len(CMutantes)>3:
				ax2.legend(loc='best',ncol=int(len(CMutantes)/2), frameon=False)
			else:
				ax2.legend(bbox_to_anchor=(1, 0.85),loc=1, frameon=False)
		else:
			ax2.legend().remove()
			
		ax2.set_title('{}h'.format(t), y=1.0, pad=-30, fontweight='bold')
		fig2.tight_layout()
		fig2.savefig(filename2, transparent=True)
		plt.close()
		
	plt.rcParams["figure.figsize"] = (10,5)
	plt.rcParams['font.size'] = 18
	
	nclustersname="{}nclusters.png".format(Ctag)
	nclushetname="{}clustohetratio.png".format(Ctag)
	fig1 = plt.figure() #figure object
	fig2 = plt.figure()
	ax1 = fig1.gca() #axis object
	ax2 = fig2.gca()
	plt.rcParams["figure.figsize"] = (7,5)	
	
	for mut in CMutantes:
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/histogram_nclusters.dat".format(sim,mut), newline = '') as _file_2:
			nclusters = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
			phet = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			
		with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_clusters.dat".format(sim,mut), newline = '') as _file_2:
			pclusters = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			
		IValues = nclusters.index.tolist()
		CValues = nclusters.iloc[:,0].tolist()
		CDeviation = nclusters.iloc[:,2]
		PHValues = phet.iloc[:,0].tolist()
		PHDeviation = phet.iloc[:,2]
		PCValues = pclusters.iloc[:,0].tolist()
		PCDeviation = pclusters.iloc[:,2]
				
		plabel=mut
		for l in lalist:
			plabel=plabel.replace(l,labels[l])
			
		ax1.plot(IValues,CValues,color=colors[mut],linewidth=2.5,label=r'{}'.format(plabel))	
		ax1.fill_between(IValues, CValues-CDeviation, CValues+CDeviation,color=colors[mut],alpha=0.2)
		
		ax2.plot(IValues,PCValues,color=colors[mut],linewidth=2.5,label=r'Clusters in {}'.format(plabel))
		ax2.fill_between(IValues, PCValues-PCDeviation, PCValues+PCDeviation,color=colors[mut],alpha=0.2)
		
		ax2.plot(IValues,PHValues,color=colors[mut],linestyle='--',linewidth=2.5,label=r'Heterocysts in {}'.format(plabel))	
		ax2.fill_between(IValues, PHValues-PHDeviation, PHValues+PHDeviation,color=colors[mut],alpha=0.2)			
	
	ax1.set_ylabel('number of clusters', fontweight='bold')
	ax1.set_yticks(np.arange(mMin-1, mMax+1, dm))
	ax1.set_ylim((mMin, mMax))
	
	ax1.set_xlabel('time (h)', fontweight='bold')
	ax1.set_xticks(np.arange(0, 100, 24))
	ax1.set_xlim((18, 80))
	
	if len(CMutantes)>3:
		ax1.legend(loc='best',ncol=int(len(CMutantes)/3), frameon=False,fontsize=plt.rcParams['font.size']*0.9)
	else:
		ax1.legend(loc='best', frameon=False)
	ax1.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.9,loc='center left', bbox_to_anchor=(1, 0.5))
		
	fig1.tight_layout()
	fig1.savefig(nclustersname, transparent=True)
	plt.close()

	ax2.set_ylabel('ratio relative to total cells', fontweight='bold')
	ax2.set_yticks(np.arange(pMin-1, pMax+1, dp))
	ax2.set_ylim((pMin, pMax))
	
	ax2.set_xlabel('time (h)', fontweight='bold')
	ax2.set_xticks(np.arange(0, 100, 24))
	ax2.set_xlim((18, 80))
	
	if len(CMutantes)>3:
		ax2.legend(loc='best',ncol=int(len(CMutantes)/3), frameon=False,fontsize=plt.rcParams['font.size']*0.9)
	else:
		ax2.legend(loc='best', frameon=False)
	ax2.legend(frameon=False,fontsize=plt.rcParams['font.size']*0.9,loc='center left', bbox_to_anchor=(1, 0.5))
	
	fig2.tight_layout()
	fig2.savefig(nclushetname, transparent=True)
	plt.close()
