import pandas
import re
import os
import os.path

with open("Histograms.dat", newline = '') as _file_1:
   histo = pandas.read_csv(_file_1, sep='\t', index_col=0)

with open("SimData.dat", newline = '') as _file_1:
	SimData = pandas.read_csv(_file_1, sep='\t', index_col=0)

with open("BorderHistograms.dat", newline = '') as _file_1:
   Bhisto = pandas.read_csv(_file_1, sep='\t', index_col=0)

sim = input("Enter the simulation name: ")

if sim!='ER':
	# new = input('Do you want to keep the previous as Old? (y/n): ')
	new='n'
	# ~ mut = input("Enter the strain (name/a): ")
	mut = 'a'
	# ~ Nsim = input("Enter number of repetitions of the simulation (nº/b): ")
	Nsim='b'
	# ~ Nexp = input("Enter number of repetitions of the experimental data (nº/b): ")
	Nexp='b'

if re.search('ER',sim):
		for mut in ['WT','patSdel','hetNdel','FpatSdel','patApatSdel','patXdel']:
			for t in [24,48,72]:
				try:
					del histo[mut+"Sim"+str(t)+"h"+"Mean"]
				except:
					print("Done")
				try:
					del histo[mut+"Sim"+str(t)+"h"+"MeanDev"]
				except:
					print("Done")
				try:
					del histo["N"+mut+"Sim"+str(t)+"h"+"Mean"]
				except:
					print("Done")
				try:
					del histo["N"+mut+"Sim"+str(t)+"h"+"MeanDev"]
				except:
					print("Done")

			try:
				del SimData[mut+"Mean"]
			except:
				print("Done")
			try:
				del SimData[mut+"MeanDev"]
			except:
				print("Done")
			try:
				del SimData[mut+"HetPer"]
			except:
				print("Done")
			try:
				del SimData[mut+"HetPerDev"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"Mean"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"MeanDev"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"HetPer"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"HetPerDev"]
			except:
				print("Done")

		for mut in ['patAdel','patAhetNdel']:
			for t in [24,48,72,96]:
				try:
					del Bhisto[mut+"Sim"+str(t)+"h"+"Mean"]
				except:
					print("Done")
				try:
					del Bhisto[mut+"Sim"+str(t)+"h"+"MeanDev"]
				except:
					print("Done")
				try:
					del Bhisto["N"+mut+"Sim"+str(t)+"h"+"Mean"]
				except:
					print("Done")
				try:
					del Bhisto["N"+mut+"Sim"+str(t)+"h"+"MeanDev"]
				except:
					print("Done")

			try:
				del SimData[mut+"HetPer"]
			except:
				print("Done")
			try:
				del SimData[mut+"HetPerDev"]
			except:
				print("Done")
			try:
				del SimData[mut+"HtRCon"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"HetPer"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"HetPerDev"]
			except:
				print("Done")
			try:
				del SimData["N"+mut+"HtRCon"]
			except:
				print("Done")
else:
	if (mut == 'a'):
        #
		for mut in ['WT','patSdel','hetNdel','patXpatSdel','patApatSdel','patXdel']:

			if new == 'y':
				SimData[mut+"Mean"] = SimData["N"+mut+"Mean"].copy()
				SimData[mut+"MeanDev"] = SimData["N"+mut+"MeanDev"].copy()
				SimData[mut+"HetPer"] = SimData["N"+mut+"HetPer"].copy()
				SimData[mut+"HetPerDev"] = SimData["N"+mut+"HetPerDev"].copy()
				SimData[mut+"HtRCon"] = SimData["N"+mut+"HtRCon"].copy()

			if not sim == 'Ori':
				nom="N"+mut
			else:
				nom=mut


			for t in [24,48,72]:
				if new == 'y':
					histo[mut+"Sim"+str(t)+"h"+"Mean"]=histo["N"+mut+"Sim"+str(t)+"h"+"Mean"].copy()
					histo[mut+"Sim"+str(t)+"h"+"MeanDev"]=histo["N"+mut+"Sim"+str(t)+"h"+"MeanDev"].copy()

				with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/histogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
					data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

				histoMean = data.iloc[:,0]

				if Nsim=='b':
					histoDev = data.iloc[:,2]
				else:
					histoDev = data.iloc[:,1]

				histo[nom+"Sim"+str(t)+"h"+"Mean"]=histoMean
				histo[nom+"Sim"+str(t)+"h"+"MeanDev"]=histoDev

			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/histogram_mean.dat".format(sim,mut), newline = '') as _file_2:
				Mean = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
				HetPer = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,mut), newline = '') as _file_2:
				HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

			MValues = Mean.iloc[:,0]
			HPValues = HetPer.iloc[:,0]
			HtCValues = HtRCon.iloc[:,0]

			if Nsim=='b':
				MDeviation = Mean.iloc[:,2]
				HPDeviation = HetPer.iloc[:,2]
			else:
				MDeviation = Mean.iloc[:,1]
				HPDeviation = HetPer.iloc[:,1]

			SimData[nom+"Mean"] = MValues
			SimData[nom+"MeanDev"] = MDeviation
			SimData[nom+"HetPer"] = HPValues
			SimData[nom+"HetPerDev"] = HPDeviation
			SimData[nom+"HtRCon"] = HtCValues

		for mut in ['patAdel','patAhetNdel']:

			if new == 'y':
				SimData[mut+"HetPer"] = SimData["N"+mut+"HetPer"].copy()
				SimData[mut+"HetPerDev"] = SimData["N"+mut+"HetPerDev"].copy()
				SimData[mut+"HtRCon"] = SimData["N"+mut+"HtRCon"].copy()

			if not sim == 'Ori':
				nom="N"+mut
			else:
				nom=mut

			for t in [24,48,72,96]:

				if new == 'y':
					Bhisto[mut+"Sim"+str(t)+"h"+"Mean"]=Bhisto["N"+mut+"Sim"+str(t)+"h"+"Mean"].copy()
					Bhisto[mut+"Sim"+str(t)+"h"+"MeanDev"]=Bhisto["N"+mut+"Sim"+str(t)+"h"+"MeanDev"].copy()

				with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/Bhistogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
					data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4,index_col=0)

				BhistoMean = data.iloc[:,0]
				if Nsim=='b':
					BhistoDev = data.iloc[:,2]
				else:
					BhistoDev = data.iloc[:,1]

				Bhisto[nom+"Sim"+str(t)+"h"+"Mean"]=BhistoMean
				Bhisto[nom+"Sim"+str(t)+"h"+"MeanDev"]=BhistoDev

			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
				HetPer = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4,index_col=0)
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,mut), newline = '') as _file_2:
				HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

			HPValues = HetPer.iloc[:,0]
			HtCValues = HtRCon.iloc[:,0]

			if Nsim=='b':
				HPDeviation = HetPer.iloc[:,2]
			else:
				HPDeviation = HetPer.iloc[:,1]

			SimData[nom+"HetPer"] = HPValues
			SimData[nom+"HetPerDev"] = HPDeviation
			SimData[nom+"HtRCon"] = HtCValues

		# if new == 'y':
		# 	SimData["hetFdelHtRCon"] = SimData["NhetFdelHtRCon"].copy()
		# with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,"hetFdel"), newline = '') as _file_2:
		# 	HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

		HtCValues = HtRCon.iloc[:,0]
		SimData["NhetFdelHtRCon"] = HtCValues

	else:
		if mut=='hetF':
			if new == 'y':
				SimData["hetFdelHtRCon"] = SimData["NhetFdelHtRCon"].copy()
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,"hetFdel"), newline = '') as _file_2:
				HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

			HtCValues = HtRCon.iloc[:,0]
			SimData["NhetFdelHtRCon"] = HtCValues

		elif not mut=='patA' and not mut=='patAhetN':
			if new == 'y':
				SimData[mut+"Mean"] = SimData["N"+mut+"Mean"].copy()
				SimData[mut+"MeanDev"] = SimData["N"+mut+"MeanDev"].copy()
				SimData[mut+"HetPer"] = SimData["N"+mut+"HetPer"].copy()
				SimData[mut+"HetPerDev"] = SimData["N"+mut+"HetPerDev"].copy()
				SimData[mut+"HtRCon"] = SimData["N"+mut+"HtRCon"].copy()

			if not sim == 'Ori':
				nom="N"+mut
			else:
				nom=mut


			for t in [24,48,72]:
				if new == 'y':
					histo[mut+"Sim"+str(t)+"h"+"Mean"]=histo["N"+mut+"Sim"+str(t)+"h"+"Mean"].copy()
					histo[mut+"Sim"+str(t)+"h"+"MeanDev"]=histo["N"+mut+"Sim"+str(t)+"h"+"MeanDev"].copy()

				with open(os.path.dirname(os.getcwd()) + "../Data_{}_{}/Histograms/histogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
					data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

				histoMean = data.iloc[:,0]

				if Nsim=='b':
					histoDev = data.iloc[:,2]
				else:
					histoDev = data.iloc[:,1]

				histo[nom+"Sim"+str(t)+"h"+"Mean"]=histoMean
				histo[nom+"Sim"+str(t)+"h"+"MeanDev"]=histoDev

			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/histogram_mean.dat".format(sim,mut), newline = '') as _file_2:
				Mean = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
				HetPer = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,mut), newline = '') as _file_2:
				HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

			MValues = Mean.iloc[:,0]
			HPValues = HetPer.iloc[:,0]
			HtCValues = HtRCon.iloc[:,0]

			if Nsim=='b':
				MDeviation = Mean.iloc[:,2]
				HPDeviation = HetPer.iloc[:,2]
			else:
				MDeviation = Mean.iloc[:,1]
				HPDeviation = HetPer.iloc[:,1]

			SimData[nom+"Mean"] = MValues
			SimData[nom+"MeanDev"] = MDeviation
			SimData[nom+"HetPer"] = HPValues
			SimData[nom+"HetPerDev"] = HPDeviation
			SimData[nom+"HtRCon"] = HtCValues

		else:

			if new == 'y':
				SimData[mut+"HetPer"] = SimData["N"+mut+"HetPer"].copy()
				SimData[mut+"HetPerDev"] = SimData["N"+mut+"HetPerDev"].copy()
				SimData[mut+"HtRCon"] = SimData["N"+mut+"HtRCon"].copy()

			if not sim == 'Ori':
				nom="N"+mut
			else:
				nom=mut

			for t in [24,48,72,96]:

				if new == 'y':
					Bhisto[mut+"Sim"+str(t)+"h"+"Mean"]=Bhisto["N"+mut+"Sim"+str(t)+"h"+"Mean"].copy()
					Bhisto[mut+"Sim"+str(t)+"h"+"MeanDev"]=Bhisto["N"+mut+"Sim"+str(t)+"h"+"MeanDev"].copy()

				with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/Histograms/Bhistogram_med_t={}.dat".format(sim,mut,t), newline = '') as _file_2:
					data = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4,index_col=0)

				BhistoMean = data.iloc[:,0]
				if Nsim=='b':
					BhistoDev = data.iloc[:,2]
				else:
					BhistoDev = data.iloc[:,1]

				Bhisto[nom+"Sim"+str(t)+"h"+"Mean"]=BhistoMean
				Bhisto[nom+"Sim"+str(t)+"h"+"MeanDev"]=BhistoDev

			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/per_cent_het.dat".format(sim,mut), newline = '') as _file_2:
				HetPer = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4,index_col=0)
			with open(os.path.dirname(os.getcwd()) + "/Data_{}_{}/concentration_HetR.dat".format(sim,mut), newline = '') as _file_2:
				HtRCon = pandas.read_csv(_file_2, sep='\t\t', header=None, engine='python', skipfooter=4, index_col=0)

			HPValues = HetPer.iloc[:,0]
			HtCValues = HtRCon.iloc[:,0]

			if Nsim=='b':
				HPDeviation = HetPer.iloc[:,2]
			else:
				HPDeviation = HetPer.iloc[:,1]

			SimData[nom+"HetPer"] = HPValues
			SimData[nom+"HetPerDev"] = HPDeviation
			SimData[nom+"HtRCon"] = HtCValues

histo.to_csv('Histograms.dat', sep='\t', float_format='%.4f')
SimData.to_csv('SimData.dat', sep='\t', float_format='%.1f')
Bhisto.to_csv('BorderHistograms.dat', sep='\t', float_format='%.4f')

if sim == 'ER':
	os.system('python erasePlots.py')
elif new == 'y':
	if Nsim=='b':
		Nsim="1"
		Nexp="1"
	os.system('python erasePlots.py')
	order="gnuplot -c 2hree.plt "+ Nsim +" "+ Nexp
	os.system(order)
else:
	order="gnuplot -c 3stupendas.plt "+ Nsim +" "+ Nexp
	os.system(order)
