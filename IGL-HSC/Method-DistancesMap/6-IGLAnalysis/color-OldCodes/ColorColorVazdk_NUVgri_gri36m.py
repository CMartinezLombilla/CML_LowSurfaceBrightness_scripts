import numpy as np
import os 
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()
import math
from astropy.stats.funcs import biweight_location
from numpy import mean, sqrt, square, arange
from matplotlib import rc, rcParams
from matplotlib import gridspec

plt.close('all')

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern'], 'size'   : 30})

galaxies=['NGC5907', 'NGC4565', 'NGC4472']
kpc_betweenProfiles = ['1.5', '3.0']
color = ['g', 'darkorange', 'purple']

plt.figure(0,figsize=(28,8),facecolor='white')
gs = gridspec.GridSpec(1,3)
subPlot = 0

for galaxy in galaxies:

	Central_gri_36m = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_Central_gri-3.6micras_BINN.txt')		
	Central_NUV_gri = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_Central_NUV-gri_BINN.txt')		

	kpc15_gri_36m = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_1.5kpc_gri-3.6micras_BINN.txt')		
	kpc15_NUV_gri = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_1.5kpc_NUV-gri_BINN.txt')		

	kpc30_gri_36m = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_3.0kpc_gri-3.6micras_BINN.txt')		
	kpc30_NUV_gri = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/Color_3.0kpc_NUV-gri_BINN.txt')		

# 	Central_gri_36m = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/3.6micras/IRAC_noS4G/Color_Central_gri-3.6micras_BINN.txt')		
# 	Central_NUV_gri = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/3.6micras/IRAC_noS4G/Color_Central_NUV-gri_BINN.txt')		
# 
# 	kpc15_gri_36m = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/3.6micras/IRAC_noS4G/Color_1.5kpc_gri-3.6micras_BINN.txt')		
# 	kpc15_NUV_gri = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/3.6micras/IRAC_noS4G/Color_1.5kpc_NUV-gri_BINN.txt')		
# 
# 	kpc30_gri_36m = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/3.6micras/IRAC_noS4G/Color_3.0kpc_gri-3.6micras_BINN.txt')		
# 	kpc30_NUV_gri = np.loadtxt('/Users/cris/Tesis/scratch/Stripe82_fitting/'+galaxy+'/3.6micras/IRAC_noS4G/Color_3.0kpc_NUV-gri_BINN.txt')		



	# MODELO DE POBLACIONES ESTELARES VAZDEKIS ET AL 2015

	modSSP=np.loadtxt('/Users/cris/Tesis/Mod_EMILES_SSP/mag_tmpRHvWmg.mag')
	metallicity = modSSP[:,0]
	age = modSSP[:,1]
	NUV_mag = modSSP[:,2]
	r_mag = modSSP[:,3]
	spitz36m_mag = modSSP[:,4]

	# COLOR NUV-gri
	NUV_gri = NUV_mag - r_mag

	NUV_gri_227 = NUV_gri[0:53]
	NUV_gri_179 = NUV_gri[53:106]
	NUV_gri_149 = NUV_gri[106:159]
	NUV_gri_126 = NUV_gri[159:212]
	NUV_gri_096 = NUV_gri[212:265]
	NUV_gri_066 = NUV_gri[265:318]
	NUV_gri_035 = NUV_gri[318:371]
	NUV_gri_025 = NUV_gri[371:424]
	NUV_gri_006 = NUV_gri[424:477]
	NUV_gri_015 = NUV_gri[477:530]
	NUV_gri_026 = NUV_gri[530:583]
	NUV_gri_040 = NUV_gri[583:636]

	# COLOR gri-3.6micras
	gri_36m = r_mag - spitz36m_mag

	gri_36m_227 = gri_36m[0:53]
	gri_36m_179 = gri_36m[53:106]
	gri_36m_149 = gri_36m[106:159]
	gri_36m_126 = gri_36m[159:212]
	gri_36m_096 = gri_36m[212:265]
	gri_36m_066 = gri_36m[265:318]
	gri_36m_035 = gri_36m[318:371]
	gri_36m_025 = gri_36m[371:424]
	gri_36m_006 = gri_36m[424:477]
	gri_36m_015 = gri_36m[477:530]
	gri_36m_026 = gri_36m[530:583]
	gri_36m_040 = gri_36m[583:636]


	# ***************************************************************
	# MODELO DE POBLACIONES ESTELARES BRUZUAL Y CHARLOT 2003

	mag_model = np.load('/Users/cris/Tesis/Mod_BC03_SSP/models_conv.npy')
	mag_model_NUV = mag_model[:,:,1]
	mag_model_r = mag_model[:,:,4]
	mag_model_36m = mag_model[:,:,10]

	# COLOR NUV-gri
	NUV_gri_BC03 = mag_model_NUV - mag_model_r
	
	NUV_gri_BC03_224 = NUV_gri_BC03[0,:]
	NUV_gri_BC03_164 = NUV_gri_BC03[1,:]
	NUV_gri_BC03_063 = NUV_gri_BC03[2,:]
	NUV_gri_BC03_033 = NUV_gri_BC03[3,:]
	NUV_gri_BC03_009 = NUV_gri_BC03[4,:]
	NUV_gri_BC03_055 = NUV_gri_BC03[5,:]

	
	
	# COLOR gri-3.6micras
	gri_36m_BC03 = mag_model_r - mag_model_36m

	gri_36m_BC03_224 = gri_36m_BC03[0,:]
	gri_36m_BC03_164 = gri_36m_BC03[1,:]
	gri_36m_BC03_063 = gri_36m_BC03[2,:]
	gri_36m_BC03_033 = gri_36m_BC03[3,:]
	gri_36m_BC03_009 = gri_36m_BC03[4,:]
	gri_36m_BC03_055 = gri_36m_BC03[5,:]


	# ***************************************************************
	# COLOR - COLOR PLOTS NUV-gri vs gri-3.6micras
	# ***************************************************************

	#size = [8.+n*2. for n in  list((range(len(Central_gri_36m[0]))))]
	size = [8.+n**2.1 for n in  list((range(len(Central_gri_36m[0]))))]

	plt.subplot(gs[subPlot])
	plt.grid()
	
	# OVERPLOT DE LOS MODELOS DE VAZDEKIS	
	cmap = plt.get_cmap('Blues')
	colors = cmap(np.linspace(0.2, 1.0, 8))

	#plt.plot(g_r_SDSS_23,r_i_SDSS_23, 'o-', lw=2, color=colors[0], label='[Fe/H] -2.3')
	plt.plot(NUV_gri_227,gri_36m_227, 'o-', lw=2, color=colors[0], mec=colors[0], label='-2.27 [Fe/H]')
	plt.plot(NUV_gri_179,gri_36m_179, 'o-', lw=2, color=colors[1], mec=colors[1], label='-1.79 [Fe/H]')
	plt.plot(NUV_gri_126,gri_36m_126 , 'o-', lw=2, color=colors[2], mec=colors[2], label='-1.26 [Fe/H]')
	plt.plot(NUV_gri_096,gri_36m_096 , 'o-', lw=2, color=colors[3], mec=colors[3], label='-0.96 [Fe/H]')
	plt.plot(NUV_gri_025,gri_36m_025 , 'o-', lw=2, color=colors[4], mec=colors[4], label='-0.25 [Fe/H]')
	plt.plot(NUV_gri_006,gri_36m_006 , 'o-', lw=2, color=colors[5], mec=colors[5], label='0.06 [Fe/H]')
	plt.plot(NUV_gri_015,gri_36m_015 , 'o-', lw=2, color=colors[6], mec=colors[6], label='0.15 [Fe/H]')
	plt.plot(NUV_gri_040,gri_36m_040 , 'o-', lw=2, color=colors[7], mec=colors[7], label='0.40 [Fe/H]')

	# OVERPLOT DATOS
	plt.errorbar(Central_NUV_gri[0]+0.05,Central_gri_36m[0]-0.01,xerr=[Central_NUV_gri[1],Central_NUV_gri[2]], yerr=[Central_gri_36m[1],Central_gri_36m[2]], marker='.',linestyle='None', color='g')
	plt.scatter(Central_NUV_gri[0]+0.05,Central_gri_36m[0]-0.01, s=size, color='g', edgecolor='None', label = 'Mid-plane')

	plt.errorbar(kpc15_NUV_gri[0]+0.05,kpc15_gri_36m[0]-0.01,xerr=[kpc15_NUV_gri[1],kpc15_NUV_gri[2]], yerr=[kpc15_gri_36m[1],kpc15_gri_36m[2]], marker='.',linestyle='None', color='darkorange')
	plt.scatter(kpc15_NUV_gri[0]+0.05,kpc15_gri_36m[0]-0.01, s=size, color='darkorange', edgecolor='None', label = '1.5 kpc')

	plt.errorbar(kpc30_NUV_gri[0]+0.05,kpc30_gri_36m[0]-0.01,xerr=[kpc30_NUV_gri[1],kpc30_NUV_gri[2]], yerr=[kpc30_gri_36m[1],kpc30_gri_36m[2]], marker='.',linestyle='None', color='purple')
	plt.scatter(kpc30_NUV_gri[0]+0.05,kpc30_gri_36m[0]-0.01, s=size, color='purple', edgecolor='None', label = '3.0 kpc')



	# OVERPLOT DE LOS MODELOS DE BC03	
	cmap = plt.get_cmap('Greys')
	colors = cmap(np.linspace(0.15, 1.0, 6))

	#plt.plot(g_r_SDSS_23,r_i_SDSS_23, 'o-', lw=2, color=colors[0], label='[Fe/H] -2.3')
# 	plt.plot(NUV_gri_BC03_224,gri_36m_BC03_224 , 'o-', lw=2, color=colors[0], mec=colors[0], label='-2.24 [Fe/H]')
# 	plt.plot(NUV_gri_BC03_164,gri_36m_BC03_164 , 'o-', lw=2, color=colors[1], mec=colors[1], label='-1.64 [Fe/H]')
# 	plt.plot(NUV_gri_BC03_063,gri_36m_BC03_063 , 'o-', lw=2, color=colors[2], mec=colors[2], label='-0.63 [Fe/H]')
# 	plt.plot(NUV_gri_BC03_033,gri_36m_BC03_033 , 'o-', lw=2, color=colors[3], mec=colors[3], label='-0.33 [Fe/H]')
# 	plt.plot(NUV_gri_BC03_009,gri_36m_BC03_009 , 'o-', lw=2, color=colors[4], mec=colors[4], label='0.09 [Fe/H]')
# 	plt.plot(NUV_gri_BC03_055,gri_36m_BC03_055 , 'o-', lw=2, color=colors[5], mec=colors[5], label='0.55 [Fe/H]')

	
	if galaxy == 'NGC5907':	
		plt.plot(2.1, -0.5 ,' *', markersize=25, mec='k', mfc='lightgreen', mew=2)
		plt.plot(2.45, -0.2, ' *', markersize=25, mec='k', mfc='orange', mew=2)
		plt.plot(2.65, -0.15, ' *', markersize=25, mec='k', mfc='blueviolet', mew=2)

		#plt.arrow( 6.6, -0.55, -0.058*14, -0.022*14, fc="k", ec="k", head_width=0.15, head_length=0.3 ) 
		#plt.text(6.1, -0.9, 'A$_{\lambda}$')

	if galaxy == 'NGC4565':
		plt.plot(2.3, -0.25, ' *', markersize=25, mec='k', mfc='lightgreen', mew=2)
		plt.plot(2.8, -0.3, ' *', markersize=25, mec='k', mfc='orange', mew=2)
		plt.plot(2.6, -0.2, ' *', markersize=25, mec='k', mfc='blueviolet', mew=2)
				
				
	plt.plot(10, 10, ' *', markersize=25, mec='k', mfc='grey', mew=2, label='Truncation')

	
		#plt.arrow(7.5, -0.52, -0.084*10, -0.033*10, fc="k", ec="k", head_width=0.15, head_length=0.3 ) 
		#plt.text(7.1, -0.9, 'A$_{\lambda}$')
	
	
	
	# LABELS INDICANDO LA EDAD EN GYR
	plt.text(NUV_gri_179[20],gri_36m_179[20],'1 Gyr',color='k', fontweight='heavy', fontsize=18)
	plt.text(NUV_gri_006[20]+0.05,gri_36m_006[20]-0.01,'1 Gyr',color='k', weight='heavy', fontsize=18)
	#plt.text(NUV_gri_040[20],gri_36m_040[20]-0.01,'1 Gyr', color='k', weight='heavy', fontsize=18)

	plt.text(NUV_gri_179[34],gri_36m_179[34],'5 Gyr', color='k', weight='heavy', fontsize=18)
	plt.text(NUV_gri_006[34],gri_36m_006[34]-0.02,'5 Gyr', color='k', weight='heavy', fontsize=18)
	#plt.text(NUV_gri_040[34],gri_36m_040[34]-0.01,'5 Gyr', color='k', weight='heavy', fontsize=18)

	plt.text(NUV_gri_179[44],gri_36m_179[44],'10 Gyr', color='k', weight='heavy', fontsize=18)
	plt.text(NUV_gri_006[44],gri_36m_006[44]-0.02,'10 Gyr', color='k', weight='heavy', fontsize=18)
	#plt.text(NUV_gri_040[44],gri_36m_040[44]-0.008,'10 Gyr', color='k', weight='heavy', fontsize=18)
	
	subPlot += 1

plt.subplot(gs[0])
plt.ylabel('gri - 3.6$\mu m$')
plt.xlabel("NUV-gri")
plt.ylim(-1.3, 1.82)
plt.xlim(1.4, 7.75)
plt.title(''+galaxy[:3]+' 5907')	

plt.subplot(gs[1])
plt.ylim(-1.3, 1.82)
plt.xlim(1.4, 7.85)
plt.xlabel("NUV-gri")
plt.title(''+galaxy[:3]+' 4565')	

plt.subplot(gs[2])
plt.xlabel("NUV-gri")
plt.ylim(-1.3, 1.82)
plt.xlim(1.4, 7.75)
plt.title(''+galaxy[:3]+' 4472')

plt.legend(loc=(1.02, 0.38),fontsize=15, numpoints=1)
plt.tight_layout() 
plt.subplots_adjust( wspace=0.22, hspace=0.15,left=0.1, right=0.79)

#plt.savefig('/Users/cris/Tesis/paper_truncationsNGC4565_5907/Fig_ColourColourVazdk.eps')
#plt.savefig('/Users/cris/Tesis/paper_truncationsNGC4565_5907/Fig_ColourColourVazdk.png', dpi=100)

