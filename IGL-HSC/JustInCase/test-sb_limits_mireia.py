from astropy.io      import fits, ascii
import numpy         as np
from astropy.stats   import sigma_clipped_stats
from astropy.table   import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#############################################################################################
filters     = ['g', 'r', 'i']
files       = ['../../masking/test_mask_bcg_g.fits', '../../masking/test_mask_bcg_i.fits']  #### Con las imagenes enmascaradas.
sb_limits   = []
sigma_sky   = []
median_sky  = []
pixscale    = 0.168

for ii in np.arange(0, len(files)):  ### Loop para todas la imagenes

      image, hdr   = fits.getdata(files[ii], header = True)
      imsize       = np.shape(image)
      ZP           = np.float(hdr['ZP'])   ### Pillo el ZP del header pq lo tengo en el header.

      ## La idea aqui es colocar  diferentes regiones por toda la imagen para que sea random. Es posible que no sea necesario pero asi me aseguro que si hay algo que no esta bien enmascarado cuenta menos
      len_sim      = 1000 # number of regions! N (el num q sea. )

      print('...Estimate the SB limits ' + filters[ii] + '...')
      x_centers   = np.random.randint(0, imsize[1], len_sim)  # crear N numero de coordenadas random en la imagen.
      y_centers   = np.random.randint(0, imsize[0], len_sim)
      box_size    = 5./pixscale  ### El tamaño de la caja, pero que vamos no importa mucho pq vas a hacerlo por pixel. Pero lo hice asi para comparar.
      counts_box  = []

      for jj in np.arange(0, len(x_centers)):

          box     = image[int(y_centers[jj] - box_size): int(y_centers[jj] + box_size), int(x_centers[jj] - box_size): int(x_centers[jj] + box_size)]
          counts_box.append(box.flatten())

      counts_box  = np.concatenate(counts_box).ravel()
      plot_name   = 'counts_' + filters[ii] + '.pdf'
      sky         = sigma_clipped_stats(counts_box, sigma = 3., maxiters = 10, cenfunc = np.nanmedian, stdfunc = np.ma.std)
      med_sky     = sky[0]


      print('...Histogram ' + filters[ii] + '...')

      plt.clf()
      plt.hist(counts_box[~np.isnan(counts_box)], bins = 100, range = [-25, 25])
      plt.xlabel(r'Counts')
      plt.axvline(x = med_sky, c = 'k')
      plt.savefig(plot_name, dpi = 200)
      plt.clf()

      sb_limits.append(-2.5*np.log10(3*sky[2]/(10/0.168)) + ZP + 5*np.log10(pixscale)) ### Esto basicamente convierte la sigma por pixel en la sigma para una cajita de 10"x10"
      sigma_sky.append(sky[2])
      median_sky.append(sky[0])


### A escribirlo en una tablita!

sb_table = Table(np.transpose(sb_limits), names = ['g', 'i'])
ascii.write(sb_table, 'sb_limits.txt', format = 'commented_header', overwrite = True)
ascii.write(np.array(median_sky), 'median_sky.txt', names = ['g', 'i'], format = 'commented_header', overwrite = True)
ascii.write(np.array(sigma_sky), 'sigma_sky.txt', names = ['g', 'i'], format = 'commented_header', overwrite = True)
