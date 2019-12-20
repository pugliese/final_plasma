import argparse
import numpy as np
import matplotlib.pylab as plt

parser = argparse.ArgumentParser(description='Analysis script for the two stream simulation')
group = parser.add_mutually_exclusive_group()
group.add_argument('--kinetic', action='store_true',
                    help='analize kinetic energy vs time')
group.add_argument('--mean', action='store_true',
                    help='plot mean density and potential after termalization')
group.add_argument('--debye', action='store_true',
                    help='plot the "debye length" vs temperature')
parser.add_argument('--zoom', action='store_true',
                    help='zoom potential plot into relevant region')
parser.add_argument('-dt', type=float, default=0.001,
                    help='timestep of evolution (default: 0.001)')
parser.add_argument('-T', type=float, default=10,
                    help='system temperature (default: 10)')
"""
group_display = parser.add_mutually_exclusive_group()
group_display.add_argument('--snapshot', type=int, default=[], nargs='*',
                    help='save snapshots at given frames')
group_display.add_argument('--video_full', action='store_true',
                    help='save full video of evolution')
group_display.add_argument('--video', type=int, default=[], nargs=2,
                    help='save full video of evolution between given frames')
"""
args = parser.parse_args()

dt = args.dt
folder = "data/point_charge/"

M = 256
N = M
cells = np.arange(0, M)*N/M
density = np.loadtxt(folder + "densidad.txt")
time = density[:, 0]*dt
steps = len(time)

def inter_cuartil(counts, bins):
  acum = np.zeros_like(bins)
  acum[0] = 0
  for i, c in enumerate(counts):
    acum[i+1] = acum[i] + c
  total = np.sum(counts)
  idx_1cuartil = np.where(np.abs(acum-total/4) == np.min(np.abs(acum-total/4)))[0][0]
  idx_3cuartil = np.where(np.abs(acum-3*total/4) == np.min(np.abs(acum-3*total/4)))[0][0]
  return (bins[idx_1cuartil]+bins[idx_1cuartil+1])/2, (bins[idx_3cuartil]+bins[idx_3cuartil+1])/2

def width(data, alpha=0.2):
  mr = M//2
  ml = M//2
  while data[mr] > alpha:
    mr += 1
  while data[ml] > alpha:
    ml -= 1
  r = mr - (alpha - data[mr])/(data[mr-1] - data[mr])
  l = ml + (alpha - data[ml])/(data[ml+1] - data[ml])
  return l, r

if args.mean:
  mean_density = np.loadtxt(folder+"densidad_medio_T=%.1f.txt" %args.T)[1:]
  mean_potential = np.loadtxt(folder+"potencial_medio_T=%.1f.txt" %args.T)[1:]
  plt.plot(cells, mean_potential/np.max(mean_potential), "b-")
  plt.plot(cells, mean_density/np.max(mean_density), "g-")
  a, b = width(mean_potential/np.max(mean_potential))
  print(a, b, 0.5*(b-a), 0.5*(a+b))
  if args.zoom:
    plt.xlim(M//2-30, M//2+30)
    plt.text(M//2-30, 0.82, r"$T=%.1f$" %args.T, fontsize=24)
  else:
    plt.xlim(0, M-1)
    plt.text(0, 0.82, r"$T=%.1f$" %args.T, fontsize=20)
  plt.ylim(min([np.min(mean_potential/np.max(mean_potential)), np.min(mean_density/np.max(mean_density))]), 1)  
  plt.grid()
  plt.xlabel(r"Posicion [$l$]")
  plt.legend([r"Potencial [$V_o$]", "Densidad [$l^{-1}$]"])

if args.debye:
  Ts = np.arange(4, 101, 3)
  Ds = []
  std_T = []
  intercuart_T = []
  for T in Ts:
    mean_potential = np.loadtxt(folder+"potencial_medio_T=%.1f.txt" %T)[1:]
    a, b = width(mean_potential/np.max(mean_potential))
    temps = 2*np.loadtxt(folder+"cinetica_T=%.1f.txt" %T)[:, 1]
    Ds.append(0.5*(b-a))
    temporal = np.histogram(temps, bins=30)
    cuartil1, cuartil3 = inter_cuartil(temporal[0], temporal[1])
    std_T.append(np.std(temps))
    intercuart_T.append(cuartil3-cuartil1)
  Ds = np.array(Ds)
  params_fit = np.polyfit(np.log(Ts), np.log(Ds), 1)
  Ts_fit = np.linspace(min(Ts), max(Ts), 10001)
  Ds_fit = np.exp(np.polyval(params_fit, np.log(Ts_fit)))
  """
  plt.plot(np.log(Ts), np.log(Ds), "ro")
  plt.plot(np.log(Ts_fit), np.log(Ds_fit), "r-")
  """
  plt.figure()
  plt.loglog(Ts, Ds, "ro")
  plt.loglog(Ts_fit, Ds_fit, "r-")
  plt.xlim(Ts[0], Ts[-1])
  plt.xticks([5, 7, 10, 20, 40, 70, 100], [5, 7, 10, 20, 40, 70, 100])
  plt.yticks([3, 4, 5, 6, 8, 10, 13], [3, 4, 5, 6, 8, 10, 13])
  print(params_fit)
  plt.xlabel(r"Temperatura [$ml^2/\tau^2$]")
  plt.ylabel(r"Longitud de decaimiento [$l$]")
  plt.grid()
  plt.figure()
  plt.plot(Ts, std_T, "ro")
  plt.plot(Ts, intercuart_T, "bo")
  plt.legend(["Desvio estandar", "Dist. intercuartil"])
  params_fit_std = np.polyfit(Ts, std_T, 1)
  params_fit_intercuart = np.polyfit(Ts, intercuart_T, 1)
  plt.plot(Ts, np.polyval(params_fit_std, Ts), "r-")
  plt.plot(Ts, np.polyval(params_fit_intercuart, Ts), "b-")
  plt.xlim(0, 100)
  plt.ylim(0, max([np.max(std_T), np.max(intercuart_T)]))
  plt.xlabel(r"Temperatura [$ml^2/\tau^2$]")
  plt.ylabel(r"Fluctuacion Temperatura [$ml^2/\tau^2$]")
if args.kinetic:
  ekin = np.loadtxt(folder + "cinetica_T=%.1f.txt" %args.T)
  time = ekin[:, 0]*dt
  epot = ekin[:, 2]
  ekin = ekin[:, 1]
  plt.figure()
  plt.plot(time, 2*ekin, "r-")
  plt.plot([time[0], time[-1]], [args.T, args.T], "k--")
  plt.xlim(time[0], time[-1])
  plt.xlabel(r"Tiempo [$\tau$]")
  plt.ylabel(r"Temperatura [$ml^2/\tau^2$]")
  plt.grid()
  plt.figure()
  temporal = plt.hist(2*ekin, bins=25)
  cuartil1, cuartil3 = inter_cuartil(temporal[0], temporal[1])
  plt.xlabel(r"Temperatura [$ml^2/\tau^2$]")
  plt.ylabel("Ocurrencia")
  height = np.max(temporal[0])*1.1
  plt.xlim(temporal[1][0], temporal[1][-1])
  plt.ylim(0, height)
  plt.plot(cuartil1*np.ones(2), [0, height], "k-")
  plt.plot(cuartil3*np.ones(2), [0, height], "k-")
  print("Media:", np.mean(2*ekin))
  print("Desvio estandar:", np.std(2*ekin))
  print("Distancia inter-cuartil:", cuartil3-cuartil1)
  print("Primer cuartil:", cuartil1)
  print("Tercer cuartil:", cuartil3)
plt.show()
