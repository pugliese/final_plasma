import argparse
import numpy as np
import matplotlib.pylab as plt
#import matplotlib.animation as animation
import scipy.optimize as sco

parser = argparse.ArgumentParser(description='Analysis script for cold plasma oscilations')

group = parser.add_mutually_exclusive_group()
group.add_argument('-d', '--dispersion', action='store_true',
                    help='calculate the dispersion relation for k in [1, 200]')
group.add_argument('-v', '--video', action='store_true',
                    help='make video for given modes')
group.add_argument('-s', '--snapshots', type=int, default=0,
                    help='show n time-equidistant density-potential snapshots of given mode')
group.add_argument('-k', '--kinetic', action='store_true',
                    help='show kinetic energy oscilations')
parser.add_argument('--mode', type=int, default=1,
                    help='modes to analyze (default: 1)')
parser.add_argument('--dt', type=float, default=0.001,
                    help='filename of loaded model (default: 0.001)')

args = parser.parse_args()

N = 2048
_dt = args.dt
k = args.mode

def calculate_frequency_fft(signal, dt):
  Nn = len(signal)
  signal_fft = np.abs(np.real(np.fft.fft(signal)))[0:(Nn+1)//2]
  F = np.fft.fftfreq(Nn, dt)[0:(Nn+1)//2]
  max_idx = np.argmax(signal_fft[1:])
  return F[max_idx+1]

func_fit = lambda t, A, f, t0, B: A*np.sin(2*np.pi*f*(t-t0)) + B

def calculate_frequency_fit(signal, dt):
  fo = calculate_frequency_fft(signal, dt)
  Ao = np.max(signal)
  Bo = Ao
  t0o = 0
  Ts = np.arange(len(signal))*dt
  paramms, _ = sco.curve_fit(func_fit, Ts, signal, p0 = [Ao, fo, t0o, Bo], maxfev=10000)
  return paramms[1]

if args.video:
  def video(vector, folder, dt, prefix=""):
    pot_inf = np.min(vector[:, 1:])
    pot_sup = np.max(vector[:, 1:])
    frames = len(vector)
    for f in range(frames):
      print("\rProgreso: %.2f%%" %(100*f/frames))
      plt.figure()
      plt.plot(vector[f, 1:])
      plt.title(r"T = %.2f $\tau$" %(f*dt))
      plt.ylim(pot_inf, pot_sup)
      plt.savefig(folder+prefix+f"{f:09}.png")
      plt.close()

  data = np.loadtxt("data/perturbacion/potencial_k=%d.txt" %k)
  video(data, "data/perturbacion/pot_vid/", _dt, "k="+str(k)+"_")
  data = np.loadtxt("data/perturbacion/densidad_k=%d.txt" %k)
  video(data, "data/perturbacion/dens_vid/", _dt, "k="+str(k)+"_")

elif args.dispersion:
  file_format = "data/perturbacion/cinetica_k=%d.txt"
  k_max = 200
  ks = np.arange(1, k_max +1)
  freqs = np.zeros(k_max)
  amps = np.zeros(k_max)
  for k in ks:
    ekin = np.loadtxt(file_format %k)[:,0]
    amps[k-1] = np.max(ekin) - np.min(ekin)
    freqs[k-1] = calculate_frequency_fit(ekin, _dt)
  plt.figure()
  plt.plot(ks, np.pi*freqs, "b-")
  #plt.plot(ks, np.sinc(.5*np.pi*ks/1024)**2, "g-")
  plt.xlabel(r"Numero de onda ($k$)")
  plt.ylabel(r"$\omega$ [$\tau^{-1}$]")
  plt.plot(ks, ks/ks, "k--")
  plt.grid()
  plt.yticks(np.arange(0.8, 1.01, 0.05))
  plt.xlim(0, k_max)
  plt.ylim(0.8, np.max(freqs)*3.2)
  plt.show()
elif args.kinetic:
  energy = np.loadtxt("data/perturbacion/cinetica_k=%d.txt" %k)[:, 0]
  time = np.arange(len(energy))*_dt
  plt.plot(time, energy, "r-")
  plt.plot([np.pi, np.pi], [0, np.max(energy)*0.45], "k--")
  plt.plot([np.pi, np.pi], [np.max(energy)*0.55, np.max(energy)], "k--")
  plt.text(np.pi, np.max(energy)*0.5, r"$\pi$", fontsize=18, ha = "center", va = "center")
  plt.plot([2*np.pi, 2*np.pi], [0, np.max(energy)*0.45], "k--")
  plt.plot([2*np.pi, 2*np.pi], [np.max(energy)*0.55, np.max(energy)], "k--")
  plt.text(2*np.pi, np.max(energy)*0.5, r"$2\pi$", fontsize=18, ha = "center", va = "center")
  plt.plot([3*np.pi, 3*np.pi], [0, np.max(energy)*0.45], "k--")
  plt.plot([3*np.pi, 3*np.pi], [np.max(energy)*0.55, np.max(energy)], "k--")
  plt.text(3*np.pi, np.max(energy)*0.5, r"$3\pi$", fontsize=18, ha = "center", va = "center")
  plt.xlim(0, 10)
  plt.ylim(0, np.max(energy))
  plt.xlabel(r"Tiempo [$\tau$]")
  plt.ylabel(r"$E_{kin}$ [$ml^2\tau^2$]")
  plt.show()
elif args.snapshots > 0:
  ekin = np.loadtxt("data/perturbacion/cinetica_k=%d.txt" %k)[:, 0]
  T = 0.5/calculate_frequency_fit(ekin, _dt)

  density = np.loadtxt("data/perturbacion/densidad_k=%d.txt" %k)
  time = density[:, 0]*_dt
  density = density[:, 1:]
  potential = np.loadtxt("data/perturbacion/potencial_k=%d.txt" %k)[:, 1:]
  delta_t = time[1]-time[0]
  steps_per_cycle = int(T/delta_t)
  idxs = np.int32(np.linspace(0, steps_per_cycle, args.snapshots))
  #idxs = np.arange(0, steps_per_cycle, steps_per_cycle//args.snapshots)

  M = density.shape[1]
  cells = np.arange(0, M)*N/M
  legend = [r"t=%.2f$\tau$" %time[i] for i in idxs]
  plt.figure()
  plt.xlabel(r"x [$l$]")
  plt.ylabel(r"Densidad [$l^{-1}$] (x1000)")
  for idx in idxs:
    plt.plot(cells, density[idx, :]*1000)
  plt.xlim(0, N-1)
  dens_lim = np.max(np.abs(density[idxs, :]))*1000
  plt.ylim(-dens_lim, dens_lim)
  plt.grid()
  plt.legend(legend)

  plt.figure()
  plt.xlabel(r"x [$l$]")
  plt.ylabel(r"Potencial [$V_o$]")
  for idx in idxs:
    plt.plot(cells, potential[idx, :])
  plt.xlim(0, N-1)
  pot_lim = np.max(np.abs(potential[idxs, :]))
  plt.ylim(-pot_lim, pot_lim)
  plt.grid()
  plt.legend(legend)
