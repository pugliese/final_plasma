import sys
import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as sco

if sys.argv[1] == "video":
  k = int(sys.argv[2])
  _dt = float(sys.argv[3])

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
elif sys.argv[1] == "freqs":
  def calculate_frequency_fft(signal, dt):
    N = len(signal)
    signal_fft = np.abs(np.real(np.fft.fft(signal)))[0:(N+1)//2]
    F = np.fft.fftfreq(N, dt)[0:(N+1)//2]
    #plt.plot(F, signal_fft, ".--")
    #plt.show()
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

  _dt = float(sys.argv[2])
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
  plt.plot(ks, np.pi*freqs, "b.--")
  plt.plot(ks, np.sinc(.5*np.pi*ks/1024)**2, "g-")
  plt.xlabel(r"Numero de onda ($k L/2\pi$)")
  plt.ylabel(r"$\omega$ [$\tau^{-1}$]")
  plt.plot(ks, ks/ks, "k--")
  #plt.plot(ks, 0.5*ks/ks, "k--")
  plt.xlim(0, k_max)
  plt.ylim(np.min(freqs)*3, np.max(freqs)*3.2)
  plt.show()
elif sys.argv[1] == "comp":
  field = np.loadtxt("data/perturbacion/campo_final_k=1.txt")[0, 1:]
  pot = np.loadtxt("data/perturbacion/potencial_final_k=1.txt")[0, 1:]
  dens = np.loadtxt("data/perturbacion/densidad_final_k=1.txt")[0, 1:]
  M = len(field)
  plt.figure(1)
  plt.plot(np.linspace(0,1, M), field)
  plt.title("Campo")
  plt.figure(2)
  plt.plot(np.linspace(0,1, M), pot)
  plt.title("Potencial")
  plt.figure(3)
  plt.plot(np.linspace(0,1, M), dens, "b-")
  dens_lap = (pot[:-2]+pot[2:]-2*pot[1:-1])/64
  #plt.plot(np.linspace(0,1, M)[1:-1], -dens_lap, "r-")
  plt.title("Densidad")
elif sys.argv[1] == "initial":
  potential = np.loadtxt("data/initial/potencial.txt")
  density = np.loadtxt("data/initial/densidad.txt")[:, 1:]
  density_fft_imag_norm = np.loadtxt("data/initial/densidad_fft_imag_norm.txt")[:, 1:]
  density_fft_imag = np.loadtxt("data/initial/densidad_fft_imag.txt")[:, 1:]
  ks = potential[:, 0] + 1
  potential = potential[:, 1:]

  amps_pot = np.max(potential, axis=1) - np.min(potential, axis=1)
  params = np.polyfit(np.log(ks[50:]), np.log(amps_pot[50:]), 1)
  amps_pot_fit = np.exp(np.polyval(params, np.log(ks)))
  amps_dens = np.max(density, axis=1) - np.min(density, axis=1)
  params = np.polyfit(ks[50:], amps_dens[50:]/ks[50:], 0)
  amps_dens_fit = np.polyval(params, ks)*ks

  plt.figure()
  plt.loglog(ks, amps_pot, "b.")
  plt.loglog(ks, amps_pot_fit*(amps_pot_fit>0), "b-")
  plt.figure()
  plt.loglog(ks, amps_dens, "g.")
  plt.loglog(ks, amps_dens_fit, "g-")

  #plt.figure()
  for k in range(50):
    #plt.plot(density_fft_imag_norm[k, :], ".-")
    print(k, ")", np.min(density_fft_imag_norm[k, :]), np.min(density_fft_imag[k, :]), np.max(density_fft_imag_norm[k, :]), np.max(density_fft_imag[k, :]))
    #plt.plot(density_fft_imag[k, :], ".-")
  plt.show()
elif sys.argv[1] == "point_charge":
  ekin = np.loadtxt("data/point_charge/cinetica.txt")
  epot = ekin[:, 1]
  ekin = ekin[:, 0]
