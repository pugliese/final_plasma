import argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation

parser = argparse.ArgumentParser(description='Analysis script for the two stream simulation')
group = parser.add_mutually_exclusive_group()
group.add_argument('--density', action='store_true',
                    help='analize data from density-potential')
group.add_argument('--energy', action='store_true',
                    help='analize data from energy (vs time)')
group.add_argument('--mean', action='store_true',
                    help='plot mean density and potential after termalization')
group.add_argument('--debye', action='store_true',
                    help='plot the "debye length" vs temperature')
parser.add_argument('-s', '--save', action='store_true',
                    help='save output image/video')
parser.add_argument('--dt', type=float, default=0.001,
                    help='timestep of evolution (default: 0.001)')
parser.add_argument('--T', type=float, default=25,
                    help='system temperature (default: 25)')
group_display = parser.add_mutually_exclusive_group()
group_display.add_argument('--snapshot', type=int, default=[], nargs='*',
                    help='save snapshots at given frames')
group_display.add_argument('--video_full', action='store_true',
                    help='save full video of evolution')
group_display.add_argument('--video', type=int, default=[], nargs=2,
                    help='save full video of evolution between given frames')

args = parser.parse_args()

dt = args.dt
folder = "data/point_charge/"

M = 256
N = M
cells = np.arange(0, M)*N/M
density = np.loadtxt(folder + "densidad.txt")
time = density[:, 0]*dt
steps = len(time)

time_template = r"Tiempo: %.2f$\tau$"

if args.density:
  density = density[:, 1:]
  potential = np.loadtxt(folder + "potencial.txt")[:, 1:]

  max_dens = np.max(density)
  min_dens = np.min(density)
  t = np.where(np.min(density, axis=1) == min_dens)[0][0]
  print("Tiempo de máximo bunching: ", time[t], r"$\tau$ (frame %d)" %t)

  fig = plt.figure(figsize=(23.85, 12.4))
  ax_dens = fig.add_subplot(211, autoscale_on=False)
  ax_dens.grid()
  ax_dens.set_xlabel(r"x [$l$]", fontsize=20)
  ax_dens.set_ylabel(r"Densidad [$V_o/l^3$]", fontsize=20)
  ax_dens.set_xlim(0, N-1)
  ax_dens.set_ylim(min_dens, max_dens)
  line_dens, = ax_dens.plot([], [], 'g-', lw=1.5)
  for tick in ax_dens.xaxis.get_major_ticks():
      tick.label.set_fontsize(16)
  for tick in ax_dens.yaxis.get_major_ticks():
      tick.label.set_fontsize(16)

  ax_pot = fig.add_subplot(212, autoscale_on=False)
  ax_pot.grid()
  ax_pot.set_xlabel(r"x [$l$]", fontsize=20)
  ax_pot.set_ylabel(r"Potencial [$V_o$]", fontsize=20)
  ax_pot.set_xlim(0, N-1)
  min_pot = np.min(potential)
  max_pot = np.max(potential)
  ax_pot.set_ylim(min_pot, max_pot)
  line_pot, = ax_pot.plot([], [], 'b-', lw=1.5)
  for tick in ax_pot.xaxis.get_major_ticks():
      tick.label.set_fontsize(16)
  for tick in ax_pot.yaxis.get_major_ticks():
      tick.label.set_fontsize(16)
  time_text = ax_pot.text(N//2, max_pot - (max_pot-min_pot)*0.2, '', fontsize=40, ha="center")

  def animate(i):
    line_dens.set_data(cells, density[i, :])
    line_pot.set_data(cells, potential[i, :])
    time_text.set_text(time_template %time[i])
    return line_dens, line_pot, time_text

  if args.snapshot != []:
    for s in args.snapshot:
      animate(s)
      if args.save:
        plt.savefig("results/two_streams/density_v" + str(args.velocity) + "_" + str(s) + ".png")
  else:
    if args.video_full:
      f_init = 0
      f_final = steps-1
    else:
      f_init = args.video[0]
      f_final = args.video[1]
    ani = animation.FuncAnimation(fig, animate, range(f_init, f_final), blit=True, interval=100)#, repeat=True)
    if args.save:
      ani.save("results/two_streams/density_v" + str(args.velocity) + ".mp4", writer="ffmpeg", fps=10)

elif args.energy:
  energy = np.loadtxt(folder + "energia.txt")
  ekin = energy[:, 0]
  epot = energy[:, 1]
  plt.figure()
  plt.plot(time, ekin, "r-")
  plt.plot(time, epot, "b-")
  plt.xlabel(r"Tiempo [$\tau$]")
  plt.ylabel(r"Energía [$l^2/\tau^2$]")
  plt.legend(["Cinética", "Potencial"])
  t = np.where(epot == np.max(epot))[0][0]
  print(time[t])
  if args.save:
    plt.savefig("results/two_streams/energy_v" + str(args.velocity) + ".png")

def width(data, a=0.2):
  mr = M//2
  ml = M//2
  while data[mr] > a:
    mr += 1
  while data[ml] > a:
    ml -= 1
  r = mr - (a - data[mr])/(data[mr-1] - data[mr])
  l = ml + (a - data[ml])/(data[ml+1] - data[ml])
  return l, r

if args.mean:
  mean_density = np.loadtxt(folder+"densidad_medio_T=%.1f.txt" %args.T)[1:]
  mean_potential = np.loadtxt(folder+"potencial_medio_T=%.1f.txt" %args.T)[1:]
  plt.plot(cells, mean_potential/np.max(mean_potential), "b-")
  plt.plot(cells, mean_density/np.max(mean_density), "g-")
  #plt.xlim(1000, 1050)
  a, b = width(mean_potential/np.max(mean_potential))
  print(a, b, 0.5*(b-a), 0.5*(a+b))
  plt.xlim(M//2-20, M//2+20)
  plt.grid()
  plt.xlabel(r"Posicion [$l$]")
  plt.legend(["Potencial", "Densidad"])

if args.debye:
  Ts = np.arange(4, 101, 3)
  Ds = []
  for T in Ts:
    mean_potential = np.loadtxt(folder+"potencial_medio_T=%.1f.txt" %T)[1:]
    a, b = width(mean_potential/np.max(mean_potential))
    Ds.append(0.5*(b-a))
  Ds = np.array(Ds)
  params_fit = np.polyfit(np.log(Ts), np.log(Ds), 1)
  Ts_fit = np.linspace(min(Ts), max(Ts), 10001)
  Ds_fit = np.exp(np.polyval(params_fit, np.log(Ts_fit)))
  """
  plt.plot(np.log(Ts), np.log(Ds), "ro")
  plt.plot(np.log(Ts_fit), np.log(Ds_fit), "r-")
  """
  plt.loglog(Ts, Ds, "ro")
  plt.loglog(Ts_fit, Ds_fit, "r-")
  plt.xlim(Ts[0], Ts[-1])
  plt.xticks([5, 7, 10, 20, 40, 70, 100], [5, 7, 10, 20, 40, 70, 100])
  plt.yticks([3, 4, 5, 6, 8, 10, 13], [3, 4, 5, 6, 8, 10, 13])
  print(params_fit)
  plt.xlabel(r"Temperatura [$ml^2/\tau^2$]")
  plt.ylabel(r"Longitud de decaimiento [$l$]")
  plt.grid()
  plt.show()
plt.show()
