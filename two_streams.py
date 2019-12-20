import argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation

parser = argparse.ArgumentParser(description='Analysis script for the two stream simulation')
group = parser.add_mutually_exclusive_group()
group.add_argument('--phase', action='store_true',
                    help='analize data from phase-space')
group.add_argument('--density', action='store_true',
                    help='analize data from density-potential')
group.add_argument('--energy', action='store_true',
                    help='analize data from energy (vs time)')
group.add_argument('--build_time', action='store_true',
                    help='plot instability building time vs stream velocity')
parser.add_argument('-v', '--velocity', type=int,
                    help='indicate velocity to analyze')
parser.add_argument('-s', '--save', action='store_true',
                    help='save output image/video')
parser.add_argument('--dt', type=float, default=0.01,
                    help='timestep of evolution (default: 0.01)')
parser.add_argument('--notime', action='store_true',
                    help='prevent time tag from displaying')
parser.add_argument('--autoscale', action='store_true',
                    help='set autoscale for graphs')
group_display = parser.add_mutually_exclusive_group()
group_display.add_argument('--snapshot', type=int, default=[], nargs='*',
                    help='save snapshots at given frames')
group_display.add_argument('--video_full', action='store_true',
                    help='save full video of evolution')
group_display.add_argument('--video', type=int, default=[], nargs=2,
                    help='save full video of evolution between given frames')

args = parser.parse_args()

dt = args.dt
folder = "data/two_streams/v" + str(args.velocity) + "/"

xs = np.loadtxt(folder + "fases_x.txt")[:, 1:]
vs = np.loadtxt(folder + "fases_v.txt")[:, 1:]
N = len(xs[0,:])

M = 1024
cells = np.arange(0, M)*N/M
density = np.loadtxt(folder + "densidad.txt")
time = density[:, 0]*dt
steps = len(time)

time_template = r"Tiempo: %.2f$\tau$"

if args.phase:
  fig_phase = plt.figure(figsize=(12.4, 12.4))
  ax = fig_phase.add_subplot(111, autoscale_on=False)
  ax.grid()
  ax.set_xlabel(r"x [$l$]", fontsize=20)
  ax.set_ylabel(r"v [$l/\tau$]", fontsize=20)
  ax.set_xlim(np.min(xs), np.max(xs))
  ax.set_ylim(np.min(vs), np.max(vs))
  line_right, = ax.plot([], [], 'b.', lw=1, markersize=1)
  line_left, = ax.plot([], [], 'r.', lw=1, markersize=1)
  for tick in ax.xaxis.get_major_ticks():
      tick.label.set_fontsize(16)
  for tick in ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(16)
  time_text_phase = ax.text(N//2, np.max(vs)-(np.max(vs)-np.min(vs))*0.2, '', fontsize=40, ha="center")

  def animate_phase(i):
    line_right.set_data(xs[i, 0:N:2], vs[i, 0:N:2])
    line_left.set_data(xs[i, 1:N:2], vs[i, 1:N:2])
    if args.notime:
      time_text_phase.set_text('')
    else:
      time_text_phase.set_text(time_template %time[i])
    return line_right, line_left, time_text_phase

  if args.snapshot != []:
    for s in args.snapshot:
      animate_phase(s)
      if args.save:
        plt.savefig("results/two_streams/phase_v" + str(args.velocity) + "_" + str(s) + ".png")
  else:
    if args.video_full:
      f_init = 0
      f_final = steps-1
    else:
      f_init = args.video[0] - 1
      f_final = args.video[1]
    ani = animation.FuncAnimation(fig_phase, animate_phase, range(f_init, f_final), blit=True, interval=100)#, repeat=True)
    if args.save:
      ani.save("results/two_streams/phase_v" + str(args.velocity) + ".mp4", writer="ffmpeg", fps=10)

elif args.density:
  density = density[:, 1:]
  potential = np.loadtxt(folder + "potencial.txt")[:, 1:]

  max_dens = np.max(density)
  min_dens = np.min(density)
  t = np.where(np.min(density, axis=1) == min_dens)[0][0]
  print("Tiempo de máximo bunching: ", time[t], r"$\tau$ (frame %d)" %t)

  fig = plt.figure(figsize=(23.85, 12.4))
  if args.autoscale:
    ax_dens = fig.add_subplot(211, autoscale_on=True)
  else:
    ax_dens = fig.add_subplot(211, autoscale_on=False)
  ax_dens.grid()
  ax_dens.set_xlabel(r"x [$l$]", fontsize=20)
  ax_dens.set_ylabel(r"Densidad [$V_o/l^3$]", fontsize=20)
  ax_dens.set_xlim(0, N-1)
  if not args.autoscale:
    ax_dens.set_ylim(min_dens, max_dens)
  line_dens, = ax_dens.plot([], [], 'g-', lw=1.5)
  for tick in ax_dens.xaxis.get_major_ticks():
      tick.label.set_fontsize(16)
  for tick in ax_dens.yaxis.get_major_ticks():
      tick.label.set_fontsize(16)

  if args.autoscale:
    ax_pot = fig.add_subplot(212, autoscale_on=True)
  else:
    ax_pot = fig.add_subplot(212, autoscale_on=False)
  ax_pot.grid()
  ax_pot.set_xlabel(r"x [$l$]", fontsize=20)
  ax_pot.set_ylabel(r"Potencial [$V_o$]", fontsize=20)
  ax_pot.set_xlim(0, N-1)
  min_pot = np.min(potential)
  max_pot = np.max(potential)
  if not args.autoscale:
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
    if args.autoscale:
      ax_dens.set_ylim(np.min(density[i, :]), np.max(density[i, :]))
      new_min_pot = np.min(potential[i,:])
      new_max_pot = np.max(potential[i,:])
      ax_pot.set_ylim(new_min_pot, new_max_pot)
      time_text.set_position((N//2, new_max_pot - (new_max_pot-new_min_pot)*0.2))
    if args.notime:
      time_text.set_text('')
    else:
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
  plt.ylabel(r"Energía [$ml^2/\tau^2$]")
  plt.legend(["Cinética", "Potencial"])
  t = np.where(epot == np.max(epot))[0][0]
  print(time[t])
  if args.save:
    plt.savefig("results/two_streams/energy_v" + str(args.velocity) + ".png")

elif args.build_time:
  N_vs = np.array([10, 20, 50, 100])
  build_times = []
  for v in range(1, 5):
    folder = "data/two_streams/v" + str(v) + "/"
    density = np.loadtxt(folder + "densidad.txt")
    time = density[:, 0]*dt
    steps = len(time)

    density = density[:, 1:]
    potential = np.loadtxt(folder + "potencial.txt")[:, 1:]

    max_dens = np.max(density)
    min_dens = np.min(density)
    t = np.where(np.min(density, axis=1) == min_dens)[0][0]
    build_times.append(time[t])
  build_times = np.array(build_times)
  #params_fit = np.polyfit(np.log(N_vs), np.log(build_times), 1)
  params_fit = np.polyfit(N_vs, build_times, 1)
  plt.figure()
  plt.plot(N_vs, build_times, "bo--")
  #plt.plot(N_vs, np.exp(np.polyval(params_fit, np.log(N_vs))), "b-")
  plt.plot(N_vs, np.polyval(params_fit, N_vs), "b-")
plt.show()
