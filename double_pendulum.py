"""
Double Pendulum Simulation

This scripts simulates a double pendulum.
Physics-related paramters such as
- the two masses in the pendulum
- the two lengths of the two rods in the pendulum
- the initial angular velocities of the pendulum
- the gravitational constant

as well as paramters to tweak the simulation (see
help message) can be provided on the command-line.

Interesting parameters are `--display` and `--show-past`.
The first one not only generates a picture of the of the
trajectory of the bottom mass but also shows the animated
simulation while the second adds a trail of all past positions
of the bottom mass to the simulation.

An image of the positions visited by the bottom mass can be 
stored using the `--output` flag.
"""
import argparse as ap
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pathlib as pl

from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap

# set the plot background color
plt.rcParams["figure.facecolor"] = '#2e3440'


def parse_args():
    parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument('--m1', type=float, default=2, help='The value of the first mass [kg]')
    parser.add_argument('--m2', type=float, default=1, help='The value of the second mass [kg]')
    parser.add_argument('--r1', type=float, default=3, help='The length of the first pendulum rod [m]')
    parser.add_argument('--r2', type=float, default=2.5, help='The length of the second pendulum rod [m]')
    parser.add_argument('--v1', type=float, default=0, help='The first initial angular velocity [rad / s]')
    parser.add_argument('--v2', type=float, default=0, help='The second initial angular velocity [rad / s]')
    parser.add_argument('--i1', type=float, default=0.85 * np.pi, help='The initial value for the first angle')
    parser.add_argument('--i2', type=float, default=0.75 * np.pi, help='The initial value for the first angle')
    parser.add_argument('--dt', type=float, default=1e-3, help='The simulation time stamp in seconds')
    parser.add_argument('-g', '--gravity', type=float, default=9.81, help='The value for the gravitational constant [N * m^2 / kg^2]')
    parser.add_argument('--steps', '-s', type=int, default=100_000, help='The number of simulation steps')
    parser.add_argument('--pause', type=float, default=1e-15, help='The pause between different renderings')
    parser.add_argument('-d', '--display', action='store_true', help='Display the simulation')
    parser.add_argument('--dpi', type=int, default=350, help='The image resolution in dpi.')
    parser.add_argument('--marker-size', type=float, default=0.03, help="The marker-size.")
    parser.add_argument('--mod', '-m', type=int, default=100, help='When the simulation is shown, only display every `mod` step')
    parser.add_argument('-p', '--show-past', type=int, default=0, help='Specify how many past positions should be shown in the tail')
    parser.add_argument('--output', '-o', type=pl.Path, required=False, help='The output path of the plot')
    parser.add_argument('--ot1', type=pl.Path, required=False, help='Output a plot of the theta1 angle over time.')
    parser.add_argument('--ot2', type=pl.Path, required=False, help='Output a plot of the theta2 angle over time.')
    parser.add_argument('--delta-t2', type=float, required=False, help='Optional small offset of theta_2. Show how trajectories are sensible to initial conditions.')
    return parser.parse_args()


class FPrime:
    """Define the system of DEs.

    We express the system of second-order order DEs
    as a coupled system of first-order DEs of the form
   
    x' = f(x)

    This class defines the function f.
    """
    def __init__(self, m1, m2, r1, r2, g):
        """Pass the parameters of the system

        Arguments
        ---------
        m1: float
            The top mass
        m2: float
            The bottom mass
        r1: float
            The length of the top rod
        r2: float
            The length of the bottom rod
        g: float
            The gravitational constant
        """
        self._m1 = m1
        self._m2 = m2
        self._r1 = r1
        self._r2 = r2
        self._g = g
    
    def __call__(self, w):
        m1 = self._m1
        m2 = self._m2
        r1 = self._r1
        r2 = self._r2
        g = self._g

        dt = w[0] - w[1]
        den = 2 * m1 + m2 * (1 - np.cos(2*dt))
        m12 = m1 + m2

        return np.array([
                w[2],
                w[3],
                (-g*(2*m1 + m2)*np.sin(w[0]) - m2 * g * np.sin(w[0] - 2*w[1]) - 2 * np.sin(dt) * m2 * (w[3]**2 * r2 + w[2]**2 *r1 * np.cos(dt))) / (r1 * den),
                (2*np.sin(dt)* (w[2]**2 * r1 * m12 + g * m12 * np.cos(w[0]) + w[3]**2 * r2 * m2 * np.cos(dt) )) / (r2 * den),
                ])


class RungeKutta4:
    """The Runge Kutta solver"""
    def __init__(self, dt, f, initial):
        """Specify the solver parameters in the constructor

        Arguments
        ---------

        dt: float
            The time step
        initial: np.array
            The initial conditions
        """
        self._dt = dt
        self._f = f
        self._prev = initial
        self._iteration = 0
    
    def __call__(self):
        """Iterate the solver

        Returns
        -------
        np.array
            The next point in the iteration of the solver.
        """
        self._iteration += 1
        x = self._prev
        k1 = self._f(x) * self._dt
        k2 = self._f(x + 0.5 * k1) * self._dt
        k3 = self._f(x + 0.5 * k2) * self._dt
        k4 = self._f(x + k3) * self._dt
        self._prev = x + 1 / 6 * (k1 + 2*k2 + 2*k3 + k4)
        return self._prev

def angles_to_cartesian(t1, t2, r1, r2):
    """Convert the angles to cartesian coordinates.

    Arguments
    ---------
    t1: float
        The theta_1 angle
    t2: float
        The theta_2 angle
    r1: float
        The length of the top pendulum rod
    r2: float
        The length of the bottom pendulum rod

    Returns
    -------
    tuple
        A tuple containing the x and y coordinates of the two 
        pendulum masses
    """
    x1 = r1 * np.sin(t1)
    y1 = -r1 * np.cos(t1)
    x2 = x1 + r2 * np.sin(t2)
    y2 = y1 - r2 * np.cos(t2)
    return x1, y1, x2, y2

def store_angle_plot(path, time, angles):
    """Store an angular plot."""
    plt.rcParams["figure.facecolor"] = '#2e3440'
    fig, ax = plt.subplots(figsize=(16, 9), dpi=args.dpi)
    ax.scatter(time, angles, marker='o', c='#eef2f6', s=args.marker_size)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(path, bbox_inches='tight', transparent=False, dpi=args.dpi, pad_inches=0)


if __name__ == '__main__':

    # parse the command-line args
    args = parse_args()

    # define the system function
    f = FPrime(args.m1, args.m2, args.r1, args.r2, args.gravity)

    # initialize the solver
    rk = RungeKutta4(args.dt, f, np.array([args.i1, args.i2, args.v1, args.v2]))

    # if we simulate a second system with small offset to the first system
    # to demonstrate sensitivity to initial conditions
    if args.delta_t2:
        rk_delta = RungeKutta4(args.dt, f, np.array([args.i1, args.i2 + args.delta_t2, args.v1, args.v2]))
        cmap_delta = mpl.colormaps['copper']

    # create the colormap
    colors = ['#eef2f6', '#5e81ab', '#81a1c0', '#88c0d1', '#81a1c0']
    cmap = LinearSegmentedColormap.from_list('blueish', colors)

    # prepare the plot
    if args.display:
        plt.ion()
        fig, ax = plt.subplots()
        line, = ax.plot([], [], lw=2)

    # simulate
    time = 0
    xs = []
    ys = []
    ts = []
    if args.ot1:
        t1l = []
    if args.ot2:
        t2l = []
    if args.delta_t2:
        xsd = []
        ysd = []

    for step in range(args.steps):
        plt.clf()
        current = rk()
        t1, t2 = current[0], current[1]
        # print(f'Step {step}: {t1}, {t2}')

        x1, y1, x2, y2 = angles_to_cartesian(t1, t2, args.r1, args.r2) 
        if args.delta_t2:
            current_delta = rk_delta()
            t1d, t2d = current_delta[0], current_delta[1]
            x1d, y1d, x2d, y2d = angles_to_cartesian(t1d, t2d, args.r1, args.r2)

        xs.append(x2)
        ys.append(y2)
        ts.append(time)
        if args.ot1:
            t1l.append(t1)
        if args.ot2:
            t2l.append(t2)
        
        if args.delta_t2:
            xsd.append(x2d)
            ysd.append(y2d)

        if args.display and step % args.mod == 0:
            plt.axis('off')
            ax.relim()  # Recompute the data limits
            ax.autoscale_view()  # Rescale the view
            if args.delta_t2:
                plt.plot([0, x1d, x2d], [0, y1d, y2d], color='red', marker='o', markerfacecolor='red')
            plt.plot([0, x1, x2], [0, y1, y2], color='#00ff00', marker='o', markerfacecolor='grey')
            plt.xlim(-1 - args.r1 - args.r2, 1 + args.r1 + args.r2)  # Set the range of the x-axis
            plt.ylim(-1 - args.r1 - args.r2, 1 + args.r1 + args.r2)  # Set the range of the y-axis
            if args.show_past > 0:
                back = min(len(xs), args.show_past)
                if args.delta_t2:
                    plt.scatter(xsd[-back:], ysd[-back:], c=ts[-back:], cmap=cmap_delta, marker='o', s=args.marker_size)
                plt.scatter(xs[-back:], ys[-back:], c=ts[-back:], cmap=cmap, marker='o', s=args.marker_size)
            plt.draw()

            plt.pause(args.pause)  # Option

        time += args.dt

    if args.display:
        plt.ioff()  # Turn off interactive mode
        plt.close(fig)

    if args.output:
        plt.rcParams["figure.facecolor"] = '#2e3440'
        fig, ax = plt.subplots(figsize=(16, 9), dpi=args.dpi)
        ax.scatter(xs, ys, c=ts, cmap=cmap, marker='o', s=args.marker_size)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(args.output, bbox_inches='tight', transparent=False, dpi=args.dpi, pad_inches=0)

    # store theta1 plot
    if args.ot1:
        store_angle_plot(args.ot1, ts, t1l)

    # store theta2 plot
    if args.ot2:
        store_angle_plot(args.ot2, ts, t2l)
