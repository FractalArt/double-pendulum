# Double Pendulum

![double-pendulum](double_pendulum.png)

This repository contains the Python code for a double-pendulum simulation.

## Implementation

The second-order differential equations for the two angles are taken from [myphysicslab.com](https://www.myphysicslab.com/pendulum/double-pendulum-en.html).
They are solved using fourth-order Runge-Kutta.

## Usage

The help message can be shown as follows:

```sh
$ python3 double_pendulum.py -h
usage: double_pendulum.py [-h] [--m1 M1] [--m2 M2] [--r1 R1] [--r2 R2] [--v1 V1] [--v2 V2] [--i1 I1] [--i2 I2] [--dt DT] [-g GRAVITY] [--steps STEPS] [--pause PAUSE] [-d] [--dpi DPI]
                          [--marker-size MARKER_SIZE] [--mod MOD] [-p SHOW_PAST] [--output OUTPUT]

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

options:
  -h, --help            show this help message and exit
  --m1 M1               The value of the first mass [kg]
  --m2 M2               The value of the second mass [kg]
  --r1 R1               The length of the first pendulum rod [m]
  --r2 R2               The length of the second pendulum rod [m]
  --v1 V1               The first initial angular velocity [rad / s]
  --v2 V2               The second initial angular velocity [rad / s]
  --i1 I1               The initial value for the first angle
  --i2 I2               The initial value for the first angle
  --dt DT               The simulation time stamp in seconds
  -g GRAVITY, --gravity GRAVITY
                        The value for the gravitational constant [N * m^2 / kg^2]
  --steps STEPS, -s STEPS
                        The number of simulation steps
  --pause PAUSE         The pause between different renderings
  -d, --display         Display the simulation
  --dpi DPI             The image resolution in dpi.
  --marker-size MARKER_SIZE
                        The marker-size.
  --mod MOD, -m MOD     When the simulation is shown, only display every `mod` step
  -p SHOW_PAST, --show-past SHOW_PAST
                        Specify how many past positions should be shown in the tail
  --output OUTPUT, -o OUTPUT
                        The output path of the plot
```

Thus, a command to display the simulated pendulum with default physical parameters
and show a trail of the past 5000 positions of the bottom pendulum as well is given by:

```sh
$ python3 double_pendulum.py --display --show-past 5000 -o double_pendulum.png
```

The `-o` flag provides the name of the output file and stores the positions of the bottom 
pendulum path in that image.

If only the simulation shall be shown without producing an output image one can run the following
command where only every hundredth step (`-m 100`) is shown to improve the rendering speed:

```sh
$ python3 double_pendulum.py -d -p 5000 -s 20000 -m 100
```
