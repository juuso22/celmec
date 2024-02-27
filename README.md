# Celmec

Celmec aims to be a rust library for celestial mechanics simulations.

Currently it contains the functionality to solve the time evolution of a two-body orbit in the reference frame of one of the bodies if either:

1. The masses of the bodies id known and the initial position and velocity of the other body relative to the "reference frame body" are known.
2. The orbital (Keplerian) elements of the other body are known.

The code attempts to use the following abbreviations:

```
a - semi-major axis
ee - vector pointing towards the perihelion
e - eccentricity (euclidean norm of ee)
f - true anomaly
h - energy integral (freely translated from Finnish 'energiaintegraali' used in Johdatus Taivaanmekaniikkaan by Karttunen, 2020)
kk - angular momentum vector
k - angular momentum scalar (euclidean norm of kk)
m - mass
n - average angular velocity
rr - position vector
r - distance (euclidean norm of rr)
vv - velocity vector
v - velocity (euclidean norm of vv)
t - time
iota - inclination
tau - perihelion time
omega - argument of periapsis
```

Note in the above, that when same letter is usually used for a scalar and a vector, the letter is doubled for the vector in the code.

## TODO

1. Make open source and publish as crate
2. Parallellize 2-body problem
3. Configurable number of Fourier series terms and iterative Kepler equation solving for true anomaly
4. First order perturbations
5. Gauss's orbit determination method

## Sources

Karttunen, Hannu: Johdatus taivaanmekaniikkaan, Ursan julkaisuja 82, 2020
