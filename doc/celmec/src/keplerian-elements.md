# Keplerian Elements

Instead of starting with a known initial position and velocity of a planet it is also possible to have collected a collection of the properties of the orbit and to calculate the time-evolution of the planet along the orbit using those. One set of such properties are called the Keplerian elements and what they are is [described in Wikipedia](https://en.wikipedia.org/wiki/Orbital_elements#Keplerian_elements) with an especially helpful picture attached at the time of the writing.

In a nutshell, the Keplerian elements as used by `celmec` consist of quantities that describe:

1. How the plane of the orbit of the planet (or other object of interest) is positioned vis-à-vis the ecliptica ie. the plane on which the erath orbits the sun. These quantities are the angle between the orbital plane of the object being studied and ecliptica: inclination (`iota`) and the angle at which the orbit rises from below above to ecliptica with respect to a fixed direction (spring equinox): longitude of the ascending node (`longitude_of_the_ascending_node`).
2. How the orbit is shaped ie. its eccentricity `e` and semi-major axis `a`.
3. The angle from longitude of the ascending node to the perihelion along the orbital plane: argument of perihelion (`omega`)
4. A time at which the object passes its perihelion `tau`. Note that if you read the Wikipedia article, this differs slightly from the "true anomaly at t<sub>0</sub>": `tau` is equal to that t<sub>0</sub> when true anomaly is 0.

Next a little simulation using the Keplerian elements. And the reason I moved from talking about the orbit of a planet to the orbit of an object is because the simulation target won't be a planet, but a comet instead.
