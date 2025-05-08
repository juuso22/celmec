# Celmec

Celmec aims to be a rust library for celestial mechanics simulations. Very much work in progress at the moment.

Currently it contains the functionality to solve the time evolution of a two-body sytem in a reference frame where one of the bodies orbits around another. It is also possible to give impulses to the "orbiting" body. See Documentation section below about how to achieve all this.

Current work item: first order perturbations.

## Usage pre-requisites

Having rust installed.

## Documentation

To access the documentation, do the following:

1. Install `mdbook`. (`cargo install mdbook`, I think)
2. Navigate to `docs/celmec`
3. `mdbook serve --open`
4. In your browser, go to `localhost:3000`

## TODO

(Not necessarily in the order they're going to be done)

1. Perturbations: equations for 2-body, some example perturbations, book content, tests
2. Patched-cone approximation
3. Make all demos work
4. Integrals
5. Clean up impulse, impulse tests, mass ejection book entry + gravity drag simulation
6. Setup website
7. Parallellize 2-body problem
8. Configurable number of Fourier series terms
9. Gauss's orbit determination method
10. Maybe publish as crate

## Sources

Karttunen, Hannu: Johdatus taivaanmekaniikkaan, Ursan julkaisuja 82, 2020
Sobotov: Orbital Mechanics
