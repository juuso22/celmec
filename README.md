# Celmec

Celmec aims to be a rust library for celestial mechanics simulations.

Currently it contains the functionality to solve the time evolution of a two-body orbit in the reference frame of one of the bodies if either:

1. The masses of the bodies id known and the initial position and velocity of the other body relative to the "reference frame body" are known.
2. The orbital (Keplerian) elements of the other body are known.

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

1. Parallellize 2-body problem
2. Configurable number of Fourier series terms
3. Publish as crate
4. First order perturbations
5. Gauss's orbit determination method

## Sources

Karttunen, Hannu: Johdatus taivaanmekaniikkaan, Ursan julkaisuja 82, 2020
