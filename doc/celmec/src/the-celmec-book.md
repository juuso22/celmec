# The Celmec Book - or learning elementary celectial mechanics via coding

This a book that want's to be 2 things at once:

1. An introduction to elementary celestial mechanics.
2. An guide to the `celmec` Rust library for celestial mechanics.

I'll try to write it so that no pre-requisite knowledge of celestial mechanics is necessary. Some knowledge of basic programming, Rust, Python (especially `matplotlib`) and vector algebra is needed, though.

After instructions on how to make the `celmec` library available for you, each chapter introduces a little simulation with visualization with Python. The chapters start with a section  about some necessary physics to understand the simulation, the simulation itself as well as its results and then some additional physics which can be seen working in the simulation results. The physics sections will mainly provide ready formulas that other people have proved/calculated from lower principles, but for those interested, there will be references to other texts.

Currently `celmec` covers only the *two-body problem* ie. how to bodies move with gravity acting between them, but I'll attempt to extend the library (and then this book) in the future.
