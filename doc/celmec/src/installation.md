# Installation

*Note*: Familiarity with managing rust project dependencies with `Cargo.toml` is assumed. To get the necessary info look [here](https://doc.rust-lang.org/book/ch01-03-hello-cargo.html) and [here](https://doc.rust-lang.org/book/ch02-00-guessing-game-tutorial.html) in "The Rust Programming Language" -book.

There is not yet a crate made of celmec. To use it at the moment you need add this to your project's `Cargo.toml` under `[dependencies]`:

```
celmec = { git = "https://github.com/juuso22/celmec.git" }
```

Celmec uses `ndarray`'s `Array` structure for its vector input and output, so you'll almost inevitably have to add `ndarray` to the dependencies as well.
