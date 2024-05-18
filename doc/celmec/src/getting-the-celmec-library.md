# Getting the Celmec Library

`Celmec` is Rust library and can be made available by adding it to your Rust project's `Cargo.toml`. To get the necessary info about what `Cargo.toml` is read [this](https://doc.rust-lang.org/book/ch01-03-hello-cargo.html) and [this](https://doc.rust-lang.org/book/ch02-00-guessing-game-tutorial.html) section from the "The Rust Programming Language" -book.

To use `celmec` you need to add this to your project's `Cargo.toml` under `[dependencies]`:

```
celmec = { git = "https://github.com/juuso22/celmec.git" }
```

Celmec uses `ndarray` for vector algebra, so you also need to add that library to your `Cargo.toml`. The examples in this book point out which version of `ndarray` was used for them. I'll try to make it so that the examples would work with the latest version, but given that this a free-time project, I will be checking if newer version work only sporadically.
