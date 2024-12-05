# Stochastic series expansion algorithm for S=1/2 AFM Heisenberg model

This project reproduces the same algorithm as that in https://physics.bu.edu/~sandvik/programs/ssebasic/ssebasic.html

## Unsafe macro
Of course you can replace all the unsafe blocks in my code with the safe version. I have tested on my MacBook Pro (M1), and found that it would lead to less than 5% performance loss. I prefer the safe version actually, but I upload the unsafe version here for those who need ultimate performance.

## Why I use Rust?

First of all, I believe that old-fashioned languages like Fortran, which lack modern features, should fade out of history rather than continue to accumulate 'legacy' code. As for C++, I was once a C++ user, but the memory safety issues often drive me crazy. 

Many people recommend Julia to me. Although I greatly appreciate the language, the elegance and rigor of Rust are highly appealing to me, as it offers both safety and performance, and now I am a big fan of Rust. 

I am sharing here a SSE Rust code for reference, hoping to encourage more people engaged in physics algorithms and numerical research to switch to Rust programming. Unlike some uncommented code on GitHub, my code contains many comments. I believe that commenting is an important habit, and open-source contributions are also very important for physics research. I also welcome valuable suggestions for improving parts of my code.
