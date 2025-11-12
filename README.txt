Introduction

GAQ is a small C++ project exploring geometric algebra computations in a way that’s meant to be readable, easy to debug, and conceptually clear.
It’s not optimized for performance — instead, it focuses on transparency, learning value, and universality, allowing the creation of arbitrary multivectors.
Internally, the project leverages standard C++ features such as std::string and std::map, which help manage memory efficiently and keep the codebase straightforward to understand and extend.
If you’re curious about how high dimensional geometric algebra can be represented and manipulated in code, this is a solid starting point.

Overview

This project was originally created as a learning exercise in C++ and geometric algebra.
The main objective was to make geometric algebra computations readable, debuggable, and intuitive, even at the cost of raw performance.
While the architecture is far from perfect, it works reliably for its intended purpose.

Structure and Usage
All calculations are currently written directly inside the main() function.
A few example usages are included to demonstrate the syntax and workflow.

In QCGA.h, you’ll find definitions for:
Basis elements 1, e1, e2, ..., e15
Null basis elements ei1, eo1, ..., ei6, eo6
These serve as the building blocks for constructing multivectors used in calculations.

⚠️ Note:
Some functions can throw exceptions if used with incompatible multivectors — a side effect of the early design choices. Users are expected to understand the structure of the algebra before experimenting.
