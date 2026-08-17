# Computational Physics & Numerical Simulations

A comprehensive collection of numerical physics models, dynamical system simulations, and mathematical algorithms developed in C++ during my Computational Physics BSc studies at Eötvös Loránd University. 

While originating from coursework, this repository demonstrates my ability to implement complex mathematical models from scratch, solve differential equations numerically, and structure reusable code.

## Core Highlight: Custom Numerical Library
* **[`/cpl`](./cpl)** — A shared C++ utility library built to support various simulations. It includes implementations for Fast Fourier Transform, ODE integration, mathematical optimization, matrix/vector operations, interpolation, and random number generation.

## Simulation Projects

### Complex & Statistical Systems
* **[`/Molecular dynamics`](./Molecular_dynamics)** — Particle interaction simulations demonstrating time integration and system evolution, including OpenGL-based visual variants. We received codes for this by default.
* **[`/Cellular automatons`](./Cellular_automatons)** — Discrete-state grid evolution models illustrating rule-based emergent behaviors.
* **[`/Population dynamics`](./Population_dynamics)** — Time evolution modeling of populations under different assumptions and parameterizations.

### Dynamical Systems & Mechanics
* **[`/Solar system simulation`](./Solar_system_simulation)** — An N-body orbital simulation solving ordinary differential equations for planetary motion.
* **[`/Pendulum`](./Pendulum)** — Solvers for nonlinear ODE systems representing simple, coupled, and double pendulum dynamics. We received codes for this by default.
* **[`/String`](./String)** — Wave and string dynamics simulations for discretized physical systems.

### Foundations & Algorithms
* **[`/Numeric integration`](./Numeric_integration)** — Implementation and comparison of various numerical quadrature techniques.
* **[`/Vectors`](./Vectors)** — Foundational linear algebra operations, memory management, and OOP in C++.

## Tech Stack, Tools & Structure
* **Core Language:** C++
* **AI-Assisted Development:** Actively leveraged GitHub Copilot and Gemini Pro (3.0 and 3.1) for pair programming, rapid prototyping, and algorithmic debugging.
* **Visualization:** OpenGL (for real-time physics rendering), Python (`.ipynb` notebooks for exploratory data analysis of the simulation outputs).
* **Structure:** Each project contains its own source files (`.cpp`), alongside `datas/` (numeric outputs) and `figures/` (generated plots).

## How to Run
Projects are designed to be compiled individually. Navigate to any project directory and compile the target `.cpp` file using a standard C++ compiler (e.g., `g++`). For OpenGL variants (`*-gl.cpp`), ensure your environment has the necessary graphics dependencies linked.
