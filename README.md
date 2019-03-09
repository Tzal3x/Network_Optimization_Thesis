# Network_Optimization_Thesis
This is my work regarding my thesis at Computer Science M.Sc. of Athen's University of Economics and Business.
It is important to record every step of progress done in the span of those months.
Therefore, this README contains two-in-one: a record of progress and a description of my work's contents.

## First part: Understanding how optimization works.
Before diving into a specific subject, it is required to understand the basic concepts of network optimization theory and to experiment with some hypothetical/artificial problem(s).

### Satelite-Stations problem:
Let N satelites orbiting earth and M stations on the ground.
xs denotes rate of information transmited by satelite s. +/-
for every satelite: maximize U(Si) \\
\t\t\t\t\t\t\t\t\t\t subject to: Kirchhoff law \\
\t\t\t\t\t\t\t\t\t\t |xs| <= Capacity \\
\t\t\t\t\t\t\t\t\t\t ++some_other_constraints \\
                                
Formulate and solve this problem using MATLAB's optimization package. \\ 
