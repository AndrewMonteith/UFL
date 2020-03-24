# Project Plan

Rough sketch of points to talk about in my report.
General Questions are highlighted in bold.

## Introduction 2-3 Pages

* Introduce workflow of a physical simulation ie: Give structure & Mathematical process (in equation), Simulate Over Time 
* Introduce symbolic computing as area that deals with representing and manipulating the equation
* Introduce notion of Domain Specific Language
  * Talk about different flavours and which best applies for symbolic computing
  * Introduce Unified Form Language as example of DSL relating to symbolic computing
* Simplification Passes on DSL's
  * Simple examples like common subexpression elimination
* Introduce Julia Programming Language
  * Cover basics of Type System & Multi Dispatch with small code snippets example
* Introduce goals of the project - Julia better implementation Language than Python for UFL
  * First class abstractions offered in Julia better than Python for making UFL?
  * More performant? - Relative to constructing and manipulating trees
    * Mention not implementing entire UFL surface API but recreating core parts

## Related Work 2-4 Pages     : Prob gonna be nearer 2 than 4 pages

* Firedrake Project
  * 2 Languages in project: C & Python.
    * Python for the expressive & ease for the developer 
      * Term rewriting and simplification is written in python
    * C for the Needs C for the speed of simulation
    * Single host&target language would make things easier
  * Uses UFL for DSL for DE
  * **Any other points to talk about here?**
  
* Current Julia Ecosystem
  * Often maths processes are DE's, this is Julia's current defacto package for representing and solving them
  * ModellingToolkit -> DSL to represent maths stuff
    * Compatible with all DE solvers
    * The metainformation we need for our solver would be hard to shoehorn into the ModelingToolkit
  * Simplify/Rewrite -> Transformation Passes
    * Could be used for our domain but would require changing our current datstructure to match their macro system

## Solution 4-7 Pages

* How I used Julia's type System
  * Using Abstract Types and multidispatch to generalise type behaviour
  * Traits as a maintainble way to describe behaviour on types, can add complexity for describing multiple 

* Julia Metaprogramming reducing boilerplate compared to pythons decorators
* Not including .T syntax because can't dispatch off data dependent types
* Design choices for hashing (eager vs lazy)
  * Added complexity by doing eager (Having to use @sig macro)
* Tree Traversal

## Results 2-3 Pages

**Which ones should I include**:

* Hasing Benchmarks
* Tree Construction
* Tree Traversal Benchmarks
  * Would I need define this pass here?
* Differentiation Benchmarks
  * Would I need define this pass here or earlier?
* Algebraic Lowering Benchmarks
* **Is this the right place to introduce the example form**
* Example Form Benchmarks

## Evaluation 1-2 Pages

* Why is this a bloody separate section?
* Maintability results
  * Small API surface due to traits
  * Reducer boilerplate because of metaprogramming means we can easily change general behaviour of type tree
  * **Any other easy points**

## Conclusion <= 1 Page

Julia is faster
Julia is less verbose
Julia has some good first class abstractions