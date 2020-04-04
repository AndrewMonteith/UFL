## Introduction 2-3 Pages

* Workflow of a physical simulation ie: Give structure & Mathematical process (in equation), Simulate Over Time
* Symbolic Languages / Domain Specific Languages
  * Different flavours and what applied to symbolic languages
* Term Rewriting -> Series of transformation passes
  * Give some example of transformation passes
* Unified Form Language in a form compiler setting
* Julia Programming Language
  * Type System & Mutli-Dispatch & Metaprogramming with small code snippets
* Project Goals:
  * First class abstractions offered in Julia better than Python?
  * More performant? - Relative to constructing and manipulating trees
    * Mention not implementing entire UFL surface API but recreating core parts

## Related Work

SciPy / Haskell / Racket

* Current Julia Ecosystem
  * Often maths processes are DE's, this is Julia's current defacto package for representing and solving them
  * ModellingToolkit -> DSL to represent maths stuff
    * Compatible with all DE solvers
    * The metainformation we need for our solver would be hard to shoehorn into the ModelingToolkit
  * Simplify/Rewrite -> Transformation Passes
    * Could be used for our domain but would require changing our current datstructure to match their macro system

## Solution

* How I used Julia's type System
  * Using Abstract Types and multidispatch to generalise type behaviour
  * Traits as a maintainble way to describe behaviour on types, can add complexity for describing multiple
  * Trees need be homogenous types else we get a combinator type explosion

* Julia Metaprogramming reducing boilerplate compared to pythons decorators
* Not including .T syntax because can't dispatch off data dependent types
* Design choices for hashing (eager vs lazy)
  * Added complexity by doing eager (Having to use @sig macro)
* Tree Traversal
  * Treating trees a DAGs to avoid blowup in size

## Results

* Discuss how we create input instances to benchmark on
  * Number of nodes in the tree.
  * Number of types in the tree (Which is an indirect measure of number of dispatches)

* Define and Benchmark
  * Hashing Benchmark
  * Tree Construction
  * Tree Traversal
  * Algebraic Lowering
  * Differentiation
  * Example Form Benchmark
    * Introduce example form here

## Evaluation

* Why is this a bloody separate section?
* Maintability results
  * Small API surface due to traits
  * Reducer boilerplate because of metaprogramming means we can easily change general behaviour of type tree
  * Discuss how effective the benchmarking is
    * How input instances reflect true forms
    * Whether results show only great gain for obnoxiously large forms 

## Conclusion

Julia is faster
Julia is less verbose
Julia has some good first class abstractions
