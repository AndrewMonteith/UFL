# Benchmarking Results

We give julia times followed by python times

## Constructing Tree

* n = 100, 26.7μs / 142μs
* n = 200, 43.6μs / 239μs
* n = 300, 64.1μs / 368μs
* n = 500, 106.2μs / 607μs
* n = 1000, 200μs / 1.17ms
* n = 5000, 987.3μs /  5.86ms
* n = 10000, 1.96ms / 11.6ms
* n = 20000, 3.91ms / 23.5ms

## Traversal

### Pre-Order Traversal

* n = 100, 12.3μs / 19.8μs
* n = 200, 22.1μs / 40.4μs
* n = 300, 34.6μs / 63.4μs
* n = 500, 54.3μs / 112μs
* n = 1000, 107μs / 233μs
* n = 5000, 522.6μs / 1.18ms
* n = 10000, 1ms / 2.42ms
* n = 20000, 2.1ms / 4.87ms

### Post-Order Traversal

* n = 100, 18.8μs / 79.7μs
* n = 200, 35.4μs / 161μs
* n = 300, 50.8μs / 250μs
* n = 500, 88.2μs / 449μs
* n = 1000, 171μs / 915μs
* n = 5000, 820μs / 4.75ms
* n = 10000, 1.65ms / 9.61ms
* n = 20000, 3.26ms / 20.3ms

## Differentiation

* Form 1 - 458.7μs / 1.21ms
* Form 2 - 1.82ms / 3.51ms
* Form 3 - 1.85ms / 3.45ms
* Form 4 - 2.96ms / 5.05ms
* Form 5 - 2.98ms / 5.74ms

## Example Form

* Algebraic Lowering - 364μs / 680μs
* Differentiation 1 - 697μs / 1300μs
* Function Pullback - 609μs / 913μs
* Differentiation 2 - 775μs / 1250μs
* Total - 2445μs / 4143μs
