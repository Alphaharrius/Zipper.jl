# JuliaZipper
Implementation of Zipper Entanglement Renormalization on Julia platform.

## Using multithreading to enhance performance
`Zipper.jl` by default enables support for multithreading, it will utilize `Threads.nthreads()` number of threads available to the Julia environment to enhance the performance of many functions related to the generation or operation of type `CrystalFock` and `CrystalFockMap`.

### Parallel APIs
`Zipper.ParallelTasks` defines a structure that stores the information of a set of parallel tasks.
`Zipper.paralleltasks` is the main entry of using parallel computing in `Zipper.jl`, it creates an instance of `Zipper.ParallelTasks`.
`Zipper.parallel` is to execute the set of parallel tasks encapsulated by `Zipper.ParallelTasks`.
`Zipper.showtaskmeter` is to control whether the progress bar will appear when parallel computing is being performed.

Progress bars are available by default for any methods that uses the parallel computing APIs so that user can track and estimate the time required to complete the operation.

### Enable parallel computing in Visual Studio Code
If the Julia environment in Visual Studio Code does not provide multithreading capabilities, `Zipper.jl` cannot take advantage of parallel computing thus no performance enhancements will be made.

To allow Julia REPL in Visual Studio Code to use multiple threads, head to Julia settings in the `settings.json`, search for the setting `julia.additionalArgs` and add `--threads=auto` to the list if you want Julia to decide how many threads it will use in the REPL, which most likely being the number of logical processors your native computing environment provides, or `--threads=N` to manually set the `N` number of threads that Julia can use.
