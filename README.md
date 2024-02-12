# Zipper.jl
Implementation of Zipper Entanglement Renormalization on Julia platform.

## Using multithreading to enhance performance
`Zipper.jl` by default disabled support for multithreading, by setting up multithreading environment by `Zipper.setmaxthreads` it will utilize all threads available to the Julia environment to enhance the performance of many functions related to the generation or operation of type `CrystalFock` and `CrystalFockMap`.

### Parallel APIs
- `ParallelTasks` defines a structure that stores the information of a set of parallel tasks.
- `paralleltasks` is the main entry of using parallel computing in `Zipper.jl`, it creates an instance of `Zipper.ParallelTasks`.
- `parallel` is to execute the set of parallel tasks encapsulated by `Zipper.ParallelTasks`.
- `setmaxthreads` is to control how many threads Zipper will use during computations.
- `showtaskmeter` is to control whether the progress bar will appear when parallel computing is being performed.

Progress bars are available by default for any methods that uses the parallel computing APIs so that user can track and estimate the time required to complete the operation.

### Enable parallel computing in Visual Studio Code
If the Julia environment in Visual Studio Code does not provide multithreading capabilities, `Zipper.jl` cannot take advantage of parallel computing thus no performance enhancements will be made.

To allow Julia REPL in Visual Studio Code to use multiple threads, head to Julia settings in the `settings.json`, search for the setting `julia.additionalArgs` and add `--threads=auto` to the list if you want Julia to decide how many threads it will use in the REPL, which most likely being the number of logical processors your native computing environment provides, or `--threads=N` to manually set the `N` number of threads that Julia can use. Then adjust the number of threads Zipper can use by `setmaxthreads(N)`.

## Save/Load Zipper data objects
`Zipper.jl` supports saving/loading data objects from JSON files, so that data that requires long period of computation time can be loaded with ease. This API also support defining custom 
serialization rules for special data types. Noted that *normal* Julia types are also supported by `fioload` and `fiosave`. 

### Fio APIs
- `fiodir` Set the current project directory of the session, this is the directory which all saving/loading data will be performed on.
- `fiolower` Define a special lower function for the serialization of the object attributes, the lower function will 
be applied to the object before `JSON` serializes the object into a JSON string.
- `fioconstructor` Define a special constructor function for the deserialization of the object, the constructor function will 
be called after all object attributes are deserialized.
- `fioparser` Define a special parser function for the deserialization of the object attributes, the parser function will
be applied to the attribute value after the object is deserialized into a `Dict`.
- `fiostoragetype` Register a storage type for the given Julia type, the storage type will be serialized instead of the actual 
type when saving the object to a file. During deserialization the storage type will be deserialized and converted 
back to the actual Julia type. Please be noted that the storage type must be a subtype of `Element` and have 
the function `convert(::Type{StorageType}, ::JuliaType)` and `convert(::Type{JuliaType}, ::StorageType)` defined.
- `fiosave` Save the given `object` to a file with the given `name` in the current project directory defined in `fiodir()`, 
the file name will be in the format of `{name}.json`.
- `fioload` Load the object from the file with the given `name` from the current project directory defined in `fiodir()`, 
the file to be loaded will be in the path `{project directory}/{name}.json`. Depands on the type some object might 
require extra data file to be loaded.
