mutable struct ParallelSettings
    showmeter::Bool
    maxthreads::Integer
    taskmeterlock::ReentrantLock
    threadmeters::Vector{ProgressUnknown}
end

global parallelsettings = ParallelSettings(
    true, # Show meter by default.
    1, # Zipper.jl will use 1 thread by default.
    ReentrantLock(),
    # I have to set this value to 128 since Threads.nthreads() at initialization of Julia env 
    # seems to be 1 for some reason...
    Vector(undef, 128)) # 128 threads should be enough for everyone.

mutable struct ParallelTasks
    tasks
    meter
end

function showtaskmeter(bool::Bool)
    @warn("Task meter visibility is set to $bool")
    parallelsettings.showmeter = bool
end
export showtaskmeter

function setmaxthreads(count::Integer)
    @warn("Max thread count is set to $count")
    parallelsettings.maxthreads = count
end
export setmaxthreads

getmaxthreads() = parallelsettings.maxthreads
export getmaxthreads

function paralleltasks(; name::String, tasks, count::Integer)
    actualcorecount::Integer = max(1, min(getmaxthreads(), Threads.nthreads()))
    batchsize::Integer = (count / actualcorecount)|>ceil
    taskpartitions = Iterators.partition(tasks, batchsize)
    paralleltasks = ParallelTasks(undef, Progress(count, desc="#threads($actualcorecount) $name"))
    @debug begin
        "paralleltasks: $name"
        "actualcorecount: $actualcorecount"
        "batchsize: $batchsize"
    end

    function runnable(tasks)::Vector
        rets = []
        for task in tasks
            ret = task()
            if typeof(ret) != Nothing
                push!(rets, ret)
            end
            if parallelsettings.showmeter
                # TODO: Might have race conditions.
                lock(()->next!(paralleltasks.meter), parallelsettings.taskmeterlock)
            end
        end
        return rets
    end

    watchprogress(desc="paralleltasks ($name)")
    partitioned::Vector = []
    for partition in taskpartitions
        push!(partitioned, (runnable, partition))
        updateprogress()
    end
    paralleltasks.tasks = partitioned
    unwatchprogress()
    return paralleltasks
end

function parallel(tasks::ParallelTasks)
    # Wait for the previous task meter to be cleared.
    lock(()->(), parallelsettings.taskmeterlock)
    threads = map(tasks.tasks) do (runnable, partition)
        Threads.@spawn runnable(partition)
    end
    buckets = fetch.(threads)
    finish!(tasks.meter)
    return (item for bucket in buckets for item in bucket)
end

function watchprogress(; desc::String)
    if !parallelsettings.showmeter
        return
    end
    progress = ProgressUnknown(desc=desc, spinner=true)
    parallelsettings.threadmeters[Threads.threadid()] = progress
    return
end
export watchprogress

function updateprogress()
    if !parallelsettings.showmeter
        return
    end
    try
        tid = Threads.threadid() # Always give 1 -> N
        progress = parallelsettings.threadmeters[tid]
        if typeof(progress) == UndefInitializer
            # This branch is for the case where the current thread have called watchprogress before.
            return
        end
        ProgressMeter.update!(progress)
    catch _ # UndefRefError
        # This branch is for the case where the current thread have not called watchprogress before.
    end
    return
end
export updateprogress

function unwatchprogress()
    if !parallelsettings.showmeter
        return
    end
    try
        tid = Threads.threadid() # Always give 1 -> N
        progress = parallelsettings.threadmeters[tid]
        if typeof(progress) == UndefInitializer
            # This branch is for the case where the current thread have called watchprogress before.
            return
        end
        ProgressMeter.finish!(progress)
        parallelsettings.threadmeters[tid] = undef
    catch _ # UndefRefError
        # This branch is for the case where the current thread have not called watchprogress before.
    end
    return
end
export unwatchprogress
