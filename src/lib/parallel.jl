mutable struct ParallelSettings
    showmeter::Bool
    maxthreads::Integer
    taskmeterlock::ReentrantLock
    mainthreadmeter::Union{UndefInitializer, ProgressUnknown}
end

global parallelsettings = ParallelSettings(
    true, # Show meter by default.
    1, # Zipper.jl will use 1 thread by default.
    ReentrantLock(),
    # I have to set this value to 128 since Threads.nthreads() at initialization of Julia env 
    # seems to be 1 for some reason...
    undef) # 128 threads should be enough for everyone.

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

function singlestageparalleldivideconquer(f::Function, iter, count::Integer)
    actualcorecount::Integer = max(1, min(getmaxthreads(), Threads.nthreads()))
    countmetric::Integer = count/actualcorecount|>ceil
    batchsize::Integer = countmetric > 1 ? count/actualcorecount|>ceil : 2
    itembatches = Iterators.partition(iter, batchsize)
    batchcount::Integer = count/batchsize|>ceil
    return batchcount, paralleltasks(
        name="paralleldivideconquer count=$count",
        tasks=(()->f(batch) for batch in itembatches),
        count=batchcount)|>parallel
end

function paralleldivideconquer(f::Function, iter; count::Integer=iter|>length)
    batchcount, current = singlestageparalleldivideconquer(f, iter, count)
    while batchcount > 2
        batchcount, current = singlestageparalleldivideconquer(f, current, batchcount)
    end
    return f(current)
end
export paralleldivideconquer
#\end
function watchprogress(; desc::String)
    if Threads.threadid() != 1 || !parallelsettings.showmeter
        return
    end
    parallelsettings.mainthreadmeter = ProgressUnknown(desc=desc, spinner=true)
    return
end
export watchprogress

function updateprogress()
    if Threads.threadid() != 1 || !parallelsettings.showmeter
        return
    end
    try
        progress = parallelsettings.mainthreadmeter
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
    if Threads.threadid() != 1 || !parallelsettings.showmeter
        return
    end
    try
        progress = parallelsettings.mainthreadmeter
        if typeof(progress) == UndefInitializer
            # This branch is for the case where the current thread have called watchprogress before.
            return
        end
        ProgressMeter.finish!(progress)
    catch _ # UndefRefError
        # This branch is for the case where the current thread have not called watchprogress before.
    end
    return
end
export unwatchprogress
