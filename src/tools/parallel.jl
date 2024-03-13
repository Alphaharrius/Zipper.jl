# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆  Global settings definition ◆
mutable struct ParallelState
    showmeter::Bool
    maxthreads::Integer
    mainthreadmeter::Union{UndefInitializer, ProgressUnknown}
    divideconquermeter::Union{UndefInitializer, ProgressUnknown}
end

global parallelstate = ParallelState(
    true, # Show meter by default.
    1, # Zipper.jl will use 1 thread by default.
    # I have to set this value to 128 since Threads.nthreads() at initialization of Julia env 
    # seems to be 1 for some reason...
    undef, undef)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ ParallelTasks definition ◆
mutable struct ParallelTasks
    taskchannel::Channel
    resultchannel::Channel
    corecount::Integer
    meter
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Parallel computing configuration APIs ◆
function showtaskmeter(bool::Bool)
    @warn("Task meter visibility is set to $bool")
    parallelstate.showmeter = bool
end
export showtaskmeter

function setmaxthreads(count::Integer)
    @warn("Max thread count is set to $count")
    parallelstate.maxthreads = count
end
export setmaxthreads

getmaxthreads() = parallelstate.maxthreads
export getmaxthreads
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Parallel computing APIs ◆
function paralleltasks(;
    tasks, 
    name::String, count::Integer, 
    showmeter::Bool = parallelstate.showmeter, corecount::Integer = getmaxthreads())

    function producer(ch::Channel)
        for task in tasks
            put!(ch, task)
        end
    end

    taskchannel = Channel(producer)
    resultchannel = Channel(count)
    actualcorecount::Integer = max(1, min(corecount, Threads.nthreads()))
    meter = showmeter ? Progress(count, desc="#threads($actualcorecount) $name", dt=0.2) : undef

    return ParallelTasks(taskchannel, resultchannel, actualcorecount, meter)
end
export paralleltasks

function parallel(tasks::ParallelTasks)
    counter = Threads.Atomic{Int}(0)
    monitorthreadid = rand(1:tasks.corecount)

    function producer()
        results = Queue{Any}()
        for task in tasks.taskchannel
            enqueue!(results, task())
            Threads.atomic_add!(counter, 1)
        end
        return results
    end

    # To prevent checking the monitor thread condition everytime, we will introduce a specific 
    # producer method for the monitor thread so it can update the progress bar directly.
    function monitorproducer()
        results = Queue{Any}()
        for task in tasks.taskchannel
            enqueue!(results, task())
            Threads.atomic_add!(counter, 1)
            ProgressMeter.update!(tasks.meter, counter[])
        end
        return results
    end

    threads = [
        Threads.@spawn (tasks.meter != undef && i == monitorthreadid ? monitorproducer : producer)() 
        for i in 1:tasks.corecount]
    resultbatches = fetch.(threads)
    tasks.meter == undef || ProgressMeter.finish!(tasks.meter)
    return (v for batch in resultbatches for v in batch)
end
export parallel

function updatedivideconquer()
    if parallelstate.divideconquermeter == undef || !parallelstate.showmeter
        return
    end
    if Threads.threadid() != 1
        return
    end
    next!(parallelstate.divideconquermeter)
end
export updatedivideconquer

getconquerer(f::Function) = function(iter)
    ret = iter|>first
    for item in Iterators.drop(iter, 1)
        ret = f(ret, item)
        updatedivideconquer()
    end
    return ret
end
export getconquerer

function paralleldivideconquer(f::Function, iter, count::Integer, desc::String)
    usablecores::Integer = max(1, min(getmaxthreads(), Threads.nthreads()))
    countmetric::Integer = count/usablecores|>ceil
    batchsize::Integer = countmetric > 1 ? count/usablecores|>ceil : 2
    corecount::Integer = count/batchsize|>ceil

    iterchannel = (function (ch::Channel)
        for el in iter
            put!(ch, el)
        end
    end)|>Channel

    conquer() = [v for v in iterchannel]|>f
    
    parallelstate.divideconquermeter = (
        ProgressUnknown(desc="divideconquer $desc count=$count batch=$batchsize", spinner=true))
    threads = [Threads.@spawn conquer() for _ in 1:corecount]
    results = fetch.(threads)
    ProgressMeter.finish!(parallelstate.divideconquermeter)
    if corecount > 3
        return paralleldivideconquer(f, results, corecount, desc)
    end

    parallelstate.divideconquermeter = ProgressUnknown(
        desc="divideconquer $desc merge", spinner=true)
    results = f(results)
    ProgressMeter.finish!(parallelstate.divideconquermeter)
    return results
end

paralleldivideconquer(f::Function, iter; count::Integer=iter|>length, desc::String="any") = (
    paralleldivideconquer(f, iter, count, desc))

export paralleldivideconquer
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Progress meter APIs ◆
function watchprogress(; desc::String)
    if Threads.threadid() != 1 || !parallelstate.showmeter
        return
    end
    parallelstate.mainthreadmeter = ProgressUnknown(desc=desc, spinner=true)
    return
end
export watchprogress

function updateprogress()
    if Threads.threadid() != 1 || !parallelstate.showmeter
        return
    end
    try
        progress = parallelstate.mainthreadmeter
        if typeof(progress) == UndefInitializer
            # This branch is for the case where the current 
            # thread have called watchprogress before.
            return
        end
        ProgressMeter.next!(progress)
    catch _ # UndefRefError
        # This branch is for the case where the current thread 
        # have not called watchprogress before.
    end
    return
end
export updateprogress

function unwatchprogress()
    if Threads.threadid() != 1 || !parallelstate.showmeter
        return
    end
    try
        progress = parallelstate.mainthreadmeter
        if typeof(progress) == UndefInitializer
            # This branch is for the case where the current thread 
            # have called watchprogress before.
            return
        end
        ProgressMeter.finish!(progress)
    catch _ # UndefRefError
        # This branch is for the case where the current thread have 
        # not called watchprogress before.
    end
    return
end
export unwatchprogress
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
