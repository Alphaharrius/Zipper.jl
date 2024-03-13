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
    parallelstate.divideconquermeter = (
        ProgressUnknown(desc="divideconquer $desc count=$count batch=$batchsize", spinner=true))
    result = tasks|>parallel|>collect
    ProgressMeter.finish!(parallelsettings.divideconquermeter)
    # Since the issue is still unknown after checking, we suspect that some variables 
    # passing through the iterative approach might be marked as garbage while it should 
    # not be. Therefore we will take a more functional approach (since Julia is functional) 
    # in hopes that the vm will lift this issue automatically.
    if batchcount > 3
        return paralleldivideconquer(f, result, batchcount, desc)
    end
    parallelstate.divideconquermeter = ProgressUnknown(
        desc="divideconquer $desc merge", spinner=true)
    result = f(result)
    ProgressMeter.finish!(parallelstate.divideconquermeter)
    return result
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
