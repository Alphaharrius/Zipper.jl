mutable struct ParallelSettings
    showmeter::Bool
    maxthreads::Integer
    taskmeterlock::ReentrantLock
end

global parallelsettings = ParallelSettings(true, 1, ReentrantLock())

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

    partitioned::Vector = []
    for partition in taskpartitions
        push!(partitioned, (runnable, partition))
    end
    paralleltasks.tasks = partitioned
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
