mutable struct ParallelSettings
    showmeter::Bool
    # TODO: Implement!
    maxthreads::Integer
    taskmeterlock::ReentrantLock
end

global parallelsettings = ParallelSettings(true, Threads.nthreads(), ReentrantLock())

mutable struct ParallelTasks
    name::String
    tasks
    meter
end

function showtaskmeter(bool::Bool)
    @warn("Task meter visibility is set to $bool")
    parallelsettings.showmeter = bool
end
export showtaskmeter

function paralleltasks(; name::String, tasks, batchsize::Integer)
    taskpartitions = Iterators.partition(tasks, batchsize)
    paralleltasks = ParallelTasks(name, undef, undef)

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

    count::Integer = 0
    partitioned::Vector = []
    for partition in taskpartitions
        count += length(partition)
        push!(partitioned, (runnable, partition))
    end
    paralleltasks.tasks = partitioned
    paralleltasks.meter = Progress(count, desc=name)
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
