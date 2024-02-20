function notimplemented()
    trace = stacktrace()
    callermethod = trace[2].func
    error("$callermethod is not implemented.")
end
export notimplemented

function notimplemented(T::Type)
    trace = stacktrace()
    callermethod = trace[2].func
    error("$callermethod is not implemented for $T.")
end
