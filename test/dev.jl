using Base.Iterators
using Zipper, ProgressMeter

setmaxthreads(Threads.nthreads())
function heavytask()
    sleep(0.01)
    return 1
end

tasks = paralleltasks(name="Test", tasks=(()->heavytask() for _ in 1:5000), count=5000)
ret = tasks|>parallel
vec = ret|>collect

using ProgressMeter
prog = ProgressUnknown(desc="Total length of characters read:")
total_length_characters = 0
for _ in 1:1000
    global total_length_characters += 1
    update!(prog, total_length_characters/1000)
    sleep(0.01)
end
