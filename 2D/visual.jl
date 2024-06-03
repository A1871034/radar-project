#println("Loading GLMakie...")
#using GLMakie

println("Loading CairoMakie...")
using CairoMakie

println("Loading Sim Data...")
f = open(raw"C:\Users\Xander\Documents\Coding\radar-project\2D\sim.debug", "r")
data_header = readline(f)

prev = position(f)
seekend(f)

SIZE_X, SIZE_Y, PRECISION_BYTES = split(data_header, ",")
SIZE_X = parse(Int64, SIZE_X)
SIZE_Y = parse(Int64, SIZE_Y)
PRECISION_BYTES = parse(Int64, PRECISION_BYTES)
SIZE_Z = Int64((position(f) - prev)/(SIZE_X*SIZE_Y*PRECISION_BYTES))

seek(f, prev)
ez = Array{Float32, 3}(undef, SIZE_X, SIZE_Y, SIZE_Z)
read!(f, ez)
close(f)

println("Loading Height Data...")
f = open(raw"C:\Users\Xander\Documents\Coding\radar-project\2D\heights.data", "r")
seekend(f)
SIZE_X_HEIGHTS = Int32(position(f)/4)
seekstart(f)
terrain_heights = Array{Int32, 1}(undef, SIZE_X_HEIGHTS)
read!(f, terrain_heights)
close(f)

println("Initializing Animation...")

SCALE = 2000
function colorscale(val)
    if val == 0
        return 0
    end

    if val < 0
        return -log10(-val*SCALE + 1)/SCALE
    end

    return log10(val*SCALE + 1)/SCALE
end

time = Observable(1)

function get_function(time)
	fn(x,y) = ez[(SIZE_Y*SIZE_X) * @lift($time - 1) + SIZE_Y * (x-1) + (y-1) + 1]
	fn
end

funct = @lift(get_function($time))
x = range(1,SIZE_X)
y = range(1,SIZE_Y)
fig = Figure(resolution=(SIZE_X, SIZE_Y))
ax = Axis(fig[1,1], aspect=DataAspect(), title="2D Terrain Demo (Julia)", limits=(1, SIZE_X, 1, SIZE_Y))
hm = heatmap!(ax, x, y, funct, colormap=:coolwarm, colorscale=colorscale, colorrange=(-0.01, 0.01))

Colorbar(fig[:, end+1], hm, width=10)
rowsize!(fig.layout, 1, size(ax.scene)[2])
resize_to_layout!(fig)
lines!(ax, terrain_heights, color="black", label=false)

time[] = 1
#display(fig)

println("Animating...")
record(fig, "time_animation.mp4", 1:SIZE_Z;
        framerate = 12) do t
    if t % 5 == 0
        println("t = $t")
    end
    time[] = t
end

println("Done!")

