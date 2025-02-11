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

SAMPLE_PERIOD = 8

funct = @lift(get_function($time))
x = range(1,1001)
y = range(1,501)
fig = Figure(resolution=(1001, 501))
ax = Axis(fig[1,1], aspect=DataAspect(), title="2D FDTD with Terrain", limits=(1, 1001, 1, 501), xlabel=@lift string("Step ", ($time[1]-1)*SAMPLE_PERIOD))
hm = heatmap!(ax, x, y, funct, colormap=:coolwarm, colorscale=colorscale, colorrange=(-0.01, 0.01))

Colorbar(fig[:, end+1], hm, width=10)
rowsize!(fig.layout, 1, size(ax.scene)[2])
resize_to_layout!(fig)
ln = lines!(ax, terrain_heights, color="black")

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