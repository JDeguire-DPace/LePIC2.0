# rotate.py — rotation 3D animation in VisIt

# 1) Load your dataset
OpenDatabase("../Outputs/state.h5")  # change to your file

AddPlot("Pseudocolor", "phi")  # change "density" to your variable
DrawPlots()

# 2) Camera + plot options (optional)
p = PseudocolorAttributes()
p.opacityType = p.Constant
p.opacity = 0.7
SetPlotOptions(p)

View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0, 0, 1)
View3DAtts.viewUp = (0, 1, 0)
View3DAtts.imageZoom = 1.2
SetView3D(View3DAtts)

# 3) Rotation loop
nframes = 180
for i in range(nframes):
    # Rotate around Z axis
    angle = 2 * 3.1415926 * float(i) / nframes
    View3DAtts.viewNormal = (math.sin(angle), math.cos(angle), 0.3)
    SetView3D(View3DAtts)

    SaveWindowAttributes().family = 0
    s = SaveWindowAttributes()
    s.format = s.PNG
    s.fileName = f"frame_{i:03d}"
    SetSaveWindowAttributes(s)
    SaveWindow()

# 4) After script runs, assemble frames
print("\n✅ Frames saved — build the video with:")
print("ffmpeg -r 30 -i frame_%03d.png -vcodec libx264 rotation.mp4")
