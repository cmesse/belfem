from corc import Cable, StraightLine, Lame, Trefoil

C = Cable()

C.tapeWidth = 4

# pitch = z-advance per radian. 3 gives a 35.6 deg lay angle from the
# cable axis on layer 0 (tan = r0/pitch, r0 = 2.149 from the packing
# perimeter) and 18.8 mm pitch length per turn -- the shallow end of
# the commercial CORC winding range. The realistic 45 deg lay
# (pitch 2.15) is NOT currently meshable: gmsh's volume pass diverged
# past the 1 h assembler timeout even in the proven 200 um / 0.2 mm
# regime (131k+ nodes against 74k for pitch 3, worst-tet-radius spikes;
# the counter-wound layers then cross at 90 deg instead of 71 deg).
# Steepening the lay needs its own meshing campaign.
C.pitch = 3
# 4/3 turns: numTurns * numTapesPerLayer = 4 is whole, so the caps close
# with zero rotation and every tape's far end is the NEXT tape's near end
# -- the periodic twist (sideset 4 -> 5 -> 6 -> 4, inner layer likewise).
# The model is then one period of an infinite cable whose three tapes per
# layer are one tape of three times the length; a defect's periodic images
# sit 3 L = 75 mm apart along the tape. 1/3 would also close but puts the
# images 19 mm apart
C.numTurns = 4.0 / 3.0

C.gap = 0.5

# in mu. tapeThickness is the radial slot per layer; solderThickness the
# solder bed between layers. 200 um is generous against a real soldered
# cable, but it is the meshable floor with this gmsh pipeline: thinner
# beds (50 um at res 0.15, and 100 um at res 0.15 = 0.75 x gap) both ran
# gmsh's volume pass past the 1 h assembler timeout with decaying
# insertion rates -- shrinking the annulus needs its own meshing
# campaign, not a parameter tweak.
# stackThickness is the physical tape stack (the BELFEM layers sum
# 20 + 2 + 1.6 + 0.2 + 50 + 1.8 + 20) used for the density correction.
C.tapeThickness = 100
C.solderThickness = 120
C.stackThickness = 95.6

# in mm
C.domainRadius = 60

C.numTapesPerLayer = 3
C.numLayers = 2
C.delta = 4

# resolutions in mm. tapeResolution tracks the inter-layer gap
# (tapeThickness + solderThickness = 300 um): Spike D wants it within
# [0.5, 2] x the gap; 0.67 x the gap (0.2 / 0.3) is the proven regime.
C.tapeResolution = 0.22
C.innerResolution = 0.25
C.domainResolution = 4

# centerline: straight (default), or e.g.
#C.centerline = Lame(76.1, 25, 1.55)      # coil end (use domainRadius < 26)
#C.centerline = Trefoil(8.718490620126541)  # show-off, shells only
C.centerline = StraightLine()

C.build()
C.print()

# Relative, so both land in the directory the script is run from -- which is
# where ../input.conf looks for them (`mesh { file : corc.msh ; }`). Only the
# gmsh scratch files stay in Cable.workdir.
C.save("corc.msh")
C.write_topology("corc_topology.conf")
