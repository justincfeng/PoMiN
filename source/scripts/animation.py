############################################################
# animation.py - create interactive 3d animation of pmcode(1)'s output
############################################################
# This file creates an interactive 3D animation (or a static plot) of
# any given number of <output-file>'s from pmcode(1); run
# "animation.py -h" for usage. Plotting a great number of objects
# should only be limited to the user's screen size.
############################################################
# for-loops are over f=file & p=particle & a=axis (x,y,z)
#!/usr/bin/python

############################################################
# IMPORTS
############################################################
import numpy as np                           # for data calculations
import pandas as pd                          # for manipulating csv files
import matplotlib.pyplot as plt              # for plotting data
import matplotlib.animation as ani           # for adding animation
import mpl_toolkits.mplot3d.axes3d as p3     # for creating a 3d plot
import argparse                              # for commandline arguments

############################################################
# FUNCTIONS
############################################################
def update_lines(num, dataLines, lines):
    if not pause:
        for f in range(len(args.inputfile)):
            for line, data in zip(lines[f], dataLines[f]):
                # NOTE: there is no .set_data() for 3 dim data...
                line.set_data(data[0:2, :num])
                line.set_3d_properties(data[2, :num])
                time_text[f].set_text("time ("+str(f+1)+") = "+str(inputfile[f]['time'][num]))
        return lines

def on_click_pause(event):      # TODO: make function actually pause animation.
    # Toggle pause effect if the spacebar is pressed.
    if event.key ==" ":
        global pause
        pause ^= True

############################################################
# COMMAND LINE ARGUMENTS
############################################################
# Create an object to accept arguments.
parser = argparse.ArgumentParser(description='contruct an animated 3D plot from the output csv files of pmcode(1)')

parser.add_argument("inputfile", nargs='*', default=["../output.csv"],
                    help="input csv file which is known as <output-file>.csv to pmcode(1) [default: ../output.csv]")
parser.add_argument("-m", "--momentum", action="store_true",
                    help="plot momentum instead of position")
parser.add_argument("-s", "--skip",action="store_true",
                    help="skip animated plotting sequence; ignores {-l,-i} flags")
parser.add_argument("-p", "--printfig", nargs=1,
                    help="print plot to file given by argument name (with extension!); assumes {-s} flag")
parser.add_argument("-l", "--loop", action="store_true",
                    help="loop repeatedly through the <inputfile>")
parser.add_argument("-c", "--cubic", action="store_true",
                    help="force cubic plot ranges")
parser.add_argument("-e", "--every", type=int, default=1,
                    help="plot every x points [default: 1]")
parser.add_argument("-i", "--interval", type=float, default=50,
                    help="interval between frames in the animation in milliseconds [default: 50]")
parser.add_argument("-t", "--title", default="Post-Minkowsikian Few-Body Code",
                    help="give plot a specific title [default: Post-Minkowskian Few-Body Code]")
parser.add_argument("-r", "--range", type=float, nargs=6, default=[0.0,0.0,0.0,0.0,0.0,0.0],
                    help="{BROKEN}set range of x,y,z axes as \"xmin xmax ymin ymax zmin zmax\" [default: 0 0 0 0 0 0] (automatic)")
#parser.add_argument("-t", "--time", action="store_true",
#                    help="plot different files with syncronized time")

# TODO: make flags for user-defined plot ranges

# Define a shorthand way to call the arguments.
args = parser.parse_args()

############################################################
# SETUP
############################################################
# Create figure with three dimensions.
fig = plt.figure()
ax = p3.Axes3D(fig)

# Assign title for the plot.
ax.set_title(args.title)

# Pause when a key is pressed on figure.
fig.canvas.mpl_connect('key_press_event', on_click_pause)
# Start in a non-paused state.
pause = False

# Read <inputfile> into array for manipulation. TODO: get rid of pandas dependeny?
inputfile = [pd.read_csv(args.inputfile[f], sep=',', header=0) for f in range(len(args.inputfile))]

# Derive number of particles from the floor divisor ("//").
# Shape is an array that gives the number of [1=columns, 2=rows].
numparticles = [(inputfile[f].shape[1] // 7) for f in range(len(args.inputfile))]

# Create array of either position or momentum from csv.
# FORMAT: data[file][particle][x,y,z][time]
if (args.momentum):
    data = [np.array([[inputfile[f]['px_'+str(p)][::args.every].values, inputfile[f]['py_'+str(p)][::args.every].values, \
                       inputfile[f]['pz_'+str(p)][::args.every].values]
                      for p in range(1,numparticles[f]+1)]) for f in range(len(args.inputfile))]
    ax.set_xlabel('px'); ax.set_ylabel('py'); ax.set_zlabel('pz')
else:
    data = [np.array([[inputfile[f]['qx_'+str(p)][::args.every].values, inputfile[f]['qy_'+str(p)][::args.every].values, \
                       inputfile[f]['qz_'+str(p)][::args.every].values]
                      for p in range(1,numparticles[f]+1)]) for f in range(len(args.inputfile))]
    ax.set_xlabel('qx'); ax.set_ylabel('qy'); ax.set_zlabel('qz')

############################################################
# PLOT
############################################################
# Initiate range variables.
# FORMAT: args.range=[xmin,xmax,ymin,ymax,zmin,zmax]
for f in range(len(args.inputfile)):
    # Find the best range to plot over.
    for p in range(numparticles[f]):
        for i,a in zip([0,2,4],[0,1,2]):
            if args.range[i] is not None: args.range[i] = min(min(data[f][p][a]),args.range[i])
            if args.range[i+1] is not None: args.range[i+1] = max(max(data[f][p][a]),args.range[i+1])
            
            # Define object to print the time from each file.
            time_text = [fig.text(0.805, fi*0.02, "") for fi in reversed(range(numparticles[f]))]
                
# Set fixed values of plot ranges based on previous calculation.
if args.cubic:
    range_max = max(abs(args.range[0]),abs(args.range[1]),abs(args.range[2]),\
                    abs(args.range[3]),abs(args.range[4]),abs(args.range[5]))
    ax.set_xlim3d([-range_max, range_max])
    ax.set_ylim3d([-range_max, range_max])
    ax.set_zlim3d([-range_max, range_max])
else:
    ax.set_xlim3d([args.range[0], args.range[1]])
    ax.set_ylim3d([args.range[2], args.range[3]])
    ax.set_zlim3d([args.range[4], args.range[5]])

# Create line segments to plot
lines = [[ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1], label=str(idx))[0]
          for idx, dat in enumerate(da, start=1)] for da in data]

# TODO: test for data not equal in time (i.e., the amount of input rows)
if args.skip:
    # Skip the animation sequence altogether.
    [[ax.plot(data[f][p][0], data[f][p][1], data[f][p][2])[0]
      for p in range(numparticles[f])] for f in range(len(args.inputfile))]
elif args.loop:
    # Repeat the animation ad infinitum.
    line_ani = [ani.FuncAnimation(fig, update_lines, data[0].shape[2], fargs=(data, lines),
                                  interval=args.interval, blit=False)]
else:
    # Otherwise, don't repeat or skip the animation (by default).
    line_ani = [ani.FuncAnimation(fig, update_lines, data[0].shape[2], fargs=(data, lines),
                                  interval=args.interval, blit=False, repeat=False)]


# Place a legend to the right of the plot.
[plt.gca().add_artist(ax.legend(handles=lines[f], loc=2, ncol=4, mode="expand",
                                title=(("..."+str(args.inputfile[f])[-35:])
                                       # Truncate long file names with leading ellipses.
                                       if len(str(args.inputfile[f])) > 35 else str(args.inputfile[f])),
                                bbox_to_anchor=(1.0, 1.2-f*(0.04*numparticles[f-1]-0.006), 0.25, -0.2)))
 for f in range(len(args.inputfile))]

# TODO: make plot legend turn gray when complete?
#ax.legend().get_frame().set_facecolor("gray")

# Shrink plot horizontally by 20%.
ax.set_position([ax.get_position().x0, ax.get_position().y0,
                 ax.get_position().width * 0.8, ax.get_position().height])

if (args.printfig): 
    plt.savefig(args.printfig[0])
else:
    plt.show()

############################################################
# END animation.py
############################################################
