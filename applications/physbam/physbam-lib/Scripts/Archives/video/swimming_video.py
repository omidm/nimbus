from video import *
import sys

header_color='#ffaa00'
slide("title.png",48,[(50,"Creature Control in",header_color),
                      (50,"a Fluid Environment",header_color),
                      " ",
                      (36,"paperid 0143"),
                      (36,"SIGGRAPH 2009")])


slide("wind_blocks.png",36,
      [(54,"Block Examples",header_color),
       " ",
       "Two blocks connected with",
       "a single joint optimizing",
       "various objective functions",
       "(simple fluids)"])

slide("fluid_blocks.png",36,
      [(54,"Block Examples",header_color),
       " ",
       "Two blocks connected with",
       "a single joint optimizing",
       "various objective functions",
       "(Navier-Stokes fluids)"])

slide("human-back-effort.png",36,
      [(54,"Human",header_color),
       " ",
       "A human with three controlled",
       "joints optimizes effort",
       "on the lower back"])

slide("human-effort.png",36,
      [(54,"Human",header_color),
       " ",
       "A human with three controlled",
       "joints optimizes effort",
       "on the lower back and arms"])

slide("human-wind-back-effort.png",36,
      [(54,"Human",header_color),
       " ",
       "A human with three controlled",
       "joints optimizes effort",
       "on the lower back"])

slide("human-wind-effort.png",36,
      [(54,"Human",header_color),
       " ",
       "A human with three controlled",
       "joints optimizes effort",
       "on the lower back and arms"])

slide("human-wind-drag.png",36,
      [(54,"Human",header_color),
       " ",
       "A human with three controlled",
       "joints optimizes effort",
       "on the lower back",
       "and drag on the arms"])

slide("bird-block.png",36,
      [(54,"Driven Flyer",header_color),
       " ",
       "A driven flyer with two",
       "controlled joints alternates",
       "maximizing and minimizing",
       "drag"])

slide("bird-mesh.png",36,
      [(54,"Driven Flyer",header_color),
       " ",
       "A driven flyer with two",
       "controlled joints alternates",
       "maximizing and minimizing",
       "drag (with flesh)"])

slide("octosquid-close.png",36,
      [(54,"Swimmer",header_color),
       " ",
       "A swimmer with five controlled",
       "joints alternates maximizing",
       "and minimizing drag, close"])

slide("octosquid-farther.png",36,
      [(54,"Swimmer",header_color),
       " ",
       "A swimmer with five controlled",
       "joints alternates maximizing",
       "and minimizing drag, far"])

slide("half.png",36,
      ["half speed"])

slide("end.png",54,["The End"])

final_render_base="/solver/vol0/swimming4/final_video"

video_dir="swimming_video_tmp"
if not os.path.exists(video_dir):
    os.mkdir(video_dir)
#if os.path.isdir(video_dir):
#    print "%s already exists... delete? (y/n) "%video_dir,
#    c=sys.stdin.readline()
#    if c=="y\n": shutil.rmtree(video_dir)
#    else: sys.exit(1)
shutil.rmtree(video_dir)    
video=VIDEO(video_dir)

using_lengths=True
testing_titles=True
skip_interval=1


def caption_length(time=3.5):
    if not using_lengths and testing_titles:
        return 1
    else:
        return int(video.fps*time)

video.add_frame("title.png",caption_length(4./skip_interval))

# blocks demos
video.add_frame("wind_blocks.png",caption_length(8./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/wind_final_render",step=skip_interval)

video.add_frame("fluid_blocks.png",caption_length(8./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/fluids_final_render",step=skip_interval)

#humans
video.add_frame("human-back-effort.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/human_no_smoke_back_render",0,274,step=skip_interval)

video.add_frame("human-effort.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/human_no_smoke_min_effort_col_render",0,274,step=skip_interval)

video.add_frame("human-wind-back-effort.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/wind_back_final_render",124,399,step=skip_interval)

video.add_frame("human-wind-effort.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/wind_effort_final_render",124,399,step=skip_interval)

video.add_frame("human-wind-drag.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming2/final_input_videos/wind_drag_final_render",124,399,step=skip_interval)

video.add_frame("bird-block.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming9/render_backup/bird/bird_bones_mike/bird_flying_bones_Smoke_Postprocess_Test_3_Resolution_15_final_render",step=skip_interval)

video.add_frame("bird-mesh.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming9/render_backup/bird/bird_mesh_mike/bird_flying_mesh_Smoke_Postprocess_Test_2_Resolution_15_final_high_render",step=2*skip_interval)
    #video.add_frame("half.png",caption_length(2./skip_interval))
    #video.add_directory("/solver/vol0/swimming9/render_backup/bird/bird_mesh_mike/bird_flying_mesh_Smoke_Postprocess_Test_2_Resolution_15_final_high_render",step=skip_interval,duplicate=2)

video.add_frame("octosquid-close.png",caption_length(6./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming9/render_backup/octosquid/octosquid_mike_2/octosquid_final_close_render",step=2*skip_interval)
    #video.add_frame("half.png",caption_length(2./skip_interval))
    #video.add_directory("/solver/vol0/swimming9/render_backup/octosquid/octosquid_mike_2/octosquid_final_close_render",step=skip_interval,duplicate=2)

video.add_frame("octosquid-farther.png",caption_length(6./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol0/swimming9/render_backup/octosquid/octosquid_mike/octosquid_final_far_render",step=2*skip_interval)
    #video.add_frame("half.png",caption_length(2./skip_interval))
    #video.add_directory("/solver/vol0/swimming9/render_backup/octosquid/octosquid_mike/octosquid_final_far_render",step=skip_interval,duplicate=2)


video.add_frame("end.png",caption_length(2./skip_interval))

video.make_movie('swimming')


