from video import *
import sys

header_color='#ffaa00'
slide("title.png",48,[(50,"A Novel Mass Spring",header_color),
                      (50,"Model for Simulating",header_color),
                      (50,"Full Hair Geometry",header_color),
                      " ",
                      (36,"paperid 0384"),
                      (36,"SIGGRAPH 2008")])


slide("tet_springs.png",36,
      [(54,"Altitude Springs",header_color),
       " ",
       "Our novel altitude",
       "springs allow a volumetric",
       "torus to recover from",
       "full collapse"])

slide("tet_springs_shell.png",36,
      [(54,"Altitude Springs",header_color),
       " ",
       "Our novel altitude",
       "springs also allow a",
       "thin shell to recover from",
       "full collapse"])

slide("model.png",36,
      [(54,"Hair Model",header_color),
       " ",
       "Our new hair model",
       "reproduces real hair",
       "phenomena"])

slide("curly.png",36,
      [(54,"Hair Model",header_color),
       " ",
       "Curliness can be varied by",
       "by adjusting spring stiffness"])

slide("cloth_frame_rate.png",36,
      [(54,"Time Integration",header_color),
       " ",
       "Our linear implicit springs",
       "allow simulating cloth using",
       "only one time step per frame"])

slide("cloth_cfl_10.png",36,
      [(54,"Time Integration",header_color),
       " ",
       "...of course smaller",
       "time steps yield",
       "more accuracy"])

slide("levelset_interpolation.png",36,
      [(54,"Time Integration",header_color),
       " ",
       "Our new level set",
       "interpolation scheme",
       "allows more accurate",
       "collision bodies"])

slide("tuft_old.png",36,
      [(54,"Hair Tuft",header_color),
       " ",
       "Basic mass-spring model",
       "with dynamic head motion",
       " ",
       "(too much stretching)"])

slide("tuft_new.png",36,
      [(54,"Hair Tuft",header_color),
       " ",
       "Our mass-spring model",
       "better handles the",
       "complex motion",
       " ",
       "(better length preservation)"])

slide("tuft_new_interact.png",36,
      [(54,"Hair Tuft",header_color),
       " ",
       "Our model with stiction",
       "and self-collisions"])

slide("shaking5k.png",36,
      [(54,"Full head",header_color),
       " ",
      "5,000 long straight hairs",
       " ",
       "(coverage too thin)"])

slide("shaking.png",36,
      [(54,"Full head",header_color),
       " ",
      "10,000 long straight hairs",
       " ",
       "(better coverage)"])

slide("spinning.png",36,
      [(54,"Full head",header_color),
       " ",
      "5,000 long curly hairs"])

slide("wind.png",36,
      [(54,"Full head",header_color),
       " ",
      "10,000 medium length",
       "straight hairs",
      "undergoing wind forces"])

slide("shrek.png",36,
      [(54,"Full head",header_color),
       " ",
      "10,000 medium length",
       "straight hairs"])


slide("half.png",36,
      ["half speed"])


slide("end.png",54,["The End"])

final_render_base="/solver/vol3/hair7/final_render_new"
#final_render_base_20080121="/solver/vol3/hair7/render/Projects/renderman/output_final_20080121/generic_hair"
final_render_base_20080122="/solver/vol3/hair7/render/Projects/renderman/output_final_20080122/generic_hair"


video_dir="hair_video_tmp"
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
testing_titles=False
skip_interval=1


def caption_length(time=3.5):
    if not using_lengths and testing_titles:
        return 1
    else:
        return int(video.fps*time)

video.add_frame("title.png",caption_length(4./skip_interval))

# tet springs demos
video.add_frame("tet_springs.png",caption_length(9./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/torus_volume",0,170,step=skip_interval)

video.add_frame("tet_springs_shell.png",caption_length(9./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/torus_shell",0,240,step=skip_interval)

video.add_frame("model.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/strand/merge",step=skip_interval)

video.add_frame("curly.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/strand/curly_drop_composite",step=skip_interval)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory("/solver/vol3/hair7/final_render/strand/curly_drop_composite",step=skip_interval,duplicate=2)


video.add_frame("cloth_frame_rate.png",caption_length(11./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/cloth_cfl10000",0,158,step=skip_interval)

# implicit spring demos
video.add_frame("cloth_cfl_10.png",caption_length(6./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/cloth_cfl10",0,158,step=skip_interval)

video.add_frame("levelset_interpolation.png",caption_length(9./skip_interval))
if not testing_titles:
    video.add_directory("/solver/vol3/hair7/final_render/interpolation/merged",1,60,step=skip_interval)
    video.add_directory("/solver/vol3/hair7/final_render/interpolation/merged",1,60,step=skip_interval)

# tuft examples
video.add_frame("tuft_old.png",caption_length(8./skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_tuft_back.old-blur_render"),0,188,step=1*skip_interval)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_tuft_back.old_render"),0,379,step=1*skip_interval)
video.add_frame("tuft_new.png",caption_length(8./skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_tuft_back.nocoladh-blur_render"),0,188,step=1*skip_interval)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_tuft_back.nocoladh_render"),0,379,step=1*skip_interval)
video.add_frame("tuft_new_interact.png",caption_length(6./skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_tuft_back.blur_render"),0,188,step=1*skip_interval)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_tuft_back_render"),0,379,step=1*skip_interval)

# full head examples
video.add_frame("shaking5k.png",caption_length(8./skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_head_front.5k-blur_render"),0,109,step=1*skip_interval) # TODO: make this whatever we get to on long (119 natural) (99 old)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_head_front.5k_render"),0,220,step=1*skip_interval)  # TODO: make this whatever we get to on long (239 natural) (199 old)

video.add_frame("shaking.png",caption_length(5.5/skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_head_front.blur_render"),0,109,step=1*skip_interval) # TODO: make this whatever we get to on long
    video.add_frame("half.png",caption_length(2.5/skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shaking_head_front_render"),0,220,step=1*skip_interval) # TODO: make this whatever we get to on long

video.add_frame("spinning.png",caption_length(5./skip_interval))
if not testing_titles:
   video.add_directory(os.path.join(final_render_base_20080122,"spinning_head_camera_lower.new-blur_render"),0,43,step=1*skip_interval)
   video.add_frame("half.png",caption_length(2./skip_interval))
   video.add_directory(os.path.join(final_render_base_20080122,"spinning_head_camera_lower.new_render"),0,86,step=1*skip_interval)

video.add_frame("wind.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shrek_head_back.blur_render"),0,144,step=1*skip_interval)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shrek_head_back_render"),0,289,step=1*skip_interval)

video.add_frame("shrek.png",caption_length(5./skip_interval))
if not testing_titles:
    video.add_directory(os.path.join(final_render_base_20080122,"shrek_head.blur_render"),0,102,step=1*skip_interval)
    video.add_frame("half.png",caption_length(2./skip_interval))
    video.add_directory(os.path.join(final_render_base_20080122,"shrek_head_render"),0,207,step=1*skip_interval)

video.add_frame("end.png",caption_length(2./skip_interval))

video.make_movie('hair')


