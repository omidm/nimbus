from video import *
import sys

header_color='#A9A9A9'
slide("title.png",48,[(50,"Mass and Momentum",header_color),
                      (50,"Conservation for",header_color),
                      (50,"Fluid Simulation",header_color),
                      " ",
                      (36,"paperid 0045"),
                      (36,"SIGGRAPH 2011")])


slide("smoke.png",36,
      [(54,"Smoke",header_color),
       " ",
       "Smoke injected from below"])

slide("smoke_low_cfl.png",36,
      [(54,"Smoke",header_color),
       " ",
       "Comparison using a typical",
       "time step size",
       "Resolution: 128x256x128"])

slide("smoke_high_cfl.png",36,
      [(54,"Smoke",header_color),
       " ",
       "Comparison using a large",
       "time step size",
       "Resolution: 128x256x128"])

slide("compare.png",36,
      ["  Standard                      Our",
      "Semi-Lagrangian           Method"])

slide("full.png",36,
      ["Full Comparison"])

slide("smoke_high_cfl_256.png",36,
      [(54,"Smoke",header_color),
       " ",
       "Comparison using a large",
       "time step size",
       "Resolution: 256x512x256"])

slide("smoke_sphere.png",36,
      [(54,"Smoke with Sphere",header_color),
       " ",
       "Smoke injected from below",
       "with a sphere"])

slide("smoke_low_cfl_sphere.png",36,
      [(54,"Smoke with Sphere",header_color),
       " ",
       "Comparison using a typical",
       "time step size",
       "Resolution: 128x256x128"])

slide("smoke_high_cfl_sphere.png",36,
      [(54,"Smoke with Sphere",header_color),
       " ",
       "Comparison using a large",
       "time step size",
       "Resolution: 128x256x128"])

slide("smoke_high_cfl_256_sphere.png",36,
      [(54,"Smoke with Sphere",header_color),
       " ",
       "Comparison using a large",
       "time step size",
       "Resolution: 256x512x256"])

slide("water.png",36,
      [(54,"Water",header_color),
       " ",
       "A water sphere dropped into a box",
       "Resolution: 128x256x128"])

slide("energy.png",36,
      [(54,"Energy Conservation",header_color),
       " "])

slide("fish.png",36,
      [(54,"Energy Conservation",header_color),
       " ",
       "A closed box simulating",
       "one-phase flow with passively",
       "advected fish particles"])

slide("bees.png",36,
      [(54,"Energy Conservation",header_color),
       " ",
       "A free surface flow with",
       "passively advected",
       "insect particles"])

slide("smoke_energy.png",36,
      [(54,"Energy Conservation",header_color),
       " ",
       "Smoke injected from below"])

slide("smoke_energy_sphere.png",36,
      [(54,"Energy Conservation",header_color),
       " ",
       "Smoke injected from below",
       "with a sphere"])

final_render_base="/n/field/data/render_backup_2011"

video_dir="conservative_video"
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
hi_res=1

def caption_length(time=3.5):
    if not using_lengths and testing_titles:
        return 1
    else:
        return int(video.fps*time)

video.add_frame("title.png",caption_length(5./skip_interval))

video.add_frame("smoke.png",caption_length(4./skip_interval))
video.add_frame("smoke_low_cfl.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(3./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t1_128_low_sl_low_cons_pngs",1,200,step=skip_interval)
    #video.add_directory(final_render_base+"/smoke_t1_128_low_sl_low_cons_pngs",1,140,step=skip_interval)
    #video.add_frame(final_render_base+"/smoke_t1_128_low_sl_low_cons_pngs/frame.00140.png",1);
    #video.add_blending(final_render_base+"/smoke_t1_128_low_sl_low_cons_pngs/frame.00140.png",final_render_base+"/labeled/128_1_frame.00140_final.png",int(1.*video.fps))
    #video.add_frame(final_render_base+"/labeled/128_1_frame.00140_final.png",caption_length(5./skip_interval));
    #video.add_directory(final_render_base+"/smoke_t1_128_low_sl_low_cons_pngs",141,200,step=skip_interval)

video.add_frame("smoke_high_cfl.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(3./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t1_128_high_sl_high_cons_pngs",1,103,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t1_128_high_sl_high_cons_pngs/frame.00103.png",1);
    video.add_blending(final_render_base+"/smoke_t1_128_high_sl_high_cons_pngs/frame.00103.png",final_render_base+"/labeled/128_1_frame.00103_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/128_1_frame.00103_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t1_128_high_sl_high_cons_pngs",104,200,step=skip_interval)

video.add_frame("smoke_high_cfl_256.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(3./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs",1,83,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs/frame.00083.png",1);
    video.add_blending(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs/frame.00083.png",final_render_base+"/labeled/256_2_frame.00083_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/256_2_frame.00083_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs",84,114,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs/frame.00114.png",1);
    video.add_blending(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs/frame.00114.png",final_render_base+"/labeled/256_2_frame.00114_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/256_2_frame.00114_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs",115,200,step=skip_interval)

video.add_frame("full.png",caption_length(1./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t1_256_high_sl_high_cons_pngs",1,200,step=skip_interval)

video.add_frame("smoke_sphere.png",caption_length(4./skip_interval))
video.add_frame("smoke_low_cfl_sphere.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(3./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t2_128_low_sl_low_cons_pngs",1,143,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t2_128_low_sl_low_cons_pngs/frame.00143.png",1);
    video.add_blending(final_render_base+"/smoke_t2_128_low_sl_low_cons_pngs/frame.00143.png",final_render_base+"/labeled/128_1_frame.00143_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/128_1_frame.00143_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t2_128_low_sl_low_cons_pngs",144,200,step=skip_interval)

video.add_frame("smoke_high_cfl_sphere.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(3./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t2_128_high_sl_high_cons_pngs",1,152,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t2_128_high_sl_high_cons_pngs/frame.00152.png",1);
    video.add_blending(final_render_base+"/smoke_t2_128_high_sl_high_cons_pngs/frame.00152.png",final_render_base+"/labeled/128_1_frame.00152_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/128_1_frame.00152_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t2_128_high_sl_high_cons_pngs",153,200,step=skip_interval)

video.add_frame("smoke_high_cfl_256_sphere.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(3./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs",1,83,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs/frame.00083.png",1);
    video.add_blending(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs/frame.00083.png",final_render_base+"/labeled/256_1_frame.00083_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/256_1_frame.00083_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs",83,123,step=skip_interval)
    video.add_frame(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs/frame.00123.png",1);
    video.add_blending(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs/frame.00123.png",final_render_base+"/labeled/256_1_frame.00123_final.png",int(1.*video.fps))
    video.add_frame(final_render_base+"/labeled/256_1_frame.00123_final.png",caption_length(5./skip_interval));
    video.add_directory(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs",124,150,step=skip_interval)

video.add_frame("full.png",caption_length(1./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke_t2_256_high_sl_high_cons_pngs",1,150,step=skip_interval)

video.add_frame("water.png",caption_length(5./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/pls/Test_3_128_3_2.70_0.00_conservative_10_1_new/pics",1,177,step=skip_interval)

video.add_frame("energy.png",caption_length(3./skip_interval))
video.add_frame("fish.png",caption_length(7./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/fish/old_data/J8/pics",1,1000,step=skip_interval)

video.add_frame("bees.png",caption_length(6./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/bees/pics",1,1000,step=skip_interval)

video.add_frame("smoke_energy.png",caption_length(4./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke9/Test_1_128_1_vc_0.00_clamped_3d_3.00_0.90_9_conservative_energy/pics/",1,40,step=skip_interval)

video.add_frame("smoke_energy_sphere.png",caption_length(4./skip_interval))
if not testing_titles:
    video.add_directory(final_render_base+"/smoke9/Test_2_128_1_vc_0.00_clamped_3d_3.00_0.90_9_conservative_energy/pics/",1,54,step=skip_interval)

#video.add_frame("end.png",caption_length(2./skip_interval))
video.make_movie('conservative')


